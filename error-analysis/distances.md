# Error analysis of `dist_sqr_r_s_t`

This bounds the rounding error of each `dist_sqr_r_s_t` method in
[`src/distances.jl`](../src/distances.jl). The bounds compare the computed
results with the exact results for the *same* floating point inputs. The
inputs are arbitrary: no assumption is made that `x[1]` is zero, or that
the coordinates are relative to anything.

Point-point and point-line are derived by hand. Point-triangle and line-line
are traced with [`ErrorBounds.jl`](ErrorBounds.jl), which runs the kernels
themselves on symbolic error bounds, described in
[Tracing the kernels](#tracing-the-kernels). [`analyze.jl`](analyze.jl)
prints the traced bounds and checks them against `BigFloat`.

## Setting

- `T` is `Float32` or `Float64`, rounding to nearest even. The unit
  roundoff is `u = eps(T)/2`, `2^-24` for `Float32` and `2^-53` for `Float64`.
- Standard model: for `+`, `-`, `*` and `fma`,
  `fl(a op b) = (a op b)(1 + δ)` with `|δ| ≤ u`, as long as the result
  does not overflow. If the exact result is in the subnormal range, `+` and
  `-` are exact, and `*` and `fma` have an absolute error of at most
  `η = nextfloat(zero(T))/2` instead.
- `γₖ = k*u/(1 - k*u)`, so `(1 + u)^k - 1 ≤ γₖ`.
- Julia's `fma` always rounds once, in scalar code and in each lane of a
  `SIMD.Vec`. The analysis relies on this; `muladd` may or may not fuse, so
  a kernel that used it would have two different error bounds.
- `dist_sqr` returns the `d²` of `dist_sqr_r_s_t` bit for bit, and the SIMD
  methods are bit identical to the scalar ones in each lane, so every bound
  here applies to all of them.
- When a result can be bounded relative to its own size, it is. Otherwise
  the bound is in units of `u*ℓ`, where `ℓ` is the largest absolute input
  coordinate.

## Point-point

```julia
d = y[1] - x[1]
dot3(d, d), d, SVector{0,T}(), SVector{0,T}()
```

with `dot3(a, b) = fma(a[3], b[3], fma(a[2], b[2], a[1]*b[1]))`.

Write `r = y[1] - x[1]` for the exact separation and `D = r⋅r` for the exact
squared distance. Hats mark computed values.

### The separation `r`

Each component is one correctly rounded subtraction:

    r̂ᵢ = rᵢ(1 + δᵢ),    |δᵢ| ≤ u.

So `|r̂ᵢ - rᵢ| ≤ u|rᵢ|` and `‖r̂ - r‖ ≤ u‖r‖`, a relative bound however
large the coordinates are. By Sterbenz's lemma the subtraction is exact
(`δᵢ = 0`) when `yᵢ/2 ≤ xᵢ ≤ 2yᵢ`, which covers points close together far
from the origin.

### The squared distance `d²`

`dot3` rounds `r̂₁²` once, then each `fma` once:

    p₁ = r̂₁²(1 + ε₁)
    p₂ = (r̂₂² + p₁)(1 + ε₂)
    d̂² = (r̂₃² + p₂)(1 + ε₃)
       = r̂₃²(1 + ε₃) + r̂₂²(1 + ε₂)(1 + ε₃) + r̂₁²(1 + ε₁)(1 + ε₂)(1 + ε₃)

with `|εⱼ| ≤ u`. Every term is nonnegative, so nothing cancels, and each
`r̂ᵢ² = rᵢ²(1 + δᵢ)²` adds two more factors. The first term has the most,
five, so

    (1 - u)⁵ D ≤ d̂² ≤ (1 + u)⁵ D,    |d̂² - D| ≤ γ₅ D ≈ 5u D.

In units of `eps(T)` that is `2.5 eps(T)` relative. For the distance,
`sqrt(d̂²)` is within about `2.5u` of `sqrt(D)` relative, plus `u` for the
rounding of the square root itself.

The bound is relative, so it holds at any `ℓ` and implies any bound in
units of `u*ℓ`.

### Backward error

`r̂` is the exact separation from `x[1]` to a point `y′` with
`|y′ᵢ - yᵢ| ≤ u|rᵢ|`. Then `d̂² = ‖r̂‖²(1 + θ)` with `|θ| ≤ γ₃`, which is the
exact squared distance after scaling `r̂` by `sqrt(1 + θ)`. So `d̂²` is the
exact squared distance from `x[1]` to a point within about `2.5u‖r‖` of
`y[1]`.

### Underflow and overflow

- **Overflow** needs `|yᵢ - xᵢ| > floatmax(T)`, or `‖r‖` above about
  `sqrt(floatmax(T))`: `1.8e19` for `Float32`, `1.3e154` for `Float64`.
- **Underflow** happens when a product falls below `floatmin(T)`, for
  example `|rᵢ| < 1.1e-19` in `Float32`. Then a product or `fma` has an
  absolute error of up to `η` instead of a relative one, and the bound
  becomes `|d̂² - D| ≤ γ₅ D + 3(1 + u)²η`. In `Float32`, `3η` is about
  `4.2e-45`.

### Measured

Maximum error against `BigFloat` on the same inputs, 200,000 pairs per
generator:

| `T` | generator | `d²`, relative | `r`, componentwise relative |
|---|---|---|---|
| `Float32` | uniform in `[0, 100]³` | 4.06u | 1.0u |
| `Float32` | within 0.01 of each other, up to 1e4 from the origin | 1.39u | 0 |
| `Float32` | coordinates from 1e-3 to 1e3, either sign | 4.31u | 1.0u |
| `Float64` | uniform in `[0, 100]³` | 4.06u | 1.0u |
| `Float64` | within 0.01 of each other, up to 1e4 from the origin | 2.51u | 0 |
| `Float64` | coordinates from 1e-3 to 1e3, either sign | 4.19u | 1.0u |

The bound of `5u` on `d²` is nearly reached. The `r` bound of `u` is
reached, and for points close together the subtraction is exact, as
Sterbenz's lemma predicts.

## Point-line

```julia
r = y[2] - y[1]
p = y[1] - x[1]
t = clamp01nan(-dot3(p, r)/dot3(r, r))
d = p + scale3(t, r)
dot3(d, d), d, SVector{0,T}(), SVector{1,T}(t)
```

`scale3` multiplies without fusing, so `dᵢ = fl(pᵢ + fl(t*rᵢ))`.

Exact, on the same inputs:

- `R = y[2] - y[1]`, `P = y[1] - x[1]`, `N = P⋅R`, `C = R⋅R`.
- `T* = -N/C`, the closest point of the line, and `T = clamp(T*, 0, 1)`.
- `d = ‖P + T R‖`, the exact distance.
- `M = Σ|PᵢRᵢ|`. By Cauchy-Schwarz, `|N| ≤ M ≤ ‖P‖‖R‖`.

### Result

To first order in `u`,

    -u(3‖P‖ + 4.5d) ≤ sqrt(d̂²) - d ≤ u(3‖P‖ + 4.5d) + τ,
    τ = min(11u‖P‖, (11u‖P‖)²/(2d)),
    |t̂ - T| ≤ 11u‖P‖/‖R‖.

Unlike point-point, the error is not relative to `d`: `p` and `t*r` cancel
when the point is near the line. It scales with `‖P‖`, the distance from
the point to `y[1]`, not with the distance from the origin. `τ` is the
error of `t̂`, which only makes the distance larger, and only to second
order unless `d` is below about `11u‖P‖`.

### The parameter `t`

The subtractions give `r̂ = R(1 + δ)` and `p̂ = P(1 + δ′)` componentwise.

In `n̂ = dot3(p̂, r̂)`, each term `PᵢRᵢ` picks up at most five factors
`(1 + ε)`: two from the subtractions and up to three in the `fma` chain.
Measured against `M` instead of `N`, cancellation between the terms
doesn't matter, so `|n̂ - N| ≤ γ₅ M`. `ĉ = dot3(r̂, r̂)` is `C(1 + θ)` with
`|θ| ≤ γ₅`, as for point-point.

The division rounds once, `t̃ = (-n̂/ĉ)(1 + ρ)`, so

    t̃ - T* = -((n̂ - N)/C)(1 + ρ)/(1 + θ) + T*((1 + ρ)/(1 + θ) - 1),
    |t̃ - T*| ≤ 5u M/C + 6u|T*| ≤ 11u M/C ≤ 11u‖P‖/‖R‖,

using `|T*| = |N|/C ≤ M/C`. Clamping can only bring `t̂` closer to `T`, so
`|t̂ - T| ≤ |t̃ - T*|`. If `n̂` is zero, `t̃` is exact zero and the same
bound holds. If `R` is zero, `ĉ` is zero, `t̂` is `0` from the `NaN`, and
the result is the point-point distance to `y[1]`, which is exact here.

### How the error in `t` moves the distance

Let `η = |t̃ - T*|`. Then

    d² ≤ ‖P + t̂R‖² ≤ d² + η²‖R‖².

The lower bound holds because `T` minimizes `‖P + sR‖` over `s ∈ [0, 1]`,
and `t̂` is in `[0, 1]`. For the upper bound, expand
`‖P + sR‖² = d² + 2(s - T)(P + TR)⋅R + (s - T)²‖R‖²`:

- **`T` inside:** `(P + TR)⋅R = 0`, and `|t̂ - T| ≤ η`.
- **`T = 0`:** then `T* ≤ 0`, so the cross term is `2t̂|T*|‖R‖²`. Since `t̂`
  is at most `max(t̃, 0)`, `t̂ + |T*| ≤ η`, and so
  `2t̂|T*| + t̂² ≤ (t̂ + |T*|)² ≤ η²`.
- **`T = 1`:** the same, mirrored.

So `‖P + t̂R‖ - d ≤ min(η‖R‖, η²‖R‖²/(2d)) ≤ τ`. The error in `t` costs a
first order error only when the point is within about `11u‖P‖` of the
segment.

### Forming the separation `d`

With `mᵢ = t̂r̂ᵢ(1 + αᵢ)` and `d̂ᵢ = (p̂ᵢ + mᵢ)(1 + βᵢ)`,

    d̂ = P + t̂R + E,
    |Eᵢ| ≤ u|Pᵢ| + 2u t̂|Rᵢ| + u|Pᵢ + t̂Rᵢ| + O(u²),
    ‖E‖ ≤ u(‖P‖ + 2t̂‖R‖ + d) + O(u²) ≤ u(3‖P‖ + 3d) + O(u²),

using `t̂‖R‖ ≤ ‖P‖ + ‖P + t̂R‖`.

The last `dot3` adds a relative error of `γ₃` to `d̂²`, about `1.5u d` on
the distance. Adding `‖E‖`, `τ` and this gives the result above.

The returned separation `r̂ = d̂` is within `‖E‖ + η‖R‖ ≤ u(14‖P‖ + 3d)` of
the exact `P + TR`. Its direction can be off by about `u‖P‖/d`, along the
segment, when the point is close to the segment.

### Measured

Maximum of the error over the bound, against `BigFloat` on the same inputs,
100,000 cases per generator:

| generator | `Float32` distance | `Float32` `t` | `Float64` distance | `Float64` `t` | `Float32` distance in `u·ℓ` |
|---|---|---|---|---|---|
| uniform in `[0, 100]³` | 0.36 | 0.33 | 0.43 | 0.36 | 1.97 |
| 0.01 from a segment 50 long, 1e4 from the origin | 0.46 | 0.40 | 0.36 | 0.38 | 0.008 |
| projection within 1e-6 of an end | 0.60 | 0.44 | 0.54 | 0.36 | 2.4 |
| segment 0.01 long, point 1000 away | 0.23 | 0 | 0.23 | 0 | 2.57 |
| point on the segment | 0.31 | 0.40 | 0.32 | 0.45 | 7.07 |
| vertices from a pool of three | 0.38 | 0.38 | 0.40 | 0.33 | 1.78 |

The bound holds with a margin of about two. The last column shows why it is
written in `‖P‖` rather than `ℓ`, the largest coordinate: far from the
origin, the error is much smaller than `u·ℓ`.

## Tracing the kernels

[`ErrorBounds.jl`](ErrorBounds.jl) defines a number type `EB` and runs
`dist_sqr_r_s_t` on it, so the analysis follows the code as written, with
nothing transcribed. An `EB` carries two bounds, as polynomials with
nonnegative coefficients in size symbols such as `|Pᵢ|`, `‖P‖` and `‖R‖⁻²`:
`mag` on the absolute value of the exact result, and `err` on the rounding
error, in units of `u`. Each operation applies a rule that holds for all
inputs, so the result is a bound for all inputs in terms of the sizes. Terms
of order `u²` are dropped.

The rules:

- **Sums and differences.** A sum or difference of two exact values, each a
  single size symbol, is a new size symbol, rounded once, so its error is
  relative to itself. For example `P = y[1] - x[1]` gives `|Pᵢ|`, and the
  bounds are in the sizes of the differences, not of the coordinates.
- **Vertex moves.** A rounded difference of two input vertices is treated
  as the exact difference after moving the first vertex. Following the
  differences back gives how far each vertex moved, and the exact distance
  changes by at most the largest move on each simplex. The rest of the
  kernel is then analysed on the moved vertices. This is why the kernels
  reuse the same rounded `E0` and `E1`: they compute the distance of a
  nearby triangle accurately, even a sliver.
- **Nonnegative sums**, such as `dot3(r, r)`, keep a relative error: the
  largest of their terms', plus one rounding.
- **Signed sums** are bounded by the sum of the magnitudes of their terms,
  as in the hand analysis. Correlations survive because the magnitudes stay
  symbolic. A sum over components that then scales one component, like `t`
  in `p + t*r`, is reduced to norms with Cauchy-Schwarz, and vector errors
  with Minkowski.
- **Division** needs an exact nonnegative divisor with a relative error,
  such as `dot3(R, R) = ‖R‖²`.
- **`cross3_accurate`** of exact vectors has each component within `2u` of
  the exact one (Kahan's 2×2 determinant, Jeannerod, Louvet and Muller 2013),
  and its components become size symbols, such as `|Nᵢ|`.
- **Plain `cross3(a, b)`** has component `i` within
  `(rel_a + rel_b + 2)u·|Xᵢ|`, where `|Xᵢ| = |aⱼ||bₖ| + |aₖ||bⱼ|` also bounds the
  component, and `‖X‖ ≤ √2‖a‖‖b‖`.
- **Clamps and `min`** are 1-Lipschitz, so they don't add error. The `NaN`
  branch of `clamp01nan` is assumed not to be taken: the divisor is not zero.
- **Squared distances** also carry `|√d̂² - ρ| ≤ u(abs + rel·ρ)` for their
  exact root `ρ`, from `dot3(r, r)` and from `h*h/‖N‖²`.

The kernels are branch free, so every select has both arms traced. A select
joins them, taking the larger bound of each term, and adds a jump in the
distance where the computed and exact conditions can disagree:

- **Running minimum**, `cand[1] < best[1]`: no jump, because
  `|min(â, b̂) - min(a, b)| ≤ max(|â - a|, |b̂ - b|)`. Either candidate that
  can be chosen has a distance within the error of the minimum, so the
  `rel·ρ` terms stay in terms of the final distance.
- **`inside` test**: a jump of `Σ err(v)·L(v)` over its distinct condition
  variables `v`. `L(v)` is how far the boundary is, in distance in the plane,
  per unit of `v`. It is the one hand derived fact for each kernel, given
  below.

For these selects, `r`, `s` and `t` are compared with the exact values of
the arm or candidate that was chosen, not of the exact choice: at a tie, or
across the `inside` boundary, the closest points can jump while the distance
does not.

## Point-triangle

With `B = y[1]`, `E0 = y[2] - B`, `E1 = y[3] - B`, `BP = B - x[1]`,
`E2 = E1 - E0`, `BP2 = BP + E0 = y[2] - x[1]` and `N = E0×E1`.

### The jump at the inside test

The test is `0 ≤ sbar`, `0 ≤ tbar` and `0 ≤ wbar`, the barycentric weights of
`y[2]`, `y[3]` and `B`. Let `Y` be the projection of `x[1]` onto the plane, at
height `h`. If `Y` is outside the triangle, at distance `g` in the plane, the
closest point is on an edge, and `d² = h² + g²`. So choosing the wrong arm
costs at most `g`. A weight `v` is the distance of `Y` from the line of the
opposite edge over the height of its vertex over that edge, so if the
computed and exact signs of `v` disagree, `Y` is within `err(v)·L(v)` of
that line, with

    L(sbar) = ‖N‖/‖E1‖,    L(tbar) = ‖N‖/‖E0‖,    L(wbar) = ‖N‖/‖E2‖.

These are the heights of the triangle. The edge lengths would also be
bounds, but they would give a spurious `1/sin θ` for slivers.

### The needle fix

The kernel used to test `sbar + tbar ≤ 1` instead of `0 ≤ wbar`. The traced
bound then had a term proportional to `u‖BP‖(‖E0‖ + ‖E1‖)/‖E2‖`. The errors of `sbar` and
`tbar` are about `u‖BP‖‖E0‖/‖N‖` and `u‖BP‖‖E1‖/‖N‖`, and their sum was
tested against the height over the short edge `E2` of a needle. It is
real: for a needle 100 long with a short edge of 0.001, a point 0.01 above
the plane that projects 0.2 past the short edge was judged inside in
`Float32`, with a distance of 0.01 instead of 0.2. `wbar` is now computed from
its own edge, `N⋅(BP2×E2)/‖N‖²`, and its error over `L(wbar)` is about
`u‖BP2‖`.

Inside, `sbar` and `sbar + tbar` can round above one, so the returned weights
are capped: `s = min(sbar, 1)` and `t = min(tbar, 1 - s)`.

### Result

    |√d̂² - d| ≤ u(62‖BP‖ + 40‖BP2‖ + 3‖E0‖ + 3‖E1‖ + 3‖E2‖ + 5d)

including the vertex moves, `u(‖BP‖ + ‖E0‖ + ‖E1‖)`. Most of it is the jump:
`(16 + 16)·(3/2)‖BP‖` and `18·(3/2)‖BP2‖` from the errors of the three weights.

Largest error over the bound, 20,000 cases per generator:

| generator | `Float32` | `Float64` |
|---|---|---|
| uniform in `[0, 100]³` | 0.022 | 0.023 |
| 0.01 above, projecting within 1e-6 of an edge | 0.011 | 0.011 |
| sliver | 0.010 | 0.011 |
| needle, near the short edge | 0.012 | 0.012 |
| needle, 0.2 past the short edge | 0.012 | 0.014 |
| vertices from a pool of three | 0.017 | 0.018 |

The bound holds everywhere but is loose by a factor of about 50: it adds the
worst error of every weight at once.

## Line-line

With `P = x[2] - x[1]`, `Q = y[2] - y[1]`, `W = y[1] - x[1]`,
`WQ = W + Q = y[2] - x[1]`, `WP = W - P = y[1] - x[2]` and `n = P×Q`.

### The jump at the inside test

The test is `0 ≤ sbar ≤ 1` and `0 ≤ tbar ≤ 1`. Project both line segments onto
a plane normal to `n`. They cross there exactly when the closest points of
the lines are inside both segments, and then `d = h`, the distance along
`n`. Otherwise the closest points are an end point and its closest point on
the other segment, and `d² = h² + g²`, with `g` the distance in the plane
between the projected segments. If `sbar` is off by `err(sbar)`, the crossing
moves by `err(sbar)‖P‖` along `P`, which is `err(sbar)‖P‖ sin θ = err(sbar)‖n‖/‖Q‖`
from the projection of `y`. So

    L(sbar) = ‖n‖/‖Q‖,    L(tbar) = ‖n‖/‖P‖.

With `‖P‖` and `‖Q‖` instead, nearly parallel line segments would get a
spurious `1/sin θ`.

### Result

    |√d̂² - d| ≤ u(10‖P‖ + 10‖Q‖ + 63‖W‖ + 2‖WP‖ + 2‖WQ‖ + 5d)

including the vertex moves, `u(‖P‖ + ‖Q‖ + ‖W‖)`. Most of it is the jump,
`(16 + 16)·(3/2)‖W‖`.

Largest error over the bound, 20,000 cases per generator:

| generator | `Float32` | `Float64` |
|---|---|---|
| uniform in `[0, 100]³` | 0.028 | 0.028 |
| crossing 0.01 apart, 1e4 from the origin | 0.015 | 0.022 |
| nearly parallel, `θ ≈ 1e-4` | 0.031 | 0.022 |
| crossing within 1e-6 of an end | 0.007 | 0.006 |
| vertices from a pool of three | 0.039 | 0.039 |

The bound has no `1/sin θ`, and nothing that grows with the distance from the
origin, but like point-triangle it is loose by a factor of about 30.

## Point-line, traced

For comparison with the hand analysis, the tracer gives

    |√d̂² - d| ≤ u(9‖P‖ + 3‖R‖ + 1.5d)

including the vertex moves, `u(‖P‖ + ‖R‖)`. The largest error over this
bound was 0.35, for a point on the segment.
