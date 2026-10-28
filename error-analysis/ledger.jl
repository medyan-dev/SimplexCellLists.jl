# Per-domain error-claim ledgers for the distance kernels, with focused sampling.
#
# Inputs live in the box [-1, 1]^(3 vertices), so overflow is impossible and the
# claims are absolute multiples of u = eps(T)/2. Each branch combination of a
# kernel is one domain of input space; on it the exact output is a single
# rational function of the inputs, and the computed output adds one (1+δ) per
# rounding. A ledger has one entry per domain and per boundary:
#
#   - claim:  |√d̂² - d| ≤ u·bound(exact quantities), first order in u, on the
#     samples the entry matches.
#   - status: where the claim comes from and its fine print (the O(u²)
#     remainder, and underflow of a denominator, are excluded from every
#     first-order claim and get their own rows).
#   - points: candidate worst-case inputs, where the per-domain analysis says
#     the claim is tightest, with the predicted sharpness.
#
# The sampler perturbs every coordinate of each worst-case point with Gaussian
# noise at scales σ from ~30 ulps (holds the input still, scrambles the
# roundings) up to 1e-2 (wanders across nearby domain boundaries). Each sample
# is routed to every entry whose (exact branch pattern, computed branch
# pattern) pair it matches: domain entries take pe == pc == their pattern,
# boundary entries take the samples where rounding flipped the branch. The
# measured max error/bound must stay ≤ 1 everywhere, and should approach the
# predicted sharpness at the sharp points.
#
# Run from the test environment, which has SimplexCellLists and StaticArrays:
#     include("error-analysis/ledger.jl")

using SimplexCellLists, StaticArrays, Random, LinearAlgebra
using SimplexCellLists: dist_sqr_r_s_t, dot3, cross3, cross3_accurate, clamp01nan, scale3

setprecision(BigFloat, 512)

# Entry matchers: a domain entry owns the samples where exact and computed
# arithmetic take the same branch pattern; a boundary entry owns the flips.
domain(pat) = (pe, pc) -> pe === pat && pc === pat
flip(a, b) = (pe, pc) -> (pe === a && pc === b) || (pe === b && pc === a)

# ---- point-point ----
# d = y[1] - x[1]; dot3(d, d). One domain, no branches. Per coordinate the
# rounding path is 2δ(sub) + δ(mul/fma) + the fma accumulations: 5u, 4u, 3u on
# d² for coordinates 1, 2, 3, so the claim 5u/2·d is attained only when one
# coordinate carries the whole separation.

function pp_exact(xs, ys)
    p = big.(ys[1]) - big.(xs[1])
    (; d = norm(p)), :all
end
pp_computed_pattern(xs, ys) = :all

const POINT_POINT = (exact = pp_exact, computed = pp_computed_pattern, entries = [
    (name = "whole box", claim = "5/2·d",
     status = "proved (distances.md; tracer: 3/2·d plus the ‖R‖ vertex move)",
     bound = q -> 5//2*q.d, match = (pe, pc) -> true,
     points = [
        (x = ((-1.0, 0.03, -0.02),), y = ((1.0, 0.06, 0.03),),
         expect = "supremum 1: the first coordinate carries the whole 5u path; sampled max ~0.8-0.9"),
        (x = ((-0.95, -1.0, 0.9),), y = ((0.9, 0.85, -1.0),),
         expect = "supremum 0.8: equal coordinates average the 5u, 4u, 3u paths to 4u; sampled max ~0.45 (more independent roundings must align)"),
     ]),
])

# ---- point-line ----
# r = y[2] - y[1]; p = y[1] - x[1]; t = clamp01nan(-dot3(p,r)/dot3(r,r));
# d = p + scale3(t, r); dot3(d, d).
# Exact domains in t* = -(p·r)/(r·r): t* ≤ 0 (c0 arm, d = p), 0 < t* < 1
# (interior, d² = (|p|²|r|² - (p·r)²)/|r|², the Gram ratio), t* ≥ 1 (c1 arm,
# d = p + r). clamp01nan sends NaN and t ≤ 0 to the 0 arm, t ≥ 1 to the 1 arm.

function pl_exact(xs, ys)
    p = big.(ys[1]) - big.(xs[1])
    r = big.(ys[2]) - big.(ys[1])
    rr = dot(r, r)
    pr = dot(p, r)
    pattern = iszero(rr) || -pr <= 0 ? :c0 : (-pr >= rr ? :c1 : :interior)
    tstar = pattern === :c0 ? big(0.0) : (pattern === :c1 ? big(1.0) : -pr/rr)
    (; d = norm(p + tstar*r), np = norm(p), nr = sqrt(rr)), pattern
end

# The branch pattern the floating point kernel takes, replicating its ops.
function pl_computed_pattern(xs, ys)
    r = ys[2] - ys[1]
    p = ys[1] - xs[1]
    traw = -dot3(p, r)/dot3(r, r)
    traw > 0 ? (traw < 1 ? :interior : :c1) : :c0
end

const POINT_LINE = (exact = pl_exact, computed = pl_computed_pattern, entries = [
    (name = "t clamped to 0 (t* < 0)", claim = "5/2·d",
     status = "proved: scale3(0,r) = 0 exactly, so d = p̂ and the point-point argument applies",
     bound = q -> 5//2*q.d, match = domain(:c0),
     points = [
        (x = ((1.0, 0.12, -0.03),), y = ((-1.0, 0.0, 0.0), (-1.0, -1.0, 0.0)),
         expect = "supremum 1: separation on the first coordinate, as in point-point; sampled max ~0.9"),
     ]),
    (name = "interior (0 < t* < 1)", claim = "9‖p‖ + 3‖r‖ + 3/2·d",
     status = "proved by the tracer (analyze.jl, all domains joined); slack expected: it charges the t rounding linearly, but d ⊥ r at the exact minimum makes it second order",
     bound = q -> 9*q.np + 3*q.nr + 3//2*q.d, match = domain(:interior),
     points = [
        (x = ((1.0, 0.0, 0.05),), y = ((-1.0, 1.0, 0.0), (-1.0, -1.0, 0.0)),
         expect = "≈0.2-0.3: generic midpoint, the 9‖p‖ term is mostly unspent"),
        (x = ((1.0, 0.0, 0.0),), y = ((-1.0, 0.1, 0.02), (-1.0, -0.1, 0.02)),
         expect = "low despite t error ~110u (‖p‖/‖r‖ = 10): gradient-zero kills it"),
     ]),
    (name = "t clamped to 1 (t* > 1)", claim = "‖p‖ + ‖r‖ + 5/2·d",
     status = "derived here: d̂ᵢ = fl(p̂ᵢ + r̂ᵢ) has error ≤ u(|pᵢ| + |rᵢ| + |pᵢ+rᵢ|), then Minkowski and the 3u/2 dot-self root",
     bound = q -> q.np + q.nr + 5//2*q.d, match = domain(:c1),
     points = [
        (x = ((1.0, 1.0, 0.04),), y = ((-1.0, -1.0, 0.0), (0.98, 0.98, 0.0)),
         expect = "supremum ≈1: p + r cancels to 0.03, so the input roundings u(‖p‖+‖r‖) dwarf d; sampled max ~0.6-0.8"),
     ]),
    (name = "boundary t* ≈ 0 (flips)", claim = "max of both arms",
     status = "proved modulo O(u²): the clamp boundary is C¹, a flip at |t*| ≤ t-err costs t*²‖r‖² ≤ (11u‖p‖)² on d²",
     bound = q -> max(5//2*q.d, 9*q.np + 3*q.nr + 3//2*q.d), match = flip(:c0, :interior),
     nmul = 10,   # flips are ~0.1% of samples, oversample to get hundreds
     points = [
        # p·r = 0 through cancelling ±1 products, so the computed condition
        # carries full u-scale error and actually flips. (With componentwise
        # orthogonality, e.g. p ∥ e1 and r ∥ e2, the dot3 error is only
        # u·Σ|pᵢrᵢ| ~ u·σ and flips essentially never happen.)
        (x = ((1.0, 1.0, 0.0),), y = ((-1.0, 0.0, 0.0), (-0.5, -1.0, 0.0)),
         expect = "flips only at small σ; ratio like the adjacent domains"),
     ]),
    (name = "boundary t* ≈ 1 (flips)", claim = "max of both arms",
     status = "proved modulo O(u²): same clamp-quadratic argument at t* = 1",
     bound = q -> max(q.np + q.nr + 5//2*q.d, 9*q.np + 3*q.nr + 3//2*q.d), match = flip(:c1, :interior),
     nmul = 10,
     points = [
        # -p·r = ‖r‖² = 2 with order-1 terms in both dots: a ~4u condition band.
        (x = ((0.5, 0.5, 0.0),), y = ((-1.0, 1.0, 0.0), (0.0, 0.0, 0.0)),
         expect = "flips only at small σ; ratio like the adjacent domains"),
     ]),
    (name = "degenerate (‖r‖² ≤ realmin)", claim = "not sampled",
     status = "argued: r = 0 gives t = NaN → 0 arm, d = p̂, exact-correct; ‖r‖² subnormal can misclamp t, costing ≤ ‖r‖ ≤ √realmin absolute, below any u-scale budget on the box. Gaussian sampling cannot reach sin-type slivers this deep",
     bound = q -> big(1.0), match = (pe, pc) -> false, points = NamedTuple[]),
])

# ---- point-triangle ----
# Exact domains are the Voronoi regions of the triangle: :face (projection
# inside), :e0/:e1/:e2 (edge interiors), :va/:vb/:vc (vertices). The inside
# arm is d = |BP·N|/‖N‖ with N from cross3_accurate (2u per component), so
# its triple-product error is ≤ 5u·Σ|NᵢBPᵢ| ≤ 5u‖N‖‖BP‖ (Cauchy-Schwarz),
# Δ = dot3(N,N) is relative 7u, and the vertex moves add u(‖BP‖+‖E0‖+‖E1‖).
# Each edge candidate is literally the point-line kernel, so its traced bound
# applies, with extra u‖BP2‖, u‖E2‖ charges for edge 2's rounded sums.

function pt_exact(xs, ys)
    A, B, C = big.(ys[1]), big.(ys[2]), big.(ys[3])
    E0 = B - A; E1 = C - A; E2 = C - B
    BP = A - big.(xs[1]); BP2 = B - big.(xs[1])
    N = cross(E0, E1)
    Δ = dot(N, N)
    q = (; d = big(0.0), nBP = norm(BP), nBP2 = norm(BP2),
        nE0 = norm(E0), nE1 = norm(E1), nE2 = norm(E2))
    # The inside test on the exact signs of the barycentric numerators.
    if !iszero(Δ) && dot(N, cross(E1, BP)) >= 0 && dot(N, cross(BP, E0)) >= 0 &&
            dot(N, cross(BP2, E2)) >= 0
        return merge(q, (; d = abs(dot(BP, N))/sqrt(Δ))), :face
    end
    f = big(Inf); pattern = :none
    for (p, r, pat0, pat1, patmid) in
            ((BP, E0, :va, :vb, :e0), (BP, E1, :va, :vc, :e1), (BP2, E2, :vb, :vc, :e2))
        rr = dot(r, r)
        u = iszero(rr) ? big(0.0) : clamp(-dot(p, r)/rr, big(0.0), big(1.0))
        v = p + u*r
        fv = dot(v, v)
        # Ties keep the earlier candidate, as in _min_d2_r_s_t.
        if fv < f
            f = fv
            pattern = iszero(u) ? pat0 : (isone(u) ? pat1 : patmid)
        end
    end
    merge(q, (; d = sqrt(f))), pattern
end

# The branch pattern the kernel takes, mirroring its ops (same helpers, same
# expression shapes, so the values are bit identical to the kernel's).
function pt_computed_pattern(xs, ys)
    B = ys[1]
    E0 = ys[2] - B
    E1 = ys[3] - B
    BP = B - xs[1]
    N = cross3_accurate(E0, E1)
    Δ = dot3(N, N)
    invΔ = inv(Δ)
    sbar = dot3(N, cross3(E1, BP))*invΔ
    tbar = dot3(N, cross3(BP, E0))*invΔ
    E2 = E1 - E0
    BP2 = BP + E0
    wbar = dot3(N, cross3(BP2, E2))*invΔ
    T = eltype(eltype(xs))
    (zero(T) <= sbar) & (zero(T) <= tbar) & (zero(T) <= wbar) && return :face
    f = T(Inf); pattern = :none
    for (p, r, pat0, pat1, patmid) in
            ((BP, E0, :va, :vb, :e0), (BP, E1, :va, :vc, :e1), (BP2, E2, :vb, :vc, :e2))
        u = clamp01nan(-dot3(r, p)/dot3(r, r))
        v = p + scale3(u, r)
        fv = dot3(v, v)
        if fv < f
            f = fv
            pattern = u == 0 ? pat0 : (u == 1 ? pat1 : patmid)
        end
    end
    pattern
end

pt_face_b(q) = 6*q.nBP + q.nE0 + q.nE1 + 8*q.d
pt_traced_b(q) = 62*q.nBP + 40*q.nBP2 + 3*q.nE0 + 3*q.nE1 + 3*q.nE2 + 5*q.d

const POINT_TRIANGLE = (exact = pt_exact, computed = pt_computed_pattern, entries = [
    (name = "face (inside)", claim = "6‖BP‖ + ‖E0‖ + ‖E1‖ + 8d",
     status = "derived here: hN err ≤ 5u·Σ|NᵢBPᵢ| ≤ 5u‖N‖‖BP‖ (Kahan 2u/comp + CS), Δ rel 7u, moves; d coefficient padded",
     bound = pt_face_b, match = domain(:face),
     points = [
        # Tilted triangle so N is generic and the hN terms can cancel.
        (x = ((-0.01, 0.01, 0.29),), y = ((-0.9, -0.6, 0.2), (0.8, -0.2, -0.3), (-0.2, 0.9, 0.1)),
         expect = "well below 1: CS equality needs BP aligned with N term by term"),
        (x = ((-0.09, 0.03, 0.02),), y = ((-0.9, -0.6, 0.2), (0.8, -0.2, -0.3), (-0.2, 0.9, 0.1)),
         expect = "the live point: d ≈ 0.02 and the hN terms (~0.5 each) cancel, so the ‖BP‖ term carries; CS is ~2x loose here and the moves/edge terms pad the denominator, sampled max ~0.1"),
     ]),
    (name = "edge 0 interior", claim = "9‖BP‖ + 3‖E0‖ + 3/2·d",
     status = "the candidate is the point-line kernel verbatim, so its traced interior bound plus moves applies",
     bound = q -> 9*q.nBP + 3*q.nE0 + 3//2*q.d, match = domain(:e0),
     points = [
        (x = ((0.2, -0.9, 0.4),), y = ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0)),
         expect = "supremum ~0.17 (bound padded for this geometry), sampled max ~0.07"),
     ]),
    (name = "edge 1 interior", claim = "9‖BP‖ + 3‖E1‖ + 3/2·d",
     status = "by symmetry with edge 0; not sampled directly, catches strays",
     bound = q -> 9*q.nBP + 3*q.nE1 + 3//2*q.d, match = domain(:e1), points = NamedTuple[]),
    (name = "edge 2 interior", claim = "9‖BP2‖ + 3‖E2‖ + ‖BP‖ + ‖E0‖ + ‖E1‖ + 3/2·d",
     status = "point-line bound on (BP2, E2) plus u‖BP2‖, u‖E2‖ for their rounded sums and the vertex moves",
     bound = q -> 9*q.nBP2 + 3*q.nE2 + q.nBP + q.nE0 + q.nE1 + 3//2*q.d, match = domain(:e2),
     points = [
        (x = ((0.56, 0.59, 0.3),), y = ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0)),
         expect = "like edge 0 with the extra BP2/E2 charges, sampled max ~0.14"),
     ]),
    (name = "vertex a", claim = "5/2·d",
     status = "proved: the winning candidate is BP + scale3(0, E) = BP exactly, point-point argument",
     bound = q -> 5//2*q.d, match = domain(:va),
     points = [
        (x = ((-1.0, -0.95, 0.35),), y = ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0)),
         expect = "like point-point with balanced coordinates, sampled max ~0.4"),
     ]),
    (name = "vertex b", claim = "‖BP‖ + ‖E0‖ + 5/2·d",
     status = "derived: both reachable candidates compute fl(BP + E0), error ≤ u(|BPᵢ|+|E0ᵢ|+|BP2ᵢ|) per component",
     bound = q -> q.nBP + q.nE0 + 5//2*q.d, match = domain(:vb),
     points = [
        (x = ((1.0, -0.45, 0.3),), y = ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0)),
         expect = "sampled max well below 1 unless BP and E0 errors align with the small d"),
     ]),
    (name = "vertex c", claim = "‖BP‖ + ‖E1‖ + ‖BP2‖ + ‖E2‖ + ‖E0‖ + 5/2·d",
     status = "covers whichever tied candidate wins (fl(BP+E1), or fl(fl(BP+E0)+fl(E1-E0))); not sampled directly",
     bound = q -> q.nBP + q.nE1 + q.nBP2 + q.nE2 + q.nE0 + 5//2*q.d, match = domain(:vc),
     points = NamedTuple[]),
    (name = "inside-test boundary (flips)", claim = "tracer all-domain bound",
     status = "proved by the tracer including the inside-select jump terms (heights ‖N‖/‖E‖); expect small ratios, the joined jump dominates the bound",
     bound = pt_traced_b, match = (pe, pc) -> (pe === :face) ⊻ (pc === :face),
     nmul = 10,
     points = [
        # x over the midpoint of edge 0, so the tbar numerator is exactly zero.
        (x = ((-0.05, -0.4, 0.6),), y = ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0)),
         expect = "flips at small σ; ratio small against the joined bound"),
     ]),
    (name = "clamp and tie flips", claim = "tracer all-domain bound",
     status = "clamp boundaries are C¹ (O(u²) flip cost) and min ties are 1-Lipschitz; the traced join covers every arm",
     bound = pt_traced_b, match = (pe, pc) -> pe !== pc && pe !== :face && pc !== :face,
     nmul = 10,
     points = [
        # x exactly on the u0 = 1 plane through b: (x-a)·E0 = ‖E0‖² by construction.
        (x = ((0.88, -0.54, 0.3),), y = ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0)),
         expect = "edge0↔vertex b flips at small σ; O(u²) cost, tiny ratios"),
     ]),
    (name = "degenerate (‖N‖² ≤ realmin)", claim = "not sampled",
     status = "argued: cross3_accurate is exactly zero for an exactly degenerate triangle, so inside is false and the edge candidates (their own ledger rows) answer; the subnormal-Δ sliver needs the inside test to fail safe, which zero/‖N‖² ≥ 0 and NaN→false give",
     bound = q -> big(1.0), match = (pe, pc) -> false, points = NamedTuple[]),
])

# ---- line-line ----
# Exact domains by the clamp state of (s*, t*), one letter each for s then t
# (z = clamped to 0, o = clamped to 1, i = interior): :ii (both interior, the
# crossing arm), :iz/:io/:zi/:oi (one endpoint against the other's interior,
# point-line candidates), :zz/:oz/:zo/:oo (corners). The crossing arm is
# d = |W·n|/‖n‖ with n = cross3_accurate(P, Q): same error shape as the
# triangle face, and no 1/sinθ anywhere — the nearly parallel point stresses
# exactly the cancellation this kernel was rewritten to avoid.

function ll_exact(xs, ys)
    P0, P1 = big.(xs[1]), big.(xs[2])
    Q0, Q1 = big.(ys[1]), big.(ys[2])
    P = P1 - P0; Q = Q1 - Q0; W = Q0 - P0
    q = (; d = big(0.0), nP = norm(P), nQ = norm(Q), nW = norm(W),
        nWP = norm(W - P), nWQ = norm(W + Q))
    n = cross(P, Q)
    Δ = dot(n, n)
    if !iszero(Δ)
        sbar = dot(n, cross(W, Q))/Δ
        tbar = dot(n, cross(W, P))/Δ
        if 0 <= sbar <= 1 && 0 <= tbar <= 1
            return merge(q, (; d = abs(dot(W, n))/sqrt(Δ))), :ii
        end
    end
    a = dot(P, P); c = dot(Q, Q)
    f = big(Inf); pattern = :none
    for (den, num, vof, pat0, pat1, patmid) in (
            (a, dot(P, W), u -> W - u*P, :zz, :oz, :iz),
            (a, dot(P, W + Q), u -> (W + Q) - u*P, :zo, :oo, :io),
            (c, -dot(Q, W), u -> W + u*Q, :zz, :zo, :zi),
            (c, dot(Q, P - W), u -> (W - P) + u*Q, :oz, :oo, :oi))
        u = iszero(den) ? big(0.0) : clamp(num/den, big(0.0), big(1.0))
        v = vof(u)
        fv = dot(v, v)
        if fv < f
            f = fv
            pattern = iszero(u) ? pat0 : (isone(u) ? pat1 : patmid)
        end
    end
    merge(q, (; d = sqrt(f))), pattern
end

function ll_computed_pattern(xs, ys)
    P0 = xs[1]
    Q0 = ys[1]
    P = xs[2] - P0
    Q = ys[2] - Q0
    Q0mP0 = Q0 - P0
    a = dot3(P, P)
    b = dot3(P, Q)
    c = dot3(Q, Q)
    d = dot3(P, Q0mP0)
    e = dot3(Q, Q0mP0)
    n = cross3_accurate(P, Q)
    Δ = dot3(n, n)
    invΔ = inv(Δ)
    sbar = dot3(n, cross3(Q0mP0, Q))*invΔ
    tbar = dot3(n, cross3(Q0mP0, P))*invΔ
    T = eltype(eltype(xs))
    ((zero(T) <= sbar) & (sbar <= one(T)) & (zero(T) <= tbar) & (tbar <= one(T))) && return :ii
    inva = inv(a)
    invc = inv(c)
    f = T(Inf); pattern = :none
    for (u, v, pat0, pat1, patmid) in (
            (s0 = clamp01nan(d*inva); (s0, Q0mP0 - scale3(s0, P), :zz, :oz, :iz)),
            (s1 = clamp01nan((b+d)*inva); (s1, (Q0mP0 + Q) - scale3(s1, P), :zo, :oo, :io)),
            (t0 = clamp01nan(-e*invc); (t0, Q0mP0 + scale3(t0, Q), :zz, :zo, :zi)),
            (t1 = clamp01nan((b-e)*invc); (t1, (Q0mP0 - P) + scale3(t1, Q), :oz, :oo, :oi)))
        fv = dot3(v, v)
        if fv < f
            f = fv
            pattern = u == 0 ? pat0 : (u == 1 ? pat1 : patmid)
        end
    end
    pattern
end

ll_ii_b(q) = 6*q.nW + q.nP + q.nQ + 8*q.d
ll_traced_b(q) = 63*q.nW + 10*q.nP + 10*q.nQ + 2*q.nWP + 2*q.nWQ + 5*q.d

const LINE_LINE = (exact = ll_exact, computed = ll_computed_pattern, entries = [
    (name = "crossing (interior-interior)", claim = "6‖W‖ + ‖P‖ + ‖Q‖ + 8d",
     status = "derived here: hn err ≤ 5u·Σ|nᵢWᵢ| ≤ 5u‖n‖‖W‖ (Kahan 2u/comp + CS), Δ rel 7u, moves; d coefficient padded",
     bound = ll_ii_b, match = domain(:ii),
     points = [
        (x = ((-1.0, -0.2, 0.0), (1.0, 0.2, 0.0)), y = ((-0.3, -1.0, 0.25), (0.3, 1.0, 0.25)),
         expect = "generic crossing, well below 1"),
        (x = ((-1.0, -0.004, 0.0), (1.0, 0.004, 0.0)), y = ((-1.0, 0.004, 0.15), (1.0, -0.004, 0.15)),
         expect = "the stress point: sinθ ≈ 0.008, where the Gram form loses ~125×; Kahan cross keeps the claim"),
     ]),
    (name = "s interior, t = 0", claim = "9‖W‖ + 3‖P‖ + ‖Q‖ + 3/2·d",
     status = "point-line traced bound on (W, P) plus moves",
     bound = q -> 9*q.nW + 3*q.nP + q.nQ + 3//2*q.d, match = domain(:iz),
     points = [
        (x = ((-1.0, -0.2, 0.0), (1.0, 0.2, 0.0)), y = ((0.1, 0.3, 0.3), (0.9, 0.9, 0.5)),
         expect = "sampled max ~0.06: the 9‖W‖ bound is padded for this geometry"),
     ]),
    (name = "s = 0, t interior", claim = "9‖W‖ + 3‖Q‖ + ‖P‖ + 3/2·d",
     status = "point-line traced bound on (W, Q) plus moves",
     bound = q -> 9*q.nW + 3*q.nQ + q.nP + 3//2*q.d, match = domain(:zi),
     points = [
        (x = ((-0.4, -0.3, 0.0), (0.56, 0.42, 0.0)), y = ((-0.31, -0.67, 0.25), (-0.73, -0.11, 0.25)),
         expect = "sampled max ~0.08"),
     ]),
    (name = "s interior, t = 1", claim = "10‖WQ‖ + 3‖P‖ + ‖W‖ + ‖Q‖ + 3/2·d",
     status = "point-line bound on (W+Q, P) plus the rounded sum and numerator charges; by symmetry, not sampled directly",
     bound = q -> 10*q.nWQ + 3*q.nP + q.nW + q.nQ + 3//2*q.d, match = domain(:io), points = NamedTuple[]),
    (name = "s = 1, t interior", claim = "10‖WP‖ + 3‖Q‖ + ‖W‖ + ‖P‖ + 3/2·d",
     status = "point-line bound on (W-P, Q) plus the rounded sum and numerator charges; by symmetry, not sampled directly",
     bound = q -> 10*q.nWP + 3*q.nQ + q.nW + q.nP + 3//2*q.d, match = domain(:oi), points = NamedTuple[]),
    (name = "corner (0,0)", claim = "5/2·d",
     status = "proved: the winning separation is W exactly, point-point argument",
     bound = q -> 5//2*q.d, match = domain(:zz),
     points = [
        (x = ((-0.4, -0.3, 0.0), (0.56, 0.42, 0.0)), y = ((-0.7, -0.8, 0.3), (-1.0, -0.9, 0.5)),
         expect = "like point-point, sampled max ~0.4-0.9"),
     ]),
    (name = "corner (1,0)", claim = "‖W‖ + ‖P‖ + 5/2·d",
     status = "derived: fl(W-P), error ≤ u(|Wᵢ|+|Pᵢ|+|WPᵢ|) per component; by symmetry, not sampled directly",
     bound = q -> q.nW + q.nP + 5//2*q.d, match = domain(:oz), points = NamedTuple[]),
    (name = "corner (0,1)", claim = "‖W‖ + ‖Q‖ + 5/2·d",
     status = "derived: fl(W+Q); by symmetry, not sampled directly",
     bound = q -> q.nW + q.nQ + 5//2*q.d, match = domain(:zo), points = NamedTuple[]),
    (name = "corner (1,1)", claim = "‖W‖ + ‖P‖ + ‖Q‖ + ‖WP‖ + ‖WQ‖ + 5/2·d",
     status = "covers whichever tied candidate wins (fl(fl(W+Q)-P) or fl(fl(W-P)+Q)); not sampled directly",
     bound = q -> q.nW + q.nP + q.nQ + q.nWP + q.nWQ + 5//2*q.d, match = domain(:oo), points = NamedTuple[]),
    (name = "inside-test boundary (flips)", claim = "tracer all-domain bound",
     status = "proved by the tracer including the inside-select jump terms (heights ‖n‖/‖P‖, ‖n‖/‖Q‖)",
     bound = ll_traced_b, match = (pe, pc) -> (pe === :ii) ⊻ (pc === :ii),
     nmul = 10,
     points = [
        # Q crosses exactly over P0: the sbar numerator is exactly zero.
        (x = ((-0.8, -0.6, 0.0), (0.8, 0.6, 0.0)), y = ((-0.68, -0.76, 0.25), (-0.92, -0.44, 0.25)),
         expect = "flips at small σ; ratio small against the joined bound"),
     ]),
    (name = "clamp and tie flips", claim = "tracer all-domain bound",
     status = "clamp boundaries are C¹ (O(u²) flip cost) and min ties are 1-Lipschitz; the traced join covers every arm",
     bound = ll_traced_b, match = (pe, pc) -> pe !== pc && pe !== :ii && pc !== :ii,
     nmul = 10,
     points = [
        # (Q0-P0) ⊥ Q exactly: the t* = 0 clamp boundary in the s* = 0 region.
        (x = ((-0.4, -0.3, 0.0), (0.56, 0.42, 0.0)), y = ((-0.52, -0.39, 0.25), (-0.82, 0.01, 0.25)),
         expect = "corner↔edge flips at small σ; O(u²) cost, tiny ratios"),
     ]),
    (name = "degenerate (‖n‖² ≤ realmin)", claim = "not sampled",
     status = "argued: cross3_accurate is exactly zero for exactly parallel segments, inside is false (NaN→false), and the edge candidates answer; a wrong clamp from subnormal Δ costs ≤ the parallel gap, covered by the edge rows",
     bound = q -> big(1.0), match = (pe, pc) -> false, points = NamedTuple[]),
])

function perturb(T, v, σ)
    SVector{3, T}(T(v[1] + σ*randn()), T(v[2] + σ*randn()), T(v[3] + σ*randn()))
end

function run_ledger(name, kernel, T; n=1000)
    u = big(eps(T))/2
    σs = 10.0 .^ range(log10(4*eps(T)), -2; length=7)
    entries = kernel.entries
    counts = zeros(Int, length(entries))
    maxratio = zeros(length(entries))
    atσ = fill(NaN, length(entries))
    # Sharpness at each entry's own points: max ratio among samples that both
    # originate from the point and match the entry.
    ownmax = [zeros(length(e.points)) for e in entries]
    unmatched = Dict{Tuple{Symbol, Symbol}, Int}()
    for (oj, oe) in enumerate(entries), (opi, pt) in enumerate(oe.points), σ in σs, _ in 1:(n*get(oe, :nmul, 1))
        xs = map(v -> perturb(T, v, σ), pt.x)
        ys = map(v -> perturb(T, v, σ), pt.y)
        q, pe = kernel.exact(xs, ys)
        pc = kernel.computed(xs, ys)
        err = abs(sqrt(big(dist_sqr_r_s_t(SVector(xs), SVector(ys))[1])) - q.d)
        routed = false
        for (j, e) in enumerate(entries)
            e.match(pe, pc) || continue
            routed = true
            b = u*e.bound(q)
            iszero(b) && continue
            ratio = Float64(err/b)
            counts[j] += 1
            ratio > maxratio[j] && (maxratio[j] = ratio; atσ[j] = σ)
            j == oj && ratio > ownmax[j][opi] && (ownmax[j][opi] = ratio)
        end
        routed || (unmatched[(pe, pc)] = get(unmatched, (pe, pc), 0) + 1)
    end
    println("=== $name ($T) ===")
    for (j, e) in enumerate(entries)
        if counts[j] > 0
            println("  ", rpad(e.name, 32), lpad(counts[j], 7), " samples  max err/bound ",
                rpad(round(maxratio[j], sigdigits=3), 9), " @σ=", round(atσ[j], sigdigits=2),
                maxratio[j] > 1 ? "  *** CLAIM VIOLATED ***" : "")
        elseif isempty(e.points)
            println("  ", rpad(e.name, 32), " not sampled")
        else
            println("  ", rpad(e.name, 32), " NO SAMPLES MATCHED — coverage gap")
        end
        println("      claim u·(", e.claim, ") — ", e.status)
        for (opi, pt) in enumerate(e.points)
            println("      point ", opi, " max ", rpad(round(ownmax[j][opi], sigdigits=3), 9),
                " expect ", pt.expect)
        end
    end
    isempty(unmatched) || println("  UNROUTED samples: ", sort!(collect(unmatched); by=last, rev=true))
end

Random.seed!(7)
for T in (Float32, Float64)
    run_ledger("point-point", POINT_POINT, T)
    run_ledger("point-line", POINT_LINE, T)
    run_ledger("point-triangle", POINT_TRIANGLE, T)
    run_ledger("line-line", LINE_LINE, T)
end
