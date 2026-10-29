"""
    dist_sqr_r_s_t(x::Simplex{N,T}, y::Simplex{M,T})::Tuple{T, Vec3{T}, SVector{N-1,T}, SVector{M-1,T}}

Return a tuple `(d², r, s, t)` of:
1. `d²`: the minimum squared distance between the two simplexes.
2. `r`: the minimum separation vector, pointing from a closest point on `x` to a closest point on `y`.
3. `s`: the barycentric coordinates of a closest point on `x`.
4. `t`: the barycentric coordinates of a closest point on `y`.

`s` holds the barycentric weights of vertices `x[2:N]`; the weight of `x[1]`
is implied as `1 - sum(s)`, so a closest point on `x` is:

    (1 - sum(s))*x[1] + sum(s[i]*x[i+1] for i in 1:N-1)

Each weight is in `[0, 1]` with `sum(s) ≤ 1`. For a `Point`, `s` is empty;
for a `LineSeg`, a closest point is `(1 - s[1])*x[1] + s[1]*x[2]`.
`t` works the same way for `y`.
"""
function dist_sqr_r_s_t end

#=
What a `dist_sqr_r_s_t` method must do to be correct, for an element type `T`. Here `d` is the exact distance between the
simplexes as given, and `ℓ` is the largest absolute coordinate.

For finite coordinates below `(floatmax(T)^(1/4))/4` (1e9 in `Float32`,
3e76 in `Float64`):

1. The result is finite: never `NaN` or `Inf`, even for degenerate simplexes
   (repeated vertices, zero length edges, collinear triangles).
2. Each weight in `s` and `t` is in `[0, 1]`, and the computed `sum(s)` and
   `sum(t)` are at most one.
3. `d²` is computed from `r` alone, the same way in every method:
   `fma(r[3], r[3], fma(r[2], r[2], r[1]*r[1]))`.
4. `|√d² - d| ≤ 8*eps(T)*ℓ + floatmin(T)^(1/4)`.
5. With `X(s)` and `Y(t)` the points the weights give, computed exactly,
   `‖r - (Y(t) - X(s))‖` is within the same bound.

The `floatmin(T)^(1/4)` term (3e-10 in `Float32`, 1e-77 in `Float64`) is only
for when something underflows. Then accuracy may be lost, but only by that
much.

Above `(floatmax(T)^(1/4))/4` anything goes, except that for any
coordinates, if one is `NaN`, `d²` is `NaN`.

The methods must also accept `SIMD.Vec` elements and be bit identical to the scalar method.

`tiny_normal_cases` and `hard_points` in test/test-mindistance.jl partially check
this at the edges of the range and at hard points.
=#

"""
    dist_sqr(x::Simplex{N,T}, y::Simplex{M,T})::T

Return the minimum squared distance between two simplexes.

This is the `d²` of [`dist_sqr_r_s_t`](@ref), bit for bit.
"""
@inline function dist_sqr(x::Simplex{N, T}, y::Simplex{M, T}) where {T, N, M}
    d2, r, s, t = dist_sqr_r_s_t(x, y)
    d2
end

# Dispatch N > M to N ≤ M
@inline function dist_sqr_r_s_t(x::Simplex{3, T}, y::Simplex{2, T}) where T
    d2, r, s, t = dist_sqr_r_s_t(y, x)
    d2, -r, t, s
end
@inline function dist_sqr_r_s_t(x::Simplex{3, T}, y::Simplex{1, T}) where T
    d2, r, s, t = dist_sqr_r_s_t(y, x)
    d2, -r, t, s
end
@inline function dist_sqr_r_s_t(x::Simplex{2, T}, y::Simplex{1, T}) where T
    d2, r, s, t = dist_sqr_r_s_t(y, x)
    d2, -r, t, s
end

# Branch free select, elementwise for an `SVector`. Works on `SIMD.Vec` lanes.
@inline _vifelse(c, a, b) = SIMD.vifelse(c, a, b)
@inline _vifelse(c, a::SVector{N}, b::SVector{N}) where N = SVector(ntuple(i -> SIMD.vifelse(c, a[i], b[i]), Val(N)))

# Branch free selection of `(d2, r, s, t)` tuples, with `r` a `Vec3` and `s`
# and `t` scalars or `SVector`s. Works on `SIMD.Vec` lanes.
@inline function _select_d2_r_s_t(c, (d2a, ra, sa, ta), (d2b, rb, sb, tb))
    _vifelse(c, d2a, d2b), _vifelse(c, ra, rb), _vifelse(c, sa, sb), _vifelse(c, ta, tb)
end

# Branch free update of the running minimum `best` with `cand`.
# Ties keep the earlier minimum. A NaN `d2` is kept, so a NaN coordinate gives a NaN `d2`.
@inline _min_d2_r_s_t(best, cand) = _select_d2_r_s_t((cand[1] < best[1]) | isnan(cand[1]), cand, best)

# The methods

@inline function dist_sqr_r_s_t(x::Point{T}, y::Point{T}) where T
    d = y[1] - x[1]
    dot3(d, d), d, SVector{0,T}(), SVector{0,T}()
end

@inline function dist_sqr_r_s_t(x::Point{T}, y::LineSeg{T}) where T
    r = y[2] - y[1]
    # Written as `y - x[1]` so it folds to `y[1]` when `x[1]` is positive zero.
    p = y[1] - x[1]
    t = clamp01nan(-dot3(p, r)/dot3(r, r))
    d = fma3(t, r, p)
    dot3(d, d), d, SVector{0,T}(), SVector{1,T}(t)
end

#=
For the triangle `A`, `B`, `C`, the closest point is the projection `X'` of
the point onto the plane if that is inside the triangle, and otherwise on an
edge. Rather than test which, this takes the closest of three points of the
triangle: the closest point on edges `AC` and `BC`, and `X'` with its weights
clamped into the triangle, which also covers edge `AB`. Each is a point of
the triangle, so none is closer than the true distance, and `d2` is
`dot(r, r)` of its own separation `r`, which does not cancel.

The weights `(u, v)` of `X'` are a least squares solve. With
`E1 = c1*E0 + e` and `e ⟂ E0`, `X - A = (u + v*c1)*E0 + v*e` plus a normal
part, so `v = e⋅(X - A)/(e⋅e)` and `u = E0⋅(X - A)/(E0⋅E0) - v*c1`. Then
`u*E0 + v*E1 = (u + v*c1)*E0 + v*e`, so the errors of `u` and `v` move the
point by about a rounding of `‖X - A‖`, not over the sine of an angle.

`e` itself is never formed. `E1 - c1*E0` cancels for a sliver, but with
`N = E0×E1`, `N×E0 = (E0⋅E0)*e`, so `v = N⋅(E0×(X - A))/(N⋅N)`. The rounding
of `E0×(X - A)` moves `X'` by about a rounding of `‖X - A‖`. The rounding of
`N` is about a rounding of `‖E0‖*‖E1‖`. Along `N` it only scales `v`, which
moves `X'` by at most about a rounding of `‖E1‖`. In the plane it tilts the
normal, which moves `X'` in the plane by about a rounding of the height over
the sine of the sliver's angle. That adds to the distance in second order,
and when it is large the sliver is thin and an edge is close: the smaller of
the two is at most about a rounding of `‖E1‖`. The two cross products don't
wait on each other or on a division, which keeps the chain of dependent
operations short.

Edge `AB` is `v = 0`, and `u` is solved after `v` is clamped. So when `X'`
is past `AB`, `v` clamps to zero and `X'` becomes exactly the closest point
on `AB`, with an exact zero weight. `v = min(v, 1 - u)` then keeps `X'` in
the triangle. Past `AC` or `BC` it is just some point of the triangle, and
the edge candidates give the closest point.

`N⋅N` and `N⋅(E0×(X - A))` are at most fourth powers of the coordinates, so
they are in range. When `N⋅N` underflows, every point of the triangle is
within about `floatmin^(1/4)` of edges `AC` and `BC`: within `‖E0‖/2` and
within the height over `AB`, whose product is `‖N‖/2`. A degenerate triangle
or edge gives NaN or infinite weights, which `clamp01nan` makes a point of
the triangle, and then the edges are exact or within the length of the edge.
=#
@inline function dist_sqr_r_s_t(x::Point{T}, y::Triangle{T}) where T
    A = y[1]
    E0 = y[2] - A
    E1 = y[3] - A
    # Written as `y - x[1]` so it folds to `A` when `x[1]` is positive zero.
    AX = A - x[1]
    d2in, rin, u, v = _tri_face_candidate(AX, E0, E1, cross3(E0, E1))
    # The closest points on the edges `AC` and `BC`.
    u1 = clamp01nan(-dot3(E1, AX)/dot3(E1, E1))
    r1 = fma3(u1, E1, AX)
    best = (dot3(r1, r1), r1, zero(T), u1)
    E2 = E1 - E0
    BX = AX + E0
    u2 = clamp01nan(-dot3(BX, E2)/dot3(E2, E2))
    r2 = fma3(u2, E2, BX)
    best = _min_d2_r_s_t(best, (dot3(r2, r2), r2, one(T) - u2, u2))
    # Last, so a tie with an edge keeps the edge's exact zero weight.
    d2, r, s, t = _min_d2_r_s_t(best, (d2in, rin, u, v))
    d2, r, SVector{0,T}(), SVector{2,T}(s, t)
end

# The point `X'` of the triangle `A`, `B`, `C` for the point `X`, as above,
# with `AX = A - X`, `E0 = B - A`, `E1 = C - A`, and `N` about `E0×E1`.
# This is the point and triangle kernel without its edges, for the line
# segment and triangle kernels, which have closer candidates on the edges
# already and share `N` between calls. Return `(d2, r, u, v)`, with `r` from `X` to `X'`, `d2 = dot(r, r)`, and `u`
# and `v` the weights of `B` and `C`.
@inline function _tri_face_candidate(AX, E0, E1, N)
    ia = inv(dot3(E0, E0))
    c1 = dot3(E0, E1)*ia
    # The weights of `B` and `C` at `X'`, clamped into the triangle, `v` first.
    v = clamp01nan(-dot3(N, cross3(E0, AX))/dot3(N, N))
    u = clamp01nan(fma(-v, c1, -dot3(E0, AX)*ia))
    # A select, as `Base.min` also checks for NaN, which can't happen here.
    w = one(u) - u
    v = SIMD.vifelse(v < w, v, w)
    r = fma3(v, E1, fma3(u, E0, AX))
    dot3(r, r), r, u, v
end

#=
The separation between the points at `s` on `x` and `t` on `y` is
`Q0mP0 + t*Q - s*P`. Rather than test where the closest points are, this
takes the closest of three pairs of points: the closest points on `y` of the
end points of `x`, and the closest point on `x` of the point at `t` on `y`,
with `t` from the closest points of the lines through the line segments,
clamped. Each is a pair of points of the line segments, so none is closer
than the true distance, `d2` is `dot(r, r)` of its own separation `r`, which
does not cancel, and `r` is the separation of the points the weights give.

With `t` clamped to zero or one, `s` is the closest point on `x` of an end
point of `y`, so those pairs need no candidates of their own. So if the
closest points have `s` at an end of `x`, an end point candidate has them,
and otherwise the pair from `t` does.

`t = -n⋅(P×Q0mP0)/(n⋅n)` with `n = P×Q` is a triple product over `‖n‖²`
(Lagrange's identity), as in Kahan, "Computing Cross-Products and Rotations
in 2- and 3-Dimensional Euclidean Spaces" (2016), §9. Unlike the Gram
determinant `a*c - b*b` of Lumelsky (1985) and Ericson, "Real-Time Collision
Detection" (2005), §5.1.9, it doesn't cancel for nearly parallel line
segments. The rounding of `P×Q0mP0` moves `t` by about a rounding of
`‖Q0mP0‖` over `‖Q‖*sinθ`, which moves the point on `y` by about a rounding
of `‖Q0mP0‖` across the line through `x`. If that moves its closest point on
`x` off an end, the end point there is within as much of `y`. The rounding
of `n` is about a rounding of `‖P‖*‖Q‖`. Along `n` it only scales `t`, which
moves the point on `y` across the line through `x` by at most about a
rounding of `‖Q‖`. In the plane it moves that point by about a rounding of
the height over `sinθ`, which adds to the distance in second order. That is
at most `‖Q‖*sinθ` too, as the point stays on `y`, and the smaller of the
two adds at most about a rounding of `‖Q‖`.

`n⋅n` and `n⋅(P×Q0mP0)` are at most fourth powers of the coordinates, so
they are in range. When `n⋅n` underflows, `‖P‖*‖Q‖*sinθ` is below about
`floatmin^(1/2)`, so `y` is shorter than about `floatmin^(1/4)`, or `x` is
within that of a line segment parallel to `y`. For parallel line segments
any `t` will do: if an end point of `x` is closest, an end point candidate
has it, and otherwise all of `y` is beside `x`. Zero length or parallel line
segments give NaN or infinite parameters, which `clamp01nan` makes points of
the line segments.

For a different robust method see Eberly, "Robust Computation of Distance
Between Line Segments" (Geometric Tools, 2018, revised 2023).
=#
@inline function dist_sqr_r_s_t(x::LineSeg{T}, y::LineSeg{T}) where T
    P0 = x[1]
    Q0 = y[1]
    P = x[2] - P0
    Q = y[2] - Q0
    # Written as `y - x[1]` so it folds to `Q0` when `x[1]` is positive zero.
    Q0mP0 = Q0 - P0
    b = dot3(P, Q)
    c = dot3(Q, Q)
    e = dot3(Q, Q0mP0)
    n = cross3(P, Q)
    # `t` of the closest points of the lines through the line segments,
    # clamped, and the closest point on `x` of the point at `t` on `y`.
    t = clamp01nan(-dot3(n, cross3(P, Q0mP0))/dot3(n, n))
    s = clamp01nan(fma(t, b, dot3(P, Q0mP0))/dot3(P, P))
    rin = fma3(-s, P, fma3(t, Q, Q0mP0))
    # The closest points on `y` of the end points of `x`.
    t0 = clamp01nan(-e/c)
    r0 = fma3(t0, Q, Q0mP0)
    best = (dot3(r0, r0), r0, zero(T), t0)
    t1 = clamp01nan((b - e)/c)
    r1 = fma3(t1, Q, Q0mP0 - P)
    best = _min_d2_r_s_t(best, (dot3(r1, r1), r1, one(T), t1))
    # Last, so a tie with an end point keeps its exact weight.
    d2, r, s, t = _min_d2_r_s_t(best, (dot3(rin, rin), rin, s, t))
    d2, r, SVector{1,T}(s), SVector{1,T}(t)
end

#=
The closest points of a line segment `x` and a triangle are on an edge of
the triangle, or are an end point of `x` over the face, or are where `x`
crosses the face, at distance zero. If `x` is parallel to the face, sliding
along it keeps the distance until an end point or an edge. So this takes the
closest of six pairs of points: the closest points of `x` and each edge, the
end points of `x` with their `X'` as for a point and a triangle, and the
point of `x` where its line meets the plane of the triangle, with `s`
clamped, with its `X'`. Each is a pair of points of the simplexes, so none is
closer than the true distance, and `d2` is `dot(r, r)` of its own separation
`r`, which does not cancel. An `X'` past an edge is some point of the
triangle, and the edge pair is at least as close.

Testing whether `x` crosses the triangle instead, as Möller and Trumbore,
"Fast, minimum storage ray-triangle intersection" (1997), decides the
distance is zero from the signs of triple products, which are rounding
errors when `x` is in the plane of the triangle. Then a line segment beside
a triangle can get a distance of zero.

`s = N⋅(A - P0)/(N⋅P)`. The rounding of the dot products moves the point
along `x` by about a rounding of `‖A - P0‖` or `s*‖P‖` over the sine of the
angle between `x` and the plane, which moves it off the plane by just a
rounding of those. A tilt of `N` by an angle `τ` moves the point off the
plane by about `τ*‖X - A‖`, and so the distance at first order, unlike for a
point. `cross3` tilts `N` by a rounding over the sine of the triangle's
smallest angle, which for a line segment through a sliver gave errors near
`√eps*ℓ`, so `N` is `cross3_accurate`, which tilts it by about a rounding.

When `N⋅N` underflows the edges are within about `floatmin^(1/4)` of every
point of the triangle, as for a point and a triangle. A line segment
parallel to the plane gives NaN or infinite `s`, which `clamp01nan` makes an
end point.

For two triangles, the closest points are on an edge of each, or are a
vertex of one over the face of the other, or are where an edge of one
crosses the other. So that takes the closest of the pairs of edges, and of
the vertices and the edges of each with the face of the other, in the same
way.
=#
@inline function dist_sqr_r_s_t(x::LineSeg{T}, y::Triangle{T}) where T
    A = y[1]
    E0 = y[2] - A
    E1 = y[3] - A
    N = cross3_accurate(E0, E1)
    z = zero(T)
    # The closest points of `x` and the edges.
    d2, r, s, t = dist_sqr_r_s_t(x, SA[y[1], y[2]])
    best = (d2, r, s, SVector(t[1], z))
    d2, r, s, t = dist_sqr_r_s_t(x, SA[y[1], y[3]])
    best = _min_d2_r_s_t(best, (d2, r, s, SVector(z, t[1])))
    d2, r, s, t = dist_sqr_r_s_t(x, SA[y[2], y[3]])
    best = _min_d2_r_s_t(best, (d2, r, s, SVector(one(T) - t[1], t[1])))
    # The end points of `x` with the face.
    for i in 1:2
        d2, r, u, v = _tri_face_candidate(A - x[i], E0, E1, N)
        best = _min_d2_r_s_t(best, (d2, r, SVector(ifelse(i == 1, z, one(T))), SVector(u, v)))
    end
    # Last, so a tie keeps an exact weight.
    d2, r, s, u, v = _seg_tri_crossing_candidate(x[1], x[2], A, E0, E1, N)
    _min_d2_r_s_t(best, (d2, r, SVector(s), SVector(u, v)))
end

# The point at `s` of the line segment from `P0` to `P1` where its line meets
# the plane of the triangle with `y[1] = A`, `E0`, `E1`, and `N` about
# `E0×E1`, with `s` clamped, and its `X'`. Return `(d2, r, s, u, v)`.
@inline function _seg_tri_crossing_candidate(P0, P1, A, E0, E1, N)
    P = P1 - P0
    AP0 = A - P0
    s = clamp01nan(dot3(N, AP0)/dot3(N, P))
    d2, r, u, v = _tri_face_candidate(fma3(-s, P, AP0), E0, E1, N)
    d2, r, s, u, v
end

# The weights of vertices 2 and 3 of a triangle at the point
# `(1 - w)*V[i] + w*V[j]` on an edge, with `i < j`, or at `V[i]` with
# `i == j` and `w = 0`.
@inline function _tri_weights(i, j, w)
    z = zero(w)
    wi = one(w) - w
    SVector(ifelse(i == 2, wi, ifelse(j == 2, w, z)), ifelse(i == 3, wi, ifelse(j == 3, w, z)))
end

@inline function dist_sqr_r_s_t(x::Triangle{T}, y::Triangle{T}) where T
    Ax = x[1]
    E0x = x[2] - Ax
    E1x = x[3] - Ax
    Nx = cross3_accurate(E0x, E1x)
    Ay = y[1]
    E0y = y[2] - Ay
    E1y = y[3] - Ay
    Ny = cross3_accurate(E0y, E1y)
    z = zero(T)
    edges = ((1, 2), (1, 3), (2, 3))
    # A pair of vertices, to start.
    r = Ay - Ax
    best = (dot3(r, r), r, SVector(z, z), SVector(z, z))
    # The closest points of the pairs of edges.
    for (i, j) in edges, (k, l) in edges
        d2, r, s, t = dist_sqr_r_s_t(SA[x[i], x[j]], SA[y[k], y[l]])
        best = _min_d2_r_s_t(best, (d2, r, _tri_weights(i, j, s[1]), _tri_weights(k, l, t[1])))
    end
    # The vertices of each with the face of the other.
    for k in 1:3
        d2, r, u, v = _tri_face_candidate(Ay - x[k], E0y, E1y, Ny)
        best = _min_d2_r_s_t(best, (d2, r, _tri_weights(k, k, z), SVector(u, v)))
        d2, r, u, v = _tri_face_candidate(Ax - y[k], E0x, E1x, Nx)
        best = _min_d2_r_s_t(best, (d2, -r, SVector(u, v), _tri_weights(k, k, z)))
    end
    # Last, so a tie keeps an exact weight. The edges of each with the face
    # of the other.
    for (i, j) in edges
        d2, r, s, u, v = _seg_tri_crossing_candidate(x[i], x[j], Ay, E0y, E1y, Ny)
        best = _min_d2_r_s_t(best, (d2, r, _tri_weights(i, j, s), SVector(u, v)))
        d2, r, s, u, v = _seg_tri_crossing_candidate(y[i], y[j], Ax, E0x, E1x, Nx)
        best = _min_d2_r_s_t(best, (d2, -r, SVector(u, v), _tri_weights(i, j, s)))
    end
    best
end

