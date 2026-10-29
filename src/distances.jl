"""
    dist_sqr_r_s_t(x::Simplex{N,T}, y::Simplex{M,T})::Tuple{T, Vec3{T}, SVector{N-1,T}, SVector{M-1,T}}

Return a tuple `(d², r, s, t)` of:
1. `d²`: the minimum squared distance between the two simplexes.
2. `r`: the minimum separation vector, pointing from the closest point on `x` to the closest point on `y`.
3. `s`: the barycentric coordinates of the closest point on `x`.
4. `t`: the barycentric coordinates of the closest point on `y`.

`s` holds the barycentric weights of vertices `x[2:N]`; the weight of `x[1]`
is implied as `1 - sum(s)`, so the closest point on `x` is:

    (1 - sum(s))*x[1] + sum(s[i]*x[i+1] for i in 1:N-1)

Each weight is in `[0, 1]` with `sum(s) ≤ 1`. For a `Point`, `s` is empty;
for a `LineSeg`, the closest point is `(1 - s[1])*x[1] + s[1]*x[2]`.
`t` works the same way for `y`.
"""
function dist_sqr_r_s_t end

#=
What a `dist_sqr_r_s_t` method must do to be correct, for an element type `T`
(or `SIMD.Vec` lanes of it). Here `d` is the exact distance between the
simplexes as given, in `T`, and `ℓ` is the largest absolute coordinate.

For finite coordinates below about `floatmax(T)^(1/4)` (4e9 in `Float32`,
1e77 in `Float64`):

1. The result is finite: never `NaN` or `Inf`, even for degenerate simplexes
   (repeated vertices, zero length edges, collinear triangles).
2. Each weight in `s` and `t` is in `[0, 1]`, and the computed `sum(s)` and
   `sum(t)` are at most one.
3. `r` is about the separation between the points the weights give, and `d²`
   is about `dot(r, r)`. Neither needs to be exact.
4. `|√d² - d| ≤ c₁*eps(T)*ℓ + c₂*floatmin(T)^(1/4)`, with `c₁` a small
   constant (the tests allow 64 right now).
   This holds for any shape, including slivers and nearly parallel line
   segments: the error must not grow like one over an angle.

The `floatmin(T)^(1/4)` term (3e-10 in `Float32`, 1e-77 in `Float64`) is only
for when something underflows. Then accuracy may be lost, but only by that
much.

Above `floatmax(T)^(1/4)` anything goes. The methods must also be branch
free, with each `SIMD.Vec` lane bit identical to the scalar method.

`tiny_normal_cases` and `worst_case_points` in test/test-mindistance.jl check
this at the edges of the range and at the worst case points.
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

# Branch free selection of `(d2, r, s, t)` tuples, with `r` a `Vec3` and `s`
# and `t` scalars. Works on `SIMD.Vec` lanes.
@inline function _select_d2_r_s_t(c, (d2a, ra, sa, ta), (d2b, rb, sb, tb))
    (
        SIMD.vifelse(c, d2a, d2b),
        SVector(SIMD.vifelse(c, ra[1], rb[1]), SIMD.vifelse(c, ra[2], rb[2]), SIMD.vifelse(c, ra[3], rb[3])),
        SIMD.vifelse(c, sa, sb),
        SIMD.vifelse(c, ta, tb),
    )
end

# Branch free update of the running minimum `best` with `cand`.
# Ties keep the earlier minimum.
@inline _min_d2_r_s_t(best, cand) = _select_d2_r_s_t(cand[1] < best[1], cand, best)

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
    ia = inv(dot3(E0, E0))
    c1 = dot3(E0, E1)*ia
    N = cross3(E0, E1)
    # The weights of `B` and `C` at `X'`, clamped into the triangle, `v` first.
    v = clamp01nan(-dot3(N, cross3(E0, AX))/dot3(N, N))
    u = clamp01nan(fma(-v, c1, -dot3(E0, AX)*ia))
    # A select, as `Base.min` also checks for NaN, which can't happen here.
    w = one(T) - u
    v = SIMD.vifelse(v < w, v, w)
    rin = fma3(v, E1, fma3(u, E0, AX))
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
    d2, r, s, t = _min_d2_r_s_t(best, (dot3(rin, rin), rin, u, v))
    d2, r, SVector{0,T}(), SVector{2,T}(s, t)
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
Given point O, ray vector D, and a triangle with points (A,B,C), test whether the ray intersects the triangle.
If there is an intersection P. Then P can be represented as

    P = A + u(B-A) + v(C-A) = O + tD

If there is not an intersection, `t`, `u` and `v` are undefined.
=#
struct SegTriangleIntersectResult{Float}
    intersect::Bool
    t::Float
    u::Float
    v::Float
end

#=
A fast ray-triangle intersection algorithm by Moller and Trumbore.
Ref: Tomas Möller and Ben Trumbore, "Fast, minimum storage ray-triangle intersection" (1997) Journal of Graphics Tools.
=#
function moller_trumbore_intersect(o, d, a, b, c)
    fzero = zero(eltype(o))
    seg_triangle_no_intersect = SegTriangleIntersectResult(false, fzero, fzero, fzero)

    rab = b - a
    rac = c - a
    cross_d_rac = cross3(d, rac)
    det = dot3(cross_d_rac, rab)

    if det == zero(det)
        return seg_triangle_no_intersect
    end

    invdet = inv(det)
    rao = o - a
    u = dot3(cross_d_rac, rao) * invdet
    if u < 0 || u > 1
        return seg_triangle_no_intersect
    end

    cross_rao_rab = cross3(rao, rab)
    v = dot3(cross_rao_rab, d) * invdet
    if v < 0 || u + v > 1
        return seg_triangle_no_intersect
    end

    t = dot3(cross_rao_rab, rac) * invdet
    if t < 0 || t > 1
        return seg_triangle_no_intersect
    end

    SegTriangleIntersectResult(true, t, u, v)
end

@inline function dist_sqr_r_s_t(x::LineSeg{T}, y::Triangle{T}) where T
    fzero = zero(T)
    o = x[1]
    d = x[2] - x[1]
    a = y[1]
    b = y[2]
    c = y[3]
    result = moller_trumbore_intersect(o, d, a, b, c)
    if result.intersect == true
        fzero, zero(Vec3{T}), SA[result.t,], SA[result.u, result.v]
    else
        # no intersection or seg and triangle are parallel.
        # min distance is on an edge of the triangle.
        d2min, minr, mins, _tmin = dist_sqr_r_s_t(x,SA[a,b])
        mint = SA[_tmin[1], fzero]
        d2, _rmin, _smin, _tmin = dist_sqr_r_s_t(x,SA[b,c])
        if  d2 < d2min
            d2min = d2
            minr = _rmin
            mins = _smin
            mint = SA[one(T) - _tmin[1], _tmin[1]]
        end
        d2, _rmin, _smin, _tmin = dist_sqr_r_s_t(x,SA[c,a])
        if  d2 < d2min
            d2min = d2
            minr = _rmin
            mins = _smin
            mint = SA[fzero, one(T) - _tmin[1]]
        end
        d2, _rmin, _, _tri_tmin = dist_sqr_r_s_t(SA[x[1],],y)
        if  d2 < d2min
            d2min = d2
            minr = _rmin
            mins = SA[fzero]
            mint = _tri_tmin
        end
        d2, _rmin, _, _tri_tmin = dist_sqr_r_s_t(SA[x[2],],y)
        if  d2 < d2min
            d2min = d2
            minr = _rmin
            mins = SA[one(T)]
            mint = _tri_tmin
        end
        d2min, minr, mins, mint
    end
end

@inline function dist_sqr_r_s_t(x::Triangle{T}, y::Triangle{T}) where T
    fzero = zero(T)
    fone = one(T)
    d2, minr, _s, mint = @inline dist_sqr_r_s_t(SA[x[1],x[2]], y)
    _s_all = SA[fzero, _s[1], fone - _s[1]]
    mins = SA[_s_all[2], _s_all[1]]
    iszero(d2) && return d2, minr, mins, mint
    ai = 2
    bi = 3
    ci = 1
    for i in 2:3
        this_d2, _r, _s, _t = @inline dist_sqr_r_s_t(SA[x[ai],x[bi]], y)
        if this_d2 < d2
            d2 = this_d2
            minr = _r
            _s_all = SA[fzero, _s[1], fone - _s[1]]
            mins = SA[_s_all[bi], _s_all[ai]]
            mint = _t
        end
        iszero(d2) && return d2, minr, mins, mint
        ai, bi, ci = bi, ci, ai
    end
    for i in 1:3
        this_d2, _r, _s, _t = @inline dist_sqr_r_s_t(SA[y[ai],y[bi]], x)
        if this_d2 < d2
            d2 = this_d2
            minr = -_r
            _s_all = SA[fzero, _s[1], fone - _s[1]]
            mint = SA[_s_all[bi], _s_all[ai]]
            mins = _t
        end
        iszero(d2) && return d2, minr, mins, mint
        ai, bi, ci = bi, ci, ai
    end
    return d2, minr, mins, mint
end
