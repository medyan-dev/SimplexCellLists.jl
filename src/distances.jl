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

"""
    dist_sqr(x::Simplex{N,T}, y::Simplex{M,T})::T

Return the minimum squared distance between two simplexes.

This is the `d²` of [`dist_sqr_r_s_t`](@ref), bit for bit.
"""
@inline function dist_sqr(x::Simplex{N, T}, y::Simplex{M, T}) where {T, N, M}
    # Inlining even the large methods lets LLVM drop the unused `r, s, t`.
    d2, r, s, t = @inline dist_sqr_r_s_t(x, y)
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
    d = p + scale3(t, r)
    dot3(d, d), d, SVector{0,T}(), SVector{1,T}(t)
end

#=
The closest point on the triangle is either the projection of the point onto
the plane of the triangle, if that is inside the triangle, or on an edge.

This avoids expanding the squared distance as a quadratic in the barycentric
coordinates, which cancels: its absolute error scales with the squared size
of the triangle, not with the distance. Here the distance inside is along the
normal, and each edge candidate is `dot(r, r)` of its own separation `r`, so
the error of the distance is about a rounding of the coordinates.

`N` is from `cross3_accurate`, so it is exactly zero for an exactly
degenerate triangle, which then always uses the edges. With `cross3`, `N` of
a triangle with two equal edges is a rounding error, and the inside test can
pass by mistake. Unlike `(b*e - c*d)/(a*c - b*b)`, the cross product form of
the barycentric coordinates is also reliable for slivers.

Each barycentric coordinate is computed from its own edge, with an error of
about a rounding of `‖BP‖` times the length of that edge over `‖N‖`, so the
inside test is off by about a rounding of `‖BP‖` in distance from that edge.
Testing `sbar + tbar ≤ 1` instead puts the errors of `sbar` and `tbar`, which
scale with the long edges, against the short edge `E2` of a needle, and in
`Float32` a point 0.2 past the short edge of a needle 100 long could be
judged inside.

The projection onto the plane is from Kahan, "Computing Cross-Products and
Rotations in 2- and 3-Dimensional Euclidean Spaces" (2016), §9, and the
edges as in Ericson, "Real-Time Collision Detection" (2005), §5.1.5. See the
`LineSeg`-`LineSeg` method for the references on `cross3_accurate`.
=#
@inline function dist_sqr_r_s_t(x::Point{T}, y::Triangle{T}) where T
    B = y[1]
    E0 = y[2] - B
    E1 = y[3] - B
    # Written as `y - x[1]` so it folds to `B` when `x[1]` is positive zero.
    BP = B - x[1]
    N = cross3_accurate(E0, E1)
    Δ = dot3(N, N)
    invΔ = inv(Δ)
    # The barycentric coordinates of the projection onto the plane, by
    # Lagrange's identity, without the cancellation in `(b*e - c*d)/(a*c - b*b)`.
    sbar = dot3(N, cross3(E1, BP))*invΔ
    tbar = dot3(N, cross3(BP, E0))*invΔ
    # The edge from `B + E0` to `B + E1`.
    E2 = E1 - E0
    BP2 = BP + E0
    # The weight of `B`, `1 - sbar - tbar`, from its own edge. The sum
    # `sbar + tbar` cancels when `E2` is much shorter than `E0` and `E1`.
    wbar = dot3(N, cross3(BP2, E2))*invΔ
    # False for NaN, from a degenerate triangle.
    inside = (zero(T) ≤ sbar) & (zero(T) ≤ tbar) & (zero(T) ≤ wbar)
    hN = dot3(BP, N)
    # `hN*invΔ` first: `hN ≈ d*‖N‖` can be near the square root of floatmax
    # when `hN*hN` would overflow, while `hN*invΔ ≈ d/‖N‖` stays in range.
    hNinvΔ = hN*invΔ
    rin = scale3(hNinvΔ, N)
    d2in = hNinvΔ*hN
    # The closest points on the edges from `B` along `E0` and `E1`, and from
    # `B + E0` along `E2`.
    u0 = clamp01nan(-dot3(E0, BP)/dot3(E0, E0))
    r0 = BP + scale3(u0, E0)
    best = (dot3(r0, r0), r0, u0, zero(T))
    u1 = clamp01nan(-dot3(E1, BP)/dot3(E1, E1))
    r1 = BP + scale3(u1, E1)
    best = _min_d2_r_s_t(best, (dot3(r1, r1), r1, zero(T), u1))
    u2 = clamp01nan(-dot3(BP2, E2)/dot3(E2, E2))
    r2 = BP2 + scale3(u2, E2)
    best = _min_d2_r_s_t(best, (dot3(r2, r2), r2, one(T) - u2, u2))
    # Inside, `sbar` and `sbar + tbar` can round above one, so cap them to
    # keep the weights a convex combination.
    s_in = Base.FastMath.min_fast(sbar, one(T))
    t_in = Base.FastMath.min_fast(tbar, one(T) - s_in)
    d2, r, s, t = _select_d2_r_s_t(inside, (d2in, rin, s_in, t_in), best)
    d2, r, SVector{0,T}(), SVector{2,T}(s, t)
end

#=
The line segments lie in parallel planes with normal `n = P×Q`. If their
projections onto one of these planes cross, the closest points are at the
crossing and the separation is along `n`. Otherwise the closest points are
an end point of one line segment and its closest point on the other.

This is accurate for the same reason as the `Point`-`Triangle` method. `n` is
from `cross3_accurate`: with `cross3`, the error of the distance between
nearly parallel line segments grows like one over the angle between them.

The nearest points of the lines through the line segments are from Kahan,
"Computing Cross-Products and Rotations in 2- and 3-Dimensional Euclidean
Spaces" (2016), §9: the separation is along `n = P×Q`, and the parameters
are triple products over `|n|²` (Lagrange's identity), not the Gram
determinant `a*c - b*b` of Lumelsky (1985) and Ericson, "Real-Time
Collision Detection" (2005), §5.1.9, which cancels for nearly parallel line
segments. `n` is computed with Kahan's algorithm for 2×2 determinants
(Higham, "Accuracy and Stability of Numerical Algorithms", 2002, Problem
2.27; bound 2u from Jeannerod, Louvet, Muller, Math. Comp. 82, 2013).
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
    a = dot3(P, P)
    b = dot3(P, Q)
    c = dot3(Q, Q)
    d = dot3(P, Q0mP0)
    e = dot3(Q, Q0mP0)
    n = cross3_accurate(P, Q)
    Δ = dot3(n, n)
    invΔ = inv(Δ)
    # The closest points of the lines through the line segments, by
    # Lagrange's identity, without the cancellation in `a*c - b*b`.
    sbar = dot3(n, cross3(Q0mP0, Q))*invΔ
    tbar = dot3(n, cross3(Q0mP0, P))*invΔ
    # False for NaN, from parallel or zero length line segments.
    inside = (zero(T) ≤ sbar) & (sbar ≤ one(T)) & (zero(T) ≤ tbar) & (tbar ≤ one(T))
    hn = dot3(Q0mP0, n)
    # `hn*invΔ` first: `hn ≈ d*‖n‖` can be near the square root of floatmax
    # when `hn*hn` would overflow, while `hn*invΔ ≈ d/‖n‖` stays in range.
    hninvΔ = hn*invΔ
    rin = scale3(hninvΔ, n)
    d2in = hninvΔ*hn
    inva = inv(a)
    invc = inv(c)
    # The closest points on `x` of the end points of `y`, and on `y` of the
    # end points of `x`. The separation between the points at `s` on `x` and
    # `t` on `y` is `Q0mP0 + t*Q - s*P`.
    s0 = clamp01nan(d*inva)
    s1 = clamp01nan((b+d)*inva)
    t0 = clamp01nan(-e*invc)
    t1 = clamp01nan((b-e)*invc)
    rs0 = Q0mP0 - scale3(s0, P)
    best = (dot3(rs0, rs0), rs0, s0, zero(T))
    rs1 = (Q0mP0 + Q) - scale3(s1, P)
    best = _min_d2_r_s_t(best, (dot3(rs1, rs1), rs1, s1, one(T)))
    rt0 = Q0mP0 + scale3(t0, Q)
    best = _min_d2_r_s_t(best, (dot3(rt0, rt0), rt0, zero(T), t0))
    rt1 = (Q0mP0 - P) + scale3(t1, Q)
    best = _min_d2_r_s_t(best, (dot3(rt1, rt1), rt1, one(T), t1))
    d2, r, s, t = _select_d2_r_s_t(inside, (d2in, rin, sbar, tbar), best)
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

function dist_sqr_r_s_t(x::LineSeg{T}, y::Triangle{T}) where T
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

function dist_sqr_r_s_t(x::Triangle{T}, y::Triangle{T}) where T
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
