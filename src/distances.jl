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
"""
function dist_sqr(x::Simplex{N, T}, y::Simplex{M, T}) where {T, N, M}
    d2, r, s, t = @inline dist_sqr_r_s_t(y, x)
    d2
end

# Dispatch N > M to N ≤ M
function dist_sqr_r_s_t(x::Simplex{3, T}, y::Simplex{2, T}) where T
    d2, r, s, t = @inline dist_sqr_r_s_t(y, x)
    d2, -r, t, s
end
function dist_sqr_r_s_t(x::Simplex{3, T}, y::Simplex{1, T}) where T
    d2, r, s, t = @inline dist_sqr_r_s_t(y, x)
    d2, -r, t, s
end
function dist_sqr_r_s_t(x::Simplex{2, T}, y::Simplex{1, T}) where T
    d2, r, s, t = @inline dist_sqr_r_s_t(y, x)
    d2, -r, t, s
end
function dist_sqr(x::Simplex{3, T}, y::Simplex{2, T}) where T
    @inline dist_sqr(y, x)
end
function dist_sqr(x::Simplex{3, T}, y::Simplex{1, T}) where T
    @inline dist_sqr(y, x)
end
function dist_sqr(x::Simplex{2, T}, y::Simplex{1, T}) where T
    @inline dist_sqr(y, x)
end

# The methods

function dist_sqr_r_s_t(x::Point{T}, y::Point{T}) where T
    d = y[1] - x[1]
    d ⋅ d, d, SVector{0,T}(), SVector{0,T}()
end

function dist_sqr_r_s_t(x::Point{T}, y::LineSeg{T}) where T
    r = y[2] - y[1]
    p = x[1] - y[1]
    t = clamp01nan((p ⋅ r)/(r ⋅ r))
    d = t*r - p
    d ⋅ d, d, SVector{0,T}(), SVector{1,T}(t)
end

#=
Based on https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
But not fully optimized.
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06
=#
function dist_sqr_r_s_t(x::Point{T}, y::Triangle{T}) where T
    P = x[1]
    B = y[1]
    E0 = y[2] - B
    E1 = y[3] - B
    BP = B - P
    a = E0 ⋅ E0
    b = E0 ⋅ E1
    c = E1 ⋅ E1
    d = E0 ⋅ BP
    e = E1 ⋅ BP
    f = BP ⋅ BP
    Δ = a*c - b^2
    invΔ = inv(Δ)
    sbar = clamp01nan((b*e - c*d)*invΔ)
    tbar = clamp01nan((-a*e + b*d)*invΔ)
    
    # this ensures sbar and tbar are in the domain
    # if they are out of the domain, the min distance is on a boundary
    outside = sbar+tbar > one(sbar)
    sbar = ifelse(outside,zero(sbar),sbar)
    tbar = ifelse(outside,zero(tbar),tbar)
    s0 = clamp01nan(-d/a)
    t0 = clamp01nan(-e/c)
    sd = clamp01nan(-(b-c+d-e)/(a-2b+c))
    td = one(sd) - sd
    r(s,t) = a*s^2 + 2b*s*t + c*t^2 + 2d*s + 2e*t

    mins = sbar
    mint = tbar
    minr = r(sbar,tbar)
    if  r(s0,false) < minr
        minr = r(s0,false)
        mins = s0
        mint = zero(a)
    end
    if  r(false,t0) < minr
        minr = r(false,t0)
        mins = zero(a)
        mint = t0
    end
    if  r(sd,td) < minr
        minr = r(sd,td)
        mins = sd
        mint = td
    end
    d2cl = Base.FastMath.max_fast(zero(a),f+minr)
    rcl = mins*E0 + mint*E1 + BP
    d2cl, rcl, SVector{0,T}(), SVector{2,T}(mins, mint)
end

#=
Based on https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
But not fully optimized.
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06
=#
function dist_sqr(x::Point{T}, y::Triangle{T}) where T
    P = x[1]
    B = y[1]
    E0 = y[3] - B
    E1 = y[2] - B
    BP = B - P
    a = E0 ⋅ E0
    b = E0 ⋅ E1
    c = E1 ⋅ E1
    d = E0 ⋅ BP
    e = E1 ⋅ BP
    f = BP ⋅ BP
    Δ = a*c - b^2
    invΔ = inv(Δ)
    sbar = clamp01nan((b*e - c*d)*invΔ)
    tbar = clamp01nan((-a*e + b*d)*invΔ)
    
    # this ensures sbar and tbar are in the domain
    # if they are out of the domain, the min distance is on a boundary
    outside = sbar+tbar > one(sbar)
    sbar = ifelse(outside,zero(sbar),sbar)
    tbar = ifelse(outside,zero(tbar),tbar)
    s0 = clamp01nan(-d/a)
    t0 = clamp01nan(-e/c)
    sd = clamp01nan(-(b-c+d-e)/(a-2b+c))
    td = one(sd) - sd
    q(s,t) = a*s^2 + 2b*s*t + c*t^2 + 2d*s + 2e*t
    
    Base.FastMath.max_fast(zero(sd),f + Base.FastMath.min_fast(
        q(s0,false),
        q(false,t0),
        q(sd,td),
        q(sbar,tbar),
    ))
end

#=
Using the simple algorithm and some comments from
https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06
Ignores the case of degenerate line segments.
=#
function dist_sqr_r_s_t(x::LineSeg{T}, y::LineSeg{T}) where T
    P0 = x[1]
    P1 = x[2]
    Q0 = y[1]
    Q1 = y[2]
    P = P1-P0
    Q = Q1-Q0
    P0mQ0 = P0 - Q0
    a = P ⋅ P
    b = P ⋅ Q
    c = Q ⋅ Q
    d = P ⋅ P0mQ0
    e = Q ⋅ P0mQ0
    f = P0mQ0 ⋅ P0mQ0
    Δ = a*c - b^2
    #assuming both segments are not zero length
    #critical points
    s0 = clamp01nan(-d/a)
    s1 = clamp01nan((b-d)/a)
    t0 = clamp01nan(e/c)
    t1 = clamp01nan((b+e)/c)
    sbar = clamp01nan((b*e - c*d)/Δ)
    tbar = clamp01nan((a*e - b*d)/Δ)
    r(s,t) = a*s^2 - 2b*s*t + c*t^2 + 2d*s - 2e*t
    mins = sbar
    mint = tbar
    minr = r(sbar,tbar)
    if  r(s0,false) < minr
        minr = r(s0,false)
        mins = s0
        mint = zero(a)
    end
    if  r(s1,true) < minr
        minr = r(s1,true)
        mins = s1
        mint = one(a)
    end
    if  r(false,t0) < minr
        minr = r(false,t0)
        mins = zero(a)
        mint = t0
    end
    if  r(true,t1) < minr
        minr = r(true,t1)
        mins = one(a)
        mint = t1
    end
    d2cl = Base.FastMath.max_fast(zero(a),f+minr)
    rcl = (mint*Q - mins*P) - P0mQ0
    d2cl, rcl, SVector{1,T}(mins), SVector{1,T}(mint)
end

#=
Using the simple algorithm and some comments from
https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06
Ignores the case of degenerate line segments.
=#
function dist_sqr(x::LineSeg{T}, y::LineSeg{T}) where T
    P0 = x[1]
    P1 = x[2]
    Q0 = y[1]
    Q1 = y[2]
    P = P1-P0
    Q = Q1-Q0
    P0mQ0 = P0 - Q0
    a = P ⋅ P
    b = P ⋅ Q
    c = Q ⋅ Q
    d = P ⋅ P0mQ0
    e = Q ⋅ P0mQ0
    f = P0mQ0 ⋅ P0mQ0
    Δ = a*c - b^2
    invΔ = inv(Δ)
    #assuming both segments are not zero length
    #critical points
    s0 = clamp01nan(-d/a)
    s1 = clamp01nan((b-d)/a)
    t0 = clamp01nan(e/c)
    t1 = clamp01nan((b+e)/c)
    sbar = clamp01nan((b*e - c*d)*invΔ)
    tbar = clamp01nan((a*e - b*d)*invΔ)
    r(s,t) = a*s^2 - 2b*s*t + c*t^2 + 2d*s - 2e*t
    Base.FastMath.max_fast(zero(s0),f+Base.FastMath.min_fast(
        r(s0,false),
        r(s1,true),
        r(false,t0),
        r(true,t1),
        r(sbar,tbar),
    ))
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
    cross_d_rac = cross(d, rac)
    det = dot(cross_d_rac, rab)

    if det == zero(det)
        return seg_triangle_no_intersect
    end

    invdet = inv(det)
    rao = o - a
    u = dot(cross_d_rac, rao) * invdet
    if u < 0 || u > 1
        return seg_triangle_no_intersect
    end

    cross_rao_rab = cross(rao, rab)
    v = dot(cross_rao_rab, d) * invdet
    if v < 0 || u + v > 1
        return seg_triangle_no_intersect
    end

    t = dot(cross_rao_rab, rac) * invdet
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

function dist_sqr(x::LineSeg{T}, y::Triangle{T})::T where T
    fzero = zero(T)
    o = x[1]
    d = x[2] - x[1]
    a = y[1]
    b = y[2]
    c = y[3]
    result = moller_trumbore_intersect(o, d, a, b, c)
    if result.intersect == true
        fzero
    else
        # no intersection or seg and triangle are parallel.
        # min distance is on an edge of the triangle.
        dab = dist_sqr(x,SA[a,b])
        dbc = dist_sqr(x,SA[b,c])
        dca = dist_sqr(x,SA[c,a])
        dx1y = dist_sqr(SA[x[1],],y)
        dx2y = dist_sqr(SA[x[2],],y)
        Base.FastMath.min_fast(
            dab,
            dbc,
            dca,
            dx1y,
            dx2y,
        )
    end
end

function dist_sqr_r_s_t(x::Triangle{T}, y::Triangle{T}) where T
    fzero = zero(T)
    fone = one(T)
    d2, minr, _s, mint = dist_sqr_r_s_t(SA[x[1],x[2]], y)
    _s_all = SA[fzero, _s[1], fone - _s[1]]
    mins = SA[_s_all[2], _s_all[1]]
    iszero(d2) && return d2, minr, mins, mint
    ai = 2
    bi = 3
    ci = 1
    for i in 2:3
        this_d2, _r, _s, _t = dist_sqr_r_s_t(SA[x[ai],x[bi]], y)
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
        this_d2, _r, _s, _t = dist_sqr_r_s_t(SA[y[ai],y[bi]], x)
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

function dist_sqr(x::Triangle{T}, y::Triangle{T}) where T
    d2 = dist_sqr(SA[x[1],x[2]], y)
    iszero(d2) && return d2
    ai = 2
    bi = 3
    ci = 1
    for i in 2:3
        this_d2 = dist_sqr(SA[x[ai],x[bi]], y)
        d2 = Base.FastMath.min_fast(d2, this_d2)
        iszero(d2) && return d2
        ai, bi, ci = bi, ci, ai
    end
    for i in 1:3
        this_d2 = dist_sqr(SA[y[ai],y[bi]], x)
        d2 = Base.FastMath.min_fast(d2, this_d2)
        iszero(d2) && return d2
        ai, bi, ci = bi, ci, ai
    end
    return d2
end
