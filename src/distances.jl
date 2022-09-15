using LinearAlgebra
using StaticArrays

clamp01nan(x) = ifelse(isnan(x), zero(x), clamp(x,zero(x),one(x)))

dist2PointPoint(a,b) = (a[1] - b[1]) ⋅ (a[1] - b[1])


"Return the minimum distance squared and where on the line it occurs, from 0 to 1"
function dist2PointLine_t(x,y)
    r= y[2]-y[1]
    p= x[1]-y[1]
    # c= a[1]-b[2]
    t = clamp01nan(p⋅r/(r⋅r))
    d = t*r - p
    d⋅d, t
end

function dist2PointLine(x,y)
    @inline dist2PointLine_t(x,y)[1]
end

"""
    dist2LineLine(x,y)

Return the squared distance between two line segments.
Using the simple algorithm and some comments from
https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06
Ignores the case of degenerate line segments.
"""
function dist2LineLine(x,y)
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
    return Base.FastMath.max_fast(zero(s0),f+Base.FastMath.min_fast(
        r(s0,0),
        r(s1,1),
        r(0,t0),
        r(1,t1),
        r(sbar,tbar),
    ))
end

"""
    dist2LineLine_s_t(x,y)

Return the squared distance between two line segments.
And where on the line segments it occurs, from 0 to 1.
Using the simple algorithm and some comments from
https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06
Ignores the case of degenerate line segments.
"""
function dist2LineLine_s_t(x,y)
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
    if  r(s0,0) < minr
        minr = r(s0,0)
        mins = s0
        mint = zero(a)
    end
    if  r(s1,1) < minr
        minr = r(s1,1)
        mins = s1
        mint = one(a)
    end
    if  r(0,t0) < minr
        minr = r(0,t0)
        mins = zero(a)
        mint = t0
    end
    if  r(1,t1) < minr
        minr = r(1,t1)
        mins = one(a)
        mint = t1
    end
    d2cl = max(0,f+minr)
    return d2cl, mins, mint
end


function dist2PointTriangle(x,y)
    P = x[1]
    B = y[1]
    E0= y[3] - B
    E1= y[2] - B
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
    
    # # this ensures sbar and tbar are in the domain
    # # if they are out of the domain, the min distance is on a boundary
    outside = sbar+tbar > one(sbar)
    sbar = ifelse(outside,zero(sbar),sbar)
    tbar = ifelse(outside,zero(tbar),tbar)
    s0 = clamp01nan(-d/a)
    t0 = clamp01nan(-e/c)
    sd = clamp01nan(-(b-c+d-e)/(a-2b+c))
    td = one(sd) - sd
    q(s,t) = a*s^2 + 2b*s*t + c*t^2 + 2d*s + 2e*t
    
    return Base.FastMath.max_fast(zero(sd),f + Base.FastMath.min_fast(
        q(s0,0),
        q(0,t0),
        q(sd,td),
        q(sbar,tbar),
    ))
end


"""
Given point O, ray vector D, and a triangle with points (A,B,C), test whether the ray intersects the triangle.
If there is an intersection P. Then P can be represented as

    P = A + u(B-A) + v(C-A) = O + tD

If there is not an intersection, `t`, `u` and `v` are undefined.
"""
struct SegTriangleIntersectResult{Float}
    intersect::Bool
    t::Float
    u::Float
    v::Float
end

"""
A fast ray-triangle intersection algorithm by Moller and Trumbore.
Ref: Tomas Möller and Ben Trumbore, "Fast, minimum storage ray-triangle intersection" (1997) Journal of Graphics Tools.
"""
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


function dist2LineTriangle(x,y)
    fzero = zero(eltype(eltype(x)))
    o = x[1]
    d = x[2] - x[1]
    a = y[1]
    b = y[2]
    c = y[3]
    result = moller_trumbore_intersect(o, d, a, b, c)
    if result.intersect == true
        return fzero
    else
        # no intersection or seg and triangle are parallel.
        # min distance is on an edge of the triangle.
        @inline dab = dist2LineLine(x,SA[a,b])
        @inline dbc = dist2LineLine(x,SA[b,c])
        @inline dcd = dist2LineLine(x,SA[c,a])
        @inline dx1y = dist2PointTriangle(SA[x[1],],y)
        @inline dx2y = dist2PointTriangle(SA[x[2],],y)
        return Base.FastMath.min_fast(
            dab,
            dbc,
            dcd,
            dx1y,
            dx2y,
        )
    end
end



function dist2TriangleTriangle(x,y)
    T = eltype(eltype(x))
    fzero = zero(T)
    d2 = typemax(T)
    ai = 1
    bi = 2
    ci = 3
    for i in 1:3
        @inline this_d2 = dist2LineTriangle(SA[x[ai],x[bi]],y)
        d2 = Base.FastMath.min_fast(d2, this_d2)
        d2 == fzero && return d2
        ai, bi, ci = bi, ci, ai
    end
    for i in 1:3
        @inline this_d2 = dist2LineTriangle(SA[y[ai],y[bi]],x)
        d2 = Base.FastMath.min_fast(d2, this_d2)
        d2 == fzero && return d2
        ai, bi, ci = bi, ci, ai
    end
    return d2
end



