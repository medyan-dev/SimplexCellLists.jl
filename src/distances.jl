using LinearAlgebra

clamp01nan(x) = ifelse(isnan(x), zero(x), clamp(x,zero(x),one(x)))

dist2PointPoint(a,b) = (a[1] - b[1]) ⋅ (a[1] - b[1])

function dist2PointLine(a,b)
    r= b[2]-b[1]
    p= a[1]-b[1]
    # c= a[1]-b[2]
    s = clamp01nan(p⋅r/(r⋅r))
    d = s*r - p
    d ⋅ d
end


"""
    dist2LineLine(a,b)

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
function dist2LineLine(a,b)
    P0 = a[1]
    P1 = a[2]
    Q0 = b[1]
    Q1 = b[2]
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
    return max(0,f+min(
        r(s0,0),
        r(s1,1),
        r(0,t0),
        r(1,t1),
        r(sbar,tbar),
    ))
end

