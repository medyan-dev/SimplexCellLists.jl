"""
    nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge, policy::CollidePolicy, calc_type)::Nothing

Add the collide force and energy of a single neighbor-list edge to `force_energy`.

Methods dispatch on the index part types of `edge` and on `policy`, so a custom
`CollidePolicy` with its own `PairParams` can define its own force laws.
"""
function nl_edge_forces! end

Base.@propagate_inbounds function nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge{PointIdxPart, PointIdxPart, DefaultPairParams}, policy::DefaultCollidePolicy, calc_type::T) where T
    x1 = pos[edge.a.i]
    x2 = pos[edge.b.i]
    d = map(calc_type, x1 - x2)
    d2cl = sum(abs2, d)
    L0 = calc_type(edge.L)
    if d2cl < L0^2
        L = _sqrt_fast(d2cl)
        k = calc_type(edge.params.k)
        ΔL = L0 - L
        E = k*calc_type(1//2)*ΔL^2
        F = (k*ΔL*inv(L))*d
        add_bead_force!(force_energy, edge.a.i, F)
        add_bead_force!(force_energy, edge.b.i, -F)
        add_energy!(force_energy, E)
    end
    nothing
end

Base.@propagate_inbounds function nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge{PointIdxPart, CLineIdxPart, DefaultPairParams}, policy::DefaultCollidePolicy, calc_type::T) where T
    x1 = pos[edge.a.i]
    y1 = pos[edge.b.i]
    y2 = pos[edge.b.i + UInt32(1)]
    r = map(calc_type, y2 - y1)
    p = map(calc_type, x1 - y1)
    t = clamp01nan((p⋅r)*inv(r⋅r))
    d = t*r - p
    d2cl = sum(abs2, d)
    L0 = calc_type(edge.L)
    if d2cl < L0^2
        L = _sqrt_fast(d2cl)
        k = calc_type(edge.params.k)
        ΔL = L0 - L
        E = k*calc_type(1//2)*ΔL^2
        F = (k*ΔL*inv(L))*d
        add_bead_force!(force_energy, edge.a.i, -F)
        add_bead_force!(force_energy, edge.b.i, (oneunit(t)-t)*F)
        add_bead_force!(force_energy, edge.b.i + UInt32(1), t*F)
        add_energy!(force_energy, E)
    end
    nothing
end

Base.@propagate_inbounds function nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge{PointIdxPart, LineIdxPart, DefaultPairParams}, policy::DefaultCollidePolicy, calc_type::T) where T
    x1 = pos[edge.a.i]
    y1 = pos[edge.b.i]
    y2 = pos[edge.b.j]
    r = map(calc_type, y2 - y1)
    p = map(calc_type, x1 - y1)
    t = clamp01nan((p⋅r)*inv(r⋅r))
    d = t*r - p
    d2cl = sum(abs2, d)
    L0 = calc_type(edge.L)
    if d2cl < L0^2
        L = _sqrt_fast(d2cl)
        k = calc_type(edge.params.k)
        ΔL = L0 - L
        E = k*calc_type(1//2)*ΔL^2
        F = (k*ΔL*inv(L))*d
        add_bead_force!(force_energy, edge.a.i, -F)
        add_bead_force!(force_energy, edge.b.i, (oneunit(t)-t)*F)
        add_bead_force!(force_energy, edge.b.j, t*F)
        add_energy!(force_energy, E)
    end
    nothing
end

"""
Core helper function for line-line collide force calculation.
Returns (do_force::Bool, u, fp1, fq0, fq1) where:
- do_force: whether forces should be applied
- u: energy (without k factor)
- fp1: force on P1 (without k factor)
- fq0: force on Q0 (without k factor)
- fq1: force on Q1 (without k factor)
Force on P0 can be computed as -(fp1+fq0+fq1)
"""
@inline function line_line_force_core(P0, P1, Q0, Q1, calc_type::T, switchover_scale_unitless, L0) where T
    P = map(calc_type, P1-P0)
    Q = map(calc_type, Q1-Q0)
    P0mQ0 = map(calc_type, P0-Q0)
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
    d2cl = relu(f+minr)
    if d2cl < L0^2
        #return the energy and force on each point give s and t
        # doesn't have the e.k factor
        function u_f(s,t)
            local Pcl = P0mQ0 + s*P
            local Qcl = t*Q
            local Rcl = Qcl - Pcl
            local L = norm_fast(Rcl)
            local ΔL = relu(L0 - L)
            local u = calc_type(1//2)*ΔL^2
            local FP = (-ΔL*inv(L))*Rcl
            (u, s*FP, -(1-t)*FP, -t*FP)
        end
        switchover_scale = calc_type(switchover_scale_unitless)*inv(L0)^2 # units of 1/length^2
        switchover = tanh(Δ/sqrt(a*c)*switchover_scale)
        #P0 is set as origin for now, subtract forces later to get final fp0
        ucl, fp1cl, fq0cl, fq1cl = u_f(mins,mint)
        us0, fp1s0, fq0s0, fq1s0 = u_f(s0,0)
        us1, fp1s1, fq0s1, fq1s1 = u_f(s1,1)
        ut0, fp1t0, fq0t0, fq1t0 = u_f(0,t0)
        ut1, fp1t1, fq0t1, fq1t1 = u_f(1,t1)
        uother = +(
            us0,
            us1,
            ut0,
            ut1,
        )
        fp1other = +(
            fp1s0,
            fp1s1,
            fp1t0,
            fp1t1,
        )
        fq0other = +(
            fq0s0,
            fq0s1,
            fq0t0,
            fq0t1,
        )
        fq1other = +(
            fq1s0,
            fq1s1,
            fq1t0,
            fq1t1,
        )
        u = switchover*ucl + calc_type(1//2)*(1-switchover)*uother
        p = sqrt(a)
        q = sqrt(c)
        g = b/(p*q)
        uswitched = (ucl + -calc_type(1//2)*uother)
        fp1 = switchover*fp1cl + calc_type(1//2)*(1-switchover)*fp1other
        d_switchover_dp1 = -switchover_scale*(1-switchover^2)*((q/p+b*g/a)*P - 2*g*Q)
        fp1 += d_switchover_dp1*uswitched

        fq0 = switchover*fq0cl + calc_type(1//2)*(1-switchover)*fq0other
        fq1 = switchover*fq1cl + calc_type(1//2)*(1-switchover)*fq1other
        d_switchover_dq0 = switchover_scale*(1-switchover^2)*((p/q+b*g/c)*Q - 2*g*P)
        d_switchover_dq1 = - d_switchover_dq0
        fq0 += d_switchover_dq0*uswitched
        fq1 += d_switchover_dq1*uswitched

        return true, u, fp1, fq0, fq1
    else
        return false, zero(L0), zero(P), zero(P), zero(P)
    end
end

Base.@propagate_inbounds function nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge{CLineIdxPart, CLineIdxPart, DefaultPairParams}, policy::DefaultCollidePolicy, calc_type::T) where T
    P0 = pos[edge.a.i]
    P1 = pos[edge.a.i + UInt32(1)]
    Q0 = pos[edge.b.i]
    Q1 = pos[edge.b.i + UInt32(1)]
    L0 = calc_type(edge.L)
    do_force, u, fp1, fq0, fq1 = line_line_force_core(P0, P1, Q0, Q1, calc_type, policy.switchover_scale_unitless, L0)
    if do_force
        ke = calc_type(edge.params.k)
        add_bead_force!(force_energy, edge.a.i, -ke*(fp1+fq0+fq1))
        add_bead_force!(force_energy, edge.a.i + UInt32(1), ke*fp1)

        add_bead_force!(force_energy, edge.b.i, ke*fq0)
        add_bead_force!(force_energy, edge.b.i + UInt32(1), ke*fq1)

        add_energy!(force_energy, u*ke)
    end
    nothing
end

Base.@propagate_inbounds function nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge{CLineIdxPart, LineIdxPart, DefaultPairParams}, policy::DefaultCollidePolicy, calc_type::T) where T
    P0 = pos[edge.a.i]
    P1 = pos[edge.a.i + UInt32(1)]
    Q0 = pos[edge.b.i]
    Q1 = pos[edge.b.j]
    L0 = calc_type(edge.L)
    do_force, u, fp1, fq0, fq1 = line_line_force_core(P0, P1, Q0, Q1, calc_type, policy.switchover_scale_unitless, L0)
    if do_force
        ke = calc_type(edge.params.k)
        add_bead_force!(force_energy, edge.a.i, -ke*(fp1+fq0+fq1))
        add_bead_force!(force_energy, edge.a.i + UInt32(1), ke*fp1)

        add_bead_force!(force_energy, edge.b.i, ke*fq0)
        add_bead_force!(force_energy, edge.b.j, ke*fq1)

        add_energy!(force_energy, u*ke)
    end
    nothing
end

Base.@propagate_inbounds function nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge{LineIdxPart, LineIdxPart, DefaultPairParams}, policy::DefaultCollidePolicy, calc_type::T) where T
    P0 = pos[edge.a.i]
    P1 = pos[edge.a.j]
    Q0 = pos[edge.b.i]
    Q1 = pos[edge.b.j]
    L0 = calc_type(edge.L)
    do_force, u, fp1, fq0, fq1 = line_line_force_core(P0, P1, Q0, Q1, calc_type, policy.switchover_scale_unitless, L0)
    if do_force
        ke = calc_type(edge.params.k)
        add_bead_force!(force_energy, edge.a.i, -ke*(fp1+fq0+fq1))
        add_bead_force!(force_energy, edge.a.j, ke*fp1)

        add_bead_force!(force_energy, edge.b.i, ke*fq0)
        add_bead_force!(force_energy, edge.b.j, ke*fq1)

        add_energy!(force_energy, u*ke)
    end
    nothing
end

#=
Based on https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06
=#
Base.@propagate_inbounds function nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge{PointIdxPart, TriangleIdxPart, DefaultPairParams}, policy::DefaultCollidePolicy, calc_type::T) where T
    P = pos[edge.a.i]
    B = pos[edge.b.i]
    E0 = map(calc_type, pos[edge.b.j] - B)
    E1 = map(calc_type, pos[edge.b.k] - B)
    BP = map(calc_type, B - P)
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
    if  r(s0,0) < minr
        minr = r(s0,0)
        mins = s0
        mint = zero(a)
    end
    if  r(0,t0) < minr
        minr = r(0,t0)
        mins = zero(a)
        mint = t0
    end
    if  r(sd,td) < minr
        minr = r(sd,td)
        mins = sd
        mint = td
    end
    d2cl = relu(f+minr)
    L0 = calc_type(edge.L)
    if d2cl < L0^2
        # dvec is the vector from P to closest point on triangle
        dvec = BP + mins*E0 + mint*E1
        L = _sqrt_fast(d2cl)
        k = calc_type(edge.params.k)
        ΔL = L0 - L
        E = k*calc_type(1//2)*ΔL^2
        F = (k*ΔL*inv(L))*dvec
        add_bead_force!(force_energy, edge.a.i, -F)
        # Forces on triangle vertices
        # Closest point on triangle is B + mins*E0 + mint*E1 = (1-mins-mint)*B + mins*(B+E0) + mint*(B+E1)
        # = (1-mins-mint)*y[1] + mins*y[2] + mint*y[3]
        w0 = oneunit(mins) - mins - mint
        add_bead_force!(force_energy, edge.b.i, w0*F)
        add_bead_force!(force_energy, edge.b.j, mins*F)
        add_bead_force!(force_energy, edge.b.k, mint*F)
        add_energy!(force_energy, E)
    end
    nothing
end

function nl_forces!(force_energy::ForceEnergy, pos, nl, policy::CollidePolicy, calc_type::T, chunk, nthreads) where T
    # Split 1:length(nl) into nthreads contiguous chunks with sizes differing
    # by at most one; chunks past the end are empty.
    q, r = divrem(length(nl), nthreads)
    start = (chunk - 1)*q + min(chunk - 1, r) + 1
    stop = chunk*q + min(chunk, r)
    for i in start:stop
        nl_edge_forces!(force_energy, pos, nl[i], policy, calc_type)
    end
    nothing
end

"""
    collide_forces!(force_energy::ForceEnergy, pos, s::NeighborLists, calc_type; chunk=1, nthreads=1)::Nothing

Add the collide forces and energy of every neighbor-list edge in `s` to
`force_energy`, doing the calculations in the floating-point type `calc_type`.

Each edge's contribution is computed by [`nl_edge_forces!`](@ref), so a custom
`CollidePolicy` can change the force laws.

For multithreading, split the work into `nthreads` chunks: thread `t` calls
`collide_forces!` with `chunk=t` and the same `nthreads`, accumulating into
its own `force_energy`, and the per-thread accumulators are merged afterwards
with [`combine_force_energy!`](@ref). Every edge is visited exactly once
across the chunks.
"""
function collide_forces!(force_energy::ForceEnergy, pos, s::NeighborLists, calc_type::T; chunk=1, nthreads=1) where T
    nl_forces!(force_energy, pos, s.PPNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.PCNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.PLNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.PTNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.CCNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.CLNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.LLNL, s.policy, calc_type, chunk, nthreads)
    nothing
end
