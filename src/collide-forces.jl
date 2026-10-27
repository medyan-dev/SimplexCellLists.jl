"""
    nl_edge_force(policy::CollidePolicy, x::Simplex{N,T}, y::Simplex{M,T}, (d2, r, s, t), L::T, params)::Tuple{T, Simplex{N,T}, Simplex{M,T}}

Return `(energy, fx, fy)`, the collide energy of one neighbor-list edge and
the forces on the vertices of `x` and `y`.

- `x` and `y` are the vertex positions of `edge.a` and `edge.b`, relative to
  the first vertex of `x`, converted to the calculation type `T`. `x[1]` is
  always `zero(Vec3{T})`, positive zero.
- `(d2, r, s, t)` is `dist_sqr_r_s_t(x, y)`, see [`dist_sqr_r_s_t`](@ref).
- `L` is `edge.L` converted to `T`, and `params` is `edge.params`.

The hook is only called for edges in range, `d2 < L^2`.

The kind of edge is given by `N` and `M`: 1 for a point, 2 for a line
segment, 3 for a triangle.

Methods dispatch on `policy`, so a custom `CollidePolicy` with its own
`PairParams` can define its own force laws. Methods should be `@inline`,
branch free, and do no memory access, so the forces of many edges can be
vectorized.
"""
function nl_edge_force end

@inline function nl_edge_force(policy::DefaultCollidePolicy, x::Simplex{N, T}, y::Simplex{M, T}, (d2, r, s, t), L::T, params::DefaultPairParams) where {N, M, T}
    k = T(params.k)
    d = _sqrt_fast(d2)
    ΔL = L - d
    # The force on the closest point of `y`, and its negative on `x`, spread
    # over the vertices by their barycentric weights.
    F = scale3(k*ΔL*inv(d), r)
    wx = SVector(one(T) - sum(s), s...)
    wy = SVector(one(T) - sum(t), t...)
    k*T(1//2)*ΔL^2, map(w -> scale3(-w, F), wx), map(w -> scale3(w, F), wy)
end

# Energy, and forces on x[2], y[1] and y[2], without the `k` factor, of the
# spring of rest length `L` between the points at `s` on `x` and `t` on `y`,
# separated by `r`.
@inline function _line_line_u_f(L::T, s::T, t::T, r) where T
    d = norm_fast(r)
    # A boundary critical point can be out of range, where the spring must
    # not pull, and the closest point can be by rounding.
    ΔL = relu(L - d)
    F = scale3(ΔL*inv(d), r)
    T(1//2)*ΔL^2, scale3(-s, F), scale3(one(T) - t, F), scale3(t, F)
end

# The closest point alone gives a discontinuous force between nearly parallel
# line segments, so it is blended with the four boundary critical points.
@inline function nl_edge_force(policy::DefaultCollidePolicy, x::LineSeg{T}, y::LineSeg{T}, (d2, r, s, t), L::T, params::DefaultPairParams) where T
    P = x[2] - x[1]
    Q = y[2] - y[1]
    Q0mP0 = y[1] - x[1]
    a = dot3(P, P)
    b = dot3(P, Q)
    c = dot3(Q, Q)
    d = dot3(P, Q0mP0)
    e = dot3(Q, Q0mP0)
    Δ = a*c - b*b
    inva = inv(a)
    invc = inv(c)
    # The boundary critical points, as in `dist_sqr_r_s_t`
    s0 = clamp01nan(d*inva)
    s1 = clamp01nan((b+d)*inva)
    t0 = clamp01nan(-e*invc)
    t1 = clamp01nan((b-e)*invc)
    ucl, fp1cl, fq0cl, fq1cl = _line_line_u_f(L, s[1], t[1], r)
    # The separation between the points at `s` on `x` and `t` on `y` is
    # `Q0mP0 + t*Q - s*P`.
    us0, fp1s0, fq0s0, fq1s0 = _line_line_u_f(L, s0, zero(T), Q0mP0 - scale3(s0, P))
    us1, fp1s1, fq0s1, fq1s1 = _line_line_u_f(L, s1, one(T), (y[2] - x[1]) - scale3(s1, P))
    ut0, fp1t0, fq0t0, fq1t0 = _line_line_u_f(L, zero(T), t0, Q0mP0 + scale3(t0, Q))
    ut1, fp1t1, fq0t1, fq1t1 = _line_line_u_f(L, one(T), t1, (Q0mP0 - P) + scale3(t1, Q))
    uother = us0 + us1 + ut0 + ut1
    fp1other = fp1s0 + fp1s1 + fp1t0 + fp1t1
    fq0other = fq0s0 + fq0s1 + fq0t0 + fq0t1
    fq1other = fq1s0 + fq1s1 + fq1t0 + fq1t1

    switchover_scale = T(policy.switchover_scale_unitless)*inv(L)^2 # units of 1/length^2
    # Smoothly switch from 0 for parallel lines to exactly 1 for x ≥ 3, where
    # x = Δ/sqrt(a*c)*switchover_scale, with continuous first and second
    # derivatives: switchover = 1 - (1 - x/3)^3 for x < 3.
    w = relu(1 - Δ/_sqrt_fast(a*c)*switchover_scale*T(1//3))
    switchover = 1 - w^3
    d_switchover = w^2 # derivative of switchover with respect to x
    u = switchover*ucl + T(1//2)*(1-switchover)*uother
    p = _sqrt_fast(a)
    q = _sqrt_fast(c)
    g = b/(p*q)
    uswitched = (ucl + -T(1//2)*uother)
    fp1 = switchover*fp1cl + T(1//2)*(1-switchover)*fp1other
    d_switchover_dp1 = -switchover_scale*d_switchover*((q/p+b*g/a)*P - 2*g*Q)
    fp1 += d_switchover_dp1*uswitched

    fq0 = switchover*fq0cl + T(1//2)*(1-switchover)*fq0other
    fq1 = switchover*fq1cl + T(1//2)*(1-switchover)*fq1other
    d_switchover_dq0 = switchover_scale*d_switchover*((p/q+b*g/c)*Q - 2*g*P)
    d_switchover_dq1 = - d_switchover_dq0
    fq0 += d_switchover_dq0*uswitched
    fq1 += d_switchover_dq1*uswitched

    k = T(params.k)
    # The force on x[1] balances the others, by translation invariance.
    k*u, SVector(-k*(fp1+fq0+fq1), k*fp1), SVector(k*fq0, k*fq1)
end

_vertex_idxs(x::PointIdxPart) = (x.i,)
_vertex_idxs(x::CLineIdxPart) = (x.i, x.i + UInt32(1))
_vertex_idxs(x::LineIdxPart) = (x.i, x.j)
_vertex_idxs(x::TriangleIdxPart) = (x.i, x.j, x.k)

# `x` relative to its own first vertex, which is exactly `zero(Vec3{T})`.
@inline function _relative_x(::Type{T}, v::SVector{N}) where {T, N}
    SVector{N, Vec3{T}}(ntuple(i -> i == 1 ? zero(Vec3{T}) : map(T, v[i] - v[1]), Val(N)))
end

"""
    nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge, policy::CollidePolicy, calc_type::Type{T})::Nothing where T

Add the collide force and energy of a single neighbor-list edge to
`force_energy`, doing the calculations in the floating-point type `calc_type`.

This is the scalar path of [`collide_forces!`](@ref) for one edge: the range
check `dist_sqr < L^2`, then the force law [`nl_edge_force`](@ref). To change
the force law, define `nl_edge_force` methods, not methods of this function.
"""
Base.@propagate_inbounds function nl_edge_forces!(force_energy::ForceEnergy, pos, edge::NeighborListEdge, policy::CollidePolicy, calc_type::Type{T}) where T
    xa = load_positions(pos, edge.a)
    x = _relative_x(T, xa)
    # Subtracted in the position type, then converted to `T`.
    y = map(p -> map(T, p - xa[1]), load_positions(pos, edge.b))
    drst = dist_sqr_r_s_t(x, y)
    L = T(edge.L)
    if drst[1] < L^2
        E, fx, fy = nl_edge_force(policy, x, y, drst, L, edge.params)
        for (n, idx) in enumerate(_vertex_idxs(edge.a))
            add_bead_force!(force_energy, idx, fx[n])
        end
        for (n, idx) in enumerate(_vertex_idxs(edge.b))
            add_bead_force!(force_energy, idx, fy[n])
        end
        add_energy!(force_energy, E)
    end
    nothing
end

# The scalar reference path, one edge at a time: works for any `calc_type`.
function nl_forces_scalar!(force_energy::ForceEnergy, pos, nl, policy::CollidePolicy, calc_type::Type{T}, chunk, nthreads) where T
    # Split 1:length(nl) into nthreads contiguous chunks with sizes differing
    # by at most one; chunks past the end are empty.
    q, rem = divrem(length(nl), nthreads)
    start = (chunk - 1)*q + min(chunk - 1, rem) + 1
    stop = chunk*q + min(chunk, rem)
    for i in start:stop
        nl_edge_forces!(force_energy, pos, nl[i], policy, calc_type)
    end
    nothing
end

function nl_forces!(force_energy::ForceEnergy, pos, nl, policy::CollidePolicy, calc_type::Type{T}, chunk, nthreads) where T
    nl_forces_scalar!(force_energy, pos, nl, policy, calc_type, chunk, nthreads)
end

"""
    collide_forces!(force_energy::ForceEnergy, pos, s::NeighborLists, calc_type::Type{T}; chunk=1, nthreads=1)::Nothing where T

Add the collide forces and energy of every neighbor-list edge in `s` to
`force_energy`, doing the calculations in the floating-point type `calc_type`.

Each edge in range, with `dist_sqr < L^2`, gets its force and energy from
[`nl_edge_force`](@ref), so a custom `CollidePolicy` can change the force laws.

For multithreading, split the work into `nthreads` chunks: thread `t` calls
`collide_forces!` with `chunk=t` and the same `nthreads`, accumulating into
its own `force_energy`, and the per-thread accumulators are merged afterwards
with [`combine_force_energy!`](@ref). Every edge is visited exactly once
across the chunks.
"""
function collide_forces!(force_energy::ForceEnergy, pos, s::NeighborLists, calc_type::Type{T}; chunk=1, nthreads=1) where T
    nl_forces!(force_energy, pos, s.PPNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.PCNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.PLNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.PTNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.CCNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.CLNL, s.policy, calc_type, chunk, nthreads)
    nl_forces!(force_energy, pos, s.LLNL, s.policy, calc_type, chunk, nthreads)
    nothing
end
