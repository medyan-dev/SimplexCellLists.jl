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
    # `‖P×Q‖^2`, which doesn't cancel for nearly parallel line segments like
    # `a*c - b*b` does, and is the `n` of `dist_sqr_r_s_t`.
    n = cross3(P, Q)
    Δ = dot3(n, n)
    # The closest points on `y` of the end points of `x`, the end point
    # candidates of `dist_sqr_r_s_t`, and the closest points on `x` of the end
    # points of `y`, its `s` at `t = 0` and `t = 1`.
    s0 = clamp01nan(d/a)
    s1 = clamp01nan((b + d)/a)
    t0 = clamp01nan(-e/c)
    t1 = clamp01nan((b - e)/c)
    ucl, fp1cl, fq0cl, fq1cl = _line_line_u_f(L, s[1], t[1], r)
    # The separation between the points at `s` on `x` and `t` on `y` is
    # `Q0mP0 + t*Q - s*P`.
    us0, fp1s0, fq0s0, fq1s0 = _line_line_u_f(L, s0, zero(T), fma3(-s0, P, Q0mP0))
    us1, fp1s1, fq0s1, fq1s1 = _line_line_u_f(L, s1, one(T), fma3(-s1, P, y[2] - x[1]))
    ut0, fp1t0, fq0t0, fq1t0 = _line_line_u_f(L, zero(T), t0, fma3(t0, Q, Q0mP0))
    ut1, fp1t1, fq0t1, fq1t1 = _line_line_u_f(L, one(T), t1, fma3(t1, Q, Q0mP0 - P))
    uother = us0 + us1 + ut0 + ut1
    fp1other = fp1s0 + fp1s1 + fp1t0 + fp1t1
    fq0other = fq0s0 + fq0s1 + fq0t0 + fq0t1
    fq1other = fq1s0 + fq1s1 + fq1t0 + fq1t1

    switchover_scale = T(policy.switchover_scale_unitless)*inv(L)^2 # units of 1/length^2
    p = _sqrt_fast(a)
    q = _sqrt_fast(c)
    pq = p*q
    # `1/(p*q)`, or zero when a line segment has zero length, or is so short
    # that the inverse would overflow. Such a line segment has no direction,
    # so it is taken as parallel to the other, the limit as its length goes
    # to zero, with no force from the switchover.
    ipq = _vifelse(pq > floatmin(T), inv(pq), zero(T))
    # Smoothly switch from 0 for parallel lines to exactly 1 for x ≥ 3, where
    # x = Δ/(p*q)*switchover_scale, with continuous first and second
    # derivatives: switchover = 1 - (1 - x/3)^3 for x < 3.
    w = relu(1 - Δ*ipq*switchover_scale*T(1//3))
    switchover = 1 - w^3
    d_switchover = w^2 # derivative of switchover with respect to x
    u = switchover*ucl + T(1//2)*(1-switchover)*uother
    g = b*ipq
    # `1 + g^2` gives `q/p + b*g/a = c/(p*q)*(1 + g^2)` and
    # `p/q + b*g/c = a/(p*q)*(1 + g^2)` with no other divisions.
    h = 1 + g^2
    uswitched = (ucl + -T(1//2)*uother)
    fp1 = switchover*fp1cl + T(1//2)*(1-switchover)*fp1other
    d_switchover_dp1 = -switchover_scale*d_switchover*((c*ipq*h)*P - 2*g*Q)
    fp1 += d_switchover_dp1*uswitched

    fq0 = switchover*fq0cl + T(1//2)*(1-switchover)*fq0other
    fq1 = switchover*fq1cl + T(1//2)*(1-switchover)*fq1other
    d_switchover_dq0 = switchover_scale*d_switchover*((a*ipq*h)*Q - 2*g*P)
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

# The vertex positions `(x, y)` of `edge.a` and `edge.b`, relative to the first
# vertex of `edge.a`: subtracted in the position type, then converted to `T`.
# `x[1]` is exactly `zero(Vec3{T})`.
Base.@propagate_inbounds function _edge_simplexes(pos, edge::NeighborListEdge, ::Type{T}) where T
    xa = load_positions(pos, edge.a)
    N = length(xa)
    x = SVector{N, Vec3{T}}(ntuple(i -> i == 1 ? zero(Vec3{T}) : map(T, xa[i] - xa[1]), Val(N)))
    y = map(p -> map(T, p - xa[1]), load_positions(pos, edge.b))
    x, y
end

# Add the forces `fx` and `fy` on the vertices of the index parts `a` and
# `b`, and the energy `E`, to `force_energy`.
Base.@propagate_inbounds function _add_edge_force_energy!(force_energy::ForceEnergy, a, b, E, fx, fy)
    for (n, idx) in enumerate(_vertex_idxs(a))
        add_bead_force!(force_energy, idx, fx[n])
    end
    for (n, idx) in enumerate(_vertex_idxs(b))
        add_bead_force!(force_energy, idx, fy[n])
    end
    add_energy!(force_energy, E)
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
    x, y = _edge_simplexes(pos, edge, T)
    drst = dist_sqr_r_s_t(x, y)
    L = T(edge.L)
    if drst[1] < L^2
        E, fx, fy = nl_edge_force(policy, x, y, drst, L, edge.params)
        _add_edge_force_energy!(force_energy, edge.a, edge.b, E, fx, fy)
    end
    nothing
end

# Split 1:n into nthreads contiguous chunks with sizes differing by at most
# one, and return chunk `chunk`. Chunks past the end are empty.
function _chunk_range(n, chunk, nthreads)
    q, rem = divrem(n, nthreads)
    start = (chunk - 1)*q + min(chunk - 1, rem) + 1
    stop = chunk*q + min(chunk, rem)
    start:stop
end

# The scalar reference path, one edge at a time: works for any `calc_type`.
function nl_forces_scalar!(force_energy::ForceEnergy, pos, nl, policy::CollidePolicy, calc_type::Type{T}, chunk, nthreads) where T
    for i in _chunk_range(length(nl), chunk, nthreads)
        nl_edge_forces!(force_energy, pos, nl[i], policy, calc_type)
    end
    nothing
end

function nl_forces!(force_energy::ForceEnergy, pos, nl, policy::CollidePolicy, calc_type::Type{T}, chunk, nthreads) where T
    nl_forces_scalar!(force_energy, pos, nl, policy, calc_type, chunk, nthreads)
end

const COLLIDE_SIMD_WIDTH = 8

# Return the simplices as one simplex whose coordinates have a lane per
# simplex.
@inline function _simplex_lanes(simplices::Tuple{SVector{N}, Vararg{SVector{N}, M}}) where {N, M}
    SVector(ntuple(Val(N)) do v
        SVector(ntuple(Val(3)) do axis
            SIMD.Vec(ntuple(l -> simplices[l][v][axis], Val(M + 1)))
        end)
    end)
end

_nvertices(::Type{PointIdxPart}) = 1
_nvertices(::Type{CLineIdxPart}) = 2
_nvertices(::Type{LineIdxPart}) = 2
_nvertices(::Type{TriangleIdxPart}) = 3

#=
A structure of arrays buffer of `W` neighbor-list edges with index parts of
types `A` and `B`, and params of type `P`, in the calculation type `T`. Lane
`l` holds an edge's positions `x[2:N]` and `y`, its `L`, its index parts, and
its params. Its forces `fx[2:N]` and `fy` are stored over `x[2:N]` and `y`,
so only `fx[1]` and the energy need columns of their own.

The `T` values share one `Memory`, a column of `W` lanes per coordinate, so
LLVM can tell from the indexes that they don't overlap. With separate
buffers it needs runtime checks, and then finds vectorizing the force loop
not worth it. The getters and setters are `@inline` and unchecked.
=#
struct EdgeBuffer{T, W, A, B, P}
    data::Memory{T}
    a::Memory{A}
    b::Memory{B}
    params::Memory{P}
end

# The first columns of the positions `x[2:N]` and `y`, `L`, the force
# `fx[1]`, and the energy, and the number of columns.
function _columns(::Type{<:EdgeBuffer{T, W, A, B}}) where {T, W, A, B}
    N = _nvertices(A)
    M = _nvertices(B)
    x = 0
    y = x + 3*(N - 1)
    L = y + 3*M
    fx1 = L + 1
    E = fx1 + 3
    (; x, y, L, fx1, E, total = E + 1)
end

function EdgeBuffer{T, W, A, B, P}() where {T, W, A, B, P}
    EdgeBuffer{T, W, A, B, P}(
        Memory{T}(undef, W*_columns(EdgeBuffer{T, W, A, B, P}).total),
        Memory{A}(undef, W),
        Memory{B}(undef, W),
        Memory{P}(undef, W),
    )
end

# Store `x[2:N]` in lane `l`. `x[1]` is not stored.
@inline function _set_x!(b::EdgeBuffer{T, W, A}, l, x) where {T, W, A}
    c = _columns(typeof(b))
    for v in 2:_nvertices(A), axis in 1:3
        @inbounds b.data[l + W*(c.x + 3*(v - 2) + axis - 1)] = x[v][axis]
    end
end

@inline function _set_y!(b::EdgeBuffer{T, W, A, B}, l, y) where {T, W, A, B}
    c = _columns(typeof(b))
    for v in 1:_nvertices(B), axis in 1:3
        @inbounds b.data[l + W*(c.y + 3*(v - 1) + axis - 1)] = y[v][axis]
    end
end

# Load `x` from lane `l`, with `x[1]` exactly `zero(Vec3{T})`.
@inline function _get_x(b::EdgeBuffer{T, W, A}, l) where {T, W, A}
    c = _columns(typeof(b))
    SVector{_nvertices(A), Vec3{T}}(ntuple(Val(_nvertices(A))) do v
        v == 1 && return zero(Vec3{T})
        SVector(ntuple(axis -> @inbounds(b.data[l + W*(c.x + 3*(v - 2) + axis - 1)]), Val(3)))
    end)
end

@inline function _get_y(b::EdgeBuffer{T, W, A, B}, l) where {T, W, A, B}
    c = _columns(typeof(b))
    SVector{_nvertices(B), Vec3{T}}(ntuple(Val(_nvertices(B))) do v
        SVector(ntuple(axis -> @inbounds(b.data[l + W*(c.y + 3*(v - 1) + axis - 1)]), Val(3)))
    end)
end

# Store `edge` in lane `l`, with its relative positions `x` and `y`, where
# `x[1]` is zero and not stored.
@inline function _set_edge!(b::EdgeBuffer{T, W}, l, x, y, edge::NeighborListEdge) where {T, W}
    _set_x!(b, l, x)
    _set_y!(b, l, y)
    @inbounds b.data[l + W*_columns(typeof(b)).L] = T(edge.L)
    @inbounds b.a[l] = edge.a
    @inbounds b.b[l] = edge.b
    @inbounds b.params[l] = edge.params
end

@inline _get_L(b::EdgeBuffer{T, W}, l) where {T, W} = @inbounds b.data[l + W*_columns(typeof(b)).L]
@inline _get_a(b::EdgeBuffer, l) = @inbounds b.a[l]
@inline _get_b(b::EdgeBuffer, l) = @inbounds b.b[l]
@inline _get_params(b::EdgeBuffer, l) = @inbounds b.params[l]

# Copy the edge in lane `src` to lane `dst`.
@inline function _copy_edge!(b::EdgeBuffer{T, W}, dst, src) where {T, W}
    for j in 0:_columns(typeof(b)).L
        @inbounds b.data[dst + W*j] = b.data[src + W*j]
    end
    @inbounds b.a[dst] = b.a[src]
    @inbounds b.b[dst] = b.b[src]
    @inbounds b.params[dst] = b.params[src]
end

# Store the energy and forces of the edge in lane `l`, over its positions.
@inline function _set_energy_forces!(b::EdgeBuffer{T, W}, l, E, fx, fy) where {T, W}
    c = _columns(typeof(b))
    _set_x!(b, l, fx)
    _set_y!(b, l, fy)
    for axis in 1:3
        @inbounds b.data[l + W*(c.fx1 + axis - 1)] = fx[1][axis]
    end
    @inbounds b.data[l + W*c.E] = E
end

@inline function _get_energy_forces(b::EdgeBuffer{T, W}, l) where {T, W}
    c = _columns(typeof(b))
    fx1 = SVector(ntuple(axis -> @inbounds(b.data[l + W*(c.fx1 + axis - 1)]), Val(3)))
    (@inbounds b.data[l + W*c.E]), Base.setindex(_get_x(b, l), fx1, 1), _get_y(b, l)
end

# Compute the forces of the `W` buffered edges, and add the first `n` to
# `force_energy`.
@inline function _flush_forces!(force_energy, policy, b::EdgeBuffer{T, W}, n) where {T, W}
    # Loaded to local variables first, so the force loop only loads from `b.data`.
    params = ntuple(l -> _get_params(b, l), Val(W))
    # A plain loop over the lanes that LLVM vectorizes.
    for l in 1:W
        x = _get_x(b, l)
        y = _get_y(b, l)
        E, fx, fy = nl_edge_force(policy, x, y, dist_sqr_r_s_t(x, y), _get_L(b, l), params[l])
        _set_energy_forces!(b, l, E, fx, fy)
    end
    for l in 1:n
        E, fx, fy = _get_energy_forces(b, l)
        @inbounds _add_edge_force_energy!(force_energy, _get_a(b, l), _get_b(b, l), E, fx, fy)
    end
end

#=
The vectorized path for `Float32` and `Float64`, with the same results as the
scalar path, as it calls the same functions on the same values.

The range check is done `W` edges at a time, with `dist_sqr` on `SIMD.Vec`
lanes, which is bit identical to the scalar `dist_sqr_r_s_t` by its contract.
Each edge in range is copied to the next lane of a structure of arrays
buffer, without `x[1]`, which is zero. When `W` are buffered, their forces
are computed in a plain loop over the lanes, then added to `force_energy` in
edge order. The last edges, fewer than `W`, take the scalar path.

The `Memory` buffers never leave the function, so LLVM keeps them on the
stack, as long as every function using them is inlined and they are under
about 2000 bytes.
=#
function nl_forces!(force_energy::ForceEnergy, pos, nl::AbstractVector{NeighborListEdge{EA, EB, P}}, policy::CollidePolicy, calc_type::Type{T}, chunk, nthreads) where {EA, EB, P, T <: Union{Float32, Float64}}
    W = COLLIDE_SIMD_WIDTH
    buf = EdgeBuffer{T, W, EA, EB, P}()
    n = 0
    range = _chunk_range(length(nl), chunk, nthreads)
    ngroups = length(range) ÷ W
    for k in first(range):W:(first(range) + W*ngroups - 1)
        lanes = ntuple(l -> @inbounds(nl[k + l - 1]), Val(W))
        simplexes = map(edge -> @inbounds(_edge_simplexes(pos, edge, T)), lanes)
        x = _simplex_lanes(map(first, simplexes))
        y = _simplex_lanes(map(last, simplexes))
        # Converted to `T` as a vector. Converting each lane made LLVM load
        # them with a slow gather instruction.
        L = convert(SIMD.Vec{W, T}, SIMD.Vec(map(edge -> edge.L, lanes)))
        hits = SIMD.bitmask(dist_sqr(x, y) < L^2)
        while !iszero(hits)
            i = k + trailing_zeros(hits)
            hits &= hits - one(hits)
            edge = @inbounds nl[i]
            xi, yi = @inbounds _edge_simplexes(pos, edge, T)
            n += 1
            _set_edge!(buf, n, xi, yi, edge)
            if n == W
                _flush_forces!(force_energy, policy, buf, W)
                n = 0
            end
        end
    end
    if n > 0
        # Fill the unused lanes with copies of the last edge, so they hold
        # valid values.
        for l in (n + 1):W
            _copy_edge!(buf, l, n)
        end
        _flush_forces!(force_energy, policy, buf, n)
    end
    # Padding the last group instead, with an index clamp, made LLVM
    # vectorize the address computations and load the positions with slow
    # gather instructions.
    for i in (first(range) + W*ngroups):last(range)
        @inbounds nl_edge_forces!(force_energy, pos, nl[i], policy, calc_type)
    end
    nothing
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
