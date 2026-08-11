using Test
using StaticArrays
using Rotations: RotMatrix
using SimplexCellLists
using SimplexCellLists: dist_sqr, nl_edge_forces!, NeighborListEdge

flat_force(force_energy) = reinterpret(Float64, get_force(force_energy))

"""
Error if the collide force and energy violates a property
`a` and `b` are vectors of three vectors
"""
function check_collide_props(a, b)
    policy = DefaultCollisionPolicy()
    positions = SVector{3, Float64}[]
    a_i, b_i = map((a,b)) do y
        local idxs = UInt32[]
        if rand(Bool)
            # Add padding
            push!(positions, SA[NaN,NaN,NaN])
        end
        for x in y
            while rand(Bool)
                # Add padding
                push!(positions, SA[NaN,NaN,NaN])
            end
            push!(positions, x)
            push!(idxs, UInt32(length(positions)))
        end
        if rand(Bool)
            # Add padding
            push!(positions, SA[NaN,NaN,NaN])
        end
        if length(y) == 1
            PointIdxPart(idxs...)
        elseif length(idxs) == 2
            if rand(Bool) || idxs[2] != idxs[1] + 1
                LineIdxPart(idxs...)
            else
                CLineIdxPart(idxs[1])
            end
        elseif length(idxs) == 3
            TriangleIdxPart(idxs...)
        else
            error("unsupported shape")
        end
    end
    # flip order if needed
    if a_i isa LineIdxPart && b_i isa CLineIdxPart
        a_i, b_i = b_i, a_i
    end
    k = 10f0
    dcl = sqrt(dist_sqr(SVector((a...,)), SVector((b...,))))
    # Out of range should have zero energy
    force_energy = ForceEnergyFloat64(length(positions))
    nl_edge_forces!(force_energy, positions, NeighborListEdge(a_i, b_i, Float32(dcl*0.9), DefaultPairParams(k)), policy, Float64)
    @test iszero(get_energy(force_energy))
    @test all(iszero, get_force(force_energy))

    # in range should have non zero energy proportional to k
    force_energy = ForceEnergyFloat64(length(positions))
    nl_edge_forces!(force_energy, positions, NeighborListEdge(a_i, b_i, Float32(dcl*1.5), DefaultPairParams(k)), policy, Float64)
    e1 = get_energy(force_energy)
    force_energy = ForceEnergyFloat64(length(positions))
    nl_edge_forces!(force_energy, positions, NeighborListEdge(a_i, b_i, Float32(dcl*1.5), DefaultPairParams(10*k)), policy, Float64)
    e2 = get_energy(force_energy)
    @test e1 > 0
    @test 10*e1 ≈ e2

    # Finite difference test
    L = Float32(dcl*1.5)
    force_energy = ForceEnergyFloat64(length(positions))
    nl_edge_forces!(force_energy, positions, NeighborListEdge(a_i, b_i, L, DefaultPairParams(k)), policy, Float64)
    e0 = get_energy(force_energy)
    f = flat_force(force_energy)
    for i in 1:length(positions)*3
        local new_pos = copy(positions)
        local Δx = 1E-5
        reinterpret(Float64, new_pos)[i] += Δx
        local force_energy = ForceEnergyFloat64(length(positions))
        nl_edge_forces!(force_energy, new_pos, NeighborListEdge(a_i, b_i, L, DefaultPairParams(k)), policy, Float64)
        local Δe = get_energy(force_energy) - e0
        @test abs(f[i]*Δx + Δe) < 1.5E-7
    end

    # Translation invariance
    L = Float32(dcl*1.5)
    force_energy = ForceEnergyFloat64(length(positions))
    nl_edge_forces!(force_energy, positions, NeighborListEdge(a_i, b_i, L, DefaultPairParams(k)), policy, Float64)
    e0 = get_energy(force_energy)
    f = flat_force(force_energy)
    new_pos = positions .+ (SA[3.0, 2.0, 1.0],)
    force_energy = ForceEnergyFloat64(length(positions))
    nl_edge_forces!(force_energy, new_pos, NeighborListEdge(a_i, b_i, L, DefaultPairParams(k)), policy, Float64)
    @test e0 ≈ get_energy(force_energy)
    @test f ≈ flat_force(force_energy)

    # Rotation invariance
    L = Float32(dcl*1.5)
    force_energy = ForceEnergyFloat64(length(positions))
    nl_edge_forces!(force_energy, positions, NeighborListEdge(a_i, b_i, L, DefaultPairParams(k)), policy, Float64)
    e0 = get_energy(force_energy)
    f0 = get_force(force_energy)
    R = rand(RotMatrix{3})
    new_pos = (R,) .* positions
    force_energy = ForceEnergyFloat64(length(positions))
    nl_edge_forces!(force_energy, new_pos, NeighborListEdge(a_i, b_i, L, DefaultPairParams(k)), policy, Float64)
    @test e0 ≈ get_energy(force_energy)
    @test (R,) .* f0 ≈ get_force(force_energy)
    return
end


@testset "collide forces" begin
    x = SA[1.0, 0.0, 0.0]
    y = SA[0.0, 1.0, 0.0]
    z = SA[0.0, 0.0, 1.0]
    check_collide_props(
        SA[0z],
        SA[x],
    )
    check_collide_props(
        SA[z],
        SA[x],
    )
    check_collide_props(
        SA[0z],
        SA[x, 2x],
    )
    check_collide_props(
        SA[0z],
        SA[x, 1.001x],
    )
    check_collide_props(
        SA[1.002x],
        SA[x, 1.001x],
    )
    check_collide_props(
        SA[2y],
        SA[-6x, 6x],
    )
    for i in 1:5000
        check_collide_props(
            SA[0.1*y + randn()*x + randn()*z],
            SA[randn()*x + randn()*z, randn()*x + randn()*z, randn()*x + randn()*z],
        )
    end
    check_collide_props(
        SA[0*x, -1x],
        SA[x, 2x],
    )
    check_collide_props(
        SA[0y, -1y],
        SA[x, 2x],
    )
    check_collide_props(
        SA[0.1y+x, 0.1y+2x],
        SA[x, 2x],
    )
    check_collide_props(
        SA[0.1y+2x, 0.1y+x],
        SA[x, 2x],
    )
    check_collide_props(
        SA[0.1y-2z, 0.1y+2z],
        SA[-2x, 2x],
    )
    check_collide_props(
        SA[20y-2z, 20y+2z],
        SA[-2x, 2x],
    )
    for i in 1:5000
        check_collide_props(
            SA[0.1*y + randn()*x + randn()*z, 0.1*y + randn()*x + randn()*z],
            SA[randn()*x + randn()*z, randn()*x + randn()*z],
        )
    end
end

@testset "collide_forces! integration" begin
    policy = DefaultCollisionPolicy()
    # 3×3×4 grid, spacing 0.4: with radius 0.3 overlaps are guaranteed, but
    # no two objects touch exactly, so all forces stay finite
    positions = [0.4*SA[mod(i-1, 3), mod((i-1)÷3, 3), (i-1)÷9] for i in 1:30]
    params = DefaultObjectParams(10.0f0, UInt32(1), UInt32(0))
    inputs = NeighborListInputs(policy;
        points = [PointIdxPart(i) for i in 1:10],
        p_radius = fill(0.3f0, 10),
        p_params = fill(params, 10),
        clines = [CLineIdxPart(i) for i in 11:2:17],
        c_radius = fill(0.3f0, 4),
        c_params = fill(params, 4),
        lines = [LineIdxPart(19, 20), LineIdxPart(21, 22)],
        l_radius = fill(0.3f0, 2),
        l_params = fill(params, 2),
        triangles = [TriangleIdxPart(23, 24, 25)],
        t_radius = [0.5f0],
        t_params = [params],
    )
    nl = NeighborLists(policy)
    setup_neighbors_sort_sweep!(nl, positions, inputs)
    # every edge type participates
    for list in (nl.PPNL, nl.PCNL, nl.PLNL, nl.PTNL, nl.CCNL, nl.CLNL, nl.LLNL)
        @test !isempty(list)
    end

    fe = ForceEnergyFloat64(length(positions))
    collide_forces!(fe, positions, nl, Float64)
    @test get_energy(fe) > 0
    @test all(f -> all(isfinite, f), get_force(fe))

    # matches summing nl_edge_forces! over every edge
    fe_manual = ForceEnergyFloat64(length(positions))
    for list in (nl.PPNL, nl.PCNL, nl.PLNL, nl.PTNL, nl.CCNL, nl.CLNL, nl.LLNL)
        for edge in list
            nl_edge_forces!(fe_manual, positions, edge, policy, Float64)
        end
    end
    @test get_energy(fe_manual) == get_energy(fe)
    @test get_force(fe_manual) == get_force(fe)

    # the chunked path visits every edge exactly once
    for nthreads in (1, 2, 3, 7, 100)
        local fe_chunked = ForceEnergyFloat64(length(positions))
        for chunk in 1:nthreads
            collide_forces!(fe_chunked, positions, nl, Float64; chunk, nthreads)
        end
        @test get_energy(fe_chunked) ≈ get_energy(fe)
        @test get_force(fe_chunked) ≈ get_force(fe)
    end
end
nothing
