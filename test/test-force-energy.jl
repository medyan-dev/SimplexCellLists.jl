using SimplexCellLists
using Test
using StaticArrays

@testset "ForceEnergy" begin
    @testset "ForceEnergyT{Float64}" begin
        nbeads = 4
        fe = ForceEnergyFloat64(nbeads)

        # Initially zero on every bead.
        for i in 1:nbeads
            @test get_bead_force(fe, i) === SA[0.0, 0.0, 0.0]
        end

        # Single push picks up on get_bead_force.
        add_bead_force!(fe, 2, SA[1.5, -2.0, 3.25])
        @test get_bead_force(fe, 1) === SA[0.0, 0.0, 0.0]
        @test get_bead_force(fe, 2) === SA[1.5, -2.0, 3.25]
        @test get_bead_force(fe, 3) === SA[0.0, 0.0, 0.0]

        # Accumulates across multiple pushes to the same bead.
        add_bead_force!(fe, 2, SA[0.5, 1.0, -0.25])
        @test get_bead_force(fe, 2) === SA[2.0, -1.0, 3.0]

        # Matches get_force! output for the same index.
        out = zeros(SVector{3, Float64}, nbeads)
        get_force!(fe, out)
        for i in 1:nbeads
            @test get_bead_force(fe, i) === out[i]
        end

        # Energy accumulates, and zero_energy! resets only the energy.
        add_energy!(fe, 1.5)
        add_energy!(fe, -0.25)
        @test get_energy(fe) === 1.25
        zero_energy!(fe)
        @test get_energy(fe) === 0.0
        @test get_bead_force(fe, 2) === SA[2.0, -1.0, 3.0]

        # zero_bead_force! wipes only that bead.
        add_bead_force!(fe, 3, SA[1.0, 1.0, 1.0])
        zero_bead_force!(fe, 2)
        @test get_bead_force(fe, 2) === SA[0.0, 0.0, 0.0]
        @test get_bead_force(fe, 3) === SA[1.0, 1.0, 1.0]

        # zero_force_energy! wipes everything.
        add_energy!(fe, 1.0)
        zero_force_energy!(fe)
        @test get_energy(fe) === 0.0
        for i in 1:nbeads
            @test get_bead_force(fe, i) === SA[0.0, 0.0, 0.0]
        end
    end

    @testset "ForceEnergyT{Float32}" begin
        # Returned type is Float64 even when the backing storage is Float32.
        fe = ForceEnergyFloat32(2)
        add_bead_force!(fe, 1, SA{Float32}[1.0f0, 2.0f0, -3.0f0])
        f = get_bead_force(fe, 1)
        @test f isa SVector{3, Float64}
        @test f == SA[1.0, 2.0, -3.0]
    end

    @testset "ForceEnergyFixedPoint" begin
        nbeads = 3
        fe = ForceEnergyFixedPoint{30, 30}(nbeads)

        for i in 1:nbeads
            @test get_bead_force(fe, i) === SA[0.0, 0.0, 0.0]
        end

        add_bead_force!(fe, 1, SA[0.25, -0.5, 1.0])
        add_bead_force!(fe, 3, SA[2.0, 0.0, -1.5])
        add_bead_force!(fe, 1, SA[0.125, 0.25, 0.0])

        f1 = get_bead_force(fe, 1)
        @test f1 isa SVector{3, Float64}
        @test f1 ≈ SA[0.375, -0.25, 1.0]
        @test get_bead_force(fe, 2) === SA[0.0, 0.0, 0.0]
        @test get_bead_force(fe, 3) ≈ SA[2.0, 0.0, -1.5]

        # Matches get_force! output.
        out = zeros(SVector{3, Float64}, nbeads)
        get_force!(fe, out)
        for i in 1:nbeads
            @test get_bead_force(fe, i) ≈ out[i]
        end

        # Energy accumulates in fixed point, and zero_energy! resets it.
        add_energy!(fe, 1.5)
        add_energy!(fe, -0.25)
        @test get_energy(fe) ≈ 1.25
        @test get_energy(fe, Float32) isa Float32
        zero_energy!(fe)
        @test get_energy(fe) === 0.0

        # zero_bead_force! wipes only that bead.
        zero_bead_force!(fe, 1)
        @test get_bead_force(fe, 1) === SA[0.0, 0.0, 0.0]
        @test get_bead_force(fe, 3) ≈ SA[2.0, 0.0, -1.5]

        # Zeroing wipes per-bead view too.
        add_energy!(fe, 1.0)
        zero_force_energy!(fe)
        @test get_energy(fe) === 0.0
        for i in 1:nbeads
            @test get_bead_force(fe, i) === SA[0.0, 0.0, 0.0]
        end

        # Quantization truncates toward zero symmetrically: adding a
        # contribution and its exact negation cancels exactly, even for values
        # that don't quantize exactly.
        v = SA[0.3, -1e-9, 7.123456789]
        add_bead_force!(fe, 1, v)
        add_bead_force!(fe, 1, -v)
        add_energy!(fe, 0.3)
        add_energy!(fe, -0.3)
        @test get_bead_force(fe, 1) === SA[0.0, 0.0, 0.0]
        @test get_energy(fe) === 0.0
    end

    @testset "ForceNoEnergyFixedPoint" begin
        nbeads = 3
        fe = ForceNoEnergyFixedPoint{30}(nbeads)

        @test get_nbeads(fe) == nbeads
        for i in 1:nbeads
            @test get_bead_force(fe, i) === SA[0.0, 0.0, 0.0]
        end

        add_bead_force!(fe, 1, SA[0.25, -0.5, 1.0])
        add_bead_force!(fe, 3, SA[2.0, 0.0, -1.5])
        add_bead_force!(fe, 1, SA[0.125, 0.25, 0.0])
        @test get_bead_force(fe, 1) ≈ SA[0.375, -0.25, 1.0]
        @test get_bead_force(fe, 2) === SA[0.0, 0.0, 0.0]
        @test get_bead_force(fe, 3) ≈ SA[2.0, 0.0, -1.5]

        # Energy is discarded: add_energy!/zero_energy! are no-ops, and reading
        # the energy errors instead of returning a made-up value.
        add_energy!(fe, 5.0)
        zero_energy!(fe)
        @test_throws MethodError get_energy(fe)

        # Matches get_force! output.
        out = zeros(SVector{3, Float64}, nbeads)
        get_force!(fe, out)
        for i in 1:nbeads
            @test get_bead_force(fe, i) ≈ out[i]
        end

        # zero_bead_force! wipes only that bead.
        zero_bead_force!(fe, 1)
        @test get_bead_force(fe, 1) === SA[0.0, 0.0, 0.0]
        @test get_bead_force(fe, 3) ≈ SA[2.0, 0.0, -1.5]

        zero_force_energy!(fe)
        for i in 1:nbeads
            @test get_bead_force(fe, i) === SA[0.0, 0.0, 0.0]
        end
    end

    @testset "combine_force_energy! bit-exact for matching ForceNoEnergyFixedPoint" begin
        nbeads = 4
        a = ForceNoEnergyFixedPoint{30}(nbeads)
        b = ForceNoEnergyFixedPoint{30}(nbeads)
        ref = ForceNoEnergyFixedPoint{30}(nbeads)

        contribs_a = [(1, SA[0.25, -0.5, 1.0]), (3, SA[2.0, 0.0, -1.5]), (1, SA[0.125, 0.25, 0.0])]
        contribs_b = [(2, SA[0.5, 0.5, 0.5]), (3, SA[-1.0, 1.0, 0.0]), (4, SA[0.75, -0.25, 0.125])]

        for (i, f) in contribs_a
            add_bead_force!(a, i, f)
            add_bead_force!(ref, i, f)
        end
        for (i, f) in contribs_b
            add_bead_force!(b, i, f)
            add_bead_force!(ref, i, f)
        end

        combine_force_energy!(a, b)

        # Underlying Int64 storage matches bit-for-bit.
        @test a.forces == ref.forces
    end

    @testset "combine_force_energy! bit-exact for matching FixedPoint" begin
        # Combining two ForceEnergyFixedPoint buffers with matching F,E should
        # produce results identical (bit-for-bit) to accumulating every
        # contribution into a single buffer from the start.
        nbeads = 4
        F, E = 30, 30
        a = ForceEnergyFixedPoint{F, E}(nbeads)
        b = ForceEnergyFixedPoint{F, E}(nbeads)
        ref = ForceEnergyFixedPoint{F, E}(nbeads)

        contribs_a = [(1, SA[0.25, -0.5, 1.0]), (3, SA[2.0, 0.0, -1.5]), (1, SA[0.125, 0.25, 0.0])]
        contribs_b = [(2, SA[0.5, 0.5, 0.5]), (3, SA[-1.0, 1.0, 0.0]), (4, SA[0.75, -0.25, 0.125])]

        for (i, f) in contribs_a
            add_bead_force!(a, i, f)
            add_bead_force!(ref, i, f)
        end
        add_energy!(a, 1.5)
        add_energy!(ref, 1.5)

        for (i, f) in contribs_b
            add_bead_force!(b, i, f)
            add_bead_force!(ref, i, f)
        end
        add_energy!(b, -0.75)
        add_energy!(ref, -0.75)

        combine_force_energy!(a, b)

        # Underlying Int64 storage matches bit-for-bit.
        @test a.forces == ref.forces
        @test a.energy === ref.energy
    end

    @testset "combine_force_energy! generic path on FixedPoint" begin
        # Using the generic combine_force_energy! (e.g. via differing F,E or a
        # mixed-type combine) goes through float round-tripping. Sanity-check
        # that the result still agrees on the values we use here, but is *not*
        # guaranteed bit-exact.
        nbeads = 2
        a = ForceEnergyFixedPoint{30, 30}(nbeads)
        b = ForceEnergyFixedPoint{20, 30}(nbeads)
        add_bead_force!(a, 1, SA[0.25, 0.5, -1.0])
        add_bead_force!(b, 1, SA[0.125, -0.25, 2.0])
        combine_force_energy!(a, b)
        @test get_bead_force(a, 1) ≈ SA[0.375, 0.25, 1.0]
    end

    @testset "DebugForceEnergy" begin
        nbeads = 3
        fe = DebugForceEnergy{Float64}(nbeads)

        for i in 1:nbeads
            @test get_bead_force(fe, i) === SA[0.0, 0.0, 0.0]
        end

        add_bead_force!(fe, 2, SA[1.0, 2.0, 3.0])
        add_bead_force!(fe, 2, SA[0.5, -1.0, 0.0])
        add_bead_force!(fe, 3, SA[-1.0, 0.0, 0.0])

        @test get_bead_force(fe, 1) === SA[0.0, 0.0, 0.0]
        @test get_bead_force(fe, 2) === SA[1.5, 1.0, 3.0]
        @test get_bead_force(fe, 3) === SA[-1.0, 0.0, 0.0]

        # Matches get_force! output.
        out = zeros(SVector{3, Float64}, nbeads)
        get_force!(fe, out)
        for i in 1:nbeads
            @test get_bead_force(fe, i) === out[i]
        end

        # Energy sums the recorded contributions, zero_energy! clears them.
        add_energy!(fe, 1.5)
        add_energy!(fe, -0.25)
        @test get_energy(fe) === 1.25
        zero_energy!(fe)
        @test get_energy(fe) === 0.0

        # zero_bead_force! wipes only that bead.
        zero_bead_force!(fe, 2)
        @test get_bead_force(fe, 2) === SA[0.0, 0.0, 0.0]
        @test get_bead_force(fe, 3) === SA[-1.0, 0.0, 0.0]

        # Out-of-range index is rejected.
        @test_throws Exception get_bead_force(fe, 0)
        @test_throws Exception get_bead_force(fe, nbeads + 1)
        @test_throws Exception add_bead_force!(fe, nbeads + 1, SA[1.0, 0.0, 0.0])
        @test_throws Exception zero_bead_force!(fe, nbeads + 1)

        add_energy!(fe, 1.0)
        zero_force_energy!(fe)
        @test get_energy(fe) === 0.0
        for i in 1:nbeads
            @test get_bead_force(fe, i) === SA[0.0, 0.0, 0.0]
        end
    end

    @testset "NullForceEnergy" begin
        nbeads = 4
        fe = NullForceEnergy(nbeads)

        @test get_nbeads(fe) == nbeads

        # All writes are no-ops.
        add_bead_force!(fe, 2, SA[1.0, 2.0, 3.0])
        add_energy!(fe, 1.5)
        zero_bead_force!(fe, 2)
        zero_energy!(fe)
        zero_force_energy!(fe)

        # Nothing is accumulated, so reading results errors.
        @test_throws MethodError get_energy(fe)
        @test_throws MethodError get_bead_force(fe, 1)
        @test_throws MethodError get_force!(fe, zeros(SVector{3, Float64}, nbeads))
        @test_throws MethodError get_force(fe)
    end
end
nothing
