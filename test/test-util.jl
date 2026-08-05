using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Test

@testset "ball_intersects_unit_cube" begin
    buc = SimplexCellLists.ball_intersects_unit_cube
    @testset "center at origin" begin
        # Any positive radius intersects
        @test buc(SA[0.0, 0.0, 0.0], 0.0)
        @test buc(SA[0.0, 0.0, 0.0], 0.1)
        @test buc(SA[0.0, 0.0, 0.0], 10.0)
    end
    @testset "center inside cube" begin
        @test buc(SA[0.4, 0.4, 0.4], 0.0)
        @test buc(SA[-0.3, 0.2, -0.1], 0.0)
        @test buc(SA[0.5, 0.5, 0.5], 0.0)  # on corner, still touching
        @test buc(SA[0.5, 0.0, 0.0], 0.0)  # on face, still touching
    end
    @testset "face separation" begin
        # Ball center along +x axis, only x exceeds cube
        @test buc(SA[0.7, 0.0, 0.0], 0.21)   # distance = 0.2, radius > distance
        @test !buc(SA[0.7, 0.0, 0.0], 0.19)  # distance = 0.2, radius < distance
        # Along -y axis
        @test buc(SA[0.0, -1.0, 0.0], 0.51)
        @test !buc(SA[0.0, -1.0, 0.0], 0.49)
        # Along +z axis
        @test buc(SA[0.0, 0.0, 0.8], 0.31)
        @test !buc(SA[0.0, 0.0, 0.8], 0.29)
    end
    @testset "edge separation" begin
        # Center near an edge: two coordinates exceed cube
        # e.g. (+x, +y) edge, center at (1, 1, 0)
        # nearest point on cube is (0.5, 0.5, 0), dist = sqrt(0.25 + 0.25) = sqrt(0.5) ≈ 0.7071
        @test buc(SA[1.0, 1.0, 0.0], 0.71)
        @test !buc(SA[1.0, 1.0, 0.0], 0.70)
        # (-x, +z) edge
        @test buc(SA[-1.0, 0.0, 1.0], 0.71)
        @test !buc(SA[-1.0, 0.0, 1.0], 0.70)
    end
    @testset "vertex separation" begin
        # Center near vertex: all three coordinates exceed cube
        # center at (1, 1, 1), nearest vertex is (0.5, 0.5, 0.5), dist = sqrt(0.75) ≈ 0.8660
        @test buc(SA[1.0, 1.0, 1.0], 0.87)
        @test !buc(SA[1.0, 1.0, 1.0], 0.86)
        # Negative vertex
        @test buc(SA[-1.0, -1.0, -1.0], 0.87)
        @test !buc(SA[-1.0, -1.0, -1.0], 0.86)
    end
    @testset "exact boundary" begin
        # Exactly touching: d == r, should intersect (<=)
        @test buc(SA[1.0, 0.0, 0.0], 0.5)          # face touch
        @test buc(SA[1.0, 1.0, 0.0], sqrt(0.5))     # edge touch
        # vertex touch: sqrt(0.75)^2 != 0.75 due to fp rounding,
        # so test just inside and just outside instead
        @test buc(SA[1.0, 1.0, 1.0], sqrt(0.75) + eps())
        @test !buc(SA[1.0, 1.0, 1.0], sqrt(0.75) - 1e-10)
    end
    @testset "large radius" begin
        @test buc(SA[100.0, 100.0, 100.0], 200.0)
        @test !buc(SA[100.0, 100.0, 100.0], 170.0)
    end
    @testset "symmetry" begin
        # All 8 octant corners should give same result
        for sx in (-1, 1), sy in (-1, 1), sz in (-1, 1)
            c = SA[sx * 0.8, sy * 0.9, sz * 0.7]
            @test buc(c, 0.5) == buc(abs.(c), 0.5)
        end
    end
    @testset "Float32" begin
        @test buc(SA[0.7f0, 0.0f0, 0.0f0], 0.21f0)
        @test !buc(SA[0.7f0, 0.0f0, 0.0f0], 0.19f0)
    end
end

@testset "ball_intersects_biunit_cube" begin
    bbc = SimplexCellLists.ball_intersects_biunit_cube
    @testset "center at origin" begin
        # Any positive radius intersects
        @test bbc(SA[0.0, 0.0, 0.0], 0.0)
        @test bbc(SA[0.0, 0.0, 0.0], 0.1)
        @test bbc(SA[0.0, 0.0, 0.0], 10.0)
    end
    @testset "center inside cube" begin
        @test bbc(SA[0.8, 0.8, 0.8], 0.0)
        @test bbc(SA[-0.6, 0.4, -0.2], 0.0)
        @test bbc(SA[1.0, 1.0, 1.0], 0.0)  # on corner, still touching
        @test bbc(SA[1.0, 0.0, 0.0], 0.0)  # on face, still touching
    end
    @testset "face separation" begin
        # Ball center along +x axis, only x exceeds cube
        @test bbc(SA[1.4, 0.0, 0.0], 0.41)   # distance = 0.4, radius > distance
        @test !bbc(SA[1.4, 0.0, 0.0], 0.39)  # distance = 0.4, radius < distance
        # Along -y axis
        @test bbc(SA[0.0, -2.0, 0.0], 1.01)
        @test !bbc(SA[0.0, -2.0, 0.0], 0.99)
        # Along +z axis
        @test bbc(SA[0.0, 0.0, 1.6], 0.61)
        @test !bbc(SA[0.0, 0.0, 1.6], 0.59)
    end
    @testset "edge separation" begin
        # Center near an edge: two coordinates exceed cube
        # e.g. (+x, +y) edge, center at (2, 2, 0)
        # nearest point on cube is (1, 1, 0), dist = sqrt(1 + 1) = sqrt(2) ≈ 1.4142
        @test bbc(SA[2.0, 2.0, 0.0], 1.42)
        @test !bbc(SA[2.0, 2.0, 0.0], 1.41)
        # (-x, +z) edge
        @test bbc(SA[-2.0, 0.0, 2.0], 1.42)
        @test !bbc(SA[-2.0, 0.0, 2.0], 1.41)
    end
    @testset "vertex separation" begin
        # Center near vertex: all three coordinates exceed cube
        # center at (2, 2, 2), nearest vertex is (1, 1, 1), dist = sqrt(3) ≈ 1.7321
        @test bbc(SA[2.0, 2.0, 2.0], 1.74)
        @test !bbc(SA[2.0, 2.0, 2.0], 1.73)
        # Negative vertex
        @test bbc(SA[-2.0, -2.0, -2.0], 1.74)
        @test !bbc(SA[-2.0, -2.0, -2.0], 1.73)
    end
    @testset "exact boundary" begin
        # Exactly touching: d == r, should intersect (<=)
        @test bbc(SA[2.0, 0.0, 0.0], 1.0)          # face touch
        @test bbc(SA[2.0, 2.0, 0.0], sqrt(2.0))    # edge touch
        # vertex touch: sqrt(3)^2 != 3 due to fp rounding,
        # so test just inside and just outside instead
        @test bbc(SA[2.0, 2.0, 2.0], sqrt(3.0) + 2*eps())
        @test !bbc(SA[2.0, 2.0, 2.0], sqrt(3.0) - 1e-10)
    end
    @testset "large radius" begin
        @test bbc(SA[100.0, 100.0, 100.0], 200.0)
        @test !bbc(SA[100.0, 100.0, 100.0], 170.0)
    end
    @testset "symmetry" begin
        # All 8 octant corners should give same result
        for sx in (-1, 1), sy in (-1, 1), sz in (-1, 1)
            c = SA[sx * 1.6, sy * 1.8, sz * 1.4]
            @test bbc(c, 1.0) == bbc(abs.(c), 1.0)
        end
    end
    @testset "Float32" begin
        @test bbc(SA[1.4f0, 0.0f0, 0.0f0], 0.41f0)
        @test !bbc(SA[1.4f0, 0.0f0, 0.0f0], 0.39f0)
    end
    @testset "agrees with unit cube version at half scale" begin
        for trial in 1:1000
            c = SVector{3,Float64}(randn(3)) * 2.0
            r = abs(randn()) * 1.5
            @test bbc(c, r) == SimplexCellLists.ball_intersects_unit_cube(c/2, r/2)
        end
    end
end
