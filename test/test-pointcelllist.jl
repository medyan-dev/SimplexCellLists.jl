using SimplexCellLists
using LinearAlgebra
using StaticArrays
using Test

# Helper to collect all entries from a PointCellList as (ix,iy,iz) => Vector{CellPointEntry}
function collect_entries(pcl::PointCellList{T,F}) where {T,F}
    result = Dict{NTuple{3,Int}, Vector{CellPointEntry{T,F}}}()
    for iz in 0:pcl.size[3]-1, iy in 0:pcl.size[2]-1, ix in 0:pcl.size[1]-1
        # li = SA[ix,iy,iz] ⋅ pcl.strides
        cs = reshape(pcl.cells, pcl.size[1], pcl.size[2], pcl.size[3])[ix+1,iy+1,iz+1]
        if cs.len > 0
            result[(ix,iy,iz)] = collect(view(cs.list, 1:cs.len))
        end
    end
    result
end

"""
Brute-force O(n) version of `map_nearby_points` for testing.
"""
function map_nearby_points_naive(
        f,
        positions::AbstractVector{SVector{3,Float32}},
        point_ids::AbstractVector{Int64},
        pos::SVector{3,Float32},
        cutoff::Float32,
        out,
    )
    cutoff ≤ 0 && return out
    cutoff2 = cutoff * cutoff
    for i in eachindex(positions, point_ids)
        diff = pos - positions[i]
        dist2 = diff ⋅ diff
        if dist2 ≤ cutoff2
            out, cont = f(CellPointEntry{Int64,Float32}(positions[i], point_ids[i]), out)
            cont || return out
        end
    end
    out
end

total_entries(pcl::PointCellList) = sum(pcl.cells[i].len for i in 1:prod(pcl.size))

@testset "pointcelllist" begin
    # Grid (4,4,4) with spacing 2.0 is centered at origin:
    # x: [-4, 4), y: [-4, 4), z: [-4, 4)
    # Cell (0,0,0) covers [-4,-2) × [-4,-2) × [-4,-2)
    # Cell (2,2,2) covers [0,2) × [0,2) × [0,2)
    # Cell (3,3,3) covers [2,4) × [2,4) × [2,4)

    @testset "invalid and degenerate sizes" begin
        @test_throws ArgumentError PointCellList{Int64,Float32}((-1,4,4), 2.0f0)
        # A zero-size axis makes an empty grid: queries find nothing.
        pcl = PointCellList{Int64,Float32}((0,4,4), 2.0f0)
        count = map_nearby_points(pcl, SA[0.5f0, 0.5f0, 0.5f0], 1.0f0, 0) do entry, out
            out + 1, true
        end
        @test count == 0
    end

    @testset "single point" begin
        pcl = PointCellList{Int64, Float32}((4,4,4), 2.0f0)
        cell_point_add!(pcl, SA[0.5f0, 0.5f0, 0.5f0], Int64(1))
        @test total_entries(pcl) == 1
        entries = collect_entries(pcl)
        @test haskey(entries, (2,2,2))
        e = entries[(2,2,2)][1]
        @test e.point_id == 1
        @test e.pos == SA[0.5f0, 0.5f0, 0.5f0]
    end

    @testset "point in corner cell" begin
        pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_point_add!(pcl, SA[-3.5f0, -3.5f0, -3.5f0], Int64(2))
        @test total_entries(pcl) == 1
        entries = collect_entries(pcl)
        @test haskey(entries, (0,0,0))
        e = entries[(0,0,0)][1]
        @test e.point_id == 2
    end

    @testset "point in last cell" begin
        pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_point_add!(pcl, SA[3.5f0, 3.5f0, 3.5f0], Int64(3))
        @test total_entries(pcl) == 1
        entries = collect_entries(pcl)
        @test haskey(entries, (3,3,3))
    end

    @testset "multiple points in same cell" begin
        pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_point_add!(pcl, SA[0.5f0, 0.5f0, 0.5f0], Int64(1))
        cell_point_add!(pcl, SA[1.5f0, 1.5f0, 1.5f0], Int64(2))
        @test total_entries(pcl) == 2
        entries = collect_entries(pcl)
        @test issetequal(entries[(2,2,2)], [
            CellPointEntry{Int64,Float32}(SA[0.5f0, 0.5f0, 0.5f0], Int64(1)),
            CellPointEntry{Int64,Float32}(SA[1.5f0, 1.5f0, 1.5f0], Int64(2)),
        ])
    end

    @testset "points in different cells" begin
        pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_point_add!(pcl, SA[-3.5f0, -3.5f0, -3.5f0], Int64(1))
        cell_point_add!(pcl, SA[3.5f0, 3.5f0, 3.5f0], Int64(2))
        @test total_entries(pcl) == 2
        entries = collect_entries(pcl)
        @test haskey(entries, (0,0,0))
        @test haskey(entries, (3,3,3))
        @test entries[(0,0,0)][1].point_id == 1
        @test entries[(3,3,3)][1].point_id == 2
    end

    @testset "point on cell boundary" begin
        pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
        # x=0 is boundary between cells 1 and 2; floor(0/2 + 2) = floor(2.0) = 2
        cell_point_add!(pcl, SA[0.0f0, -3.5f0, -3.5f0], Int64(1))
        @test total_entries(pcl) == 1
        entries = collect_entries(pcl)
        # Cells are half-open [lo, hi), so a point exactly on a boundary
        # deterministically belongs to the upper cell.
        @test haskey(entries, (2,0,0))
    end

    @testset "non-cubic grid" begin
        # (2,4,6) grid with spacing 1.0
        # x: [-1,1), y: [-2,2), z: [-3,3)
        pcl = PointCellList{Int64,Float32}((2,4,6), 1.0f0)
        cell_point_add!(pcl, SA[-0.5f0, -1.5f0, -2.5f0], Int64(1))  # cell (0,0,0)
        cell_point_add!(pcl, SA[0.5f0, 1.5f0, 2.5f0], Int64(2))     # cell (1,3,5)
        entries = collect_entries(pcl)
        @test haskey(entries, (0,0,0))
        @test haskey(entries, (1,3,5))
        @test total_entries(pcl) == 2
    end

    @testset "many points in one cell" begin
        pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
        n = 100
        for i in 1:n
            cell_point_add!(pcl, SA[0.5f0, 0.5f0, 0.5f0], Int64(i))
        end
        entries = collect_entries(pcl)
        es = entries[(2,2,2)]
        @test length(es) == n
        @test issetequal(map(e->e.point_id, es), 1:n)
        @test total_entries(pcl) == n
    end

    @testset "multiple points same position different IDs" begin
        pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
        pos = SA[0.5f0, 0.5f0, 0.5f0]
        cell_point_add!(pcl, pos, Int64(10))
        cell_point_add!(pcl, pos, Int64(20))
        count = map_nearby_points(pcl, pos, 0.1f0, 0) do entry, out
            out + 1, true
        end
        @test count == 2
    end

    @testset "cell_point_clear!" begin
        pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_point_add!(pcl, SA[0.5f0, 0.5f0, 0.5f0], Int64(1))
        cell_point_add!(pcl, SA[-3.5f0, -3.5f0, -3.5f0], Int64(2))
        @test total_entries(pcl) == 2
        cell_point_clear!(pcl)
        @test total_entries(pcl) == 0
    end

    @testset "points outside the grid are clamped to boundary cells" begin
        pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
        # Grid covers [-4,4)^3
        cell_point_add!(pcl, SA[100.0f0, 0.5f0, 0.5f0], Int64(1))   # beyond +x
        cell_point_add!(pcl, SA[-100.0f0, 0.5f0, 0.5f0], Int64(2))  # beyond -x
        cell_point_add!(pcl, SA[1.0f10, -1.0f10, 3.0f0], Int64(3))  # far outside in x and y
        cell_point_add!(pcl, SA[4.0f0, 0.5f0, 0.5f0], Int64(4))     # exactly on the upper bound
        cell_point_add!(pcl, SA[-4.0f0, 0.5f0, 0.5f0], Int64(5))    # exactly on the lower bound
        @test total_entries(pcl) == 5
        entries = collect_entries(pcl)
        @test map(e->e.point_id, entries[(3,2,2)]) == [1, 4]
        @test map(e->e.point_id, entries[(0,2,2)]) == [2, 5]
        @test entries[(3,0,3)][1].point_id == 3
    end

    @testset "map_nearby_points" begin
        @testset "empty" begin
            pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
            count = map_nearby_points(pcl, SA[0.5f0, 0.5f0, 0.5f0], 0.5f0, 0) do entry, out
                out + 1, true
            end
            @test count == 0
        end

        @testset "single point found" begin
            pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_point_add!(pcl, SA[0.5f0, 0.5f0, 0.5f0], Int64(1))
            count = map_nearby_points(pcl, SA[0.5f0, 0.5f0, 0.5f0], 0.5f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
        end

        @testset "no matches - outside cutoff" begin
            pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_point_add!(pcl, SA[0.5f0, 0.5f0, 0.5f0], Int64(1))
            count = map_nearby_points(pcl, SA[-3.0f0, -3.0f0, -3.0f0], 0.5f0, 0) do entry, out
                out + 1, true
            end
            @test count == 0
        end

        @testset "early exit" begin
            pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_point_add!(pcl, SA[0.5f0, 0.5f0, 0.5f0], Int64(1))
            cell_point_add!(pcl, SA[0.6f0, 0.5f0, 0.5f0], Int64(2))
            count = map_nearby_points(pcl, SA[0.5f0, 0.5f0, 0.5f0], 1.0f0, 0) do entry, out
                out + 1, false
            end
            @test count == 1
        end

        @testset "point across cell boundary within cutoff" begin
            pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
            # Point just inside cell (1,2,2), query just inside cell (2,2,2)
            cell_point_add!(pcl, SA[-0.1f0, 0.5f0, 0.5f0], Int64(1))
            count = map_nearby_points(pcl, SA[0.1f0, 0.5f0, 0.5f0], 0.5f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
        end

        @testset "query from outside grid finds nearby point" begin
            pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
            # Grid covers [-4,4). Point near edge.
            cell_point_add!(pcl, SA[3.5f0, 3.5f0, 3.5f0], Int64(1))
            # Query from outside grid, distance = 1.0
            count = map_nearby_points(pcl, SA[4.5f0, 3.5f0, 3.5f0], 2.0f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
        end

        @testset "query from outside grid finds points also outside grid" begin
            pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
            # Grid covers [-4,4). Point and query are both beyond the +x face.
            cell_point_add!(pcl, SA[10.0f0, 0.5f0, 0.5f0], Int64(1))
            count = map_nearby_points(pcl, SA[11.0f0, 0.5f0, 0.5f0], 2.0f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
            # The same query must not return an outside point beyond the cutoff.
            cell_point_add!(pcl, SA[11.0f0, 0.5f0, 30.0f0], Int64(2))
            count = map_nearby_points(pcl, SA[11.0f0, 0.5f0, 0.5f0], 2.0f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
        end

        @testset "query near grid edge with cutoff extending outside" begin
            pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
            # Points near the low-x boundary of the grid
            cell_point_add!(pcl, SA[-3.5f0, 0.5f0, 0.5f0], Int64(1))
            # Query with large cutoff that extends past grid boundary on -x side
            count = map_nearby_points(pcl, SA[-2.5f0, 0.5f0, 0.5f0], 3.0f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
        end

        @testset "zero cutoff returns out" begin
            pcl = PointCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_point_add!(pcl, SA[0.5f0, 0.5f0, 0.5f0], Int64(1))
            result = map_nearby_points(pcl, SA[0.5f0, 0.5f0, 0.5f0], 0.0f0, 42) do entry, out
                out + 1, true
            end
            @test result == 42
        end

        @testset "brute force comparison" begin
            for trial in 1:10000
                pcl = PointCellList{Int64,Float32}((10,10,10), 2.0f0)
                # Grid covers [-10,10)^3, but points and queries are sampled in
                # [-12,12)^3 so some fall outside the grid and get clamped.
                n_pts = 500
                positions = [SVector{3,Float32}(rand(Float32, 3) .* 24f0 .- 12f0) for _ in 1:n_pts]
                for i in 1:n_pts
                    cell_point_add!(pcl, positions[i], Int64(i))
                end

                function collect_id(entry, out)
                    push!(out, entry.point_id)
                    out, true
                end
                for _ in 1:20
                    cutoff = Float32(rand()) * 8f0
                    cutoff2 = cutoff * cutoff
                    q = SVector{3,Float32}(rand(Float32, 3) .* 24f0 .- 12f0)
                    cell_ids = map_nearby_points(collect_id, pcl, q, cutoff, Set{Int64}())
                    brute_ids = map_nearby_points_naive(collect_id, positions, Int64(1):n_pts, q, cutoff, Set{Int64}())
                    # Any id in one set but not the other must be right on the cutoff boundary.
                    for id in symdiff(cell_ids, brute_ids)
                        diff = q - positions[id]
                        d2 = diff ⋅ diff
                        @test d2 ≈ cutoff2 rtol=1e-6
                    end
                end
            end
        end
    end
end
nothing
