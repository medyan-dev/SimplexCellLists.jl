using SimplexCellLists
using StaticArrays
using Random
using Test

total_tri_entries(tcl::TriangleCellList) = sum(c.len for c in tcl.cells; init=0)

# All (cell coordinates, entry) pairs for one triangle id.
function tri_entries_for_id(tcl::TriangleCellList{T,F}, id::T) where {T,F}
    sz = tcl.size
    cells3 = reshape(tcl.cells, sz[1], sz[2], sz[3])
    es = Tuple{NTuple{3,Int}, CellTriangleEntry{T,F}}[]
    for iz in 0:sz[3]-1, iy in 0:sz[2]-1, ix in 0:sz[1]-1
        cs = cells3[ix+1,iy+1,iz+1]
        for j in 1:cs.len
            cs.list[j].id == id && push!(es, ((ix,iy,iz), cs.list[j]))
        end
    end
    es
end

# Check the entries of one triangle fill a box of cells from `lo` to `hi`,
# with `box_ends` marking the first and last cells along each axis.
function check_tri_box(es, lo, hi)
    @test sort(first.(es)) == sort(vec([(ix,iy,iz) for ix in lo[1]:hi[1], iy in lo[2]:hi[2], iz in lo[3]:hi[3]]))
    for (c, e) in es
        for a in 1:3
            @test (e.box_ends >> (a-1)) & 1 == (c[a] == lo[a])
            @test (e.box_ends >> (a+2)) & 1 == (c[a] == hi[a])
        end
        @test e.box_ends & 0xc0 == 0
    end
    @test allequal(last(x).tri for x in es)
end

function count_nearby_tris(tcl, q, cutoff)
    map_nearby_triangles(tcl, q, cutoff, 0) do entry, sep, out
        out + 1, true
    end
end

function nearby_tri_ids(tcl::TriangleCellList{T}, q, cutoff) where T
    map_nearby_triangles(tcl, q, cutoff, T[]) do entry, sep, out
        push!(out, entry.id)
        out, true
    end
end

@testset "trianglecelllist" begin
    # Grid (4,4,4) with spacing 2.0 is centered at origin:
    # x: [-4, 4), y: [-4, 4), z: [-4, 4)
    # Cell (0,0,0) covers [-4,-2) × [-4,-2) × [-4,-2)
    # Cell (2,2,2) covers [0,2) × [0,2) × [0,2)
    # Cell (3,3,3) covers [2,4) × [2,4) × [2,4)

    @testset "invalid and degenerate sizes" begin
        @test_throws ArgumentError TriangleCellList{Int64,Float32}((-1,4,4), 2.0f0)
        @test_throws ArgumentError TriangleCellList{Int64,Float32}((4,4,4), 0.0f0)
        # A zero-size axis makes an empty grid: queries find nothing, adds throw.
        tcl = TriangleCellList{Int64,Float32}((0,4,4), 2.0f0)
        @test count_nearby_tris(tcl, SA[0.5f0, 0.5f0, 0.5f0], 1.0f0) == 0
        @test_throws ArgumentError cell_triangle_add!(tcl,
            SA[0.5f0, 0.5f0, 0.5f0], SA[1.0f0, 0.5f0, 0.5f0], SA[0.5f0, 1.0f0, 0.5f0], Int64(1))
    end

    @testset "single triangle inside one cell" begin
        tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
        p0, p1, p2 = SA[0.3f0, 0.5f0, 0.5f0], SA[1.2f0, 0.5f0, 0.5f0], SA[0.5f0, 1.5f0, 0.7f0]
        cell_triangle_add!(tcl, p0, p1, p2, Int64(1))
        @test total_tri_entries(tcl) == 1
        es = tri_entries_for_id(tcl, Int64(1))
        @test first(only(es)) == (2,2,2)
        e = last(only(es))
        @test e.id == 1
        @test e.tri == SA[p0, p1, p2]
        @test e.box_ends == 0x3f
    end

    @testset "triangle spanning cells" begin
        tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
        # bounding box [-3, 3] × [-1, 0.5] × [0.5, 0.5]
        cell_triangle_add!(tcl, SA[-3.0f0, -1.0f0, 0.5f0], SA[3.0f0, -1.0f0, 0.5f0], SA[0.0f0, 0.5f0, 0.5f0], Int64(1))
        es = tri_entries_for_id(tcl, Int64(1))
        @test length(es) == 4*2*1
        check_tri_box(es, (0,1,2), (3,2,2))
    end

    @testset "large diagonal triangle fills its bounding box" begin
        tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_triangle_add!(tcl, SA[-3.9f0, -3.9f0, -3.9f0], SA[3.9f0, -3.9f0, 3.9f0], SA[3.9f0, 3.9f0, 3.9f0], Int64(1))
        es = tri_entries_for_id(tcl, Int64(1))
        @test length(es) == 64
        check_tri_box(es, (0,0,0), (3,3,3))
    end

    @testset "triangle on a cell boundary" begin
        tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
        # x = 0 is the boundary between x-cells 1 and 2; half-open cells put it in cell 2.
        cell_triangle_add!(tcl, SA[0.0f0, -3.0f0, 0.5f0], SA[-0.0f0, 3.0f0, 0.5f0], SA[0.0f0, 0.0f0, 1.5f0], Int64(1))
        es = tri_entries_for_id(tcl, Int64(1))
        check_tri_box(es, (2,0,2), (2,3,2))
    end

    @testset "degenerate triangles" begin
        tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
        pos = SA[0.5f0, 0.5f0, 0.5f0]
        # a single point
        cell_triangle_add!(tcl, pos, pos, pos, Int64(1))
        # collinear across cells
        cell_triangle_add!(tcl, SA[-1.0f0, 0.5f0, 0.5f0], SA[0.0f0, 0.5f0, 0.5f0], SA[1.0f0, 0.5f0, 0.5f0], Int64(2))
        @test length(tri_entries_for_id(tcl, Int64(1))) == 1
        check_tri_box(tri_entries_for_id(tcl, Int64(2)), (1,2,2), (2,2,2))
        @test sort(nearby_tri_ids(tcl, SA[0.5f0, 0.5f0, 1.0f0], 1.0f0)) == [1, 2]
        @test nearby_tri_ids(tcl, SA[-0.5f0, 0.5f0, 1.0f0], 0.6f0) == [2]
    end

    @testset "non-finite triangles error" begin
        tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
        ok = SA[0.5f0, 0.5f0, 0.5f0]
        for bad in (SA[NaN32, 0.5f0, 0.5f0], SA[0.5f0, Inf32, 0.5f0], SA[0.5f0, 0.5f0, -Inf32])
            @test_throws ArgumentError cell_triangle_add!(tcl, bad, ok, ok, Int64(1))
            @test_throws ArgumentError cell_triangle_add!(tcl, ok, bad, ok, Int64(1))
            @test_throws ArgumentError cell_triangle_add!(tcl, ok, ok, bad, Int64(1))
        end
        # nothing was added
        @test total_tri_entries(tcl) == 0
    end

    @testset "cell_triangles_clear!" begin
        tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_triangle_add!(tcl, SA[-3.0f0, 0.5f0, 0.5f0], SA[3.0f0, 0.5f0, 0.5f0], SA[0.0f0, 1.5f0, 0.5f0], Int64(1))
        @test total_tri_entries(tcl) == 4
        cell_triangles_clear!(tcl)
        @test total_tri_entries(tcl) == 0
        @test count_nearby_tris(tcl, SA[0.0f0, 0.5f0, 0.5f0], 1.0f0) == 0
    end

    @testset "triangles outside the grid are clamped to boundary cells" begin
        tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_triangle_add!(tcl, SA[10.0f0, 10.0f0, 10.0f0], SA[12.0f0, 10.0f0, 10.0f0], SA[10.0f0, 12.0f0, 10.0f0], Int64(1))
        check_tri_box(tri_entries_for_id(tcl, Int64(1)), (3,3,3), (3,3,3))
        cell_triangle_add!(tcl, SA[-100.0f0, 0.5f0, 0.5f0], SA[100.0f0, 0.5f0, 0.5f0], SA[0.0f0, 0.5f0, 1.5f0], Int64(2))
        check_tri_box(tri_entries_for_id(tcl, Int64(2)), (0,2,2), (3,2,2))
        cell_triangle_add!(tcl, SA[1.0f10, -1.0f10, 3.0f0], SA[1.0f10, -1.0f10, 3.5f0], SA[1.0f10, -1.0f10+1f4, 3.5f0], Int64(3))
        check_tri_box(tri_entries_for_id(tcl, Int64(3)), (3,0,3), (3,0,3))
    end

    @testset "map_nearby_triangles" begin
        @testset "empty" begin
            tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
            @test count_nearby_tris(tcl, SA[0.5f0, 0.5f0, 0.5f0], 0.5f0) == 0
        end

        @testset "single triangle found" begin
            tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_triangle_add!(tcl, SA[0.3f0, 0.5f0, 0.5f0], SA[1.2f0, 0.5f0, 0.5f0], SA[0.3f0, 1.5f0, 0.5f0], Int64(1))
            @test count_nearby_tris(tcl, SA[0.5f0, 0.7f0, 0.9f0], 0.5f0) == 1
            # outside cutoff
            @test count_nearby_tris(tcl, SA[0.5f0, 0.7f0, 0.9f0], 0.3f0) == 0
            @test count_nearby_tris(tcl, SA[-3.0f0, -3.0f0, -3.0f0], 0.5f0) == 0
        end

        @testset "early exit" begin
            tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_triangle_add!(tcl, SA[0.3f0, 0.5f0, 0.5f0], SA[1.2f0, 0.5f0, 0.5f0], SA[0.3f0, 1.5f0, 0.5f0], Int64(1))
            cell_triangle_add!(tcl, SA[0.3f0, 0.5f0, 0.6f0], SA[1.2f0, 0.5f0, 0.6f0], SA[0.3f0, 1.5f0, 0.6f0], Int64(2))
            count = map_nearby_triangles(tcl, SA[0.5f0, 0.7f0, 0.5f0], 1.0f0, 0) do entry, sep, out
                out + 1, false
            end
            @test count == 1
        end

        @testset "reported exactly once across cells" begin
            tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_triangle_add!(tcl, SA[-3.5f0, -3.5f0, -3.5f0], SA[3.5f0, -3.5f0, 3.5f0], SA[3.5f0, 3.5f0, 3.5f0], Int64(1))
            # From every cell, with cutoffs that do and don't cover the whole grid
            for ix in 0:3, iy in 0:3, iz in 0:3, cutoff in (0.5f0, 3.0f0, 20.0f0)
                q = SA[2ix - 3.0f0, 2iy - 3.0f0, 2iz - 3.0f0]
                d2 = dist_sqr(SA[q], tcl.cells[1].list[1].tri)
                @test nearby_tri_ids(tcl, q, cutoff) == (d2 ≤ cutoff^2 ? [1] : Int64[])
            end
        end

        @testset "triangle across cell boundary within cutoff" begin
            tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
            # Triangle just inside x-cell 1, query just inside x-cell 2
            cell_triangle_add!(tcl, SA[-0.1f0, 0.5f0, 0.5f0], SA[-0.1f0, 1.5f0, 0.5f0], SA[-0.1f0, 1.0f0, 1.5f0], Int64(1))
            @test count_nearby_tris(tcl, SA[0.1f0, 1.0f0, 0.7f0], 0.5f0) == 1
        end

        @testset "query near a corner of a large triangle's bounding box" begin
            tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
            # The triangle is in the plane x + y = 0, so the corners of its
            # bounding box at x = y = ±3.5 are far from it.
            cell_triangle_add!(tcl, SA[-3.5f0, 3.5f0, -3.5f0], SA[3.5f0, -3.5f0, -3.5f0], SA[0.0f0, 0.0f0, 3.5f0], Int64(1))
            @test count_nearby_tris(tcl, SA[3.0f0, 3.0f0, -3.0f0], 1.0f0) == 0
            @test count_nearby_tris(tcl, SA[3.0f0, 3.0f0, -3.0f0], 4.5f0) == 1
            @test count_nearby_tris(tcl, SA[5.0f0, 5.0f0, -5.0f0], 7.5f0) == 1
        end

        @testset "query from outside grid" begin
            tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_triangle_add!(tcl, SA[3.5f0, 3.0f0, 3.5f0], SA[3.5f0, 3.8f0, 3.5f0], SA[3.0f0, 3.4f0, 3.5f0], Int64(1))
            @test count_nearby_tris(tcl, SA[4.5f0, 3.5f0, 3.5f0], 2.0f0) == 1
            # Both beyond the +x face.
            cell_triangle_add!(tcl, SA[10.0f0, 0.5f0, 0.5f0], SA[12.0f0, 0.5f0, 0.5f0], SA[11.0f0, 1.5f0, 0.5f0], Int64(2))
            @test nearby_tri_ids(tcl, SA[11.0f0, 0.5f0, 1.5f0], 2.0f0) == [2]
            # Not one beyond the cutoff.
            cell_triangle_add!(tcl, SA[10.0f0, 0.5f0, 30.0f0], SA[12.0f0, 0.5f0, 30.0f0], SA[11.0f0, 1.5f0, 30.0f0], Int64(3))
            @test nearby_tri_ids(tcl, SA[11.0f0, 0.5f0, 1.5f0], 2.0f0) == [2]
        end

        @testset "zero and negative cutoff return out" begin
            tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_triangle_add!(tcl, SA[0.3f0, 0.5f0, 0.5f0], SA[1.2f0, 0.5f0, 0.5f0], SA[0.3f0, 1.5f0, 0.5f0], Int64(1))
            for cutoff in (0.0f0, -1.0f0)
                result = map_nearby_triangles(tcl, SA[0.5f0, 0.5f0, 0.5f0], cutoff, 42) do entry, sep, out
                    out + 1, true
                end
                @test result == 42
            end
        end

        @testset "sep argument" begin
            tcl = TriangleCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_triangle_add!(tcl, SA[-1.0f0, -1.0f0, 0.0f0], SA[1.0f0, -1.0f0, 0.0f0], SA[0.0f0, 1.0f0, 0.0f0], Int64(1))
            # closest point on the face
            seps = map_nearby_triangles(tcl, SA[0.0f0, 0.0f0, 0.75f0], 1.0f0, SVector{3,Float32}[]) do entry, sep, out
                push!(out, sep)
                out, true
            end
            @test seps ≈ [SA[0.0f0, 0.0f0, -0.75f0]]
            # closest point on an edge
            seps = map_nearby_triangles(tcl, SA[0.5f0, -1.5f0, 0.5f0], 1.0f0, SVector{3,Float32}[]) do entry, sep, out
                push!(out, sep)
                out, true
            end
            @test seps ≈ [SA[0.0f0, 0.5f0, -0.5f0]]
        end

        @testset "Float64" begin
            tcl = TriangleCellList{Int32,Float64}((4,4,4), 2.0)
            cell_triangle_add!(tcl, SA[-3.0, -3.0, 0.5], SA[3.0, -3.0, 0.5], SA[0.0, 3.0, 0.5], Int32(7))
            @test nearby_tri_ids(tcl, SA[0.0, 0.0, 1.0], 0.6) == Int32[7]
            @test nearby_tri_ids(tcl, SA[0.0, 0.0, 1.0], 0.4) == Int32[]
        end

        @testset "brute force comparison" begin
            rng = Xoshiro(1234)
            for trial in 1:300
                tcl = TriangleCellList{Int64,Float32}((10,10,10), 2.0f0)
                # Grid covers [-10,10)^3, but triangles and queries are sampled in
                # [-12,12)^3 so some fall outside the grid and get clamped.
                n_tris = 300
                tris = map(1:n_tris) do i
                    p0 = rand(rng, SVector{3,Float32}) .* 24f0 .- 12f0
                    if i % 20 == 0
                        SA[p0, p0, p0] # occasional single point
                    elseif i % 20 == 1
                        d = rand(rng, SVector{3,Float32}) .* 4f0 .- 2f0
                        SA[p0, p0 + d, p0 + 0.5f0*d] # occasional collinear
                    elseif i % 10 == 2
                        # occasional large triangle
                        SA[p0, rand(rng, SVector{3,Float32}) .* 24f0 .- 12f0, rand(rng, SVector{3,Float32}) .* 24f0 .- 12f0]
                    else
                        SA[p0, p0 + rand(rng, SVector{3,Float32}) .* 4f0 .- 2f0, p0 + rand(rng, SVector{3,Float32}) .* 4f0 .- 2f0]
                    end
                end
                for i in 1:n_tris
                    cell_triangle_add!(tcl, tris[i][1], tris[i][2], tris[i][3], Int64(i))
                end

                for _ in 1:20
                    cutoff = Float32(rand(rng)) * 8f0
                    cutoff2 = cutoff * cutoff
                    q = rand(rng, SVector{3,Float32}) .* 24f0 .- 12f0
                    seps = map_nearby_triangles(tcl, q, cutoff, Pair{Int64, SVector{3,Float32}}[]) do entry, sep, out
                        push!(out, entry.id => sep)
                        out, true
                    end
                    cell_ids = first.(seps)
                    # each triangle reported at most once
                    @test allunique(cell_ids)
                    # exactly the triangles within the cutoff by dist_sqr
                    brute_ids = [Int64(i) for i in 1:n_tris if dist_sqr(SA[q], tris[i]) ≤ cutoff2]
                    @test sort(cell_ids) == brute_ids
                    # with the separation of dist_sqr_r_s_t, bit for bit
                    @test all(seps) do (id, sep)
                        sep === dist_sqr_r_s_t(SA[q], tris[id])[2]
                    end
                end
            end
        end
    end

    @testset "show" begin
        tcl = TriangleCellList{Int,Float32}((4,4,4), 1f0)
        # Crosses x = -1, 0, and 1, so has entries in 4 cells.
        cell_triangle_add!(tcl, SA[-1.5f0, 0.5f0, 0.5f0], SA[1.5f0, 0.5f0, 0.5f0], SA[0.0f0, 0.7f0, 0.5f0], 1)
        # Inside one of those cells.
        cell_triangle_add!(tcl, SA[0.2f0, 0.5f0, 0.5f0], SA[0.3f0, 0.5f0, 0.5f0], SA[0.2f0, 0.6f0, 0.5f0], 2)
        @test repr(MIME"text/plain"(), tcl) == """
            TriangleCellList{Int64, Float32} with 2 triangles:
              size: 4×4×4 cells
              grid spacing: 1.0
              cell entries: 5
              most entries in one cell: 2"""
        cell_triangles_clear!(tcl)
        cell_triangle_add!(tcl, SA[0.2f0, 0.5f0, 0.5f0], SA[0.3f0, 0.5f0, 0.5f0], SA[0.2f0, 0.6f0, 0.5f0], 2)
        @test startswith(repr(MIME"text/plain"(), tcl), "TriangleCellList{Int64, Float32} with 1 triangle:\n")
    end
end
nothing
