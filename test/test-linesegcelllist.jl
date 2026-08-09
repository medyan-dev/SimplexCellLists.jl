using SimplexCellLists
using LinearAlgebra
using StaticArrays
using Test

# Helper to collect all entries from a LineSegCellList as (ix,iy,iz) => Vector{CellLineSegEntry}
function collect_seg_entries(scl::LineSegCellList{T,F}) where {T,F}
    result = Dict{NTuple{3,Int}, Vector{CellLineSegEntry{T,F}}}()
    sz = scl.size
    cells3 = reshape(scl.cells, sz[1], sz[2], sz[3])
    for iz in 0:sz[3]-1, iy in 0:sz[2]-1, ix in 0:sz[1]-1
        cs = cells3[ix+1,iy+1,iz+1]
        if cs.len > 0
            result[(ix,iy,iz)] = collect(view(cs.list, 1:cs.len))
        end
    end
    result
end

total_seg_entries(scl::LineSegCellList) = sum(scl.cells[i].len for i in 1:prod(scl.size); init=0)

# All entries for one line segment id, sorted by tmin.
function seg_entries_for_id(scl::LineSegCellList{T,F}, id::T) where {T,F}
    es = CellLineSegEntry{T,F}[]
    for c in scl.cells
        for j in 1:c.len
            c.list[j].id == id && push!(es, c.list[j])
        end
    end
    sort(es; by=e->e.tmin)
end

# Check that the entries of one line segment exactly partition [0, seg_len]:
# tmin of the first entry is 0, ranges chain with no gaps or overlaps, and
# only the last entry has is_end set, with tmax == seg_len.
function check_seg_partition(es, seg_len)
    @test !isempty(es)
    @test es[1].tmin == 0
    for k in 1:length(es)-1
        @test cell_line_seg_tmax(es[k]) == es[k+1].tmin
        @test !cell_line_seg_is_end(es[k])
    end
    @test cell_line_seg_is_end(es[end])
    # If more than one entry, the last entry cannot be empty
    if length(es) > 1
        @test cell_line_seg_tmax(es[end]) != es[end].tmin
    end
    # ≈ because recomputing seg_len here can differ from the package's
    # internal computation by an ulp (fma contraction).
    @test cell_line_seg_tmax(es[end]) ≈ seg_len rtol=4*eps(Float32)
    # all entries share the segment's geometry
    @test allequal(e.p0 for e in es)
    @test allequal(e.d_hat for e in es)
end

@testset "linesegcelllist" begin
    # Grid (4,4,4) with spacing 2.0 is centered at origin:
    # x: [-4, 4), y: [-4, 4), z: [-4, 4)
    # Cell (0,0,0) covers [-4,-2) × [-4,-2) × [-4,-2)
    # Cell (2,2,2) covers [0,2) × [0,2) × [0,2)
    # Cell (3,3,3) covers [2,4) × [2,4) × [2,4)

    @testset "invalid and degenerate sizes" begin
        @test_throws ArgumentError LineSegCellList{Int64,Float32}((-1,4,4), 2.0f0)
        # A zero-size axis makes an empty grid: queries find nothing, adds throw.
        scl = LineSegCellList{Int64,Float32}((0,4,4), 2.0f0)
        count = map_nearby_line_segs(scl, SA[0.5f0, 0.5f0, 0.5f0], 1.0f0, 0) do entry, out
            out + 1, true
        end
        @test count == 0
        @test_throws ArgumentError cell_line_seg_add!(scl, SA[0.5f0, 0.5f0, 0.5f0], SA[1.0f0, 0.5f0, 0.5f0], Int64(1))
    end

    @testset "single segment inside one cell" begin
        scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_line_seg_add!(scl, SA[0.3f0, 0.5f0, 0.5f0], SA[1.2f0, 0.5f0, 0.5f0], Int64(1))
        @test total_seg_entries(scl) == 1
        entries = collect_seg_entries(scl)
        @test haskey(entries, (2,2,2))
        e = entries[(2,2,2)][1]
        @test e.id == 1
        @test e.p0 == SA[0.3f0, 0.5f0, 0.5f0]
        @test e.d_hat == SA[1.0f0, 0.0f0, 0.0f0]
        @test e.tmin == 0
        @test cell_line_seg_tmax(e) == 1.2f0 - 0.3f0
        @test cell_line_seg_is_end(e)
    end

    @testset "axis-aligned segment spanning cells" begin
        scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
        # x from -3 to 3 crosses cell boundaries at x = -2, 0, 2
        cell_line_seg_add!(scl, SA[-3.0f0, 0.5f0, 0.5f0], SA[3.0f0, 0.5f0, 0.5f0], Int64(1))
        @test total_seg_entries(scl) == 4
        entries = collect_seg_entries(scl)
        for ix in 0:3
            @test haskey(entries, (ix,2,2))
        end
        es = seg_entries_for_id(scl, Int64(1))
        check_seg_partition(es, 6.0f0)
        @test [e.tmin for e in es] == [0.0f0, 1.0f0, 3.0f0, 5.0f0]
        @test cell_line_seg_tmax.(es) == [1.0f0, 3.0f0, 5.0f0, 6.0f0]
    end

    @testset "reversed direction segment" begin
        scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_line_seg_add!(scl, SA[3.0f0, 0.5f0, 0.5f0], SA[-3.0f0, 0.5f0, 0.5f0], Int64(1))
        @test total_seg_entries(scl) == 4
        es = seg_entries_for_id(scl, Int64(1))
        check_seg_partition(es, 6.0f0)
        @test es[1].d_hat == SA[-1.0f0, 0.0f0, 0.0f0]
    end

    @testset "diagonal segment" begin
        scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
        p0 = SA[-3.9f0, -3.9f0, -3.9f0]
        p1 = SA[3.9f0, 3.9f0, 3.9f0]
        cell_line_seg_add!(scl, p0, p1, Int64(1))
        d = p1 - p0
        seg_len = sqrt(d ⋅ d)
        es = seg_entries_for_id(scl, Int64(1))
        # crosses 3 boundaries per axis: 1 + 9 cells visited
        @test length(es) == 10
        check_seg_partition(es, seg_len)
        entries = collect_seg_entries(scl)
        @test haskey(entries, (0,0,0))
        @test haskey(entries, (3,3,3))
    end

    @testset "zero-length segment" begin
        scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
        pos = SA[0.5f0, 0.5f0, 0.5f0]
        cell_line_seg_add!(scl, pos, pos, Int64(1))
        @test total_seg_entries(scl) == 1
        e = collect_seg_entries(scl)[(2,2,2)][1]
        @test e.d_hat == zero(SVector{3,Float32})
        @test e.tmin == 0
        @test cell_line_seg_tmax(e) == 0
        @test cell_line_seg_is_end(e)
        count = map_nearby_line_segs(scl, SA[0.5f0, 0.5f0, 1.0f0], 1.0f0, 0) do entry, out
            out + 1, true
        end
        @test count == 1
    end

    @testset "non-finite segments error" begin
        scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
        ok = SA[0.5f0, 0.5f0, 0.5f0]
        nanp = SA[NaN32, 0.5f0, 0.5f0]
        infp = SA[Inf32, 0.5f0, 0.5f0]
        @test_throws ArgumentError cell_line_seg_add!(scl, nanp, ok, Int64(1))
        @test_throws ArgumentError cell_line_seg_add!(scl, ok, nanp, Int64(1))
        @test_throws ArgumentError cell_line_seg_add!(scl, infp, ok, Int64(1))
        @test_throws ArgumentError cell_line_seg_add!(scl, ok, infp, Int64(1))
        @test_throws ArgumentError cell_line_seg_add!(scl, infp, infp, Int64(1))
        # finite endpoints whose squared length overflows Float32
        big = SA[1.0f25, 0.5f0, 0.5f0]
        @test_throws ArgumentError cell_line_seg_add!(scl, -big, big, Int64(1))
        # nothing was added
        @test total_seg_entries(scl) == 0
    end

    @testset "segment along a cell boundary" begin
        scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
        # x = 0 is the boundary between x-cells 1 and 2; half-open cells put it in cell 2.
        # p1 x of -0.0 makes d_x == -0.0, exercising the 0/0 -> NaN -> Inf tMax fixup.
        cell_line_seg_add!(scl, SA[0.0f0, -3.0f0, 0.5f0], SA[-0.0f0, 3.0f0, 0.5f0], Int64(1))
        es = seg_entries_for_id(scl, Int64(1))
        @test length(es) == 4
        check_seg_partition(es, 6.0f0)
        entries = collect_seg_entries(scl)
        for iy in 0:3
            @test haskey(entries, (2,iy,2))
        end
    end

    @testset "multiple segments same position different IDs" begin
        scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
        p0 = SA[0.5f0, 0.5f0, 0.5f0]
        p1 = SA[1.5f0, 0.5f0, 0.5f0]
        cell_line_seg_add!(scl, p0, p1, Int64(10))
        cell_line_seg_add!(scl, p0, p1, Int64(20))
        count = map_nearby_line_segs(scl, p0, 0.1f0, 0) do entry, out
            out + 1, true
        end
        @test count == 2
    end

    @testset "cell_line_segs_clear!" begin
        scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
        cell_line_seg_add!(scl, SA[-3.0f0, 0.5f0, 0.5f0], SA[3.0f0, 0.5f0, 0.5f0], Int64(1))
        @test total_seg_entries(scl) == 4
        cell_line_segs_clear!(scl)
        @test total_seg_entries(scl) == 0
    end

    @testset "segments outside the grid are clamped to boundary cells" begin
        @testset "fully outside" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_line_seg_add!(scl, SA[10.0f0, 10.0f0, 10.0f0], SA[12.0f0, 10.0f0, 10.0f0], Int64(1))
            @test total_seg_entries(scl) == 1
            e = collect_seg_entries(scl)[(3,3,3)][1]
            @test e.tmin == 0
            @test cell_line_seg_tmax(e) == 2.0f0
            @test cell_line_seg_is_end(e)
        end
        @testset "crossing the whole grid" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_line_seg_add!(scl, SA[-100.0f0, 0.5f0, 0.5f0], SA[100.0f0, 0.5f0, 0.5f0], Int64(1))
            @test total_seg_entries(scl) == 4
            es = seg_entries_for_id(scl, Int64(1))
            check_seg_partition(es, 200.0f0)
            # boundary cells own the out-of-grid parts:
            # crossings at x = -2, 0, 2, i.e. t = 98, 100, 102
            @test [e.tmin for e in es] == [0.0f0, 98.0f0, 100.0f0, 102.0f0]
            @test cell_line_seg_tmax.(es) == [98.0f0, 100.0f0, 102.0f0, 200.0f0]
        end
        @testset "far outside in several axes" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_line_seg_add!(scl, SA[1.0f10, -1.0f10, 3.0f0], SA[1.0f10, -1.0f10, 3.5f0], Int64(1))
            @test total_seg_entries(scl) == 1
            @test haskey(collect_seg_entries(scl), (3,0,3))
        end
    end

    @testset "map_nearby_line_segs" begin
        @testset "empty" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            count = map_nearby_line_segs(scl, SA[0.5f0, 0.5f0, 0.5f0], 0.5f0, 0) do entry, out
                out + 1, true
            end
            @test count == 0
        end

        @testset "single segment found" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_line_seg_add!(scl, SA[0.3f0, 0.5f0, 0.5f0], SA[1.2f0, 0.5f0, 0.5f0], Int64(1))
            count = map_nearby_line_segs(scl, SA[0.7f0, 0.9f0, 0.5f0], 0.5f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
        end

        @testset "no matches - outside cutoff" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_line_seg_add!(scl, SA[0.3f0, 0.5f0, 0.5f0], SA[1.2f0, 0.5f0, 0.5f0], Int64(1))
            count = map_nearby_line_segs(scl, SA[-3.0f0, -3.0f0, -3.0f0], 0.5f0, 0) do entry, out
                out + 1, true
            end
            @test count == 0
        end

        @testset "early exit" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_line_seg_add!(scl, SA[0.3f0, 0.5f0, 0.5f0], SA[1.2f0, 0.5f0, 0.5f0], Int64(1))
            cell_line_seg_add!(scl, SA[0.3f0, 0.6f0, 0.5f0], SA[1.2f0, 0.6f0, 0.5f0], Int64(2))
            count = map_nearby_line_segs(scl, SA[0.5f0, 0.5f0, 0.5f0], 1.0f0, 0) do entry, out
                out + 1, false
            end
            @test count == 1
        end

        @testset "reported exactly once across cells" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_line_seg_add!(scl, SA[-3.5f0, -3.5f0, -3.5f0], SA[3.5f0, 3.5f0, 3.5f0], Int64(1))
            # cutoff big enough that every cell the segment touches is visited
            ids = map_nearby_line_segs(scl, SA[0.0f0, 0.0f0, 0.0f0], 20.0f0, Int64[]) do entry, out
                push!(out, entry.id)
                out, true
            end
            @test ids == [1]
        end

        @testset "segment across cell boundary within cutoff" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            # Segment just inside x-cell 1, query just inside x-cell 2
            cell_line_seg_add!(scl, SA[-0.1f0, 0.5f0, 0.5f0], SA[-0.1f0, 1.5f0, 0.5f0], Int64(1))
            count = map_nearby_line_segs(scl, SA[0.1f0, 1.0f0, 0.5f0], 0.5f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
        end

        @testset "query from outside grid finds nearby segment" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_line_seg_add!(scl, SA[3.5f0, 3.0f0, 3.5f0], SA[3.5f0, 3.8f0, 3.5f0], Int64(1))
            count = map_nearby_line_segs(scl, SA[4.5f0, 3.5f0, 3.5f0], 2.0f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
        end

        @testset "query from outside grid finds segments also outside grid" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            # Grid covers [-4,4). Segment and query are both beyond the +x face.
            cell_line_seg_add!(scl, SA[10.0f0, 0.5f0, 0.5f0], SA[12.0f0, 0.5f0, 0.5f0], Int64(1))
            count = map_nearby_line_segs(scl, SA[11.0f0, 0.5f0, 1.5f0], 2.0f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
            # The same query must not return an outside segment beyond the cutoff.
            cell_line_seg_add!(scl, SA[10.0f0, 0.5f0, 30.0f0], SA[12.0f0, 0.5f0, 30.0f0], Int64(2))
            count = map_nearby_line_segs(scl, SA[11.0f0, 0.5f0, 1.5f0], 2.0f0, 0) do entry, out
                out + 1, true
            end
            @test count == 1
        end

        @testset "zero and negative cutoff return out" begin
            scl = LineSegCellList{Int64,Float32}((4,4,4), 2.0f0)
            cell_line_seg_add!(scl, SA[0.3f0, 0.5f0, 0.5f0], SA[1.2f0, 0.5f0, 0.5f0], Int64(1))
            for cutoff in (0.0f0, -1.0f0)
                result = map_nearby_line_segs(scl, SA[0.5f0, 0.5f0, 0.5f0], cutoff, 42) do entry, out
                    out + 1, true
                end
                @test result == 42
            end
        end

        @testset "brute force comparison" begin
            for trial in 1:5000
                scl = LineSegCellList{Int64,Float32}((10,10,10), 2.0f0)
                # Grid covers [-10,10)^3, but segments and queries are sampled in
                # [-12,12)^3 so some fall outside the grid and get clamped.
                n_segs = 500
                segs = map(1:n_segs) do i
                    p0 = rand(SVector{3,Float32}) .* 24f0 .- 12f0
                    p1 = if i % 20 == 0
                        p0 # occasional zero-length segment
                    elseif i % 5 == 0
                        rand(SVector{3,Float32}) .* 24f0 .- 12f0 # long segment
                    else
                        p0 + rand(SVector{3,Float32}) .* 4f0 .- 2f0
                    end
                    SA[p0, p1]
                end
                for i in 1:n_segs
                    cell_line_seg_add!(scl, segs[i][1], segs[i][2], Int64(i))
                end

                for _ in 1:20
                    cutoff = Float32(rand()) * 8f0
                    cutoff2 = cutoff * cutoff
                    q = rand(SVector{3,Float32}) .* 24f0 .- 12f0
                    cell_ids = map_nearby_line_segs(scl, q, cutoff, Int64[]) do entry, out
                        push!(out, entry.id)
                        out, true
                    end
                    # each segment reported at most once
                    @test allunique(cell_ids)
                    brute_ids = [Int64(i) for i in 1:n_segs if dist_sqr(SA[q], segs[i]) ≤ cutoff2]
                    # Any id in one set but not the other must be right on the cutoff boundary.
                    for id in symdiff(cell_ids, brute_ids)
                        d2 = dist_sqr(SA[q], SA[segs[id][1], segs[id][2]])
                        @test d2 ≈ cutoff2 rtol=1e-4 atol=1f-4
                    end
                end
            end
        end
    end
end
nothing
