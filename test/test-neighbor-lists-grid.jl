using SimplexCellLists
using Test
using StaticArrays
using Random

@testset "neighborlists grid" begin
    policy = DefaultCollidePolicy()
    @testset "_grid_size" begin
        _grid_size = SimplexCellLists._grid_size

        # Range of allowed min_cell_width, the caller clamps it to this
        tiny = nextfloat(0.0f0)
        huge = floatmax(Float32)

        # Exact fits
        @test _grid_size(SA[10.0, 10.0, 10.0], tiny, 1000) == SA[10, 10, 10]
        @test _grid_size(SA[100f0, 100f0, 100f0], 4f0, 10^6) == SA[25, 25, 25]
        # Axes with equal cell widths merge together, leaving cells unused
        @test _grid_size(SA[10.0, 10.0, 10.0], tiny, 999) == SA[9, 9, 9]
        @test _grid_size(SA[10.0, 10.0, 10.0], 1.0, 7) == SA[1, 1, 1]
        # Narrowest cells are merged first
        @test _grid_size(SA[40.0, 20.0, 10.0], tiny, 8) == SA[4, 2, 1]
        @test _grid_size(SA[10.0, 20.0, 40.0], tiny, 8) == SA[1, 2, 4]
        # Near a merge on every axis, the bisection stops early with fewer cells
        @test _grid_size(SA[100.0, 1.9, 1.0], 0.01, 150) == SA[104, 1, 1]

        # Zero width axes get one cell
        @test _grid_size(SA[10.0, 0.0, 0.0], tiny, 100) == SA[100, 1, 1]
        @test _grid_size(SA[10.0, 10.0, 0.0], tiny, 100) == SA[10, 10, 1]

        # The minimum is a 1x1x1 grid
        @test _grid_size(SA[0.0, 0.0, 0.0], tiny, 1) == SA[1, 1, 1]
        @test _grid_size(SA[0.0, 0.0, 0.0], 1.0, 100) == SA[1, 1, 1]
        @test _grid_size(SA[10.0, 10.0, 10.0], tiny, 1) == SA[1, 1, 1]
        @test _grid_size(SA[10.0, 10.0, 10.0], huge, 1000) == SA[1, 1, 1]

        # Each axis is capped at 256 cells, even if more cells are allowed.
        # Other axes keep splitting, giving stretched cells.
        @test _grid_size(SA[1e4, 1.0, 1.0], tiny, 10^9) == SA[256, 256, 256]
        @test _grid_size(SA[1e4, 1.0, 1.0], 0.1, 10^9) == SA[256, 10, 10]
        @test _grid_size(SA[1e4, 1e4, 1e4], tiny, 10^9) == SA[256, 256, 256]

        # Huge widths don't overflow
        @test _grid_size(SA[1e300, 1e300, 1e300], tiny, 8) == SA[2, 2, 2]
        @test _grid_size(SA[3f38, 3f38, 3f38], tiny, 27) == SA[3, 3, 3]
        # Infinite widths from overflow in max_p - min_p get a 1x1x1 grid
        @test _grid_size(SA[Inf, 0.0, 0.0], tiny, 100) == SA[1, 1, 1]
        @test _grid_size(SA[Inf, 0.0, 0.0], huge, 1000) == SA[1, 1, 1]

        # Tiny widths
        @test _grid_size(SA[1f-45, 0f0, 0f0], tiny, 100) == SA[1, 1, 1]
        @test _grid_size(SA[3f-45, 0f0, 0f0], tiny, 100) == SA[2, 1, 1]

        # min_cell_width must be finite and positive
        @test_throws ArgumentError _grid_size(SA[10.0, 10.0, 10.0], 0.0, 1000)
        @test_throws ArgumentError _grid_size(SA[10.0, 10.0, 10.0], -1.0, 1000)
        @test_throws ArgumentError _grid_size(SA[10.0, 10.0, 10.0], NaN, 1000)
        @test_throws ArgumentError _grid_size(SA[10.0, 10.0, 10.0], Inf, 1000)
        @test_throws ArgumentError _grid_size(SA[10.0, 10.0, 10.0], tiny, 0)

        # Random widths satisfy the constraints
        rng = Xoshiro(1234)
        for trial in 1:10^4
            w = SVector{3}(rand(rng, (0.0, 1e-3, 1.0, 10.0, 1e3)) * rand(rng) for _ in 1:3)
            m = rand(rng, (tiny, 1e-3, 0.1, 1.0, 10.0))
            N = rand(rng, 1:10^5)
            n = _grid_size(w, m, N)
            # Cells are at least min_cell_width wide, unless the axis has one cell
            @test prod(n) ≤ N && all(1 .≤ n .≤ 256) && all((n .== 1) .| (w ./ n .≥ m))
        end

        # The bisection in _grid_size matches merging the narrowest cells one
        # step at a time, starting from n_max. Axes with equal cell widths
        # merge together.
        function grid_size_reference(width, min_cell_width, max_cells)
            all(isfinite, width) || return SA[1, 1, 1]
            n = map(w -> Int(clamp(floor(w / min_cell_width), 1, 256)), width)
            prod(n) ≤ max_cells && return n
            while prod(n) > max_cells
                h = minimum(i -> n[i] > 1 ? width[i] / n[i] : Inf, 1:3)
                n = map((w, k) -> k > 1 && w / k ≤ h ? k - 1 : k, width, n)
            end
            n
        end
        for trial in 1:10^5
            # Repeated widths give exact ties between axes
            w = SVector{3}(rand(rng, (0.0, 1e-3, 1.0, 10.0, 10.0, 10.0, 1e3, 1e300, Inf)) * rand(rng, (rand(rng), 1.0)) for _ in 1:3)
            m = rand(rng, (tiny, 1e-3, 0.1, 1.0, 10.0, huge))
            N = rand(rng, (rand(rng, 1:100), rand(rng, 1:10^5), rand(rng, 1:10^9), typemax(Int)))
            n = _grid_size(w, m, N)
            n_ref = grid_size_reference(w, m, N)
            # The bisection stops once it is within one cell on each axis of
            # the reference, and picks the smaller grid. The result must
            # still fit.
            @test all(n_ref .- 1 .≤ n .≤ n_ref) && prod(n) ≤ N
        end
    end

    @testset "Early return with fewer than two objects" begin
        pos = [SVector{3, Float32}(0, 0, 0), SVector{3, Float32}(1, 1, 1)]
        on = DefaultObjectParams(1.0f0, UInt32(1), UInt32(0))
        off = DefaultObjectParams(0.0f0, UInt32(1), UInt32(0))
        two_points = NeighborListInputs(policy;
            points = [PointIdxPart(1), PointIdxPart(2)],
            p_radius = [1.0f0, 1.0f0],
            p_params = [on, on],
        )
        cases = [
            # No objects
            NeighborListInputs(policy),
            # One object
            NeighborListInputs(policy;
                points = [PointIdxPart(1)], p_radius = [1.0f0], p_params = [on]),
            # Two objects, one filtered out
            NeighborListInputs(policy;
                points = [PointIdxPart(1), PointIdxPart(2)],
                p_radius = [1.0f0, 1.0f0],
                p_params = [on, off]),
            # Two objects of different types, both filtered out
            NeighborListInputs(policy;
                points = [PointIdxPart(1)], p_radius = [1.0f0], p_params = [off],
                triangles = [TriangleIdxPart(1, 2, 2)], t_radius = [1.0f0], t_params = [off]),
        ]
        for inputs in cases
            nl = NeighborLists(policy)
            # Stale edges from a previous build must be cleared
            setup_neighbors_naive!(nl, pos, two_points)
            @test length(nl.PPNL) == 1
            @test isnothing(setup_neighbors_grid!(nl, pos, inputs))
            @test is_neighbor_list_subset(nl, NeighborLists(policy))
        end

        # Fewer than two positions
        for pos in ([SVector{3, Float32}(0, 0, 0)], SVector{3, Float32}[])
            nl = NeighborLists(policy)
            @test isnothing(setup_neighbors_grid!(nl, pos, NeighborListInputs(policy)))
            @test is_neighbor_list_subset(nl, NeighborLists(policy))
        end
    end

    all_lists(nl) = (nl.PPNL, nl.PCNL, nl.PLNL, nl.PTNL, nl.CCNL, nl.CLNL, nl.LLNL)
    # The grid builds exactly the edges of the naive build. Subsets both ways
    # plus equal lengths also rules out duplicate edges.
    function test_grid_matches_naive(pos, inputs)
        nl_naive = NeighborLists(policy)
        nl_grid = NeighborLists(policy)
        setup_neighbors_naive!(nl_naive, pos, inputs)
        setup_neighbors_grid!(nl_grid, pos, inputs)
        @test is_neighbor_list_subset(nl_naive, nl_grid)
        @test is_neighbor_list_subset(nl_grid, nl_naive)
        @test map(length, all_lists(nl_grid)) == map(length, all_lists(nl_naive))
    end

    @testset "Naive vs grid consistency" begin
        # Random objects of every type with random layers, masks, zero
        # stiffness, and exclusions for every pair type.
        function random_inputs(rng; T = Float32, box_size = 50.0 + rand(rng) * 100.0,
                rmax = 5.0, skin = 0.0f0, geometry = :cube)
            n_pos = rand(rng, 50:200)
            pos = map(1:n_pos) do _
                p = rand(rng, 3) .* box_size
                geometry == :plane && (p[3] = 5)
                geometry == :line && (p[2] = 5; p[3] = 5)
                geometry == :point && (p .= 5)
                SVector{3, T}(p)
            end
            n = (
                p = rand(rng, 5:300), c = rand(rng, 3:150), l = rand(rng, 3:150),
                lp = rand(rng, 2:80), t = rand(rng, 2:100), tl = rand(rng, 2:80), tp = rand(rng, 2:80),
            )
            radius(k) = Float32.(rand(rng, k) .* rmax .+ 1.0)
            function params(k)
                stiffness = Float32.(rand(rng, k) .* 10.0)
                for i in rand(rng, 1:k, rand(rng, 0:3))
                    stiffness[i] = 0.0f0
                end
                DefaultObjectParams.(stiffness, UInt32.(rand(rng, 1:3, k)), UInt32.(rand(rng, (0, 0, 0, 1, 2), k)))
            end
            no_collide_pairs = empty_no_collide_pairs()
            pair_sizes = (
                (SimplexCellLists.CollidePoint_Point, n.p, n.p),
                (SimplexCellLists.CollidePoint_CLine, n.p, n.c),
                (SimplexCellLists.CollidePoint_Line, n.p, n.l),
                (SimplexCellLists.CollidePoint_Triangle, n.p, n.t),
                (SimplexCellLists.CollideLinePoint_Triangle, n.lp, n.t),
                (SimplexCellLists.CollideTrianglePoint_Triangle, n.tp, n.t),
                (SimplexCellLists.CollideCLine_CLine, n.c, n.c),
                (SimplexCellLists.CollideCLine_Line, n.c, n.l),
                (SimplexCellLists.CollideCLine_TriangleLine, n.c, n.tl),
                (SimplexCellLists.CollideLine_Line, n.l, n.l),
                (SimplexCellLists.CollideLine_TriangleLine, n.l, n.tl),
                (SimplexCellLists.CollideTriangleLine_TriangleLine, n.tl, n.tl),
            )
            for (pair, n_a, n_b) in pair_sizes, _ in 1:rand(rng, 0:20)
                push!(no_collide_pairs[Int(pair)], UInt32(rand(rng, 1:n_a)) => UInt32(rand(rng, 1:n_b)))
            end
            inputs = NeighborListInputs(policy;
                points = [PointIdxPart(rand(rng, 1:n_pos)) for _ in 1:n.p],
                p_radius = radius(n.p), p_params = params(n.p),
                clines = [CLineIdxPart(rand(rng, 1:n_pos-1)) for _ in 1:n.c],
                c_radius = radius(n.c), c_params = params(n.c),
                lines = [LineIdxPart(rand(rng, 1:n_pos), rand(rng, 1:n_pos)) for _ in 1:n.l],
                l_radius = radius(n.l), l_params = params(n.l),
                line_points = [PointIdxPart(rand(rng, 1:n_pos)) for _ in 1:n.lp],
                lp_radius = radius(n.lp), lp_params = params(n.lp),
                triangles = [TriangleIdxPart(rand(rng, 1:n_pos), rand(rng, 1:n_pos), rand(rng, 1:n_pos)) for _ in 1:n.t],
                t_radius = radius(n.t), t_params = params(n.t),
                triangle_lines = [LineIdxPart(rand(rng, 1:n_pos), rand(rng, 1:n_pos)) for _ in 1:n.tl],
                tl_radius = radius(n.tl), tl_params = params(n.tl),
                triangle_points = [PointIdxPart(rand(rng, 1:n_pos)) for _ in 1:n.tp],
                tp_radius = radius(n.tp), tp_params = params(n.tp),
                no_collide_pairs, skin,
            )
            pos, inputs
        end
        cases = (
            (;),
            (; skin = 2.0f0),
            (; T = Float64),
            # Dense, many objects per cell
            (; box_size = 10.0),
            # Sparse, cells merged to fit the object count
            (; box_size = 2000.0),
            # Boxes spanning many cells
            (; rmax = 100.0),
            (; geometry = :plane),
            (; geometry = :line),
            (; geometry = :point),
        )
        rng = Xoshiro(2024)
        for kw in cases, trial in 1:50
            test_grid_matches_naive(random_inputs(rng; kw...)...)
        end
    end

    @testset "Filaments" begin
        # Random walk filaments of points, clines, and triangles, like a
        # simulation, optionally excluding bonded neighbors.
        function filament_inputs(rng, n_beads; exclude_bonded, skin)
            L = 50
            box = (n_beads / 0.5)^(1/3)
            pos = SVector{3, Float32}[]
            points = PointIdxPart[]
            clines = CLineIdxPart[]
            triangles = TriangleIdxPart[]
            for f in 1:cld(n_beads, L)
                p = SVector{3, Float64}(rand(rng, 3) .* box)
                for b in 1:L
                    push!(pos, p)
                    i = length(pos)
                    push!(points, PointIdxPart(i))
                    b < L && push!(clines, CLineIdxPart(i))
                    b ≤ L - 2 && isodd(b) && push!(triangles, TriangleIdxPart(i, i + 1, i + 2))
                    d = randn(rng, SVector{3, Float64})
                    p += d / sqrt(sum(abs2, d))
                end
            end
            no_collide_pairs = empty_no_collide_pairs()
            if exclude_bonded
                for i in 1:length(points) - 1
                    push!(no_collide_pairs[Int(SimplexCellLists.CollidePoint_Point)], UInt32(i) => UInt32(i + 1))
                end
                for i in 1:length(clines) - 1
                    push!(no_collide_pairs[Int(SimplexCellLists.CollideCLine_CLine)], UInt32(i) => UInt32(i + 1))
                end
            end
            object_params(k) = fill(DefaultObjectParams(1.0f0, UInt32(1), UInt32(0)), k)
            inputs = NeighborListInputs(policy;
                points, p_radius = fill(0.5f0, length(points)), p_params = object_params(length(points)),
                clines, c_radius = fill(0.5f0, length(clines)), c_params = object_params(length(clines)),
                triangles, t_radius = fill(0.5f0, length(triangles)), t_params = object_params(length(triangles)),
                no_collide_pairs, skin,
            )
            pos, inputs
        end
        rng = Xoshiro(7)
        for exclude_bonded in (false, true), skin in (0.0f0, 0.1f0), trial in 1:3
            test_grid_matches_naive(filament_inputs(rng, 1000; exclude_bonded, skin)...)
        end
    end

    @testset "Degenerate geometry" begin
        points = [PointIdxPart(i) for i in 1:4]
        p_radius = fill(5.0f0, 4)
        p_params = fill(DefaultObjectParams(1.0f0, UInt32(1), UInt32(0)), 4)
        inputs = NeighborListInputs(policy; points, p_radius, p_params, skin = 7.0f0)
        nl = NeighborLists(policy)

        # Coincident positions, near and far from the origin, give a 1x1x1 grid
        for p in (SVector{3, Float32}(1, 2, 3), SVector{3, Float32}(1f7, 1f7, 1f7))
            pos = fill(p, 4)
            setup_neighbors_grid!(nl, pos, inputs)
            @test length(nl.PPNL) == 6
            test_grid_matches_naive(pos, inputs)
        end

        # The bounds width overflowing to Inf errors
        pos = [
            SVector{3, Float32}(-3f38, 0, 0), SVector{3, Float32}(3f38, 0, 0),
            SVector{3, Float32}(0, 0, 0), SVector{3, Float32}(1, 0, 0),
        ]
        @test_throws ErrorException setup_neighbors_grid!(nl, pos, inputs)

        # Non finite positions error
        for bad in (NaN32, Inf32, -Inf32)
            pos_bad = fill(SVector{3, Float32}(1, 2, 3), 4)
            pos_bad[2] = SVector{3, Float32}(1, bad, 3)
            @test_throws ErrorException setup_neighbors_grid!(nl, pos_bad, inputs)
        end
    end

    @testset "Negative or NaN radius or skin" begin
        pos = [SVector{3, Float32}(0, 0, 0), SVector{3, Float32}(1, 1, 1)]
        on = DefaultObjectParams(1.0f0, UInt32(1), UInt32(0))
        off = DefaultObjectParams(0.0f0, UInt32(1), UInt32(0))
        points = [PointIdxPart(1), PointIdxPart(2)]
        nl = NeighborLists(policy)
        for bad in (-1.0f0, -Inf32, NaN32)
            # Radius of an object that passes filter_object
            inputs = NeighborListInputs(policy; points, p_radius = [1.0f0, bad], p_params = [on, on])
            @test_throws ArgumentError setup_neighbors_grid!(nl, pos, inputs)
            # Even if it is the only object that passes filter_object
            inputs = NeighborListInputs(policy; points, p_radius = [1.0f0, bad], p_params = [off, on])
            @test_throws ArgumentError setup_neighbors_grid!(nl, pos, inputs)
            # Skin, even with fewer than two positions
            inputs = NeighborListInputs(policy; points, p_radius = [1.0f0, 1.0f0], p_params = [on, on], skin = bad)
            @test_throws ArgumentError setup_neighbors_grid!(nl, pos, inputs)
            @test_throws ArgumentError setup_neighbors_grid!(nl, pos[1:1], NeighborListInputs(policy; skin = bad))
        end
        # Radii of objects that fail filter_object are ignored
        inputs = NeighborListInputs(policy;
            points = [PointIdxPart(1), PointIdxPart(2), PointIdxPart(1)],
            p_radius = [1.0f0, 1.0f0, -1.0f0],
            p_params = [on, on, off],
        )
        setup_neighbors_grid!(nl, pos, inputs)
        @test length(nl.PPNL) == 1
        # Negative zero is allowed
        inputs = NeighborListInputs(policy; points, p_radius = [-0.0f0, 1.0f0], p_params = [on, on], skin = -0.0f0)
        test_grid_matches_naive(pos, inputs)
    end

    @testset "Repeated builds" begin
        # Building twice into the same lists gives the same edges, not more
        pos = [SVector{3, Float32}(i, 0, 0) for i in 1:10]
        inputs = NeighborListInputs(policy;
            points = [PointIdxPart(i) for i in 1:10],
            p_radius = fill(1.0f0, 10),
            p_params = fill(DefaultObjectParams(1.0f0, UInt32(1), UInt32(0)), 10),
        )
        nl = NeighborLists(policy)
        setup_neighbors_grid!(nl, pos, inputs)
        edges = copy(nl.PPNL)
        @test length(edges) == 9
        setup_neighbors_grid!(nl, pos, inputs)
        @test nl.PPNL == edges
    end
end
