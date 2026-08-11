using SimplexCellLists
using Test
using StaticArrays
using Random

@testset "neighborlists" begin
    policy = DefaultCollidePolicy()
    @testset "Naive vs Sort-Sweep consistency" begin

        for trial in 1:2000
            # Random box size
            box_size = 50.0 + rand() * 100.0

            # Generate random positions
            n_pos = rand(50:200)
            pos = [SVector{3, Float32}(rand(3) .* box_size) for _ in 1:n_pos]

            # Random number of each object type
            n_points = rand(5:300)
            n_clines = rand(3:150)
            n_lines = rand(3:150)
            n_triangles = rand(2:100)
            n_line_points = rand(2:80)
            n_triangle_lines = rand(2:80)
            n_triangle_points = rand(2:80)

            # Build random objects
            points = [PointIdxPart(rand(1:n_pos)) for _ in 1:n_points]
            p_radius = Float32.(rand(n_points) .* 5.0 .+ 1.0)
            p_stiffness = Float32.(rand(n_points) .* 10.0)
            # Make some zero stiffness
            for i in rand(1:n_points, rand(0:3))
                p_stiffness[i] = 0.0f0
            end

            clines = [CLineIdxPart(rand(1:n_pos-1)) for _ in 1:n_clines]
            c_radius = Float32.(rand(n_clines) .* 5.0 .+ 1.0)
            c_stiffness = Float32.(rand(n_clines) .* 10.0)
            for i in rand(1:n_clines, rand(0:2))
                c_stiffness[i] = 0.0f0
            end

            lines = [LineIdxPart(rand(1:n_pos), rand(1:n_pos)) for _ in 1:n_lines]
            l_radius = Float32.(rand(n_lines) .* 5.0 .+ 1.0)
            l_stiffness = Float32.(rand(n_lines) .* 10.0)
            for i in rand(1:n_lines, rand(0:2))
                l_stiffness[i] = 0.0f0
            end

            line_points = [PointIdxPart(rand(1:n_pos)) for _ in 1:n_line_points]
            lp_radius = Float32.(rand(n_line_points) .* 5.0 .+ 1.0)
            lp_stiffness = Float32.(rand(n_line_points) .* 10.0)

            triangles = [TriangleIdxPart(rand(1:n_pos), rand(1:n_pos), rand(1:n_pos)) for _ in 1:n_triangles]
            t_radius = Float32.(rand(n_triangles) .* 5.0 .+ 1.0)
            t_stiffness = Float32.(rand(n_triangles) .* 10.0)
            for i in rand(1:n_triangles, rand(0:2))
                t_stiffness[i] = 0.0f0
            end

            triangle_lines = [LineIdxPart(rand(1:n_pos), rand(1:n_pos)) for _ in 1:n_triangle_lines]
            tl_radius = Float32.(rand(n_triangle_lines) .* 5.0 .+ 1.0)
            tl_stiffness = Float32.(rand(n_triangle_lines) .* 10.0)

            triangle_points = [PointIdxPart(rand(1:n_pos)) for _ in 1:n_triangle_points]
            tp_radius = Float32.(rand(n_triangle_points) .* 5.0 .+ 1.0)
            tp_stiffness = Float32.(rand(n_triangle_points) .* 10.0)

            # Build random exclusion pairs
            no_collide_pairs = empty_no_collide_pairs()
            # Add some random exclusions
            for _ in 1:rand(0:10)
                push!(no_collide_pairs[Int(SimplexCellLists.CollidePoint_Point)], UInt32(rand(1:n_points)) => UInt32(rand(1:n_points)))
            end
            for _ in 1:rand(0:5)
                push!(no_collide_pairs[Int(SimplexCellLists.CollidePoint_CLine)], UInt32(rand(1:n_points)) => UInt32(rand(1:n_clines)))
            end
            for _ in 1:rand(0:5)
                push!(no_collide_pairs[Int(SimplexCellLists.CollidePoint_Line)], UInt32(rand(1:n_points)) => UInt32(rand(1:n_lines)))
            end
            for _ in 1:rand(0:5)
                push!(no_collide_pairs[Int(SimplexCellLists.CollidePoint_Triangle)], UInt32(rand(1:n_points)) => UInt32(rand(1:n_triangles)))
            end
            for _ in 1:rand(0:5)
                push!(no_collide_pairs[Int(SimplexCellLists.CollideCLine_CLine)], UInt32(rand(1:n_clines)) => UInt32(rand(1:n_clines)))
            end
            for _ in 1:rand(0:5)
                push!(no_collide_pairs[Int(SimplexCellLists.CollideCLine_Line)], UInt32(rand(1:n_clines)) => UInt32(rand(1:n_lines)))
            end
            for _ in 1:rand(0:5)
                push!(no_collide_pairs[Int(SimplexCellLists.CollideLine_Line)], UInt32(rand(1:n_lines)) => UInt32(rand(1:n_lines)))
            end

            # Object params (default: all on layer 1, no collisions disabled)
            p_params = DefaultObjectParams.(p_stiffness, UInt32(1), UInt32(0))
            c_params = DefaultObjectParams.(c_stiffness, UInt32(1), UInt32(0))
            l_params = DefaultObjectParams.(l_stiffness, UInt32(1), UInt32(0))
            lp_params = DefaultObjectParams.(lp_stiffness, UInt32(1), UInt32(0))
            t_params = DefaultObjectParams.(t_stiffness, UInt32(1), UInt32(0))
            tl_params = DefaultObjectParams.(tl_stiffness, UInt32(1), UInt32(0))
            tp_params = DefaultObjectParams.(tp_stiffness, UInt32(1), UInt32(0))

            inputs = NeighborListInputs(policy;
                points, p_radius, p_params,
                clines, c_radius, c_params,
                lines, l_radius, l_params,
                line_points, lp_radius, lp_params,
                triangles, t_radius, t_params,
                triangle_lines, tl_radius, tl_params,
                triangle_points, tp_radius, tp_params,
                no_collide_pairs,
                skin = 0.0f0,
            )

            nl_naive = NeighborLists(policy)
            nl_sweep = NeighborLists(policy)

            setup_neighbors_naive!(nl_naive, pos, inputs)
            setup_neighbors_sort_sweep!(nl_sweep, pos, inputs)

            # Check equivalence: both are subsets of each other
            @test is_neighbor_list_subset(nl_naive, nl_sweep)
            @test is_neighbor_list_subset(nl_sweep, nl_naive)
        end
    end
    @testset "can_collide" begin
        # can_collide(layers_a, no_collide_mask_a, layers_b, no_collide_mask_b)
        # Default: both on layer 1, blocking nothing → collide
        @test can_collide(UInt32(1), UInt32(0), UInt32(1), UInt32(0))
        # A on layer 1, B blocks layer 1 → no collide
        @test !can_collide(UInt32(1), UInt32(0), UInt32(2), UInt32(1))
        # A on layer 2 blocks layer 2, B on layer 2 blocks layer 2 → no collide (both block each other)
        @test !can_collide(UInt32(2), UInt32(2), UInt32(2), UInt32(2))
        # A on layer 1 blocks nothing, B on layer 2 blocks nothing → collide (no blocking)
        @test can_collide(UInt32(1), UInt32(0), UInt32(2), UInt32(0))
        # A on layer 2 blocks layer 1, B on layer 1 blocks nothing → no collide (A blocks B)
        @test !can_collide(UInt32(2), UInt32(1), UInt32(1), UInt32(0))
        # One-sided block: A on layer 1 blocks layer 2, B on layer 2 blocks nothing → no collide
        @test !can_collide(UInt32(1), UInt32(2), UInt32(2), UInt32(0))
        # Zero layers, zero masks → collide (nothing to block)
        @test can_collide(UInt32(0), UInt32(0), UInt32(0), UInt32(0))
        # A on layers 1&2, B blocks layer 2 → no collide (B blocks A's layer 2)
        @test !can_collide(UInt32(3), UInt32(0), UInt32(1), UInt32(2))
    end

    @testset "No-collide mask filtering" begin
        pos = [SVector{3, Float32}(0,0,0), SVector{3, Float32}(1,0,0)]

        points = [PointIdxPart(1), PointIdxPart(2)]
        p_radius = Float32[5.0, 5.0]
        p_stiffness = Float32[1.0, 1.0]

        # Default: same layer, blocking nothing → collide
        inputs_default = NeighborListInputs(policy;
            points, p_radius,
            p_params = DefaultObjectParams.(p_stiffness, UInt32[1, 1], UInt32[0, 0]),
        )
        nl = NeighborLists(policy)
        setup_neighbors_naive!(nl, pos, inputs_default)
        @test length(nl.PPNL) == 1
        nl_sweep = NeighborLists(policy)
        setup_neighbors_sort_sweep!(nl_sweep, pos, inputs_default)
        @test length(nl_sweep.PPNL) == 1

        # Both on layer 2, both block layer 2 → no collision (piston example)
        inputs_block = NeighborListInputs(policy;
            points, p_radius,
            p_params = DefaultObjectParams.(p_stiffness, UInt32[2, 2], UInt32[2, 2]),
        )
        nl2 = NeighborLists(policy)
        setup_neighbors_naive!(nl2, pos, inputs_block)
        @test length(nl2.PPNL) == 0
        nl2_sweep = NeighborLists(policy)
        setup_neighbors_sort_sweep!(nl2_sweep, pos, inputs_block)
        @test length(nl2_sweep.PPNL) == 0

        # One-sided block: A blocks layer 2, B on layer 2 → no collision
        inputs_onesided = NeighborListInputs(policy;
            points, p_radius,
            p_params = DefaultObjectParams.(p_stiffness, UInt32[1, 2], UInt32[2, 0]),
        )
        nl3 = NeighborLists(policy)
        setup_neighbors_naive!(nl3, pos, inputs_onesided)
        @test length(nl3.PPNL) == 0
        nl3_sweep = NeighborLists(policy)
        setup_neighbors_sort_sweep!(nl3_sweep, pos, inputs_onesided)
        @test length(nl3_sweep.PPNL) == 0

        # Different layers, no blocking → collide
        inputs_diff_no_block = NeighborListInputs(policy;
            points, p_radius,
            p_params = DefaultObjectParams.(p_stiffness, UInt32[1, 2], UInt32[0, 0]),
        )
        nl4 = NeighborLists(policy)
        setup_neighbors_naive!(nl4, pos, inputs_diff_no_block)
        @test length(nl4.PPNL) == 1
        nl4_sweep = NeighborLists(policy)
        setup_neighbors_sort_sweep!(nl4_sweep, pos, inputs_diff_no_block)
        @test length(nl4_sweep.PPNL) == 1
    end

    @testset "Collision layers with CLines" begin
        # Two CLines close together on different layers should not collide
        pos = [
            SVector{3, Float32}(0,0,0), SVector{3, Float32}(10,0,0),
            SVector{3, Float32}(0,1,0), SVector{3, Float32}(10,1,0),
        ]
        clines = [CLineIdxPart(1), CLineIdxPart(3)]
        c_radius = Float32[5.0, 5.0]
        c_stiffness = Float32[1.0, 1.0]

        # Same layer, no blocking → collide
        inputs_same = NeighborListInputs(policy;
            clines, c_radius,
            c_params = DefaultObjectParams.(c_stiffness, UInt32[1, 1], UInt32[0, 0]),
        )
        nl = NeighborLists(policy)
        setup_neighbors_naive!(nl, pos, inputs_same)
        @test length(nl.CCNL) == 1

        # Both on layer 2, both block layer 2 → no collide
        inputs_block = NeighborListInputs(policy;
            clines, c_radius,
            c_params = DefaultObjectParams.(c_stiffness, UInt32[2, 2], UInt32[2, 2]),
        )
        nl2 = NeighborLists(policy)
        setup_neighbors_naive!(nl2, pos, inputs_block)
        @test length(nl2.CCNL) == 0
    end

    @testset "Collision layers naive vs sweep consistency with random layers" begin
        for trial in 1:50
            box_size = 50.0 + rand() * 100.0
            n_pos = rand(50:200)
            pos = [SVector{3, Float32}(rand(3) .* box_size) for _ in 1:n_pos]

            n_points = rand(5:50)
            points = [PointIdxPart(rand(1:n_pos)) for _ in 1:n_points]
            p_radius = Float32.(rand(n_points) .* 5.0 .+ 1.0)
            p_stiffness = Float32.(rand(n_points) .* 10.0 .+ 0.1)
            # Random layers from a small set to get interesting interactions
            p_layers = UInt32.(rand(1:3, n_points))
            p_ncmask = UInt32.(rand(1:3, n_points))

            n_clines = rand(3:30)
            clines = [CLineIdxPart(rand(1:n_pos-1)) for _ in 1:n_clines]
            c_radius = Float32.(rand(n_clines) .* 5.0 .+ 1.0)
            c_stiffness = Float32.(rand(n_clines) .* 10.0 .+ 0.1)
            c_layers = UInt32.(rand(1:3, n_clines))
            c_ncmask = UInt32.(rand(1:3, n_clines))

            inputs = NeighborListInputs(policy;
                points, p_radius,
                p_params = DefaultObjectParams.(p_stiffness, p_layers, p_ncmask),
                clines, c_radius,
                c_params = DefaultObjectParams.(c_stiffness, c_layers, c_ncmask),
                skin = 0.0f0,
            )

            nl_naive = NeighborLists(policy)
            nl_sweep = NeighborLists(policy)
            setup_neighbors_naive!(nl_naive, pos, inputs)
            setup_neighbors_sort_sweep!(nl_sweep, pos, inputs)

            @test is_neighbor_list_subset(nl_naive, nl_sweep)
            @test is_neighbor_list_subset(nl_sweep, nl_naive)
        end
    end

    @testset "Degenerate tiny bounding box quantization" begin
        # All positions coincident with ordinary radii: the quantization scale
        # explodes, which used to overflow the Int32 conversion in the sweep.
        pos = [SVector{3, Float32}(1, 2, 3) for _ in 1:4]
        points = [PointIdxPart(i) for i in 1:4]
        p_radius = fill(5.0f0, 4)
        p_params = fill(DefaultObjectParams(1.0f0, UInt32(1), UInt32(0)), 4)
        inputs = NeighborListInputs(policy; points, p_radius, p_params, skin = 7.0f0)
        nl_naive = NeighborLists(policy)
        nl_sweep = NeighborLists(policy)
        setup_neighbors_naive!(nl_naive, pos, inputs)
        setup_neighbors_sort_sweep!(nl_sweep, pos, inputs)
        @test length(nl_sweep.PPNL) == 6
        @test is_neighbor_list_subset(nl_naive, nl_sweep)
        @test is_neighbor_list_subset(nl_sweep, nl_naive)

        # Same, but coincident far from the origin, where the minimum
        # bounding box width must scale with the coordinate magnitude
        pos_far = [SVector{3, Float32}(1f7, 1f7, 1f7) for _ in 1:4]
        setup_neighbors_naive!(nl_naive, pos_far, inputs)
        setup_neighbors_sort_sweep!(nl_sweep, pos_far, inputs)
        @test length(nl_sweep.PPNL) == 6
        @test is_neighbor_list_subset(nl_naive, nl_sweep)
        @test is_neighbor_list_subset(nl_sweep, nl_naive)

        # NaN geometry errors instead of building a bogus list
        pos_nan = copy(pos)
        pos_nan[1] = SVector{3, Float32}(NaN32, 2, 3)
        @test_throws ErrorException setup_neighbors_sort_sweep!(nl_sweep, pos_nan, inputs)
    end

    @testset "Planar and collinear geometry" begin
        # Zero bounding box width on non-sweep axes with a finite scale
        for trial in 1:20
            n_pos = rand(20:60)
            n_points = rand(5:40)
            p_radius = Float32.(rand(n_points) .* 5.0 .+ 1.0)
            p_params = fill(DefaultObjectParams(1.0f0, UInt32(1), UInt32(0)), n_points)
            points = [PointIdxPart(rand(1:n_pos)) for _ in 1:n_points]
            inputs = NeighborListInputs(policy; points, p_radius, p_params)
            # All positions in a z plane, and all positions on an x line
            planar = [SVector{3, Float32}(rand()*100, rand()*100, 5) for _ in 1:n_pos]
            collinear = [SVector{3, Float32}(rand()*100, 5, 5) for _ in 1:n_pos]
            for pos in (planar, collinear)
                nl_naive = NeighborLists(policy)
                nl_sweep = NeighborLists(policy)
                setup_neighbors_naive!(nl_naive, pos, inputs)
                setup_neighbors_sort_sweep!(nl_sweep, pos, inputs)
                @test is_neighbor_list_subset(nl_naive, nl_sweep)
                @test is_neighbor_list_subset(nl_sweep, nl_naive)
            end
        end
    end
end
nothing
