# Test Painter
using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Random
using Test
using ArgCheck
using Rotations

include("common.jl")


@testset "constructor" begin
    kwargs = (; grid_start=SA_F32[0,0,0], grid_size=SA[4,4,4], voxel_length=1.5)
    empty_s = SimplexCellLists.Painter(0, 0, 0; kwargs..., max_range=SA[Float64[],Float64[],Float64[]] )
    onepointgroup_s = SimplexCellLists.Painter(1, 0, 0; kwargs..., max_range=SA[[1.2],Float64[],Float64[]] )
    onelinegroup_s = SimplexCellLists.Painter(0, 1, 0; kwargs..., max_range=SA[Float64[],[1.2],Float64[]] )
    onetrianglegroup_s = SimplexCellLists.Painter(0, 0, 1; kwargs..., max_range=SA[Float64[],Float64[],[1.2]] )
    oneallgroup_s = SimplexCellLists.Painter(1, 1, 1; kwargs..., max_range=SA[[1.2],[1.2],[1.2]] )
    tenallgroup_s = SimplexCellLists.Painter(10, 10, 10; kwargs..., max_range=SA[fill(1.2,10),fill(1.2,10),fill(1.2,10)] )
end
@testset "setElements!" begin
    kwargs = (; grid_start=SA_F32[0,0,0], grid_size=SA[4,4,4], voxel_length=1.5)
    empty_s = SimplexCellLists.Painter(0, 0, 0; kwargs..., max_range=SA[Float64[],Float64[],Float64[]])
    setElements!(empty_s, [], [], [])

    onepointgroup_s = SimplexCellLists.Painter(1, 0, 0; kwargs..., max_range=SA[Float64[1.3],Float64[],Float64[]])
    setElements!(onepointgroup_s, [[[[1,2,3]],[[3,4,6]]]], [], [])

    oneallgroup_s = SimplexCellLists.Painter(1, 1, 1; kwargs..., max_range=SA[Float64[1.2],Float64[1.5],Float64[1.3]])
    setElements!(oneallgroup_s, 
        [[[[1,2,3]],[[3,4,6]]]],
        [[[[1,2,3],[7,2,5]],[[3,4,6],[3,9,6]]]],
        [[[[1,2,3],[7,2,5],[7,9,5]],[[3,4,6],[3,9,6],[1,2,3]]]],
    )

    twoallgroup_s = SimplexCellLists.Painter(2, 2, 2; kwargs..., max_range=SA[Float64[0.3,0.7],Float64[1.0,2.2],Float64[0.5,0.7]])
    setElements!(twoallgroup_s, 
        [[[[1,2,3]],[[3,4,6]]],[]],
        [[],[[[1,2,3],[7,2,5]],[[3,4,6],[3,9,6]]]],
        [
            [[[1,2,3],[7,2,5],[7,9,5]],[[3,4,6],[3,9,6],[1,2,3]]],
            [[[1,2,3],[7,2,5],[7,9,5]],[[3,4,6],[3,9,6],[1,2,3]]],
        ],
    )
end
@testset "addElement!" begin
    kwargs = (; grid_start=SA_F32[0,0,0], grid_size=SA[4,4,4], voxel_length=1.5)
    onepointgroup_s = SimplexCellLists.Painter(1, 0, 0; kwargs..., max_range=SA[Float64[0.2],Float64[],Float64[]])
    setElements!(onepointgroup_s, [[[[1,2,3]],[[3,4,6]]]], [], [])
    @test 3 == addElement!(onepointgroup_s, 1, SA[SA_F32[3,4,6]])
    @test 4 == addElement!(onepointgroup_s, 1, SA[SA_F32[3,4,6]])

    oneallgroup_s = SimplexCellLists.Painter(1, 1, 1; kwargs..., max_range=SA[Float64[0.3],Float64[0.5],Float64[0.7]])
    setElements!(oneallgroup_s, 
        [[[[1,2,3]],[[3,4,6]]]],
        [[[[1,2,3],[7,2,5]],[[3,4,6],[3,9,6]]]],
        [[[[1,2,3],[7,2,5],[7,9,5]],[[3,4,6],[3,9,6],[1,2,3]]]],
    )
    @test 3 == addElement!(oneallgroup_s, 1, SA[SA_F32[3,4,6]])
    @test 3 == addElement!(oneallgroup_s, 1, SA[SA_F32[3,4,6], SA_F32[4,4,6]])
    @test 3 == addElement!(oneallgroup_s, 1, SA[SA_F32[3,4,6], SA_F32[4,4,6], SA_F32[5,5,6]])

    twoallgroup_s = SimplexCellLists.Painter(2, 2, 2; kwargs..., max_range=SA[Float64[0.3,0.7],Float64[1.0,2.2],Float64[0.5,0.7]])
    setElements!(twoallgroup_s, 
        [[[[1,2,3]],[[3,4,6]]],[]],
        [[],[[[1,2,3],[7,2,5]],[[3,4,6],[3,9,6]]]],
        [[],[]],
    )
    @test 1 == addElement!(twoallgroup_s, 2, SA[SA_F32[3,4,6]])
    @test 3 == addElement!(twoallgroup_s, 1, SA[SA_F32[3,4,6]])
    @test 3 == addElement!(twoallgroup_s, 2, SA[SA_F32[3,4,6], SA_F32[4,4,6]])
    @test 1 == addElement!(twoallgroup_s, 1, SA[SA_F32[3,4,6], SA_F32[4,4,6]])
    @test 1 == addElement!(twoallgroup_s, 1, SA[SA_F32[3,4,6], SA_F32[4,4,6], SA_F32[5,5,6]])
end
@testset "mapSimplexElements" begin
    @testset "empty cell list" begin
        s = SimplexCellLists.Painter(2, 2, 2; grid_start=SA_F32[0,0,0], grid_size=SA[4,4,4], voxel_length=1.5, max_range=SA[Float64[0.3,0.7],Float64[10.1,2.2],Float64[0.5,0.7]])
        out = zeros(Int,2)
        mapSimplexElements(test_f_mapSimplexElements!, out, s, SA[SA_F32[3,4,6]], 1 ,SimplexCellLists.Line, 10.0f0)
        @test out == zero(out)
    end
    @testset "mapSimplexElements tests" begin
        for trial in 1:10
            scale = Float32(exp(rand()*8.0 - 4))
            rotation = rand(RotMatrix{3,Float32})
            translation = randn(SVector{3,Float32}) * scale
            function trans(in)
                pts = reinterpret(SVector{3,Float32},in)
                map(pts) do pt
                    rotation * (pt * scale) + translation 
                end
            end
            s = SimplexCellLists.Painter(2, 2, 2; grid_start=SA_F32[-5,-5,-5], grid_size=SA[10,10,10], voxel_length=scale*rand()*10, max_range=6.1*scale)
            elements = makeBasicCellList!(
                scale,
                rotation,
                translation,
                s,
            )
            simplex_types = (SimplexCellLists.Point, SimplexCellLists.Line, SimplexCellLists.Triangle)
            group_idxs = (1,2)
            i = 0
            for group_type in simplex_types
                for idx in group_idxs
                    i+=1
                    cutoff = scale*6.001f0
                    f(x...) = test_f_mapSimplexElements!(x..., elements[i], cutoff)
                    out = mapSimplexElements( f, zeros(Int,2), s, trans(SA[SA_F32[0,0,0]]), idx ,group_type, cutoff)
                    @test out == [1,1]

                    cutoff = scale*5.99f0
                    out = mapSimplexElements( f, zeros(Int,2), s, trans(SA[SA_F32[0,0,0]]), idx ,group_type, cutoff)
                    @test out == [0,0]

                    cutoff = scale*6.0f0
                    out = mapSimplexElements( f, zeros(Int,2), s, trans(SA[SA_F32[0.01,0,0]]), idx ,group_type, cutoff)
                    @test out == (idx==1 ? [0,1] : [1,0])
                end
            end
            i = 0
            for group_type in simplex_types
                for idx in group_idxs
                    i+=1
                    cutoff = scale*6.001f0
                    f(x...) = test_f_mapSimplexElements!(x..., elements[i], cutoff)
                    out = mapSimplexElements( f, zeros(Int,2), s, trans(SA[SA_F32[0,0,0],SA_F32[0,0,1]]), idx ,group_type, cutoff)
                    @test out == [1,1]

                    cutoff = scale*5.99f0
                    out = mapSimplexElements( f, zeros(Int,2), s, trans(SA[SA_F32[0,0,0],SA_F32[0,0,1]]), idx ,group_type, cutoff)
                    @test out == [0,0]

                    cutoff = scale*6.0f0
                    out = mapSimplexElements( f, zeros(Int,2), s, trans(SA[SA_F32[0.01,0,0],SA_F32[0,0,1]]), idx ,group_type, cutoff)
                    @test out == (idx==1 ? [0,1] : [1,0])
                end
            end
            i = 0
            for group_type in simplex_types
                for idx in group_idxs
                    i+=1
                    cutoff = scale*6.001f0
                    f(x...) = test_f_mapSimplexElements!(x..., elements[i], cutoff)
                    out = mapSimplexElements( f, zeros(Int,2), s, trans(SA[SA_F32[0,0,0],SA_F32[0,0,1],SA_F32[0,1,1]]), idx ,group_type, cutoff)
                    @test out == [1,1]

                    cutoff = scale*5.99f0
                    out = mapSimplexElements( f, zeros(Int,2), s, trans(SA[SA_F32[0,0,0],SA_F32[0,0,1],SA_F32[0,1,1]]), idx ,group_type, cutoff)
                    @test out == [0,0]

                    cutoff = scale*6.0f0
                    out = mapSimplexElements( f, zeros(Int,2), s, trans(SA[SA_F32[0.01,0,0],SA_F32[0,0,1],SA_F32[0,1,1]]), idx ,group_type, cutoff)
                    @test out == (idx==1 ? [0,1] : [1,0])
                end
            end
        end
    end
end