# Test Naive
using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Random
using Test
using ArgCheck
using Rotations

include("common.jl")


@testset "constructor" begin
    empty_s = SimplexCellLists.Naive(0, 0, 0)
    onepointgroup_s = SimplexCellLists.Naive(1, 0, 0)
    onelinegroup_s = SimplexCellLists.Naive(0, 1, 0)
    onetrianglegroup_s = SimplexCellLists.Naive(0, 0, 1)
    oneallgroup_s = SimplexCellLists.Naive(1, 1, 1)
    tenallgroup_s = SimplexCellLists.Naive(10, 10, 10)
end
@testset "setElements!" begin
    empty_s = SimplexCellLists.Naive(0, 0, 0)
    setElements!(empty_s, [], [], [])

    onepointgroup_s = SimplexCellLists.Naive(1, 0, 0)
    setElements!(onepointgroup_s, [[[[1,2,3]],[[3,4,6]]]], [], [])

    oneallgroup_s = SimplexCellLists.Naive(1, 1, 1)
    setElements!(oneallgroup_s, 
        [[[[1,2,3]],[[3,4,6]]]],
        [[[[1,2,3],[7,2,5]],[[3,4,6],[3,9,6]]]],
        [[[[1,2,3],[7,2,5],[7,9,5]],[[3,4,6],[3,9,6],[1,2,3]]]],
    )

    twoallgroup_s = SimplexCellLists.Naive(2, 2, 2)
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
    onepointgroup_s = SimplexCellLists.Naive(1, 0, 0)
    setElements!(onepointgroup_s, [[[[1,2,3]],[[3,4,6]]]], [], [])
    @test 3 == addElement!(onepointgroup_s, 1, SA[SA_F32[3,4,6]])
    @test 4 == addElement!(onepointgroup_s, 1, SA[SA_F32[3,4,6]])

    oneallgroup_s = SimplexCellLists.Naive(1, 1, 1)
    setElements!(oneallgroup_s, 
        [[[[1,2,3]],[[3,4,6]]]],
        [[[[1,2,3],[7,2,5]],[[3,4,6],[3,9,6]]]],
        [[[[1,2,3],[7,2,5],[7,9,5]],[[3,4,6],[3,9,6],[1,2,3]]]],
    )
    @test 3 == addElement!(oneallgroup_s, 1, SA[SA_F32[3,4,6]])
    @test 3 == addElement!(oneallgroup_s, 1, SA[SA_F32[3,4,6], SA_F32[4,4,6]])
    @test 3 == addElement!(oneallgroup_s, 1, SA[SA_F32[3,4,6], SA_F32[4,4,6], SA_F32[5,5,6]])

    twoallgroup_s = SimplexCellLists.Naive(2, 2, 2)
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
        s = SimplexCellLists.Naive(2, 2, 2)
        mapSimplexElements(test_f_mapSimplexElements!, 0, s, SA[SA_F32[3,4,6]], 1 ,SimplexCellLists.Line, 10.0f0)
    end
    @testset "basic tests" begin
        for i in 10*trials
            scale = rand()*100.0
            rotation = rand(RotMat{3,Float32})
            translation = randn(SVector{3,Float32}) * scale
            s = SimplexCellLists.Naive(2, 2, 2)
            makeBasicCellList!(
                scale = 1.0f0,
                rotation = one(RotMat{3,Float32}),
                translation = SA_F32[0,0,0],
                emptycelllist::SimplexCellList,
            )
        mapSimplexElements(test_f_mapSimplexElements!, 0, s, SA[SA_F32[3,4,6]], 1 ,SimplexCellLists.Line, 10.0f0)
    end
end