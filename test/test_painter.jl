# Test Painter

using StaticArrays
using LinearAlgebra
using MultiShapeCellLists
using Random
using Test


Random.seed!(1234)

@testset "just points" begin
    N = 1000
    points = rand(SVector{1,SVector{3,Float32}},N)
    naive = MultiShapeCellLists.Naive(1, 0)
    painter = MultiShapeCellLists.Painter(1,0;
        grid_start= SA[0.0,0.0,0.0],
        grid_size= SA[10,10,10],
        voxel_length= 1/10,
        max_range= SA[[0.1],Float64[]],
    )
    MultiShapeCellLists.setElements(naive,[points],[])
    MultiShapeCellLists.setElements(painter,[points],[])
    MultiShapeCellLists.addElement(naive,1,points[1])
    MultiShapeCellLists.addElement(painter,1,points[1])
    cutoff = 0.09f0
    f(x,y,i,j,d2,output) = output+((√(d2)-cutoff)^2)
    naive_out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0.0,
        naive,
        1,
        points[1],
        MultiShapeCellLists.Point,
        cutoff,
    )
    @show naive_out
    painter_out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0.0,
        painter,
        1,
        points[1],
        MultiShapeCellLists.Point,
        cutoff,
    )
    @show painter_out

end