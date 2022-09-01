using BenchmarkTools
using JET
using MultiShapeCellLists
using StaticArrays

function run1(celllist, points)
    cutoff = 0.09f0
    function f(x::MultiShapeCellLists.Point,y::MultiShapeCellLists.Point,i::Int32,j::Int32,d2::Float32,output::Float64)::Float64
        output+((√(d2)-0.09f0)^2)
    end
    sum(points) do point
        MultiShapeCellLists.mapSimplexElements!(
            f,
            0.0,
            celllist,
            1,
            point,
            MultiShapeCellLists.Point,
            cutoff,
        )
    end
end

function run2(celllist, point)
    cutoff = 0.09f0
    f(x,y,i,j,d2,output) = output+((√(d2)-cutoff)^2)
    MultiShapeCellLists.mapSimplexElements!(
        f,
        0.0,
        celllist,
        1,
        point,
        MultiShapeCellLists.Point,
        cutoff,
    )
end


N = 100000
M = 100
points = rand(SVector{1,SVector{3,Float32}},N)
otherpoints = rand(SVector{1,SVector{3,Float32}},100)
naive = MultiShapeCellLists.Naive(1, 0)
painter = MultiShapeCellLists.Painter(1,0;
    grid_start= SA[0.0,0.0,0.0],
    grid_size= SA[10,10,10],
    voxel_length= 1/10,
    max_range= SA[[0.1],Float64[]],
)
MultiShapeCellLists.setElements(naive,[points],[])
MultiShapeCellLists.setElements(painter,[points],[])
@btime run1(painter, otherpoints)
@btime run2(painter, otherpoints[1])