using BenchmarkTools
using JET
using SimplexCellLists
using StaticArrays

function run1(celllist, points)
    cutoff = 0.09f0
    function f(x::SimplexCellLists.Point,y::SimplexCellLists.Point,i::Int32,j::Int32,d2::Float32,output::Float64)::Float64
        output+((√(d2)-0.09f0)^2)
    end
    sum(points) do point
        SimplexCellLists.mapSimplexElements(
            f,
            0.0,
            celllist,
            1,
            point,
            SimplexCellLists.Point,
            cutoff,
        )
    end
end

function run2(celllist, point)
    cutoff = 0.09f0
    f(x,y,i,j,d2,output) = output+((√(d2)-cutoff)^2)
    SimplexCellLists.mapSimplexElements(
        f,
        0.0,
        celllist,
        1,
        point,
        SimplexCellLists.Point,
        cutoff,
    )
end


N = 100000
M = 100
points = rand(SVector{1,SVector{3,Float32}},N)
otherpoints = rand(SVector{1,SVector{3,Float32}},100)
naive = SimplexCellLists.Naive(1, 0)
painter = SimplexCellLists.Painter(1,0;
    grid_start= SA[0.0,0.0,0.0],
    grid_size= SA[10,10,10],
    voxel_length= 1/10,
    max_range= SA[[0.1],Float64[]],
)
SimplexCellLists.setElements!(naive,[points],[])
SimplexCellLists.setElements!(painter,[points],[])
@btime run1(painter, otherpoints)
@btime run2(painter, otherpoints[1])