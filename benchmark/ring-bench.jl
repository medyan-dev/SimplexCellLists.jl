using BenchmarkTools
using JET
using SimplexCellLists
using StaticArrays
using HDF5
import H5Zblosc

function setup()
    file = h5open(joinpath(@__DIR__, "ringsystemsnapshot.h5"), "r")

    nodepositions = read(file["filaments/1/nodepositions"])
    numcylinders = read(file["filaments/1/num_cylinders"])
    numfilaments = length(numcylinders)

    cylinders = SimplexCellLists.Line[]

    node_ptr = 0
    for filament_i in 1:numfilaments
        node_ptr += 1
        for cylinder_i in 1:numcylinders[filament_i]
            push!(cylinders,[nodepositions[node_ptr,:],nodepositions[node_ptr+1,:]])
            node_ptr += 1
        end
    end
    return cylinders
end


function run1(cylinders)
    output = 0
    f(x,y,i,j,d2,output) = output + 1
    painter = SimplexCellLists.Painter(0,1;
        grid_start= SA[0.0,0.0,0.0],
        grid_size= SA[101,101,11],
        voxel_length= 40.0,
        max_range= SA[Float64[],[226.0]],
    )
    naive = SimplexCellLists.Naive(0, 1)
    # SimplexCellLists.setElements!(naive,[cylinders],[])
    SimplexCellLists.setElements!(painter,[],[cylinders])
    SimplexCellLists.mapPairElements(f, output, painter, 1, SimplexCellLists.Line, 225.0f0)
end

function run1(cylinders)
    f(x,y,i,j,d2,output) = output + 1
    painter = SimplexCellLists.Painter(0,1;
        grid_start= SA[200.0,200.0,200.0],
        grid_size= SA[10,10,1],
        voxel_length= 400.0,
        max_range= SA[Float64[],[226.0]],
    )
    #naive = SimplexCellLists.Naive(0, 1)
    #SimplexCellLists.setElements!(naive,[],[cylinders])
    SimplexCellLists.setElements!(
        painter,
        Vector{SimplexCellLists.Point}[],
        Vector{SimplexCellLists.Line}[cylinders])
    SimplexCellLists.mapPairElements(f, 0, painter, 1, SimplexCellLists.Line, 225.0f0)
    #SimplexCellLists.mapPairElements(f, output, naive, 1, SimplexCellLists.Line, 225.0f0)
end