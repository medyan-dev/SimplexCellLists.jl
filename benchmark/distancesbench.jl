using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Random
using BenchmarkTools

function dist2linelinevect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = SimplexCellLists.dist2LineLine(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{2,SVector{3,Float32}}, N)
b = rand(SVector{2,SVector{3,Float32}}, N)
r = zeros(Float32, N)

@btime dist2linelinevect!($r,$a,$b)


function dist2pointlinevect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = SimplexCellLists.dist2PointLine(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{1,SVector{3,Float32}}, N)
b = rand(SVector{2,SVector{3,Float32}}, N)
r = zeros(Float32, N)

@btime dist2pointlinevect!($r,$a,$b)


function dist2pointpointvect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = SimplexCellLists.dist2PointPoint(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{1,SVector{3,Float32}}, N)
b = rand(SVector{1,SVector{3,Float32}}, N)
r = zeros(Float32, N)

@btime dist2pointpointvect!($r,$a,$b)


function dist2pointtrianglevect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = SimplexCellLists.dist2PointTriangle(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{1,SVector{3,Float32}}, N)
b = rand(SVector{3,SVector{3,Float32}}, N)
r = zeros(Float32, N)

@btime dist2pointtrianglevect!($r,$a,$b)