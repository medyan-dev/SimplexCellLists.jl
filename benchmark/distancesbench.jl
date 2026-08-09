using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Random
using Chairmarks

function dist2linelinevect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = dist_sqr(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{2,SVector{3,Float32}}, N)
b = rand(SVector{2,SVector{3,Float32}}, N)
r = zeros(Float32, N)
println("line line")
display(@b dist2linelinevect!($r,$a,$b))
println()


function dist2pointlinevect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = dist_sqr(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{1,SVector{3,Float32}}, N)
b = rand(SVector{2,SVector{3,Float32}}, N)
r = zeros(Float32, N)
println("point line")
display(@b dist2pointlinevect!($r,$a,$b))
println()

function dist2pointpointvect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = dist_sqr(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{1,SVector{3,Float32}}, N)
b = rand(SVector{1,SVector{3,Float32}}, N)
r = zeros(Float32, N)
println("point point")
display(@b dist2pointpointvect!($r,$a,$b))
println()

function dist2pointtrianglevect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = dist_sqr(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{1,SVector{3,Float32}}, N)
b = rand(SVector{3,SVector{3,Float32}}, N)
r = zeros(Float32, N)
println("point triangle")
display(@b dist2pointtrianglevect!($r,$a,$b))
println()

function dist2linetrianglevect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = dist_sqr(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{2,SVector{3,Float32}}, N)
b = rand(SVector{3,SVector{3,Float32}}, N)
r = zeros(Float32, N)
println("line triangle")
display(@b dist2linetrianglevect!($r,$a,$b))
println()

function dist2triangletrianglevect!(r,a,b)
    @inbounds for i in eachindex(r, a, b)
        @inline r[i] = dist_sqr(a[i],b[i])
    end
end

N = 10000
a = rand(SVector{3,SVector{3,Float32}}, N)
b = rand(SVector{3,SVector{3,Float32}}, N)
r = zeros(Float32, N)
println("triangle triangle")
display(@b dist2triangletrianglevect!($r,$a,$b))
println()