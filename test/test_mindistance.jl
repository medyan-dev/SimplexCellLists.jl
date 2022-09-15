# Test min distance functions

using JuMP
using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Random
using Test
using Clarabel


function refDist2(points1,points2)
    model = Model(Clarabel.Optimizer)
    set_optimizer_attribute(model, "tol_gap_rel", 1E-12)
    set_optimizer_attribute(model, "tol_gap_abs", 0.0)
    set_silent(model)
    n1 = length(points1)
    n2 = length(points2)
    dims = length(points1[1])
    @assert length(points1[1])==length(points2[1])
    @variable(model, s[i = 1:n1] ≥ 0)
    @variable(model, t[i = 1:n2] ≥ 0)
    @constraint(model, sum(s[i] for i in 1:n1) == 1)
    @constraint(model, sum(t[i] for i in 1:n2) == 1)
    @expression(model, pt1[i = 1:dims], sum(s[j]*points1[j][i] for j in 1:n1))
    @expression(model, pt2[i = 1:dims], sum(t[j]*points2[j][i] for j in 1:n2))
    @expression(model, d[i = 1:dims], pt1[i]-pt2[i])
    @objective(model, Min, sum(d[i]*d[i] for i in 1:dims))
    optimize!(model)
    objective_value(model)
end

Random.seed!(1234)

@testset "point to point" begin
    N = 1000
    a = rand(SVector{1,SVector{3,Float64}},N)
    b = rand(SVector{1,SVector{3,Float64}},N)
    for i in 1:N
        @test refDist2(a[i],b[i]) ≈ SimplexCellLists.dist2PointPoint(a[i],b[i])
    end
end

@testset "point to line" begin
    N = 1000
    a = rand(SVector{1,SVector{3,Float64}},N)
    b = rand(SVector{2,SVector{3,Float64}},N)
    for i in 1:N
        ref_d2min = refDist2(a[i],b[i])
        d2min = SimplexCellLists.dist2PointLine(a[i],b[i])
        d2min_t, tmin = SimplexCellLists.dist2PointLine_t(a[i],b[i])
        @test 0 ≤ d2min_t
        @test d2min ≈ d2min_t
        @test 0 ≤ d2min
        @test ref_d2min ≈ d2min
        @test 0 ≤ tmin ≤ 1
        s_full = SA[1.0]
        t_full = SA[1.0-tmin, tmin]
        r = (sum(a[i] .* s_full) - sum(b[i] .* t_full))
        @test d2min ≈ r⋅r
    end
end

@testset "line to line" begin
    N = 1000
    a = rand(SVector{2,SVector{3,Float64}},N)
    b = rand(SVector{2,SVector{3,Float64}},N)
    for i in 1:N
        ref_d2min = refDist2(a[i],b[i])
        d2min = SimplexCellLists.dist2LineLine(a[i],b[i])
        d2min_s_t, smin, tmin = SimplexCellLists.dist2LineLine_s_t(a[i],b[i])
        @test 0 ≤ d2min_s_t
        @test d2min ≈ d2min_s_t
        @test 0 ≤ d2min
        @test ref_d2min ≈ d2min
        @test 0 ≤ tmin ≤ 1
        s_full = SA[1.0-smin, smin]
        t_full = SA[1.0-tmin, tmin]
        r = (sum(a[i] .* s_full) - sum(b[i] .* t_full))
        @test d2min ≈ r⋅r
    end
end

@testset "point to triangle" begin
    N = 1000
    a = rand(SVector{1,SVector{3,Float64}},N)
    b = rand(SVector{3,SVector{3,Float64}},N)
    for i in 1:N
        ref_d2min = refDist2(a[i],b[i])
        d2min = SimplexCellLists.dist2PointTriangle(a[i],b[i])
        d2min_s_t, tmin = SimplexCellLists.dist2PointTriangle_t(a[i],b[i])
        @test 0 ≤ d2min_s_t
        @test d2min ≈ d2min_s_t
        @test 0 ≤ d2min
        @test ref_d2min ≈ d2min
        @test 0 ≤ tmin[1]
        @test 0 ≤ tmin[2]
        @test sum(tmin) ≤ 1
        s_full = SA[1.0]
        t_full = SA[1.0-sum(tmin), tmin...]
        r = (sum(a[i] .* s_full) - sum(b[i] .* t_full))
        @test d2min ≈ r⋅r
    end
end

@testset "line to triangle edge case" begin
    a = SVector{3, Float64}[
        [0.7962602087717146, 0.7678176789137129, 0.33365738129157096], 
        [0.8787232878884521, 0.22723122900042247, 0.29743853859367475],
    ]
    b = SVector{3, Float64}[
        [0.32244253613423346, 0.07951105154799964, 0.3048568965486481], 
        [0.48609530074630236, 0.14522371134511514, 0.7921987380603223], 
        [0.1937122261241142, 0.9991063549453244, 0.3332458906215262],
    ]
    ref_d2min = refDist2(a,b)
    d2min = SimplexCellLists.dist2LineTriangle(a,b)
    # This test is broken due to an issue with Clarabel.jl
    # See https://github.com/oxfordcontrol/Clarabel.jl/issues/103
    @test_broken ref_d2min ≈ d2min atol = 1E-11 rtol = 1E-11
end
    

@testset "line to triangle" begin
    N = 1000
    a = rand(SVector{2,SVector{3,Float64}},N)
    b = rand(SVector{3,SVector{3,Float64}},N)
    for i in 1:N
        ref_d2min = refDist2(a[i],b[i])
        d2min = SimplexCellLists.dist2LineTriangle(a[i],b[i])
        d2min_s_t, smin, tmin = SimplexCellLists.dist2LineTriangle_s_t(a[i],b[i])
        @test 0 ≤ d2min_s_t
        @test d2min ≈ d2min_s_t atol = 1E-11
        @test 0 ≤ d2min
        if abs(ref_d2min - d2min) > 1E-11
            @show abs(ref_d2min - d2min)
            @show d2min
            @show ref_d2min
            @show a[i]
            @show b[i]
            println()
        end
        @test ref_d2min ≈ d2min atol = 1E-11
        @test 0 ≤ smin[1]
        @test sum(smin) ≤ 1
        @test 0 ≤ tmin[1]
        @test 0 ≤ tmin[2]
        @test sum(tmin) ≤ 1
        s_full = SA[1.0-sum(smin), smin...]
        t_full = SA[1.0-sum(tmin), tmin...]
        r = (sum(a[i] .* s_full) - sum(b[i] .* t_full))
        @test d2min ≈ r⋅r atol = 1E-11
    end
end

@testset "triangle to triangle" begin
    N = 1000
    a = rand(SVector{3,SVector{3,Float64}},N)
    b = rand(SVector{3,SVector{3,Float64}},N)
    for i in 1:N
        ref_d2min = refDist2(a[i],b[i])
        d2min = SimplexCellLists.dist2TriangleTriangle(a[i],b[i])
        d2min_s_t, smin, tmin = SimplexCellLists.dist2TriangleTriangle_s_t(a[i],b[i])
        @test 0 ≤ d2min_s_t
        @test d2min ≈ d2min_s_t atol = 1E-11
        @test 0 ≤ d2min
        @test ref_d2min ≈ d2min atol = 1E-11
        @test 0 ≤ smin[1]
        @test 0 ≤ smin[2]
        @test sum(smin) ≤ 1
        @test 0 ≤ tmin[1]
        @test 0 ≤ tmin[2]
        @test sum(tmin) ≤ 1
        s_full = SA[1.0-sum(smin), smin...]
        t_full = SA[1.0-sum(tmin), tmin...]
        r = (sum(a[i] .* s_full) - sum(b[i] .* t_full))
        @test d2min ≈ r⋅r atol = 1E-11
    end
end