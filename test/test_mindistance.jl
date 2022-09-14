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
        # d2min_s_t, smin, tmin = SimplexCellLists.dist2LineLine_s_t(a[i],b[i])
        # @test 0 ≤ d2min_s_t
        # @test d2min ≈ d2min_s_t
        # @test 0 ≤ d2min
        @test ref_d2min ≈ d2min
        # @test 0 ≤ tmin ≤ 1
        # s_full = SA[1.0-smin, smin]
        # t_full = SA[1.0-tmin, tmin]
        # r = (sum(a[i] .* s_full) - sum(b[i] .* t_full))
        # @test d2min ≈ r⋅r
    end
end


# @testset "line to triangle" begin
#     N = 10
#     a = rand(SVector{2,SVector{3,Float64}},N)
#     b = rand(SVector{3,SVector{3,Float64}},N)
#     for i in 1:N
#         ref_d2min = refDist2(a[i],b[i])
#         d2min = SimplexCellLists.dist2LineTriangle(a[i],b[i])
#         # d2min_s_t, smin, tmin = SimplexCellLists.dist2LineLine_s_t(a[i],b[i])
#         # @test 0 ≤ d2min_s_t
#         # @test d2min ≈ d2min_s_t
#         # @test 0 ≤ d2min
#         @test ref_d2min ≈ d2min
#         # @test 0 ≤ tmin ≤ 1
#         # s_full = SA[1.0-smin, smin]
#         # t_full = SA[1.0-tmin, tmin]
#         # r = (sum(a[i] .* s_full) - sum(b[i] .* t_full))
#         # @test d2min ≈ r⋅r
#     end
# end