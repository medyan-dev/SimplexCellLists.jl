# Test min distance functions

using JuMP
using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Random
using Test
using Clarabel
using SimplexCellLists: dist_sqr, dist_sqr_r_s_t


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

function distance_tests(as, bs)
    for i in eachindex(as, bs)
        local a = as[i]
        local b = bs[i]
        local d2, r, s, t = dist_sqr_r_s_t(a, b)
        @test refDist2(a, b) ≈ d2 atol = 1E-10 rtol = 1E-11
        @test d2 ≈ dist_sqr(a, b) atol = 1E-10 rtol = 1E-11
        @test r ⋅ r ≈ d2 atol = 1E-10 rtol = 1E-11
        local s_all = [one(eltype(s)) - sum(s); s;]
        local t_all = [one(eltype(t)) - sum(t); t;]
        @test norm(r - (sum(b .* t_all) - sum(a .* s_all))) ≤ 1E-12
        @test all(≥(0), s)
        @test all(≥(0), t)
        @test sum(s) ≤ 1
        @test sum(t) ≤ 1
    end
end

@testset "distance $N to $M" for N in 1:3, M in 1:3
    trials = 1000
    as = rand(SVector{N,SVector{3,Float64}},trials)
    bs = rand(SVector{M,SVector{3,Float64}},trials)
    distance_tests(as, bs)
end

@testset "line to triangle edge case" begin
    a = SA[
        SA[0.7962602087717146, 0.7678176789137129, 0.33365738129157096], 
        SA[0.8787232878884521, 0.22723122900042247, 0.29743853859367475],
    ]
    b = SA[
        SA[0.32244253613423346, 0.07951105154799964, 0.3048568965486481], 
        SA[0.48609530074630236, 0.14522371134511514, 0.7921987380603223], 
        SA[0.1937122261241142, 0.9991063549453244, 0.3332458906215262],
    ]
    ref_d2min = refDist2(a,b)
    d2min = SimplexCellLists.dist_sqr(a,b)
    # This test was broken due to an issue with Clarabel.jl
    # See https://github.com/oxfordcontrol/Clarabel.jl/issues/103
    # This issue was resolved in https://github.com/oxfordcontrol/Clarabel.jl/pull/114
    @test ref_d2min ≈ d2min atol = 1E-11 rtol = 1E-11
end
nothing
