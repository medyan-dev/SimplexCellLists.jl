# Common functions used in tests

using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Random
using Test
using ArgCheck

Random.seed!(1234)

function test_f_mapSimplexElements!(x,y,i,j,d2,output, real_ys, real_cutoff)
    real_d2 = SimplexCellLists.distSqr(x, real_ys[j])
    @argcheck d2 ≥ 0
    @argcheck real_d2 ≥ 0
    @argcheck real_d2 ≈ d2
    @argcheck d2 - real_cutoff^2 ≤ 1E-5
    if √(d2) < real_cutoff-1E-5
        output[j] += 1
    end
    output
end

function test_f_mapPairElements!(x,y,i,j,d2,output, real_xs, real_cutoff)
    real_d2 = SimplexCellLists.distSqr(real_xs[i], real_xs[j])
    @argcheck d2 ≥ 0
    @argcheck real_d2 ≥ 0
    @argcheck real_d2 ≈ d2
    @argcheck d2 - real_cutoff^2 ≤ 1E-5
    if √(real_d2) < real_cutoff-1E-5
        output[i,j] += 1
        output[j,i] += 1
    end
    output
end

function test_f_mapElementsElements!(x,y,i,j,d2,output, real_xs, real_ys, real_cutoff)
    real_d2 = SimplexCellLists.distSqr(real_xs[i], real_ys[j])
    @argcheck d2 ≥ 0
    @argcheck real_d2 ≥ 0
    @argcheck real_d2 ≈ d2
    @argcheck d2 - real_cutoff^2 ≤ 1E-5
    if √(real_d2) < real_cutoff-1E-5
        output[i,j] += 1
    end
    output
end


