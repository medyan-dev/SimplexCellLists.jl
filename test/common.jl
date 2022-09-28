# Common functions used in tests

using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Random
using Test
using ArgCheck
using Rotations

trials = 1

Random.seed!(1234)

function test_f_mapSimplexElements!(x,y,i,j,d2,output, real_ys, real_cutoff)
    real_d2 = SimplexCellLists.distSqr(x, real_ys[j])
    @argcheck d2 ≥ 0
    @argcheck real_d2 ≥ 0
    @argcheck isapprox(real_d2, d2; atol=1E-5, rtol=1E-3)
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
    @argcheck isapprox(real_d2, d2; atol = 1E-4, rtol=1E-3)
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
    @argcheck isapprox(real_d2, d2; atol = 1E-4, rtol=1E-3)
    @argcheck d2 - real_cutoff^2 ≤ 1E-5
    if √(real_d2) < real_cutoff-1E-5
        output[i,j] += 1
    end
    output
end


function makeBasicCellList!(
        scale,
        rotation,
        translation,
        emptycelllist::SimplexCellLists.SimplexCellList,
    )
    function transform(in)
        pts = reinterpret(SVector{3,Float32},in)
        outpts = map(pts) do pt
            rotation * (pt * scale) + translation 
        end
        copy(reinterpret(eltype(in),outpts))
    end
    points1 = transform([
        SA[SA_F32[-6,0,0]],
        SA[SA_F32[6,0,0]],
    ])
    reinterpret
    points2 = reverse(points1)
    lines1 = transform([
        SA[SA_F32[-6,0,0],SA_F32[-6,2,0]],
        SA[SA_F32[6,0,0],SA_F32[6,2,0]],
    ])
    lines2 = reverse(lines1)
    triangles1 = transform([
        SA[SA_F32[-6,0,0],SA_F32[-6,2,0],SA_F32[-9,1,0]],
        SA[SA_F32[6,0,0],SA_F32[6,2,0],SA_F32[9,1,0]],
    ])
    triangles2 = reverse(triangles1)
    # now try and set some elements
    # some get added in 
    np1 = rand(0:2)
    np2 = rand(0:2)
    nl1 = rand(0:2)
    nl2 = rand(0:2)
    nt1 = rand(0:2)
    nt2 = rand(0:2)
    setElements!(emptycelllist,
        [points1[1:np1],points2[1:np2]],
        [lines1[1:nl1],lines2[1:nl2]],
        [triangles1[1:nt1],triangles2[1:nt2]],
    )
    addElement!.(Ref(emptycelllist),1,points1[np1+1:end])
    addElement!.(Ref(emptycelllist),2,points2[np2+1:end])
    addElement!.(Ref(emptycelllist),1,lines1[nl1+1:end])
    addElement!.(Ref(emptycelllist),2,lines2[nl2+1:end])
    addElement!.(Ref(emptycelllist),1,triangles1[nt1+1:end])
    addElement!.(Ref(emptycelllist),2,triangles2[nt2+1:end])
    points1, points2, lines1, lines2, triangles1, triangles2 
end

