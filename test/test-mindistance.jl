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

function distance_tests(as, bs; check_contact_points=true)
    for i in eachindex(as, bs)
        local a = as[i]
        local b = bs[i]
        local d2, r, s, t = dist_sqr_r_s_t(a, b)
        @test refDist2(a, b) ≈ d2 atol = 1E-10 rtol = 1E-11
        @test d2 === dist_sqr(a, b)
        @test r ⋅ r ≈ d2 atol = 1E-10 rtol = 1E-11
        local s_all = [one(eltype(s)) - sum(s); s;]
        local t_all = [one(eltype(t)) - sum(t); t;]
        if check_contact_points || !iszero(d2)
            @test norm(r - (sum(b .* t_all) - sum(a .* s_all))) ≤ 1E-12
        end
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

# Vertices drawn from a small pool, so simplexes share vertices and are often
# degenerate: zero length line segments, triangles with repeated vertices.
@testset "distance with shared vertices $N to $M" for N in 1:3, M in 1:3
    trials = 300
    pools = [rand(SVector{3,Float64}, 4) for _ in 1:trials]
    as = [SVector{N}(ntuple(_ -> rand(pool), N)) for pool in pools]
    bs = [SVector{M}(ntuple(_ -> rand(pool), M)) for pool in pools]
    # `moller_trumbore_intersect` reports line segments in the plane of a
    # triangle, such as its edges, as intersecting, with arbitrary `s` and `t`.
    distance_tests(as, bs; check_contact_points = min(N, M) == 1 || max(N, M) < 3)
end

@testset "Float32 distance with shared vertices $N to $M" for N in 1:3, M in 1:3
    for trial in 1:3000
        box = 50 + 100rand()
        pool = [SVector{3, Float32}(rand(3) .* box) for _ in 1:6]
        a = SVector{N}(ntuple(_ -> rand(pool), N))
        b = SVector{M}(ntuple(_ -> rand(pool), M))
        ref = sqrt(dist_sqr(map(p -> Float64.(p), a), map(p -> Float64.(p), b)))
        @test abs(sqrt(dist_sqr(a, b)) - ref) ≤ 4*eps(Float32)*box
    end
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

# Simplexes about 100 across, 0.01 apart. Expanding the squared distance, as
# Float64 does, gives a Float32 `d2` that is far off or zero here, and for
# nearly parallel line segments `a*c - b*b` cancels and moves the closest point.
@testset "Float32 large simplexes close together $name" for (name, a, b) in (
        ("point-triangle",
            SA[SA[3f0, 2f0, 0.01f0]],
            SA[SA[-23f0, -7f0, 0f0], SA[31f0, -9f0, 0f0], SA[2f0, 27f0, 0f0]]),
        ("line-line crossing",
            SA[SA[-37f0, 0f0, 0f0], SA[63f0, 0f0, 0f0]],
            SA[SA[-20.5f0, -35.507042f0, 0.01f0], SA[29.5f0, 51.0955f0, 0.01f0]]),
        ("line-line nearly parallel",
            SA[SA[-50f0, 0f0, 0f0], SA[50f0, 0f0, 0f0]],
            SA[SA[-30f0, -0.1f0, 0.01f0], SA[70f0, 0.1f0, 0.01f0]]),
        ("line-line nearly parallel closest point",
            SA[SA[-50f0, 0f0, 0f0], SA[50f0, 0f0, 0f0]],
            SA[SA[-32.51f0, -0.17f0, 0.01f0], SA[67.49f0, 0.17f0, 0.01f0]]),
        ("line-triangle",
            SA[SA[-40f0, 3f0, 0.01f0], SA[60f0, -2f0, 0.02f0]],
            SA[SA[-23f0, -7f0, 0f0], SA[31f0, -9f0, 0f0], SA[2f0, 27f0, 0f0]]),
        ("triangle-triangle",
            SA[SA[-20f0, -5f0, 0.01f0], SA[25f0, -6f0, 0.02f0], SA[0f0, 30f0, 0.015f0]],
            SA[SA[-23f0, -7f0, 0f0], SA[31f0, -9f0, 0f0], SA[2f0, 27f0, 0f0]]),
    )
    for (x, y) in ((a, b), (b, a))
        # Float64 of the same coordinates, which is accurate here.
        ref = sqrt(dist_sqr(map(p -> Float64.(p), x), map(p -> Float64.(p), y)))
        @test 0.0099 < ref < 0.013
        d2, r, s, t = dist_sqr_r_s_t(x, y)
        @test sqrt(d2) ≈ ref rtol = 1E-5
        @test norm(r) ≈ ref rtol = 1E-5
        local s_all = [one(eltype(s)) - sum(s); s;]
        local t_all = [one(eltype(t)) - sum(t); t;]
        @test norm(r - (sum(y .* t_all) - sum(x .* s_all))) ≤ 1E-4
    end
end

# A needle triangle 100 long with a short edge 0.001 long, and a point 0.01
# above its plane that projects 0.2 past the short edge. Testing
# `sbar + tbar ≤ 1` judged it inside in `Float32`, giving a distance of 0.01.
@testset "Float32 point past the short edge of a needle" begin
    Random.seed!(5)
    for trial in 1:3000
        n = normalize(randn(3))
        e1 = normalize(cross(n, randn(3)))
        e2 = cross(n, e1)
        a = 100randn(3)
        b = a .+ 100e1
        c = b .+ 0.001e2
        p = b .+ 0.001rand()*e2 .+ 0.2rand()*e1 .+ 0.01n
        x = SA[SVector{3, Float32}(p)]
        y = SA[SVector{3, Float32}(a), SVector{3, Float32}(b), SVector{3, Float32}(c)]
        ref = sqrt(dist_sqr(map(v -> Float64.(v), x), map(v -> Float64.(v), y)))
        @test abs(sqrt(dist_sqr(x, y)) - ref) ≤ 16*eps(Float32)*200
        @test abs(sqrt(dist_sqr(y, x)) - ref) ≤ 16*eps(Float32)*200
    end
end

# Worst-case accuracy points ported from error-analysis/ledger.jl. Each case
# pins one branch domain (or branch boundary) of a distance kernel at the
# ledger's sharpest known input, with the ledger's first-order claim
# |√d̂² - d| ≤ u·bound(exact sizes), u = eps(T)/2, evaluated at the point and
# rounded up ~1.3x as margin for how the sizes drift under the noise below.
# Each point is tested as given, then with Gaussian coordinate noise of a few
# ulps: enough to scramble every rounding, and at the boundary points to make
# the branch actually flip (their tolerance is the traced all-domain bound,
# which includes the branch jump terms).
const ledger_cases = [
    # point-point: 5u/2·d, attained when one coordinate carries the separation.
    ("point-point, one coordinate separation", 7.0,
        ((-1.0, 0.03, -0.02),), ((1.0, 0.06, 0.03),)),
    ("point-point, balanced coordinates", 11.0,
        ((-0.95, -1.0, 0.9),), ((0.9, 0.85, -1.0),)),
    # point-line: t clamped to 0 is exactly point-point (5u/2·d); the interior
    # arm has the traced bound u(9‖p‖+3‖r‖+3/2·d); t clamped to 1 computes
    # fl(p+r), u(‖p‖+‖r‖+5/2·d), tight when p+r cancels; the t* = 0 and t* = 1
    # boundary points have exactly zero conditions built from cancelling
    # order-1 products, so rounding picks either arm.
    ("point-line, t clamped to 0", 7.0,
        ((1.0, 0.12, -0.03),), ((-1.0, 0.0, 0.0), (-1.0, -1.0, 0.0))),
    ("point-line, interior", 38.0,
        ((1.0, 0.0, 0.05),), ((-1.0, 1.0, 0.0), (-1.0, -1.0, 0.0))),
    ("point-line, interior, aspect 10", 29.0,
        ((1.0, 0.0, 0.0),), ((-1.0, 0.1, 0.02), (-1.0, -0.1, 0.02))),
    ("point-line, t clamped to 1, p + r cancels", 8.0,
        ((1.0, 1.0, 0.04),), ((-1.0, -1.0, 0.0), (0.98, 0.98, 0.0))),
    ("point-line, boundary t* = 0", 35.0,
        ((1.0, 1.0, 0.0),), ((-1.0, 0.0, 0.0), (-0.5, -1.0, 0.0))),
    ("point-line, boundary t* = 1", 26.0,
        ((0.5, 0.5, 0.0),), ((-1.0, 1.0, 0.0), (0.0, 0.0, 0.0))),
    # point-triangle: the face arm claim is u(6‖BP‖+‖E0‖+‖E1‖+8d); the tilted
    # triangle makes N generic, so the terms of hN = dot3(BP, N) are order 0.5
    # and cancel — with a z = 0 triangle (as in the testsets above) hN is a
    # single product and this error path is never exercised. Edges and
    # vertices reuse the point-line rows on (BP, E) with rounded-sum charges.
    ("point-triangle, face of a tilted triangle", 17.0,
        ((-0.01, 0.01, 0.29),), ((-0.9, -0.6, 0.2), (0.8, -0.2, -0.3), (-0.2, 0.9, 0.1))),
    ("point-triangle, near face of a tilted triangle, hN cancels", 13.0,
        ((-0.09, 0.03, 0.02),), ((-0.9, -0.6, 0.2), (0.8, -0.2, -0.3), (-0.2, 0.9, 0.1))),
    ("point-triangle, edge 0 interior", 23.0,
        ((0.2, -0.9, 0.4),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    ("point-triangle, edge 2 interior", 24.0,
        ((0.56, 0.59, 0.3),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    ("point-triangle, vertex a", 1.7,
        ((-1.0, -0.95, 0.35),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    ("point-triangle, vertex b", 6.5,
        ((1.0, -0.45, 0.3),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    # x over the midpoint of edge 0: the tbar numerator of the inside test is
    # exactly zero, so rounding decides face vs edge.
    ("point-triangle, inside-test boundary", 164.0,
        ((-0.05, -0.4, 0.6),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    # x exactly on the u0 = 1 plane through b: edge 0 vs vertex b flips.
    ("point-triangle, clamp boundary at vertex b", 192.0,
        ((0.88, -0.54, 0.3),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    # line-line: the crossing arm claim is u(6‖W‖+‖P‖+‖Q‖+8d) with no 1/sinθ;
    # the near parallel point (sinθ ≈ 0.008, where the Gram form a*c - b*b
    # loses ~125x) must meet the same bound as the generic crossing. Edge and
    # corner rows reuse the point-line rows on (W, P), (W, Q), etc.
    ("line-line, generic crossing", 17.0,
        ((-1.0, -0.2, 0.0), (1.0, 0.2, 0.0)), ((-0.3, -1.0, 0.25), (0.3, 1.0, 0.25))),
    ("line-line, nearly parallel crossing", 8.0,
        ((-1.0, -0.004, 0.0), (1.0, 0.004, 0.0)), ((-1.0, 0.004, 0.15), (1.0, -0.004, 0.15))),
    ("line-line, s interior, t = 0", 25.0,
        ((-1.0, -0.2, 0.0), (1.0, 0.2, 0.0)), ((0.1, 0.3, 0.3), (0.9, 0.9, 0.5))),
    ("line-line, s = 0, t interior", 10.5,
        ((-0.4, -0.3, 0.0), (0.56, 0.42, 0.0)), ((-0.31, -0.67, 0.25), (-0.73, -0.11, 0.25))),
    ("line-line, corner (0,0)", 2.2,
        ((-0.4, -0.3, 0.0), (0.56, 0.42, 0.0)), ((-0.7, -0.8, 0.3), (-1.0, -0.9, 0.5))),
    # Q crosses exactly over P0: the sbar numerator of the inside test is zero.
    ("line-line, inside-test boundary", 66.0,
        ((-0.8, -0.6, 0.0), (0.8, 0.6, 0.0)), ((-0.68, -0.76, 0.25), (-0.92, -0.44, 0.25))),
    # (Q0 - P0) ⊥ Q exactly: the t* = 0 clamp boundary in the s* = 0 region.
    ("line-line, clamp boundary at corner (0,0)", 53.0,
        ((-0.4, -0.3, 0.0), (0.56, 0.42, 0.0)), ((-0.52, -0.39, 0.25), (-0.82, 0.01, 0.25))),
]

@testset "ledger worst-case point: $name" for (name, tolu, xt, yt) in ledger_cases
    rng = Random.Xoshiro(2026)
    for T in (Float32, Float64)
        tol = tolu*eps(T)/2
        for σ in (0, 8, 64, 1024) .* eps(T), _ in 1:(iszero(σ) ? 1 : 100)
            x = SVector(map(v -> SVector{3, T}(v .+ σ .* randn(rng, 3)), xt))
            y = SVector(map(v -> SVector{3, T}(v .+ σ .* randn(rng, 3)), yt))
            ref = sqrt(dist_sqr(map(v -> big.(v), x), map(v -> big.(v), y)))
            @test abs(sqrt(dist_sqr(x, y)) - ref) ≤ tol
            @test abs(sqrt(dist_sqr(y, x)) - ref) ≤ tol
        end
    end
end

# Return `W` simplices as one simplex with a lane per simplex.
function simplex_lanes(xs::NTuple{W, SVector{N, SVector{3, T}}}) where {W, N, T}
    SVector{N}(ntuple(Val(N)) do v
        SVector{3}(ntuple(Val(3)) do k
            SimplexCellLists.SIMD.Vec{W, T}(ntuple(l -> xs[l][v][k], Val(W)))
        end)
    end)
end

# Return a random simplex, half the time with coordinates in -1:1, so
# parallel, collinear, coincident, and zero length cases are common.
function rand_simplex(::Type{T}, N) where T
    if rand(Bool)
        SVector{N}(ntuple(_ -> SVector{3, T}(rand(-1:1, 3)), N))
    else
        rand(SVector{N, SVector{3, T}})
    end
end

# The grid neighbor list checks distances several pairs at a time, which
# must give bit identical results to the scalar version.
@testset "dist_sqr on SIMD lanes $T $N to $M" for T in (Float32, Float64),
        (N, M) in ((1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 2))
    W = 8
    for trial in 1:500
        as = ntuple(_ -> rand_simplex(T, N), Val(W))
        bs = ntuple(_ -> rand_simplex(T, M), Val(W))
        d2 = dist_sqr(simplex_lanes(as), simplex_lanes(bs))
        @test all(l -> d2[l] === dist_sqr(as[l], bs[l]), 1:W)
    end
end
nothing
