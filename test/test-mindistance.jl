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
        @test d2 === dist_sqr(a, b)
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

# Vertices drawn from a small pool, so simplexes share vertices and are often
# degenerate: zero length line segments, triangles with repeated vertices.
@testset "distance with shared vertices $N to $M" for N in 1:3, M in 1:3
    trials = 300
    pools = [rand(SVector{3,Float64}, 4) for _ in 1:trials]
    as = [SVector{N}(ntuple(_ -> rand(pool), N)) for pool in pools]
    bs = [SVector{M}(ntuple(_ -> rand(pool), M)) for pool in pools]
    distance_tests(as, bs)
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

# Line segments on nearly the same line, or nearly parallel. For collinear
# segments to within rounding, `n = P×Q` and `Q0mP0×Q` were both rounding
# errors, so `sbar` and `tbar` were arbitrary, the inside test could pass,
# and the distance between disjoint segments 0.56 apart came out as 3e-8.
# Nearly parallel overlapping segments got the distance right, but with
# weights giving points far apart. Integer coordinates miss this, as an
# exactly collinear `n` is zero.
@testset "$T nearly collinear and nearly parallel line segments" for T in (Float32, Float64)
    rng = Random.Xoshiro(29)
    # The case 0.56 apart, kept since the random ones below fail it only
    # about once in 700 for the old method.
    cases = Any[(
        SA[SA[-1.1168967f0, -0.4735017f0, -0.30628234f0], SA[-0.43016562f0, 0.06529092f0, -0.054257788f0]],
        SA[SA[0.96934664f0, 1.1633145f0, 0.45935148f0], SA[-0.007197912f0, 0.39714116f0, 0.10096781f0]],
    )]
    for trial in 1:10000
        q = normalize(randn(rng, 3))
        e1 = normalize(cross(q, randn(rng, 3)))
        e2 = cross(q, e1)
        o = randn(rng, 3)
        if isodd(trial)
            # Four points on a line, with offsets below a rounding or a few.
            ts = randn(rng, 4)
            pts = [o .+ t .* q .+ rand(rng, (0, 0.1, 4)) * eps(T) .* (randn(rng)*e1 .+ randn(rng)*e2) for t in ts]
            xt, yt = (pts[1], pts[2]), (pts[3], pts[4])
        else
            # Nearly parallel, offset along the line and across it.
            θ = 10.0^(-rand(rng)*12)
            dq = cos(θ)*q .+ sin(θ)*e1
            c = o .+ 2randn(rng)*q .+ 10.0^(-rand(rng)*6)*randn(rng)*e2
            xt = (o .- q, o .+ q)
            yt = (c .- exp(randn(rng)) .* dq, c .+ exp(randn(rng)) .* dq)
        end
        push!(cases, (SVector{2}(map(v -> SVector{3, T}(v), xt)), SVector{2}(map(v -> SVector{3, T}(v), yt))))
    end
    for (x, y) in cases
        x = map(v -> SVector{3, T}(v), x)
        y = map(v -> SVector{3, T}(v), y)
        tol = 8*eps(T)*maximum(v -> maximum(abs, v), (x..., y...))
        for (a, b) in ((x, y), (y, x))
            ref = sqrt(dist_sqr(map(v -> big.(v), a), map(v -> big.(v), b)))
            d2, r, s, t = dist_sqr_r_s_t(a, b)
            @test abs(sqrt(d2) - ref) ≤ tol
            # `r` is the separation of the points the weights give.
            pa = (1 - big(s[1])) .* a[1] .+ big(s[1]) .* a[2]
            pb = (1 - big(t[1])) .* b[1] .+ big(t[1]) .* b[2]
            @test norm(r .- (pb .- pa)) ≤ tol
            @test abs(sqrt(d2) - norm(r)) ≤ tol
        end
    end
end

# Line segments in, or nearly in, the plane of a triangle beside it, and line
# segments through the face of a sliver, and triangles with those as an edge.
# Testing whether the line segment crosses the triangle, from the signs of
# triple products that are rounding errors in the plane, gave a distance of
# zero for line segments 1.7 from the triangle, about one time in 300. Through
# a sliver, `N = E0×E1` from `cross3` tilted the plane enough for errors near
# `√eps*ℓ`.
@testset "$T line segments and triangles near the plane of a triangle" for T in (Float32, Float64)
    rng = Random.Xoshiro(31)
    # The case 1.7 apart.
    cases = Any[(
        SA[SA[0.27440169104786705, -1.8352000938383282, -1.1414620437379672], SA[-0.4237474034660068, -1.74076257090927, -0.005187092841748675]],
        SA[SA[-0.051006641386534685, 0.8264827967799175, 0.247674370931405], SA[0.3840271906908225, 0.44295798307158735, -0.5669754202169636], SA[-0.3934509020879434, -0.08075321479510644, 0.491917584990893]],
    )]
    for trial in 1:3000
        n = normalize(randn(rng, 3))
        e1 = normalize(cross(n, randn(rng, 3)))
        e2 = cross(n, e1)
        o = randn(rng, 3)
        P(a, b, c = 0) = o .+ a .* e1 .+ b .* e2 .+ c .* n
        θ = rand(rng, (0.0, 10.0^(-rand(rng)*12)))
        if isodd(trial)
            # A unit triangle, and a line segment beside it, tilted out of the
            # plane by about θ.
            yt = (P(0, 0), P(1, 0), P(0, 1))
            xt = (P(2 + rand(rng), -1 + 3rand(rng), -θ), P(-1 + 3rand(rng), 2 + rand(rng), θ))
        else
            # A sliver, and a line segment through its face.
            yt = (P(0, 0), P(1, 0), P(0.5 + rand(rng), 10.0^(-rand(rng)*12)))
            w = rand(rng, 3)
            X = sum(w ./ sum(w) .* collect(yt))
            d = normalize(randn(rng, 3))
            a = rand(rng)
            xt = (X .- a .* d, X .+ (1 - a) .* d)
        end
        x = SVector{2}(map(v -> SVector{3, T}(v), xt))
        y = SVector{3}(map(v -> SVector{3, T}(v), yt))
        # A triangle with `x` as an edge, its third vertex 1 or θ away.
        c = SVector{3, T}(xt[1] .+ ifelse(rand(rng, Bool), 1, θ) .* normalize(randn(rng, 3)))
        push!(cases, (x, y), (SA[x[1], x[2], c], y))
    end
    for (x, y) in cases
        x = map(v -> SVector{3, T}(v), x)
        y = map(v -> SVector{3, T}(v), y)
        tol = 8*eps(T)*maximum(v -> maximum(abs, v), (x..., y...))
        for (a, b) in ((x, y), (y, x))
            ref = sqrt(dist_sqr(map(v -> big.(v), a), map(v -> big.(v), b)))
            d2, r, s, t = dist_sqr_r_s_t(a, b)
            @test abs(sqrt(d2) - ref) ≤ tol
            @test all(≥(0), s) && all(≥(0), t) && sum(s) ≤ 1 && sum(t) ≤ 1
            # `r` is the separation of the points the weights give.
            pa = (1 - sum(big.(s))) .* a[1] .+ sum(big.(s) .* a[2:end])
            pb = (1 - sum(big.(t))) .* b[1] .+ sum(big.(t) .* b[2:end])
            @test norm(r .- (pb .- pa)) ≤ tol
        end
    end
end

# A point on the face of a triangle, within 1e-69 of it, with coordinates
# near 1e7. The distance came out as 9.96, which is 8.07*eps*ℓ, just over 8.
@testset "Float32 point on the face of a triangle with large coordinates" begin
    A = Float32.(SA[-10351314, -10323355, -10169968])
    B = Float32.(SA[ 10351314,   6399932,   9067019])
    C = Float32.(SA[  9871146,  10323356,  10169969])
    P = Float32.(SA[  9718146,  10124128,  10000623])
    x, y = SA[P], SA[A, B, C]
    tol = 64*eps(Float32)*maximum(v -> maximum(abs, v), (x..., y...))
    ref = sqrt(dist_sqr(map(v -> big.(v), x), map(v -> big.(v), y)))
    for (a, b) in ((x, y), (y, x))
        d2, r, s, t = dist_sqr_r_s_t(a, b)
        @test abs(sqrt(d2) - ref) ≤ tol
        @test all(≥(0), s) && all(≥(0), t) && sum(s) ≤ 1 && sum(t) ≤ 1
        # `r` is the separation of the points the weights give.
        pa = (1 - sum(big.(s))) .* a[1] .+ sum(big.(s) .* a[2:end]; init = zero(SVector{3, BigFloat}))
        pb = (1 - sum(big.(t))) .* b[1] .+ sum(big.(t) .* b[2:end]; init = zero(SVector{3, BigFloat}))
        @test norm(r .- (pb .- pa)) ≤ tol
    end
end

# The end point of a long line segment 0.0135 from a short one, with the long
# one across the range. With the short one first, `t` on the long one from
# the lines through them, which is 1 - 2e-8, rounded to 3 eps short of one,
# 9.7 along it. The distance came out as 9.43, which is 8.12*eps*ℓ.
@testset "Float32 end point of a long line segment near a short one" begin
    x = SA[SA[1.9877764f6, 9.505106f6, 9.729709f6], SA[-9.696928f6, -5.704686f6, -9.494322f6]]
    y = SA[SA[-9.676423f6, -5.6547575f6, -9.56445f6], SA[-9.731833f6, -5.7896785f6, -9.374942f6]]
    tol = 64*eps(Float32)*maximum(v -> maximum(abs, v), (x..., y...))
    ref = sqrt(dist_sqr(map(v -> big.(v), x), map(v -> big.(v), y)))
    for (a, b) in ((x, y), (y, x))
        d2, r, s, t = dist_sqr_r_s_t(a, b)
        @test abs(sqrt(d2) - ref) ≤ tol
        @test 0 ≤ s[1] ≤ 1 && 0 ≤ t[1] ≤ 1
        # `r` is the separation of the points the weights give.
        pa = (1 - big(s[1])) .* a[1] .+ big(s[1]) .* a[2]
        pb = (1 - big(t[1])) .* b[1] .+ big(t[1]) .* b[2]
        @test norm(r .- (pb .- pa)) ≤ tol
    end
end

# Hard points: one coordinate carrying the separation, separations that
# cancel, nearly parallel line segments, and points on the border between
# two closest features, where rounding picks either. Each is tested as given,
# then with Gaussian noise of a few ulps on the coordinates, which scrambles
# the roundings and moves the border points to either side.
const hard_points = [
    ("point-point, one coordinate separation",
        ((-1.0, 0.03, -0.02),), ((1.0, 0.06, 0.03),)),
    ("point-point, balanced coordinates",
        ((-0.95, -1.0, 0.9),), ((0.9, 0.85, -1.0),)),
    ("point-line, closest at the first end point",
        ((1.0, 0.12, -0.03),), ((-1.0, 0.0, 0.0), (-1.0, -1.0, 0.0))),
    ("point-line, closest inside",
        ((1.0, 0.0, 0.05),), ((-1.0, 1.0, 0.0), (-1.0, -1.0, 0.0))),
    ("point-line, closest inside a short line segment",
        ((1.0, 0.0, 0.0),), ((-1.0, 0.1, 0.02), (-1.0, -0.1, 0.02))),
    # The second end point is near the point, so `y[2] - x` cancels.
    ("point-line, closest at the second end point",
        ((1.0, 1.0, 0.04),), ((-1.0, -1.0, 0.0), (0.98, 0.98, 0.0))),
    ("point-line, border at the first end point",
        ((1.0, 1.0, 0.0),), ((-1.0, 0.0, 0.0), (-0.5, -1.0, 0.0))),
    ("point-line, border at the second end point",
        ((0.5, 0.5, 0.0),), ((-1.0, 1.0, 0.0), (0.0, 0.0, 0.0))),
    # A tilted triangle, so no product in `N` or `N⋅(X - A)` is exact.
    ("point-triangle, over the face of a tilted triangle",
        ((-0.01, 0.01, 0.29),), ((-0.9, -0.6, 0.2), (0.8, -0.2, -0.3), (-0.2, 0.9, 0.1))),
    ("point-triangle, near the face of a tilted triangle",
        ((-0.09, 0.03, 0.02),), ((-0.9, -0.6, 0.2), (0.8, -0.2, -0.3), (-0.2, 0.9, 0.1))),
    ("point-triangle, closest inside edge AB",
        ((0.2, -0.9, 0.4),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    ("point-triangle, closest inside edge BC",
        ((0.56, 0.59, 0.3),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    ("point-triangle, closest at vertex A",
        ((-1.0, -0.95, 0.35),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    ("point-triangle, closest at vertex B",
        ((1.0, -0.45, 0.3),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    ("point-triangle, border of the face and edge AB, over its midpoint",
        ((-0.05, -0.4, 0.6),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    ("point-triangle, border of edge AB and vertex B",
        ((0.88, -0.54, 0.3),), ((-0.9, -0.6, 0.0), (0.8, -0.2, 0.0), (-0.2, 0.9, 0.0))),
    ("line-line, crossing",
        ((-1.0, -0.2, 0.0), (1.0, 0.2, 0.0)), ((-0.3, -1.0, 0.25), (0.3, 1.0, 0.25))),
    ("line-line, nearly parallel crossing",
        ((-1.0, -0.004, 0.0), (1.0, 0.004, 0.0)), ((-1.0, 0.004, 0.15), (1.0, -0.004, 0.15))),
    ("line-line, inside x and at the first end point of y",
        ((-1.0, -0.2, 0.0), (1.0, 0.2, 0.0)), ((0.1, 0.3, 0.3), (0.9, 0.9, 0.5))),
    ("line-line, at the first end point of x and inside y",
        ((-0.4, -0.3, 0.0), (0.56, 0.42, 0.0)), ((-0.31, -0.67, 0.25), (-0.73, -0.11, 0.25))),
    ("line-line, at the first end points",
        ((-0.4, -0.3, 0.0), (0.56, 0.42, 0.0)), ((-0.7, -0.8, 0.3), (-1.0, -0.9, 0.5))),
    ("line-line, border with y crossing over the first end point of x",
        ((-0.8, -0.6, 0.0), (0.8, 0.6, 0.0)), ((-0.68, -0.76, 0.25), (-0.92, -0.44, 0.25))),
    ("line-line, border of the first end points and inside y",
        ((-0.4, -0.3, 0.0), (0.56, 0.42, 0.0)), ((-0.52, -0.39, 0.25), (-0.82, 0.01, 0.25))),
]

@testset "hard point: $name" for (name, xt, yt) in hard_points
    rng = Random.Xoshiro(2026)
    for T in (Float32, Float64)
        for σ in (0, 8, 64, 1024) .* eps(T), _ in 1:(iszero(σ) ? 1 : 100)
            x = SVector(map(v -> SVector{3, T}(v .+ σ .* randn(rng, 3)), xt))
            y = SVector(map(v -> SVector{3, T}(v .+ σ .* randn(rng, 3)), yt))
            tol = 8*eps(T)*maximum(v -> maximum(abs, v), (x..., y...))
            ref = sqrt(dist_sqr(map(v -> big.(v), x), map(v -> big.(v), y)))
            @test abs(sqrt(dist_sqr(x, y)) - ref) ≤ tol
            @test abs(sqrt(dist_sqr(y, x)) - ref) ≤ tol
        end
    end
end

# Small or thin simplexes, where a squared length, or a squared area like
# `Δ = |n|²`, is out of range. When something underflows the distance may
# lose accuracy, but it must stay finite and off by at most about
# `floatmin^(1/4)`, which bounds the inradius of a triangle with `Δ` below
# `floatmin`. The point-triangle method used to form `inv(Δ)`, and returned
# `Inf` for a point over the face when it overflowed, or `NaN` in the plane,
# even for a unit sliver. When `Δ` underflows the line-line `t` is arbitrary,
# and when `‖r‖²` is subnormal the point-line quotient for `t` keeps only a
# few bits, but both are within the floor. Nearly parallel unit line
# segments do not reach this: their end points are within about ℓ·sinθ of
# the crossing.
const tiny_normal_cases = [
    ("point-triangle, tiny, inv(Δ) overflows", Float32,
        ((0.25e-10, 0.25e-10, 1e-10),), ((0, 0, 0), (1e-10, 0, 0), (0, 1e-10, 0))),
    ("point-triangle, tiny, Δ underflows", Float32,
        ((0.25e-12, 0.25e-12, 1e-12),), ((0, 0, 0), (1e-12, 0, 0), (0, 1e-12, 0))),
    ("point-triangle, unit sliver, inv(Δ) overflows", Float32,
        ((1.5, 0.6e-20, 0.1),), ((0, 0, 0), (1, 0, 0), (2, 1e-20, 0))),
    ("point-triangle, unit sliver, in plane", Float32,
        ((1.5, 0.6e-20, 0),), ((0, 0, 0), (1, 0, 0), (2, 1e-20, 0))),
    ("line-line, tiny crossing, inv(Δ) overflows", Float32,
        ((-0.5e-10, 0, 0), (0.5e-10, 0, 0)), ((0, -0.5e-10, 1e-10), (0, 0.5e-10, 1e-10))),
    ("line-line, tiny crossing, Δ underflows", Float32,
        ((-0.5e-12, 0, 0), (0.5e-12, 0, 0)), ((0, -0.5e-12, 1e-12), (0, 0.5e-12, 1e-12))),
    ("point-triangle, tiny, inv(Δ) overflows", Float64,
        ((0.25e-80, 0.25e-80, 1e-80),), ((0, 0, 0), (1e-80, 0, 0), (0, 1e-80, 0))),
    ("point-triangle, tiny, Δ underflows", Float64,
        ((0.25e-90, 0.25e-90, 1e-90),), ((0, 0, 0), (1e-90, 0, 0), (0, 1e-90, 0))),
    ("point-triangle, unit sliver, inv(Δ) overflows", Float64,
        ((1.5, 0.6e-160, 0.1),), ((0, 0, 0), (1, 0, 0), (2, 1e-160, 0))),
    ("point-triangle, unit sliver, in plane", Float64,
        ((1.5, 0.6e-160, 0),), ((0, 0, 0), (1, 0, 0), (2, 1e-160, 0))),
    ("line-line, tiny crossing, inv(Δ) overflows", Float64,
        ((-0.5e-80, 0, 0), (0.5e-80, 0, 0)), ((0, -0.5e-80, 1e-80), (0, 0.5e-80, 1e-80))),
    ("line-line, tiny crossing, Δ underflows", Float64,
        ((-0.5e-90, 0, 0), (0.5e-90, 0, 0)), ((0, -0.5e-90, 1e-90), (0, 0.5e-90, 1e-90))),
    ("point-line, tiny, subnormal norm", Float32,
        ((0.5e-22, 1e-22, 0),), ((0, 0, 0), (1e-22, 0, 0))),
    ("point-line, tiny, subnormal norm", Float64,
        ((0.5e-160, 1e-160, 0),), ((0, 0, 0), (1e-160, 0, 0))),
    ("line-line, tiny, inv of a norm overflows", Float32,
        ((0, 0, 0), (1e-21, 0, 0)), ((0.5e-21, 1e-21, 0), (0.5e-21, 2e-21, 0))),
    ("line-line, tiny, inv of a norm overflows", Float64,
        ((0, 0, 0), (1e-155, 0, 0)), ((0.5e-155, 1e-155, 0), (0.5e-155, 2e-155, 0))),
]

@testset "tiny normal: $name $T" for (name, T, xt, yt) in tiny_normal_cases
    x = SVector(map(v -> SVector{3, T}(v), xt))
    y = SVector(map(v -> SVector{3, T}(v), yt))
    ref = sqrt(dist_sqr(map(v -> big.(v), x), map(v -> big.(v), y)))
    # A few roundings of the largest coordinate, plus the underflow floor.
    tol = 8*eps(T)*maximum(v -> maximum(abs, v), (x..., y...)) + floatmin(T)^(1/4)
    @test abs(sqrt(dist_sqr(x, y)) - ref) ≤ tol
    @test abs(sqrt(dist_sqr(y, x)) - ref) ≤ tol
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
        (N, M) in ((1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 2), (2, 3), (3, 2), (3, 3))
    W = 8
    for trial in 1:500
        as = ntuple(_ -> rand_simplex(T, N), Val(W))
        bs = ntuple(_ -> rand_simplex(T, M), Val(W))
        d2 = dist_sqr(simplex_lanes(as), simplex_lanes(bs))
        @test all(l -> d2[l] === dist_sqr(as[l], bs[l]), 1:W)
    end
end

# Coordinates just below `floatmax(T)^(1/4)/4`, the top of the range, where
# `N⋅N` and the triple products are near `floatmax`: the hard points scaled
# up, and `rand_simplex` scaled up, half of them with coordinates in
# `(-L, 0, L)`, where the edges are longest, and repeated vertices, parallel
# edges, and coplanar triangles are common.
@testset "$T near the top of the range" for T in (Float32, Float64)
    L = prevfloat(floatmax(T)^T(1/4)/4)
    cases = Any[]
    for (name, xt, yt) in hard_points
        push!(cases, (SVector(map(v -> SVector{3, T}(L .* v), xt)), SVector(map(v -> SVector{3, T}(L .* v), yt))))
    end
    for N in 1:3, M in 1:3, _ in 1:300
        push!(cases, (rand_simplex(T, N) .* L, rand_simplex(T, M) .* L))
    end
    for (x, y) in cases, (a, b) in ((x, y), (y, x))
        d2, r, s, t = dist_sqr_r_s_t(a, b)
        @test isfinite(d2) && all(isfinite, r)
        @test all(≥(0), s) && all(≥(0), t) && sum(s) ≤ 1 && sum(t) ≤ 1
        tol = 8*eps(T)*maximum(v -> maximum(abs, v), (a..., b...))
        ref = sqrt(dist_sqr(map(v -> big.(v), a), map(v -> big.(v), b)))
        @test abs(sqrt(d2) - ref) ≤ tol
        pa = (1 - sum(big.(s))) .* a[1] .+ sum(big.(s) .* a[2:end]; init = zero(SVector{3, BigFloat}))
        pb = (1 - sum(big.(t))) .* b[1] .+ sum(big.(t) .* b[2:end]; init = zero(SVector{3, BigFloat}))
        @test norm(r .- (pb .- pa)) ≤ tol
    end
end

# A NaN coordinate anywhere gives a NaN `d²`, whatever the candidate it would
# lose to, also in SIMD lanes next to lanes without NaN. The weights stay in
# range.
@testset "$T NaN coordinate $N to $M" for T in (Float32, Float64), N in 1:3, M in 1:3
    for _ in 1:20, i in 1:(N + M), k in 1:3
        x = rand_simplex(T, N)
        y = rand_simplex(T, M)
        if i ≤ N
            x = setindex(x, setindex(x[i], T(NaN), k), i)
        else
            y = setindex(y, setindex(y[i - N], T(NaN), k), i - N)
        end
        for (a, b) in ((x, y), (y, x))
            d2, r, s, t = dist_sqr_r_s_t(a, b)
            @test isnan(d2)
            @test all(w -> 0 ≤ w ≤ 1, s) && all(w -> 0 ≤ w ≤ 1, t)
            as = (a, ntuple(_ -> rand_simplex(T, length(a)), Val(7))...)
            bs = (b, ntuple(_ -> rand_simplex(T, length(b)), Val(7))...)
            d2s = dist_sqr(simplex_lanes(as), simplex_lanes(bs))
            @test isnan(d2s[1]) && all(l -> d2s[l] === dist_sqr(as[l], bs[l]), 2:8)
        end
    end
end
nothing
