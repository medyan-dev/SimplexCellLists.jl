# Trace the first order error bounds of the distance kernels with
# ErrorBounds.jl, then check them against BigFloat on adversarial inputs.
#
# Run from the test environment, which has SimplexCellLists and StaticArrays:
#     include("error-analysis/analyze.jl")

using SimplexCellLists, StaticArrays, Random, LinearAlgebra
using SimplexCellLists: dist_sqr_r_s_t
include(joinpath(@__DIR__, "ErrorBounds.jl"))
using .ErrorBounds
using .ErrorBounds: DEFS

setprecision(BigFloat, 512)

norm_sym(n) = SPoly(Sym(n, 0))
ratio_sym(a, b) = SPoly(Sym(a, 0)) * SPoly(Sym(b, 0), -1)

const KERNELS = (
    point_point = (x = (:x,), y = (:y,),
        names = Dict((:-, :y, :x) => :R)),
    point_line = (x = (:x,), y = (:a, :b),
        names = Dict((:-, :a, :x) => :P, (:-, :b, :a) => :R)),
    point_triangle = (x = (:x,), y = (:a, :b, :c),
        names = Dict((:-, :b, :a) => :E0, (:-, :c, :a) => :E1, (:-, :a, :x) => :BP,
            (:-, :E1, :E0) => :E2, (:+, :BP, :E0) => :BP2),
        cross = [:N],
        # In plane distance from the boundary per unit of `sbar`, `tbar` and
        # `sbar + tbar`: the heights of the triangle over its edges.
        jumps = [[ratio_sym(:N, :E1), ratio_sym(:N, :E0), ratio_sym(:N, :E2)]]),
    line_line = (x = (:x1, :x2), y = (:y1, :y2),
        names = Dict((:-, :x2, :x1) => :P, (:-, :y2, :y1) => :Q, (:-, :y1, :x1) => :W,
            (:+, :W, :Q) => :WQ, (:-, :W, :P) => :WP),
        cross = [:n],
        # In plane distance between the projected line segments per unit of
        # `sbar` and `tbar`: `sin θ` times the length of the other segment.
        jumps = [[ratio_sym(:n, :Q), ratio_sym(:n, :P)]]),
)

"""
    trace(k)

Run kernel `k` on traced inputs. Return the distance bound
`|√d̂² - d| ≤ u*(abs + rel*d + moves)` as `(; abs, rel, moves)`, with `abs`
still in components, and the traced outputs.
"""
function trace(k)
    reset!(; names=k.names, cross=get(k, :cross, Symbol[]), jumps=get(k, :jumps, Vector{SPoly}[]))
    d2, r, s, t = dist_sqr_r_s_t(input_simplex(k.x...), input_simplex(k.y...))
    farthest(vs) = reduce(maxpoly, [vecnormbound(vertex_moves(v)) for v in vs])
    moves = farthest(k.x) + farthest(k.y)
    # The bound in norms only, with the plain cross products bounded.
    simple = substitute_abscross(normbound(d2.root[1]))
    (; abs = d2.root[1], simple, rel = d2.root[2], moves, d2, r, s, t,
        defs = copy(DEFS), cross_defs = copy(CROSS_DEFS), abscross_defs = copy(ABSCROSS_DEFS),
        notes = unique(ErrorBounds.NOTES), assumptions = unique(ASSUMPTIONS))
end

# Values of every size symbol for the exact vertices.
function symvalues(k, tr, xs, ys)
    vecs = Dict{Symbol, Vector{BigFloat}}()
    for (n, v) in zip(k.x, xs); vecs[n] = big.(v); end
    for (n, v) in zip(k.y, ys); vecs[n] = big.(v); end
    function value(n)
        get!(vecs, n) do
            if haskey(tr.defs, n)
                op, a, b = tr.defs[n]
                op == :- ? value(a) - value(b) : value(a) + value(b)
            else
                a, b = tr.cross_defs[n]
                cross(value(a), value(b))
            end
        end
    end
    for (n, (a, b)) in tr.abscross_defs
        va, vb = abs.(value(a)), abs.(value(b))
        vecs[n] = [va[2]*vb[3] + va[3]*vb[2], va[3]*vb[1] + va[1]*vb[3], va[1]*vb[2] + va[2]*vb[1]]
    end
    vals = Dict{Sym, BigFloat}()
    for n in union(keys(vecs), keys(tr.defs), keys(tr.cross_defs))
        v = value(n)
        vals[Sym(n, 0)] = norm(v)
        for i in 1:3; vals[Sym(n, i)] = abs(v[i]); end
    end
    vals
end

"""
    check(k, tr, xs, ys)

The error of the distance over its traced bound, for vertices `xs` and `ys`.
"""
function check(k, tr, xs::NTuple{N, SVector{3, T}}, ys::NTuple{M, SVector{3, T}}) where {N, M, T}
    d2 = dist_sqr_r_s_t(SVector(xs), SVector(ys))[1]
    d = sqrt(dist_sqr_r_s_t(SVector(map(v -> big.(v), xs)), SVector(map(v -> big.(v), ys)))[1])
    vals = symvalues(k, tr, xs, ys)
    u = big(eps(T))/2
    bound = u*(boundvalue(tr.simple, vals) + tr.rel*d + boundvalue(tr.moves, vals))
    Float64(abs(sqrt(big(d2)) - d)/bound)
end

function report(name, k, gens; n=20_000)
    tr = trace(k)
    println("=== $name ===")
    println("  |√d̂² - d| ≤ u*(", tr.simple, " + ", tr.rel, "d)  in the moved vertices")
    println("            + u*(", tr.moves, ")  for the vertex moves")
    foreach(s -> println("  note: ", s), tr.notes)
    foreach(s -> println("  assumes: ", s), tr.assumptions)
    for T in (Float32, Float64), (gname, g) in gens
        ratios = filter(isfinite, [check(k, tr, g(T)...) for _ in 1:n])
        println("  ", rpad("$T $gname", 40), " max error/bound ", round(maximum(ratios), sigdigits=3))
    end
    tr
end

randdir() = normalize(randn(3))
vec3(T, v) = SVector{3, T}(v)
function inplane(n)
    e1 = normalize(cross(n, randn(3)))
    (e1, cross(n, e1))
end

const GENS = (
    point_line = (
        "random [0,100]" => T -> ((vec3(T, 100rand(3)),), (vec3(T, 100rand(3)), vec3(T, 100rand(3)))),
        "point on the segment" => T -> begin
            a = 100randn(3); b = 100randn(3); s = rand()
            ((vec3(T, a .+ s .* (b .- a)),), (vec3(T, a), vec3(T, b)))
        end,
        "near an end, 1e4 from the origin" => T -> begin
            c = 1e4randdir(); dd = randdir(); n = normalize(cross(dd, randn(3)))
            s = rand((0.0, 1.0)) + 1e-6randn()
            ((vec3(T, c .+ 0.3n),), (vec3(T, c .- 100s*dd), vec3(T, c .+ 100(1-s)*dd)))
        end,
    ),
    point_triangle = (
        "random [0,100]" => T -> ((vec3(T, 100rand(3)),), Tuple(vec3(T, 100rand(3)) for _ in 1:3)),
        "0.01 above, near an edge" => T -> begin
            n = randdir(); e1, e2 = inplane(n); a = 100randn(3)
            b = a .+ 100e1; c = a .+ 100e2; w = rand((0.0, 0.5)) + 1e-6randn()
            ((vec3(T, a .+ w*100e1 .+ (0.5 - w)*100e2 .+ 0.01n),), (vec3(T, a), vec3(T, b), vec3(T, c)))
        end,
        "sliver" => T -> begin
            dd = randdir(); n = normalize(cross(dd, randn(3))); a = 100randn(3)
            b = a .+ 100dd; c = a .+ 50dd .+ 1e-3n
            ((vec3(T, a .+ 40dd .+ 0.01randdir()),), (vec3(T, a), vec3(T, b), vec3(T, c)))
        end,
        "needle, near the short edge" => T -> begin
            n = randdir(); e1, e2 = inplane(n); a = 100randn(3)
            b = a .+ 100e1; c = b .+ 0.01e2; w = rand()
            ((vec3(T, b .+ w*0.01e2 .+ 1e-7randn()*e1 .+ 0.01n),), (vec3(T, a), vec3(T, b), vec3(T, c)))
        end,
        "needle, 0.2 past the short edge" => T -> begin
            n = randdir(); e1, e2 = inplane(n); a = 100randn(3)
            b = a .+ 100e1; c = b .+ 0.001e2
            ((vec3(T, b .+ 0.001rand()*e2 .+ 0.2rand()*e1 .+ 0.01n),), (vec3(T, a), vec3(T, b), vec3(T, c)))
        end,
        "shared vertex pool" => T -> begin
            pool = [vec3(T, 100rand(3)) for _ in 1:3]
            ((rand(pool),), Tuple(rand(pool) for _ in 1:3))
        end,
    ),
    line_line = (
        "random [0,100]" => T -> (Tuple(vec3(T, 100rand(3)) for _ in 1:2), Tuple(vec3(T, 100rand(3)) for _ in 1:2)),
        "crossing 0.01 apart, 1e4 from the origin" => T -> begin
            c = 1e4randdir(); n = randdir(); e1, e2 = inplane(n)
            d1 = normalize(randn()*e1 + randn()*e2); d2 = normalize(randn()*e1 + randn()*e2); s = rand(); t = rand()
            ((vec3(T, c .- 100s*d1), vec3(T, c .+ 100(1-s)*d1)), (vec3(T, c .- 100t*d2 .+ 0.01n), vec3(T, c .+ 100(1-t)*d2 .+ 0.01n)))
        end,
        "nearly parallel" => T -> begin
            dd = randdir(); n = normalize(cross(dd, randn(3))); e = normalize(cross(dd, n)); a = 100randn(3)
            θ = 1e-4randn()
            ((vec3(T, a), vec3(T, a .+ 100dd)), (vec3(T, a .+ 30dd .+ 0.01n .- 50θ*e), vec3(T, a .+ 130dd .+ 0.01n .+ 50θ*e)))
        end,
        "crossing near an end" => T -> begin
            n = randdir(); e1, e2 = inplane(n); a = 100randn(3)
            s = rand((0.0, 1.0)) + 1e-6randn(); d2 = normalize(e1 + e2)
            ((vec3(T, a), vec3(T, a .+ 100e1)), (vec3(T, a .+ 100s*e1 .- 50d2 .+ 0.01n), vec3(T, a .+ 100s*e1 .+ 50d2 .+ 0.01n)))
        end,
        "shared vertex pool" => T -> begin
            pool = [vec3(T, 100rand(3)) for _ in 1:3]
            (Tuple(rand(pool) for _ in 1:2), Tuple(rand(pool) for _ in 1:2))
        end,
    ),
)

Random.seed!(4)
report("point-line", KERNELS.point_line, GENS.point_line)
report("point-triangle", KERNELS.point_triangle, GENS.point_triangle)
report("line-line", KERNELS.line_line, GENS.line_line)
