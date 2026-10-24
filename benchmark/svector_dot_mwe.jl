# MWE: a non-inlined function taking SVector{3,Float32} arguments runs 2-4x slower
# depending on which two of the three floats LLVM pairs into one 8-byte load.
# `LinearAlgebra.dot` on SVector{3,Float32} happens to trigger the slow pairing.
#
# Chain of causes, as far as I can tell:
#   1. Julia passes a 12-byte SVector{3,Float32} to a non-inlined callee by pointer
#      to a stack copy, and the caller writes that copy as an 8-byte store at +0
#      plus a 4-byte store at +8 (see the caller listing printed below).
#   2. In the callee, LLVM's SLP vectorizer pairs two of the three elements into
#      one 8-byte load. If it pairs elements (2,3) the load is at offset +4 and
#      straddles both stores, so store-to-load forwarding fails on every call.
#      Pairing (1,2) gives a load at +0 that is covered by the 8-byte store.
#   3. StaticArrays' `_vecdot` seeds its `@simd` reduction with
#      `ret = product(za, zb) + product(za, zb)`, i.e. an unfoldable `+ 0.0`.
#      That is enough to flip the SLP pairing to (2,3). A hand-written dot with
#      the same `+ 0f0` seed reproduces it; seeding with `-0f0` (foldable) or
#      with the first product does not.
#
# Measured on Julia 1.13.1 / LLVM 20.1.8 / StaticArrays 1.9.23, Ryzen 7 9800X3D (znver5):
#   Part A  sum12 1.3 ns/call   sum23 5.9 ns/call      (4.5x, no dot product involved)
#   Part B  plain 2.6 ns/call   dot 7.0 ns/call  plus0 7.0  neg0 2.6  (2.7x)
#   Part C  16-byte padded elements: all ~2.3 ns; loading from the array by index: all 1.1 ns
# Float64 SVector{3} (24-byte copy as 16 + 8, 16-byte load at +8) shows the same effect.

using StaticArrays, LinearAlgebra, Chairmarks, Random, InteractiveUtils

const V3 = SVector{3, Float32}
rng = Xoshiro(1)
n = 10_000
As = [randn(rng, V3) * 50f0 for _ in 1:n]
Bs = [randn(rng, V3) * 50f0 for _ in 1:n]

function count_below(f, As, Bs, cut)
    m = 0
    @inbounds for k in eachindex(As)
        m += f(As[k], Bs[k]) < cut
    end
    m
end
per_call(f, As, Bs, cut) = round((@b count_below($f, $As, $Bs, $cut)).time * 1e9 / length(As); digits = 2)

# ---- Part A: pure load pattern, no dot product -------------------------------
@noinline sum12(a, b) = (a[1] + b[1]) * (a[2] + b[2])   # LLVM: one 8-byte load at +0 per argument
@noinline sum23(a, b) = (a[2] + b[2]) * (a[3] + b[3])   # LLVM: one 8-byte load at +4 per argument
println("Part A (ns/call): sum12 ", per_call(sum12, As, Bs, 100f0), "   sum23 ", per_call(sum23, As, Bs, 100f0))

# ---- Part B: how `dot` gets there -------------------------------------------
plaindot(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
plus0dot(a, b) = ((a[1]*b[1] + 0f0) + a[2]*b[2]) + a[3]*b[3]    # what StaticArrays' _vecdot effectively does
neg0dot(a, b)  = ((a[1]*b[1] + -0f0) + a[2]*b[2]) + a[3]*b[3]   # x + -0.0 == x, so LLVM folds it
@noinline one_plain(P0, P1, Q0, Q1) = plaindot(P1 - P0, Q1 - Q0)
@noinline one_dot(P0, P1, Q0, Q1)   = dot(P1 - P0, Q1 - Q0)
@noinline one_plus0(P0, P1, Q0, Q1) = plus0dot(P1 - P0, Q1 - Q0)
@noinline one_neg0(P0, P1, Q0, Q1)  = neg0dot(P1 - P0, Q1 - Q0)
function count_below4(f, P0s, P1s, Q0s, Q1s, cut)
    m = 0
    @inbounds for k in eachindex(P0s)
        m += f(P0s[k], P1s[k], Q0s[k], Q1s[k]) < cut
    end
    m
end
Cs = [randn(rng, V3) * 50f0 for _ in 1:n]
Ds = [randn(rng, V3) * 50f0 for _ in 1:n]
per_call4(f) = round((@b count_below4($f, $As, $Bs, $Cs, $Ds, 400f0)).time * 1e9 / n; digits = 2)
println("Part B (ns/call): plain ", per_call4(one_plain), "   dot ", per_call4(one_dot), "   plus0 ", per_call4(one_plus0), "   neg0 ", per_call4(one_neg0))

# ---- Part C: the slow pairing is harmless without the split stack copy ------
const V4 = SVector{4, Float32}
As4 = [V4(v..., 0f0) for v in As]; Bs4 = [V4(v..., 0f0) for v in Bs]
println("Part C padded to 16 bytes (ns/call): sum12 ", per_call(sum12, As4, Bs4, 100f0), "   sum23 ", per_call(sum23, As4, Bs4, 100f0))
@noinline sum23_idx(As, Bs, k) = @inbounds (As[k][2] + Bs[k][2]) * (As[k][3] + Bs[k][3])
function count_below_idx(f, As, Bs, cut)
    m = 0
    @inbounds for k in eachindex(As)
        m += f(As, Bs, k) < cut
    end
    m
end
println("Part C loading from the array, no copy (ns/call): sum23 ", round((@b count_below_idx(sum23_idx, $As, $Bs, 100f0)).time * 1e9 / n; digits = 2))

# ---- The code -----------------------------------------------------------------
println("\n--- caller: how the SVector arguments are copied to the stack (count_below with sum23) ---")
let io = IOBuffer()
    code_native(io, count_below, Tuple{typeof(sum23), Vector{V3}, Vector{V3}, Float32}; debuginfo = :none)
    for l in split(String(take!(io)), '\n')
        occursin(r"ptr \[rbp - \d+\], |call\t", l) && !occursin("Spill", l) && println(l)
    end
end
println("\n--- callee sum23: 8-byte load at +4 straddles the 8+4 byte stores ---")
code_native(sum23, Tuple{V3, V3}; debuginfo = :none)
println("\n--- callee sum12: 8-byte load at +0 is covered by the 8-byte store ---")
code_native(sum12, Tuple{V3, V3}; debuginfo = :none)
println("\n--- LLVM IR of dot(::V3, ::V3): the `fadd ... 0.0` seed from StaticArrays' _vecdot ---")
code_llvm(dot, Tuple{V3, V3}; debuginfo = :none)
println("\nJulia ", VERSION, ", LLVM ", Base.libllvm_version, ", StaticArrays ", pkgversion(StaticArrays), ", ", Sys.cpu_info()[1].model)
