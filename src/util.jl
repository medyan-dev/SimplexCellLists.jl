# Clamp `x` to `[0, 1]`, with NaN to zero. Written as selects that LLVM makes
# one max and one min, as `x > 0` is false for NaN.
@inline function clamp01nan(x)
    y = SIMD.vifelse(x > zero(x), x, zero(x))
    SIMD.vifelse(y < one(x), y, one(x))
end

"""
    dot3(a::Vec3, b::Vec3)

Return the dot product of two 3D vectors, computed with fused multiply adds.

Also works on `SVector{3, SIMD.Vec{N,T}}`, bit identical to the scalar version in each lane.
"""
@inline dot3(a::Vec3, b::Vec3) = fma(a[3], b[3], fma(a[2], b[2], a[1]*b[1]))

"""
    cross3(a::Vec3, b::Vec3)

Return the cross product of two 3D vectors, computed with fused multiply adds.

Also works on `SVector{3, SIMD.Vec{N,T}}`, bit identical to the scalar version in each lane.
"""
@inline function cross3(a::Vec3, b::Vec3)
    SVector(
        fma(a[2], b[3], -(a[3]*b[2])),
        fma(a[3], b[1], -(a[1]*b[3])),
        fma(a[1], b[2], -(a[2]*b[1])),
    )
end

"""
    scale3(t, v::Vec3)
    scale3(v::Vec3, t)

Return `v` scaled by `t`, as `t*v` or `v*t` elementwise.

Unlike `t*v`, this also works when `t` is a `SIMD.Vec{N,T}`.
"""
@inline scale3(t::Union{Number, SIMD.Vec}, v::Vec3) = SVector(t*v[1], t*v[2], t*v[3])
@inline scale3(v::Vec3, t::Union{Number, SIMD.Vec}) = SVector(v[1]*t, v[2]*t, v[3]*t)

"""
    fma3(t, a::Vec3, b::Vec3)

Return `t*a + b` elementwise, computed with fused multiply adds.

Also works when `t` is a `SIMD.Vec{N,T}`, bit identical to the scalar version in each lane.
"""
@inline fma3(t::Union{Number, SIMD.Vec}, a::Vec3, b::Vec3) = SVector(fma(t, a[1], b[1]), fma(t, a[2], b[2]), fma(t, a[3], b[3]))

function relu(x::T)::T where T
    SIMD.vifelse(x > zero(x), x, zero(x))
end

_sqrt_fast(x::Union{Float64, Float32}) = Core.Intrinsics.sqrt_llvm(x)
_sqrt_fast(x) = sqrt(x)

"""
    norm_fast(v::Vec3)

Return the Euclidean norm of a 3D vector, using [`dot3`](@ref) and an unchecked square root.

Also works on `SVector{3, SIMD.Vec{N,T}}`.
"""
@inline norm_fast(v::Vec3) = _sqrt_fast(dot3(v, v))

# Return the zero based XYZ grid index.
# `pos_norm` has an origin at the center of the grid
# and a scale of 2 units per grid spacing.
# `pos_norm` is expected to be inside the grid.
round2grid(pos_norm::Vec3, sz::Vec3{Int})::Vec3{Int} = (floor.((Int,), pos_norm) .+ sz) .>> 1
