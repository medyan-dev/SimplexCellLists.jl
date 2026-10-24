@inline function clamp01nan(x)
    SIMD.vifelse(isnan(x), zero(x), Base.FastMath.min_fast(Base.FastMath.max_fast(x, zero(x)), one(x)))
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
norm_fast(v::Vec3) = _sqrt_fast(dot3(v, v))

# Return the zero based XYZ grid index.
# `pos_norm` has an origin at the center of the grid
# and a scale of 2 units per grid spacing.
# `pos_norm` is expected to be inside the grid.
round2grid(pos_norm::Vec3, sz::Vec3{Int})::Vec3{Int} = (floor.((Int,), pos_norm) .+ sz) .>> 1
