# Clamp to [0, 1], with NaN to 0. Works on `Vec` lanes too.
@inline clamp01nan(x) = vifelse(isnan(x), zero(x), Base.FastMath.max_fast(Base.FastMath.min_fast(x, one(x)), zero(x)))

# Plain dot product of 3-vectors. `LinearAlgebra.dot` on `SVector`s reduces
# with `@simd`, which compiles to noticeably slower scalar code.
@inline _dot(a::SVector{3}, b::SVector{3}) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

function relu(x::T)::T where T
    (x > zero(x) ? x : zero(x))
end

_sqrt_fast(x::Union{Float64, Float32}) = Core.Intrinsics.sqrt_llvm(x)
_sqrt_fast(x) = sqrt(x)

norm_fast(v) = _sqrt_fast(v ⋅ v)

# Return the zero based XYZ grid index.
# `pos_norm` has an origin at the center of the grid
# and a scale of 2 units per grid spacing.
# `pos_norm` is expected to be inside the grid.
round2grid(pos_norm::Vec3, sz::Vec3{Int})::Vec3{Int} = (floor.((Int,), pos_norm) .+ sz) .>> 1
