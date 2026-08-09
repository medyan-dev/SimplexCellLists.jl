clamp01nan(x) = ifelse(isnan(x), zero(x), clamp(x,zero(x),one(x)))

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
