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

"""
    ball_intersects_unit_cube(center, radius) -> Bool

Check if a ball (sphere) intersects with the unit cube centered at the origin,
extending from -0.5 to 0.5 in each axis.
"""
function ball_intersects_unit_cube(center, radius)
    # Squared distance from sphere center to nearest point on [-0.5, 0.5]^3
    h = oftype(radius, 1//2)
    d2 = relu(abs(center[1]) - h)^2 + relu(abs(center[2]) - h)^2 + relu(abs(center[3]) - h)^2
    d2 <= radius * radius
end

"""
    ball_intersects_biunit_cube(center, radius) -> Bool

Check if a ball (sphere) intersects with the biunit cube centered at the origin,
extending from -1 to 1 in each axis.
"""
function ball_intersects_biunit_cube(center, radius)
    # Squared distance from sphere center to nearest point on [-1, 1]^3
    h = oneunit(radius)
    d2 = relu(abs(center[1]) - h)^2 + relu(abs(center[2]) - h)^2 + relu(abs(center[3]) - h)^2
    d2 <= radius * radius
end
