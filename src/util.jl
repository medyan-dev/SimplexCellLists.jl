clamp01nan(x) = ifelse(isnan(x), zero(x), clamp(x,zero(x),one(x)))

function relu(x::T)::T where T
    (x > zero(x) ? x : zero(x))
end

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
