struct PointIdxPart
    i::UInt32
end
Base.@propagate_inbounds function _load_positions(pos, x::PointIdxPart)
    SA[pos[x.i]]
end
Base.@propagate_inbounds function _load_axis_bounds(pos, x::PointIdxPart, axis::Int)
    v = pos[x.i][axis]
    (v, v)
end

struct CLineIdxPart
    i::UInt32
    # j = i + 1
end
Base.@propagate_inbounds function _load_positions(pos, x::CLineIdxPart)
    SA[pos[x.i], pos[x.i+UInt32(1)]]
end
Base.@propagate_inbounds function _load_axis_bounds(pos, x::CLineIdxPart, axis::Int)
    a = pos[x.i][axis]
    b = pos[x.i+UInt32(1)][axis]
    (Base.FastMath.min_fast(a, b), Base.FastMath.max_fast(a, b))
end

struct LineIdxPart
    i::UInt32
    j::UInt32
end
Base.@propagate_inbounds function _load_positions(pos, x::LineIdxPart)
    SA[pos[x.i], pos[x.j]]
end
Base.@propagate_inbounds function _load_axis_bounds(pos, x::LineIdxPart, axis::Int)
    local a = pos[x.i][axis]
    local b = pos[x.j][axis]
    (Base.FastMath.min_fast(a, b), Base.FastMath.max_fast(a, b))
end

struct TriangleIdxPart
    i::UInt32
    j::UInt32
    k::UInt32
end
Base.@propagate_inbounds function _load_positions(pos, x::TriangleIdxPart)
    SA[pos[x.i], pos[x.j], pos[x.k]]
end
Base.@propagate_inbounds function _load_axis_bounds(pos, x::TriangleIdxPart, axis::Int)
    local a = pos[x.i][axis]
    local b = pos[x.j][axis]
    local c = pos[x.k][axis]
    (Base.FastMath.min_fast(a, b, c), Base.FastMath.max_fast(a, b, c))
end

@enum CollidePairs begin
    CollidePoint_Point=1
    CollidePoint_CLine
    CollidePoint_Line
    CollidePoint_Triangle
    CollideLinePoint_Triangle
    CollideTrianglePoint_Triangle
    CollideCLine_CLine
    CollideCLine_Line
    CollideCLine_TriangleLine
    CollideLine_Line
    CollideLine_TriangleLine
    CollideTriangleLine_TriangleLine
end

const N_COLLIDE_PAIRS = length(instances(CollidePairs))

"""
Return an empty `no_collide_pairs` tuple with one `Set{Pair{UInt32, UInt32}}`
for each `CollidePairs` instance, indexed by `Int(::CollidePairs)`.

See [`NeighborListInputs`](@ref) for the canonical order pairs must be stored in.
"""
function empty_no_collide_pairs()
    ntuple(N_COLLIDE_PAIRS) do i
        Set{Pair{UInt32, UInt32}}()
    end
end

@enum CollideObjectTypes begin
    CollidePoint=1
    CollideCLine
    CollideLine
    CollideLinePoint
    CollideTriangle
    CollideTriangleLine
    CollideTrianglePoint
end

const N_COLLIDE_OBJECT_TYPES = length(instances(CollideObjectTypes))

"""
Dispatch token and global configuration for neighbor-list construction.

`ObjectParams` is per-object. `PairParams` is per-pair.
"""
abstract type CollisionPolicy{ObjectParams, PairParams} end

object_params_type(::Type{<:CollisionPolicy{O, P}}) where {O, P} = O
pair_params_type(::Type{<:CollisionPolicy{O, P}}) where {O, P} = P
object_params_type(policy::CollisionPolicy) = object_params_type(typeof(policy))
pair_params_type(policy::CollisionPolicy) = pair_params_type(typeof(policy))

"""
    filter_object(
        policy::CollisionPolicy{ObjectParams, PairParams},
        object_params::ObjectParams,
        object_type::CollideObjectTypes,
        object_index::UInt32,
        radius::Float32,
    )::Bool

Return whether an object participates in neighbor-list construction.

Must be a pure function. Objects that fail this filter are excluded from the
neighbor lists entirely, and are never passed to `filter_pair` or `mix_params`.
"""
function filter_object end

"""
    filter_pair(
        policy::CollisionPolicy{ObjectParams, PairParams},
        a::ObjectParams,
        b::ObjectParams,
        object_type_a::CollideObjectTypes,
        object_type_b::CollideObjectTypes,
        object_index_a::UInt32,
        object_index_b::UInt32,
        radius_a::Float32,
        radius_b::Float32,
    )::Bool

Return whether a pair of objects participates in the neighbor list.

Must be a pure function. It filters pairs from being added to the neighbor
lists. It is only checked for pairs where both objects passed `filter_object`,
but the order relative to the distance check and the `no_collide_pairs`
exclusions is unspecified, so it may be called on pairs that are out of range
or explicitly excluded.
"""
function filter_pair end

"""
    mix_params(
        policy::CollisionPolicy{ObjectParams, PairParams},
        a::ObjectParams,
        b::ObjectParams,
        object_type_a::CollideObjectTypes,
        object_type_b::CollideObjectTypes,
        object_index_a::UInt32,
        object_index_b::UInt32,
        radius_a::Float32,
        radius_b::Float32,
    )::PairParams

Return the `PairParams`. Mixing cannot reject a pair.
"""
function mix_params end

struct DefaultObjectParams
    stiffness::Float32
    layers::UInt32
    no_collide_mask::UInt32
end

struct DefaultPairParams
    k::Float32
end

Base.@kwdef struct DefaultCollisionPolicy <: CollisionPolicy{
        DefaultObjectParams,
        DefaultPairParams,
    }
    # This is a scale for ensuring the line line force is smooth
    switchover_scale_unitless::Float32 = 8.0f0
end

@inline function filter_object(
        ::DefaultCollisionPolicy,
        params::DefaultObjectParams,
        object_type::CollideObjectTypes,
        object_index::UInt32,
        radius::Float32,
    )::Bool
    !iszero(params.stiffness)
end

"""
    can_collide(layers_a::UInt32, no_collide_mask_a::UInt32, layers_b::UInt32, no_collide_mask_b::UInt32)::Bool

Check if two objects can collide based on their collision layers and no-collide masks.

Each object has `collision_layers` (which layers it is on) and `no_collide_mask`
(which layers it will not collide with). Two objects cannot collide if either one
has disabled collisions with the other's layer:

    no_collide = (A.layers & B.no_collide_mask) ≠ 0  OR  (B.layers & A.no_collide_mask) ≠ 0

Default `collision_layers = UInt32(1)` and `no_collide_mask = UInt32(0)` (collide with all layers),
so all objects collide by default.
"""
@inline function can_collide(layers_a::UInt32, no_collide_mask_a::UInt32, layers_b::UInt32, no_collide_mask_b::UInt32)::Bool
    iszero(layers_a & no_collide_mask_b) && iszero(layers_b & no_collide_mask_a)
end

@inline function filter_pair(
        ::DefaultCollisionPolicy,
        a::DefaultObjectParams,
        b::DefaultObjectParams,
        object_type_a::CollideObjectTypes,
        object_type_b::CollideObjectTypes,
        object_index_a::UInt32,
        object_index_b::UInt32,
        radius_a::Float32,
        radius_b::Float32,
    )::Bool
    can_collide(a.layers, a.no_collide_mask, b.layers, b.no_collide_mask)
end

@inline function mix_params(
        ::DefaultCollisionPolicy,
        a::DefaultObjectParams,
        b::DefaultObjectParams,
        object_type_a::CollideObjectTypes,
        object_type_b::CollideObjectTypes,
        object_index_a::UInt32,
        object_index_b::UInt32,
        radius_a::Float32,
        radius_b::Float32,
    )::DefaultPairParams
    DefaultPairParams(
        (a.stiffness * b.stiffness) / (abs(a.stiffness) + abs(b.stiffness)),
    )
end

"""
Input data for building neighbor lists in the collision detection system.

`no_collide_pairs` holds explicitly excluded pairs of object indexes, one
`Set{Pair{UInt32, UInt32}}` per `CollidePairs` instance, indexed by
`Int(::CollidePairs)`. Pairs only take effect in canonical order; pairs stored
in the reverse order are silently ignored:

- Pairs of two different object types are ordered as named by the
  `CollidePairs` instance, e.g. `CollidePoint_CLine` excludes
  `point_index => cline_index`.
- Pairs of the same object type are ordered smaller index first, e.g.
  `CollideLine_Line` excludes `line_index_a => line_index_b` where
  `line_index_a < line_index_b`.
"""
Base.@kwdef mutable struct NeighborListInputs{Policy <: CollisionPolicy, ObjectParams}
    policy::Policy

    points::Vector{PointIdxPart} = PointIdxPart[]
    p_radius::Vector{Float32} = Float32[]
    p_params::Vector{ObjectParams} = ObjectParams[]

    clines::Vector{CLineIdxPart} = CLineIdxPart[]
    c_radius::Vector{Float32} = Float32[]
    c_params::Vector{ObjectParams} = ObjectParams[]

    lines::Vector{LineIdxPart} = LineIdxPart[]
    l_radius::Vector{Float32} = Float32[]
    l_params::Vector{ObjectParams} = ObjectParams[]

    line_points::Vector{PointIdxPart} = PointIdxPart[]
    lp_radius::Vector{Float32} = Float32[]
    lp_params::Vector{ObjectParams} = ObjectParams[]

    triangles::Vector{TriangleIdxPart} = TriangleIdxPart[]
    t_radius::Vector{Float32} = Float32[]
    t_params::Vector{ObjectParams} = ObjectParams[]

    triangle_lines::Vector{LineIdxPart} = LineIdxPart[]
    tl_radius::Vector{Float32} = Float32[]
    tl_params::Vector{ObjectParams} = ObjectParams[]

    triangle_points::Vector{PointIdxPart} = PointIdxPart[]
    tp_radius::Vector{Float32} = Float32[]
    tp_params::Vector{ObjectParams} = ObjectParams[]

    no_collide_pairs::NTuple{N_COLLIDE_PAIRS, Set{Pair{UInt32, UInt32}}} = empty_no_collide_pairs()

    skin::Float32 = 0.0f0
    nthreads::Int = 1
end

function NeighborListInputs(policy::Policy; kwargs...) where {Policy <: CollisionPolicy}
    ObjectParams = object_params_type(policy)
    NeighborListInputs{Policy, ObjectParams}(; policy, kwargs...)
end

struct NeighborListEdge{A, B, PairParams}
    a::A # IdxPart
    b::B # IdxPart
    L::Float32 # sum of radii
    params::PairParams
end

mutable struct NeighborLists{Policy <: CollisionPolicy, PairParams}
    policy::Policy
    PPNL::Vector{NeighborListEdge{PointIdxPart, PointIdxPart, PairParams}}
    PCNL::Vector{NeighborListEdge{PointIdxPart, CLineIdxPart, PairParams}}
    PLNL::Vector{NeighborListEdge{PointIdxPart, LineIdxPart, PairParams}}
    PTNL::Vector{NeighborListEdge{PointIdxPart, TriangleIdxPart, PairParams}}
    CCNL::Vector{NeighborListEdge{CLineIdxPart, CLineIdxPart, PairParams}}
    CLNL::Vector{NeighborListEdge{CLineIdxPart, LineIdxPart, PairParams}}
    LLNL::Vector{NeighborListEdge{LineIdxPart, LineIdxPart, PairParams}}
end
function NeighborLists(policy::Policy) where {Policy <: CollisionPolicy}
    PairParams = pair_params_type(policy)
    NeighborLists{Policy, PairParams}(
        policy,
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
end

"""
Check if all edges in `subset` are also present in `superset`.
Used to verify that the skinned neighbor list contains all the edges
that would be found with zero skin.
"""
function is_neighbor_list_subset(subset::NeighborLists, superset::NeighborLists)::Bool
    issubset(subset.PPNL, superset.PPNL) || return false
    issubset(subset.PCNL, superset.PCNL) || return false
    issubset(subset.PLNL, superset.PLNL) || return false
    issubset(subset.PTNL, superset.PTNL) || return false
    issubset(subset.CCNL, superset.CCNL) || return false
    issubset(subset.CLNL, superset.CLNL) || return false
    issubset(subset.LLNL, superset.LLNL) || return false
    return true
end

function _prepare_neighbor_lists!(s::NeighborLists, inputs::NeighborListInputs)
    s.policy = inputs.policy
    empty!(s.PPNL)
    empty!(s.PCNL)
    empty!(s.PLNL)
    empty!(s.PTNL)
    empty!(s.CCNL)
    empty!(s.CLNL)
    empty!(s.LLNL)
    nothing
end

function _argcheck_neighbor_list_input(inputs::NeighborListInputs)
    (; policy,
       points, p_radius, p_params,
       clines, c_radius, c_params,
       lines, l_radius, l_params,
       line_points, lp_radius, lp_params,
       triangles, t_radius, t_params,
       triangle_lines, tl_radius, tl_params,
       triangle_points, tp_radius, tp_params,
       no_collide_pairs, skin, nthreads) = inputs
    # Assert matching array lengths
    @argcheck length(points) == length(p_radius) == length(p_params)
    @argcheck length(clines) == length(c_radius) == length(c_params)
    @argcheck length(lines) == length(l_radius) == length(l_params)
    @argcheck length(line_points) == length(lp_radius) == length(lp_params)
    @argcheck length(triangles) == length(t_radius) == length(t_params)
    @argcheck length(triangle_lines) == length(tl_radius) == length(tl_params)
    @argcheck length(triangle_points) == length(tp_radius) == length(tp_params)
    
    # Assert no UInt32 overflow
    @argcheck length(points) < typemax(UInt32)
    @argcheck length(clines) < typemax(UInt32)
    @argcheck length(lines) < typemax(UInt32)
    @argcheck length(line_points) < typemax(UInt32)
    @argcheck length(triangles) < typemax(UInt32)
    @argcheck length(triangle_lines) < typemax(UInt32)
    @argcheck length(triangle_points) < typemax(UInt32)
end

function setup_neighbors_naive!(s::NeighborLists, pos, inputs::NeighborListInputs)::Nothing
    _argcheck_neighbor_list_input(inputs)
    (; policy,
       points, p_radius, p_params,
       clines, c_radius, c_params,
       lines, l_radius, l_params,
       line_points, lp_radius, lp_params,
       triangles, t_radius, t_params,
       triangle_lines, tl_radius, tl_params,
       triangle_points, tp_radius, tp_params,
       no_collide_pairs, skin, nthreads) = inputs

    _prepare_neighbor_lists!(s, inputs)
    if length(pos) < 2
        # Only one or zero points
        return
    end

    extra_cutoff = 2skin

    # Calculate Interactions
    foreach((
        (s.PPNL, no_collide_pairs[Int(CollidePoint_Point)], points, points, true, CollidePoint, CollidePoint, p_radius, p_params, p_radius, p_params),
        (s.PCNL, no_collide_pairs[Int(CollidePoint_CLine)], points, clines, false, CollidePoint, CollideCLine, p_radius, p_params, c_radius, c_params),
        (s.PLNL, no_collide_pairs[Int(CollidePoint_Line)], points, lines, false, CollidePoint, CollideLine, p_radius, p_params, l_radius, l_params),
        (s.PTNL, no_collide_pairs[Int(CollidePoint_Triangle)], points, triangles, false, CollidePoint, CollideTriangle, p_radius, p_params, t_radius, t_params),
        (s.PTNL, no_collide_pairs[Int(CollideLinePoint_Triangle)], line_points, triangles, false, CollideLinePoint, CollideTriangle, lp_radius, lp_params, t_radius, t_params),
        (s.PTNL, no_collide_pairs[Int(CollideTrianglePoint_Triangle)], triangle_points, triangles, false, CollideTrianglePoint, CollideTriangle, tp_radius, tp_params, t_radius, t_params),
        (s.CCNL, no_collide_pairs[Int(CollideCLine_CLine)], clines, clines, true, CollideCLine, CollideCLine, c_radius, c_params, c_radius, c_params),
        (s.CLNL, no_collide_pairs[Int(CollideCLine_Line)], clines, lines, false, CollideCLine, CollideLine, c_radius, c_params, l_radius, l_params),
        (s.CLNL, no_collide_pairs[Int(CollideCLine_TriangleLine)], clines, triangle_lines, false, CollideCLine, CollideTriangleLine, c_radius, c_params, tl_radius, tl_params),
        (s.LLNL, no_collide_pairs[Int(CollideLine_Line)], lines, lines, true, CollideLine, CollideLine, l_radius, l_params, l_radius, l_params),
        (s.LLNL, no_collide_pairs[Int(CollideLine_TriangleLine)], lines, triangle_lines, false, CollideLine, CollideTriangleLine, l_radius, l_params, tl_radius, tl_params),
        (s.LLNL, no_collide_pairs[Int(CollideTriangleLine_TriangleLine)], triangle_lines, triangle_lines, true, CollideTriangleLine, CollideTriangleLine, tl_radius, tl_params, tl_radius, tl_params),
    )) do (nl, excl, a_objs, b_objs, self, a_type, b_type, a_r, a_params, b_r, b_params)
        for i in UInt32(1):UInt32(length(a_objs))
            local a = a_objs[i]
            local a_pos = _load_positions(pos, a)
            local r_i = a_r[i]
            local params_i = a_params[i]
            filter_object(policy, params_i, a_type, i, r_i) || continue
            for j in (self ? (i+UInt32(1):UInt32(length(a_objs))) : (UInt32(1):UInt32(length(b_objs))))
                local b = b_objs[j]
                local r_j = b_r[j]
                local params_j = b_params[j]
                filter_object(policy, params_j, b_type, j, r_j) || continue
                filter_pair(policy, params_i, params_j, a_type, b_type, i, j, r_i, r_j) || continue
                local b_pos = _load_positions(pos, b_objs[j])
                local d2 = dist_sqr(
                    a_pos,
                    b_pos,
                )
                local cutoff = r_i + r_j + extra_cutoff
                d2 < cutoff^2 || continue
                (i => j) in excl && continue
                local pair_params = mix_params(
                    policy, params_i, params_j, a_type, b_type, i, j, r_i, r_j)
                push!(nl, NeighborListEdge(a, b, r_i + r_j, pair_params))
            end
        end
    end
    nothing
end


function setup_neighbors_sort_sweep!(s::NeighborLists, pos, inputs::NeighborListInputs)::Nothing
    _argcheck_neighbor_list_input(inputs)
    (; policy,
       points, p_radius, p_params,
       clines, c_radius, c_params,
       lines, l_radius, l_params,
       line_points, lp_radius, lp_params,
       triangles, t_radius, t_params,
       triangle_lines, tl_radius, tl_params,
       triangle_points, tp_radius, tp_params,
       no_collide_pairs, skin, nthreads) = inputs

    _prepare_neighbor_lists!(s, inputs)
    if length(pos) < 2
        # Only one or zero points
        return
    end

    extra_cutoff = 2skin
    PTYPE = eltype(eltype(pos))

    # Get simulation bounds
    min_p::SVector{3, PTYPE} = SVector{3, PTYPE}(Inf, Inf, Inf)
    max_p::SVector{3, PTYPE} = SVector{3, PTYPE}(-Inf, -Inf, -Inf)
    for p in pos
        all(isfinite, p) || error("positions must be finite")
        min_p = min.(min_p, p)
        max_p = max.(max_p, p)
    end
    box_width::SVector{3, PTYPE} = max_p - min_p
    max_box_width::PTYPE, sweep_axis = findmax(box_width)
    other_axis1, other_axis2 = mod1(1+sweep_axis, 3), mod1(2+sweep_axis, 3)
    # For a zero or tiny width this overflows to Inf, which is safe: all
    # bounding boxes then clamp to the full UInt16 range and only the
    # distance checks filter.
    scale_16bit = typemax(UInt16) / max_box_width
    # Include extra epsilon padding to account for floating-point rounding errors
    skin_with_eps = skin +  8 * eps(max_box_width)

    # Convert min and max interval to a 16 bit integer interval
    # scaling with max_box_width and with r + skin added padding
    # Round conservatively, clamping to the UInt16 range so huge radii or a
    # tiny bounding box can't overflow the conversion. NaN geometry errors.
    function quantize_bounds(axis_bounds::Tuple, r, base)::NTuple{2, UInt16}
        local lo, hi = axis_bounds
        # Add radius and skin padding, then shift relative to min_p
        local padded_lo = (lo - base) - (r + skin_with_eps)
        local padded_hi = (hi - base) + (r + skin_with_eps)
        # Scale to UInt16 range. The lower bound clamps to [0, typemax-1] and
        # the upper bound to [1, typemax] so that lower < upper holds even if
        # both saturate at the same end.
        local top = oftype(scale_16bit, typemax(UInt16))
        local quantized_lo = trunc(UInt16, clamp(padded_lo * scale_16bit, zero(top), top - one(top)))
        local quantized_hi = ceil(UInt16, clamp(padded_hi * scale_16bit, one(top), top))
        (quantized_lo, quantized_hi)
    end

    # Prepare data structures
    ActiveStruct = @NamedTuple{
        index::UInt32,
        bounds1_m::UInt16,
        bounds1_p::UInt16,
        bounds2_m::UInt16,
        bounds2_p::UInt16,
    }
    # List of active elements during the sweep
    active_objs = ntuple(x->Vector{ActiveStruct}(), N_COLLIDE_OBJECT_TYPES)
    # Place to keep track of each objects position in the active_objs list
    # This is to accelerate removing.
    active_objs_idx = (
        zeros(UInt32, length(p_radius)),
        zeros(UInt32, length(c_radius)),
        zeros(UInt32, length(l_radius)),
        zeros(UInt32, length(lp_radius)),
        zeros(UInt32, length(t_radius)),
        zeros(UInt32, length(tl_radius)),
        zeros(UInt32, length(tp_radius)),
    )
    # Save quantized 16 bit bounding box edges in the sweep_axis | 1 bit 0 - end, 1 - start | 15 bit object type tag | 32 bit index
    # This is what gets sorted
    sweep_base = min_p[sweep_axis]
    n_objects = (
        length(points) + length(clines) + length(lines) + length(line_points) +
        length(triangles) + length(triangle_lines) + length(triangle_points)
    )
    endpoints = Vector{UInt64}()
    # Two endpoints per object; a slight overestimate because some objects
    # get filtered out
    sizehint!(endpoints, 2 * n_objects)
    foreach((
        (CollidePoint, points, p_radius, p_params),
        (CollideCLine, clines, c_radius, c_params),
        (CollideLine, lines, l_radius, l_params),
        (CollideLinePoint, line_points, lp_radius, lp_params),
        (CollideTriangle, triangles, t_radius, t_params),
        (CollideTriangleLine, triangle_lines, tl_radius, tl_params),
        (CollideTrianglePoint, triangle_points, tp_radius, tp_params),
    )) do (obj_type, a_objs, a_r, a_params)
        for index in UInt32(1):UInt32(length(a_objs))
            local a = a_objs[index]
            local r = a_r[index]
            local params = a_params[index]
            filter_object(policy, params, obj_type, index, r) || continue
            local start::UInt16, stop::UInt16 = quantize_bounds(_load_axis_bounds(pos, a, sweep_axis), r, sweep_base)
            local e_start = (UInt64(start)<<48) | (UInt64(1)<<47) | (UInt64(obj_type)<<32) | UInt64(index)
            local e_stop  = (UInt64(stop )<<48) | (UInt64(0)<<47) | (UInt64(obj_type)<<32) | UInt64(index)
            push!(endpoints, e_start)
            push!(endpoints, e_stop)
        end
    end
    sort!(endpoints)
    
    # Helper to check if two 2D bounding boxes overlap
    function boxes_overlap(
        a_b1_m::UInt16, a_b1_p::UInt16, a_b2_m::UInt16, a_b2_p::UInt16,
        b_b1_m::UInt16, b_b1_p::UInt16, b_b2_m::UInt16, b_b2_p::UInt16
    )::Bool
        # Check overlap in both axes
        # Thanks to conservative rounding this can be inequalities.
        (a_b1_m < b_b1_p && b_b1_m < a_b1_p) && (a_b2_m < b_b2_p && b_b2_m < a_b2_p)
    end

    for i_endpoints in 1:length(endpoints)
        local e = endpoints[i_endpoints]
        local index = e%UInt32
        local obj_type = ((e>>32) & 0xFF)%Int
        local same_type_active_objs = active_objs[obj_type]
        local active_idxs = active_objs_idx[obj_type]
        local n_same_type_active_objs = length(same_type_active_objs)
        local isend = iszero(e & (1<<47))
        if !isend
            local new_obj::ActiveStruct
            # - Load the positions, radius, and collision params
            # - Get the 16 bit bounds in other_axis1 and other_axis2
            # Then go through each active elements
            #   First check other axis bounding boxes
            #   Then check if not an excluded pair
            #   Finally do the full distance test
            #   If these checks pass push to the pair list
            # Afterwards push the element to the active objs
            # and store its index in active objs so it can be quickly removed
            function _collide_active_list(index, a, a_pos, r_a, params_a, a_type,
                    bounds1, bounds2, objs, radii, params, other_type,
                    excl, active_list, NL, should_swap)
                for other in active_list
                    boxes_overlap(bounds1..., bounds2..., other.bounds1_m, other.bounds1_p, other.bounds2_m, other.bounds2_p) || continue
                    local j = other.index
                    if should_swap(index, j)
                        local _a = objs[j]
                        local _b = a
                        local _r_a = radii[j]
                        local _r_b = r_a
                        local _params_a = params[j]
                        filter_pair(policy, _params_a, params_a,
                            other_type, a_type, j, index, _r_a, _r_b) || continue
                        local _a_pos = _load_positions(pos, _a)
                        local _b_pos = a_pos
                        local _d2 = dist_sqr(_a_pos, _b_pos)
                        local _cutoff = _r_a + _r_b + extra_cutoff
                        _d2 < _cutoff^2 || continue
                        (j => index) in excl && continue
                        local pair_params = mix_params(
                            policy, _params_a, params_a,
                            other_type, a_type, j, index, _r_a, _r_b)
                        push!(NL, NeighborListEdge(_a, _b, _r_a + _r_b, pair_params))
                    else
                        local b = objs[j]
                        local r_b = radii[j]
                        local params_b = params[j]
                        filter_pair(policy, params_a, params_b,
                            a_type, other_type, index, j, r_a, r_b) || continue
                        local b_pos = _load_positions(pos, b)
                        local d2 = dist_sqr(a_pos, b_pos)
                        local cutoff = r_a + r_b + extra_cutoff
                        d2 < cutoff^2 || continue
                        (index => j) in excl && continue
                        local pair_params = mix_params(
                            policy, params_a, params_b,
                            a_type, other_type, index, j, r_a, r_b)
                        push!(NL, NeighborListEdge(a, b, r_a + r_b, pair_params))
                    end
                end
            end
            if obj_type == Int(CollidePoint)
                let # Points interact with: Points, CLines, Lines, Triangles
                    local a = points[index]
                    local a_pos = _load_positions(pos, a)
                    local r_a = p_radius[index]
                    local params_a = p_params[index]
                    local bounds1 = quantize_bounds(_load_axis_bounds(pos, a, other_axis1), r_a, min_p[other_axis1])
                    local bounds2 = quantize_bounds(_load_axis_bounds(pos, a, other_axis2), r_a, min_p[other_axis2])
                    # Point-Point (self interaction, need index < other index)
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollidePoint, bounds1, bounds2,
                        points, p_radius, p_params, CollidePoint, no_collide_pairs[Int(CollidePoint_Point)],
                        active_objs[Int(CollidePoint)], s.PPNL, >,
                    )
                    # Point-CLine
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollidePoint, bounds1, bounds2,
                        clines, c_radius, c_params, CollideCLine, no_collide_pairs[Int(CollidePoint_CLine)],
                        active_objs[Int(CollideCLine)], s.PCNL, Returns(false),
                    )
                    # Point-Line
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollidePoint, bounds1, bounds2,
                        lines, l_radius, l_params, CollideLine, no_collide_pairs[Int(CollidePoint_Line)],
                        active_objs[Int(CollideLine)], s.PLNL, Returns(false),
                    )
                    # Point-Triangle
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollidePoint, bounds1, bounds2,
                        triangles, t_radius, t_params, CollideTriangle, no_collide_pairs[Int(CollidePoint_Triangle)],
                        active_objs[Int(CollideTriangle)], s.PTNL, Returns(false),
                    )
                    new_obj = (;index, bounds1_m=bounds1[1], bounds1_p=bounds1[2], bounds2_m=bounds2[1], bounds2_p=bounds2[2])
                end
            elseif obj_type == Int(CollideCLine)
                let # CLines interact with: CLines, Lines, TriangleLines
                    local a = clines[index]
                    local a_pos = _load_positions(pos, a)
                    local r_a = c_radius[index]
                    local params_a = c_params[index]
                    local bounds1 = quantize_bounds(_load_axis_bounds(pos, a, other_axis1), r_a, min_p[other_axis1])
                    local bounds2 = quantize_bounds(_load_axis_bounds(pos, a, other_axis2), r_a, min_p[other_axis2])
                    # CLine-CLine (self interaction)
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideCLine, bounds1, bounds2,
                        clines, c_radius, c_params, CollideCLine, no_collide_pairs[Int(CollideCLine_CLine)],
                        active_objs[Int(CollideCLine)], s.CCNL, >,
                    )
                    # CLine-Line
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideCLine, bounds1, bounds2,
                        lines, l_radius, l_params, CollideLine, no_collide_pairs[Int(CollideCLine_Line)],
                        active_objs[Int(CollideLine)], s.CLNL, Returns(false),
                    )
                    # CLine-TriangleLine
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideCLine, bounds1, bounds2,
                        triangle_lines, tl_radius, tl_params, CollideTriangleLine, no_collide_pairs[Int(CollideCLine_TriangleLine)],
                        active_objs[Int(CollideTriangleLine)], s.CLNL, Returns(false),
                    )
                    # Reverse: active Points check against this new CLine
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideCLine, bounds1, bounds2,
                        points, p_radius, p_params, CollidePoint, no_collide_pairs[Int(CollidePoint_CLine)],
                        active_objs[Int(CollidePoint)], s.PCNL, Returns(true),
                    )
                    new_obj = (;index, bounds1_m=bounds1[1], bounds1_p=bounds1[2], bounds2_m=bounds2[1], bounds2_p=bounds2[2])
                end
            elseif obj_type == Int(CollideLine)
                let # Lines interact with: Lines, TriangleLines
                    local a = lines[index]
                    local a_pos = _load_positions(pos, a)
                    local r_a = l_radius[index]
                    local params_a = l_params[index]
                    local bounds1 = quantize_bounds(_load_axis_bounds(pos, a, other_axis1), r_a, min_p[other_axis1])
                    local bounds2 = quantize_bounds(_load_axis_bounds(pos, a, other_axis2), r_a, min_p[other_axis2])
                    # Line-Line (self interaction)
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideLine, bounds1, bounds2,
                        lines, l_radius, l_params, CollideLine, no_collide_pairs[Int(CollideLine_Line)],
                        active_objs[Int(CollideLine)], s.LLNL, >,
                    )
                    # Line-TriangleLine
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideLine, bounds1, bounds2,
                        triangle_lines, tl_radius, tl_params, CollideTriangleLine, no_collide_pairs[Int(CollideLine_TriangleLine)],
                        active_objs[Int(CollideTriangleLine)], s.LLNL, Returns(false),
                    )
                    # Reverse: active Points check against this new Line
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideLine, bounds1, bounds2,
                        points, p_radius, p_params, CollidePoint, no_collide_pairs[Int(CollidePoint_Line)],
                        active_objs[Int(CollidePoint)], s.PLNL, Returns(true),
                    )
                    # Reverse: active CLines check against this new Line
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideLine, bounds1, bounds2,
                        clines, c_radius, c_params, CollideCLine, no_collide_pairs[Int(CollideCLine_Line)],
                        active_objs[Int(CollideCLine)], s.CLNL, Returns(true),
                    )
                    new_obj = (;index, bounds1_m=bounds1[1], bounds1_p=bounds1[2], bounds2_m=bounds2[1], bounds2_p=bounds2[2])
                end
            elseif obj_type == Int(CollideLinePoint)
                let # LinePoints interact with: Triangles
                    local a = line_points[index]
                    local a_pos = _load_positions(pos, a)
                    local r_a = lp_radius[index]
                    local params_a = lp_params[index]
                    local bounds1 = quantize_bounds(_load_axis_bounds(pos, a, other_axis1), r_a, min_p[other_axis1])
                    local bounds2 = quantize_bounds(_load_axis_bounds(pos, a, other_axis2), r_a, min_p[other_axis2])
                    # LinePoint-Triangle
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideLinePoint, bounds1, bounds2,
                        triangles, t_radius, t_params, CollideTriangle, no_collide_pairs[Int(CollideLinePoint_Triangle)],
                        active_objs[Int(CollideTriangle)], s.PTNL, Returns(false),
                    )
                    new_obj = (;index, bounds1_m=bounds1[1], bounds1_p=bounds1[2], bounds2_m=bounds2[1], bounds2_p=bounds2[2])
                end
            elseif obj_type == Int(CollideTriangle)
                let # Triangles don't initiate interactions (points/linepoints interact with them)
                    # But we need to check active points/linepoints/trianglepoints against this new triangle
                    local a = triangles[index]
                    local a_pos = _load_positions(pos, a)
                    local r_a = t_radius[index]
                    local params_a = t_params[index]
                    local bounds1 = quantize_bounds(_load_axis_bounds(pos, a, other_axis1), r_a, min_p[other_axis1])
                    local bounds2 = quantize_bounds(_load_axis_bounds(pos, a, other_axis2), r_a, min_p[other_axis2])
                    # Reverse: active Points check against this new Triangle
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideTriangle, bounds1, bounds2,
                        points, p_radius, p_params, CollidePoint, no_collide_pairs[Int(CollidePoint_Triangle)],
                        active_objs[Int(CollidePoint)], s.PTNL, Returns(true),
                    )
                    # Reverse: active LinePoints check against this new Triangle
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideTriangle, bounds1, bounds2,
                        line_points, lp_radius, lp_params, CollideLinePoint, no_collide_pairs[Int(CollideLinePoint_Triangle)],
                        active_objs[Int(CollideLinePoint)], s.PTNL, Returns(true),
                    )
                    # Reverse: active TrianglePoints check against this new Triangle
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideTriangle, bounds1, bounds2,
                        triangle_points, tp_radius, tp_params, CollideTrianglePoint, no_collide_pairs[Int(CollideTrianglePoint_Triangle)],
                        active_objs[Int(CollideTrianglePoint)], s.PTNL, Returns(true),
                    )
                    new_obj = (;index, bounds1_m=bounds1[1], bounds1_p=bounds1[2], bounds2_m=bounds2[1], bounds2_p=bounds2[2])
                end
            elseif obj_type == Int(CollideTriangleLine)
                let # TriangleLines interact with: TriangleLines (self)
                    local a = triangle_lines[index]
                    local a_pos = _load_positions(pos, a)
                    local r_a = tl_radius[index]
                    local params_a = tl_params[index]
                    local bounds1 = quantize_bounds(_load_axis_bounds(pos, a, other_axis1), r_a, min_p[other_axis1])
                    local bounds2 = quantize_bounds(_load_axis_bounds(pos, a, other_axis2), r_a, min_p[other_axis2])
                    # TriangleLine-TriangleLine (self interaction)
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideTriangleLine, bounds1, bounds2,
                        triangle_lines, tl_radius, tl_params, CollideTriangleLine, no_collide_pairs[Int(CollideTriangleLine_TriangleLine)],
                        active_objs[Int(CollideTriangleLine)], s.LLNL, >,
                    )
                    # Reverse: active CLines check against this new TriangleLine
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideTriangleLine, bounds1, bounds2,
                        clines, c_radius, c_params, CollideCLine, no_collide_pairs[Int(CollideCLine_TriangleLine)],
                        active_objs[Int(CollideCLine)], s.CLNL, Returns(true),
                    )
                    # Reverse: active Lines check against this new TriangleLine
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideTriangleLine, bounds1, bounds2,
                        lines, l_radius, l_params, CollideLine, no_collide_pairs[Int(CollideLine_TriangleLine)],
                        active_objs[Int(CollideLine)], s.LLNL, Returns(true),
                    )
                    new_obj = (;index, bounds1_m=bounds1[1], bounds1_p=bounds1[2], bounds2_m=bounds2[1], bounds2_p=bounds2[2])
                end
            elseif obj_type == Int(CollideTrianglePoint)
                let # TrianglePoints interact with: Triangles
                    local a = triangle_points[index]
                    local a_pos = _load_positions(pos, a)
                    local r_a = tp_radius[index]
                    local params_a = tp_params[index]
                    local bounds1 = quantize_bounds(_load_axis_bounds(pos, a, other_axis1), r_a, min_p[other_axis1])
                    local bounds2 = quantize_bounds(_load_axis_bounds(pos, a, other_axis2), r_a, min_p[other_axis2])
                    # TrianglePoint-Triangle
                    _collide_active_list(index, a, a_pos, r_a, params_a, CollideTrianglePoint, bounds1, bounds2,
                        triangles, t_radius, t_params, CollideTriangle, no_collide_pairs[Int(CollideTrianglePoint_Triangle)],
                        active_objs[Int(CollideTriangle)], s.PTNL, Returns(false),
                    )
                    new_obj = (;index, bounds1_m=bounds1[1], bounds1_p=bounds1[2], bounds2_m=bounds2[1], bounds2_p=bounds2[2])
                end
            else
                error("unreachable")
            end
            @assert new_obj.index == index
            push!(same_type_active_objs, new_obj)
            active_idxs[index] = n_same_type_active_objs + 1
        else
            # Remove object from active list doing a swap
            local active_idx = active_idxs[index]
            @assert !iszero(active_idx)
            if n_same_type_active_objs != active_idx
                local last_same_type_active_obj = last(same_type_active_objs)
                same_type_active_objs[active_idx] = last_same_type_active_obj
                active_idxs[last_same_type_active_obj.index] = active_idx
            end
            pop!(same_type_active_objs)
            active_idxs[index] = 0
        end
    end
    nothing
end
