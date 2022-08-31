const MAX_SAMPLE_SPACING = Float32(1//8)

"""
This implementation is based on the idea in https://gitlab.com/f-nedelec/cytosim/-/blob/master/src/sim/fiber_grid.h
And in figure 10 of "Francois Nedelec and Dietrich Foethke 2007 New J. Phys. 9 427"
to paint a grid with elements based on a maximum range.
"""
struct Painter
    grid_start::SVector{3,Float32}

    grid_size::SVector{3,Int32}

    voxel_length::Float32

    inv_voxel_length::Float32

    "Max ranges of each group in grid units"
    max_range::SVector{2, Vector{Float32}}

    "indexed by shape id then (groupid,x,y,z) to get a vector of elements painted"
    grids::SVector{2, Array{Vector{Tuple{Int32, Float32}}, 4}}

    "in grid units"
    data::Tuple{Vector{Vector{Point}}, Vector{Vector{Line}}}

    "true if element exists"
    exists::Tuple{Vector{Vector{Bool}}, Vector{Vector{Bool}}}
end


function Painter(numpointgroups::Integer, numlinegroups::Integer ;
        grid_start::SVector{3,AbstractFloat},
        grid_size::SVector{3,Integer},
        voxel_length::AbstractFloat,
        max_range::SVector{2, Vector{AbstractFloat}},
    )
    pointgrid = [Tuple{Int32, Float32}[] for i in 1:numpointgroups, j in 1:grid_size[1], k in 1:grid_size[2], l in 1:grid_size[3]]
    linegrid = [Tuple{Int32, Float32}[] for i in 1:numlinegroups, j in 1:grid_size[1], k in 1:grid_size[2], l in 1:grid_size[3]]
    Painter(
        grid_start,
        grid_size,
        voxel_length,
        inv(voxel_length),
        max_range= max_range .* inv(voxel_length),
        SA[pointgrid, linegrid],
        ([Point[] for i in 1:numpointgroups], [Line[] for i in 1:numlinegroups]),
        ([Bool[] for i in 1:numpointgroups],  [Bool[] for i in 1:numlinegroups]),
    )
end

"Return nearest grid point, x is in grid units"
_nearestGridPoint(m::Painter, x::SVector{3,Float32}) = clamp.(round.(Int32, bottom), Int32(0), (m.grid_size .- Int32(1)))

"""
Paint an element onto the grid.
All inputs have grid units.
Max range is extended by `MAX_SAMPLE_SPACING` because of how elements are sampled. 
Any voxel that may contain points within `max_range + MAX_SAMPLE_SPACING` of `element`
will have `i`, and the distance squared, added to that voxel's list.

In grid units, voxel (0,0,0) is a unit cube with center at SA[0.0, 0.0, 0.0]
"""
function _paintElement(m::Painter, groupid, i, element::Simplex{N}, max_range) where {N}
    grid = grids[N]
    #TODO try removing Float32(1//2*√(3)) and using an axis aligned unit cube to simplex distance function
    extended_max_range = max_range + MAX_SAMPLE_SPACING
    voxelcutoff = Float32(1//2*√(3)) + max_range + MAX_SAMPLE_SPACING
    voxelcutoff2 = voxelcutoff^2
    #TODO try cytosim's method of rasterizing instead of using an axis aligned bounding box
    # Get the axis aligned bounding box
    bottom = min.(element...)
    top = max.(element...)
    # Add max_range + MAX_SAMPLE_SPACING
    bottom = bottom .- extended_max_range
    top = top .+ extended_max_range
    bottom_int = _nearestGridPoint(m, bottom)
    top_int = _nearestGridPoint(m, top)
    # Go through all voxel center points and check distance to element and paint if in range.
    for voxel_int in Iterators.product(range.(bottom_int,top_int))
        @inline voxel_float = convert.(Float32,voxel_int)
        @inline d2 = distSqr(SA[voxel_float,], element)
        if d2 ≤ voxelcutoff2
            voxelid = voxel_int .+ Int32(1)
            voxel = grid[groupid, voxelid...]
            push!(voxel, (i,d2))
        end
    end
end


"""
Set the elements stored in the cell list.
"""
function setElements(m::Painter, points, lines)
    @argcheck length(points) == length(m.data[1])
    @argcheck length(lines) == length(m.data[2])
    for (n, in_data) in enumerate((points, lines,))
        for (groupid, in_group) in enumerate(in_data)
            group = m.data[n][groupid]
            groupexists = m.exists[n][groupid]
            resize!(group, length(in_group))
            resize!(groupexists, length(in_group))
            group .= (in_group .- m.grid_start) .* m.inv_voxel_length
            groupexists .= true
        end
    end
    for (groupid, group) in enumerate(m.data[1])
        max_range = m.max_range[1][groupid]
        for (i, element) in enumerate(group)
            _paintElement(m, groupid, i, element, max_range)
        end
    end
    for (groupid, group) in enumerate(m.data[2])
        max_range = m.max_range[2][groupid]
        for (i, element) in enumerate(group)
            _paintElement(m, groupid, i, element, max_range)
        end
    end

end

"""
add a point, return the point id.
"""
function addElement(m::Painter, groupid::Integer, element::Simplex{N}) where {N} ::Int32
    x = (element .- m.grid_start) .* m.inv_voxel_length
    group = push!(m.data[N][groupid], x)
    push!(m.exists[N][groupid], true)
    i = length(group)
    _paintElement(m, groupid, i, x, m.max_range[N][groupid])
    return i
end

"""
delete an element, the other element ids are stable.
"""
function deleteElement(m::Painter, groupid::Integer, elementid::Integer, elementtype::Type{Simplex{N}}) where {N} ::Nothing
    m.exists[N][groupid][elementid] = false
    return
end

"""
map a function to all elements in `groupid` in range of one element

The function f should have the same form as used in CellListMap.jl
Except here `x` and `y` are `SVector{N, SVector{3, Float32}}`, `SVector{M, SVector{3, Float32}}`

`x` is always `element` and i is always 0.

    function f(x,y,i,j,d2,output)
        # update output
        return output
    end
"""
function mapSimplexElements!(
        f, 
        output, 
        m::Painter, 
        groupid::Integer, 
        in_x::Simplex{N}, 
        elementstype::Type{Simplex{M}}, 
        in_cutoff::Float32
        ;
        x = (in_x .- m.grid_start) .* m.inv_voxel_length,
        i::Integer=0,
        filterj=Returns(true),
    ) where {N, M}
    @argcheck cutoff ≤ m.max_range[M][groupid]*m.voxel_length
    cutoff = in_cutoff * m.inv_voxel_length
    group = m.data[M][groupid]
    exists = m.exists[M][groupid]
    grid = m.grids[M]
    voxelcutoff::Float32 = Float32(1//2*√(3)) + cutoff + MAX_SAMPLE_SPACING
    voxelcutoff_sqr = voxelcutoff^2
    cutoff_sqr = cutoff^2
    # sample unique voxels that the element goes in based on sampling points that are at most 1/8 voxel_length away from any point in the element.
    mapSampleVoxels(x) do voxelid, numvoxels
        voxel = grid[groupid, voxelid...]
        for (j, sqrdist) in voxel
            if filterj(j)
                if sqrdist ≤ voxelcutoff_sqr
                    # load other element
                    if exists[j]
                        y = group[j]
                        @inline d2, s_min = distSqr_sMin(x, y)
                        if d2 ≤ cutoff_sqr
                            canonical_voxelid = if numvoxels > 1
                                # check that this voxel is the canonical one.
                                closest_sampled_s = getClosestSampledS(s_min, x)
                                closest_sampled_norm_x_pt = closest_sampled_s ⋅ x
                                getVoxelId(closest_sampled_norm_x_pt)
                            else
                                voxelid
                            end
                            if voxelid == canonical_voxelid
                                in_y = (y .* m.voxel_length) .+ m.grid_start
                                @inline output = f(in_x, in_y, i, j, d2, output)
                            end
                        end
                    end
                end
            end
        end
    end
    return output
end


# """
# map a function to all pairs of elements in the same group in range of each other.

# The function f should have the same form as used in CellListMap.jl
# Except here `x` and `y` are `SVector{N, SVector{3, Float32}}`, `SVector{N, SVector{3, Float32}}`

#     function f(x,y,i,j,d2,output)
#         # update output
#         return output
#     end
# """
# function mapPairElements!(f, output, m::Painter, groupid::Integer, elementstype::Type{Simplex{N}}, cutoff::Float64) where {N}
#     # just double loop through all element in groupid
#     norm_cutoff = cutoff*m.inv_voxel_length
#     group = m.data[N][groupid]
#     exists = m.exists[N][groupid]
#     n = length(group)
#     for i in 1:(n-1)
#         if exists[i]
#             x = group[i]
#             output = mapSimplexElements!(f, output, m, groupid, x, elementstype, cutoff;
#                 i,
#                 filterj = >(i),
#             )
#         end
#     end
#     return output
# end


# """
# map a function to all pairs of elements in different groups in range of each other.

# The function f should have the same form as used in CellListMap.jl
# Except here `x` and `y` are `SVector{N, SVector{3, Float32}}`, `SVector{M, SVector{3, Float32}}`

#     function f(x,y,i,j,d2,output)
#         # update output
#         return output
#     end
# """
# function mapPairElementsElements!(
#         f, 
#         output, 
#         m::Painter, 
#         x_groupid::Integer, 
#         x_elementstype::Type{Element{N}}, 
#         y_groupid::Integer, 
#         y_elementstype::Type{Element{M}}, 
#         cutoff_sqr::Float32,
#     ) where {N, M}
#     # just double loop through all element in groupid
#     x_group = m.data[N][x_groupid]
#     x_exists = m.exists[N][x_groupid]
#     y_group = m.data[M][y_groupid]
#     y_exists = m.exists[M][y_groupid]
#     xn = length(x_group)
#     yn = length(y_group)
#     for i in 1:xn
#         if x_exists[i]
#             x = x_group[i]
#             for j in 1:yn
#                 if y_exists[j]
#                     y = y_group[j]
#                     @inline d2 = distSqr(x, y)
#                     if d2 ≤ cutoff_sqr
#                         @inline output = f(x, y, i, j, d2, output)
#                     end
#                 end
#             end
#         end
#     end
#     return output
# end