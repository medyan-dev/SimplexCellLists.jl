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
_nearestGridPoint(grid_size::SVector{3,Int32}, x::SVector{3,Float32}) = clamp.(round.(Int32, x), Int32(0), (grid_size .- Int32(1)))

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
    bottom_int = _nearestGridPoint(m.grid_size, bottom)
    top_int = _nearestGridPoint(m.grid_size, top)
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
This function is to avoid any floating point weirdness
"""
@noinline function _getSampleGridPoint(inv_n::Float32, n::Int32, R::Int32, x::Simplex{2}, grid_size)
    Q::Int32 = n-R
    p1 = R*x[1] + Q*x[2]
    p2 = inv_n*p1
    _nearestGridPoint(grid_size, p2)
end

"""
Get the grid point that is sampled by x, to y.
x and y are in grid units
"""
@noinline function _canonicalGridPoint(x::Simplex{N}, y::Simplex{M}, grid_size, extra) where {N,M} ::SVector{3, Int32}
    if N == 1
        return _nearestGridPoint(grid_size, x[1])
    elseif N == 2
        if isnothing(extra)
            return _nearestGridPoint(grid_size, 1//2*(x[2] + x[1]))
        else
            halfn::Int32, inv_n::Float32 = extra
            d2, smin = distSqr_sMin(x,y)
            n = Int32(2)*halfn 
            R::Int32 = trunc(Int32, smin[1]*halfn)*Int32(2) + Int32(1)
            return _getSampleGridPoint(inv_n, n, R, x, grid_size)
        end
    else
        error("Triangle not implemented yet")
    end
end

"""
Call output = f(gridpoint, is_one_voxel, extra, output)
for each unique voxel sampled, then return the final output.
All points on x must be within MAX_SAMPLE_SPACING of a sampled voxel.
"""
@inline function _mapSampleVoxels(f, grid_size, x::Simplex{N}, output) where {N}
    if N == 1
        gridpoint = _nearestGridPoint(grid_size, x[1])
        return @inline f(gridpoint, true, nothing, output)
    elseif N == 2
        # if the line segment is short enough just us the mid point
        r = x[2] - x[1]
        len2 = r ⋅ r
        if len2 < (MAX_SAMPLE_SPACING*2)^2
            gridpoint = _nearestGridPoint(grid_size, 1//2*(x[2] + x[1]))
            return f(gridpoint, true, nothing, output)
        end
        # if both ends, moved in by MAX_SAMPLE_SPACING are in the same voxel, just sample that voxel.
        len = √len2
        small_s = MAX_SAMPLE_SPACING/len
        start_pt = (1-small_s)*x[1] + (small_s)*x[2]
        end_pt = (1-small_s)*x[2] + (small_s)*x[1]
        gridpoint1 = _nearestGridPoint(grid_size, start_pt)
        gridpoint2 = _nearestGridPoint(grid_size, end_pt)
        if gridpoint1 == gridpoint2
            return f(gridpoint1, true, nothing, output)
        end
        # Now more than one voxel needs to be sampled, yikes.
        # I need to ensure that _canonicalGridPoint can find a gridpoint that is actually sampled
        # 
        num_sub_div = ceil(inv(2*MAX_SAMPLE_SPACING)*len)
        inv_n::Float32 = inv(2num_sub_div)
        halfn = Int32(num_sub_div)
        n = Int32(2)*halfn
        prevgridpoint = SA[Int32(-1),Int32(-1),Int32(-1)]
        for i in Int32(1):Int32(2):n
            gridpoint = _getSampleGridPoint(inv_n, n, i, x, grid_size)
            if gridpoint != prevgridpoint
                output = f(gridpoint, false, (halfn, inv_n), output)
                prevgridpoint = gridpoint
            end
        end
        return output
    else
        error("Triangle not implemented yet")
    end
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
    voxel_length_sqr = m.voxel_length^2
    # sample unique voxels that the element goes in based on sampling points that are at most 1/8 voxel_length away from any point in the element.
    return _mapSampleVoxels(x, m.grid_size, output) do gridpoint, is_one_voxel, extra, output
        local voxelid = gridpoint .+ Int32(1)
        local voxel = grid[groupid, voxelid...]
        for (j, sqrdist) in voxel
            if filterj(j)
                if sqrdist ≤ voxelcutoff_sqr
                    # load other element
                    if exists[j]
                        local y = group[j]
                        # this case should only happen for long line segments or big triangles
                        # The check ensure no double counting
                        if !is_one_voxel
                            @inline local d2 = distSqr(x, y)
                            if d2 ≤ cutoff_sqr
                                if gridpoint == _canonicalGridPoint(x, y, m.grid_size, extra)
                                    local in_y = (y .* m.voxel_length) .+ m.grid_start
                                    @inline output = f(in_x, in_y, i, j, d2*voxel_length_sqr, output)
                                end
                            end
                        else
                            @inline local d2 = distSqr(x, y)
                            if d2 ≤ cutoff_sqr
                                local in_y = (y .* m.voxel_length) .+ m.grid_start
                                @inline output = f(in_x, in_y, i, j, d2*voxel_length_sqr, output)
                            end
                        end
                    end
                end
            end
        end
        return output
    end
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