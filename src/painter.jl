const MAX_SAMPLE_SPACING = Float32(1//8)

"""
This implementation is based on the idea in https://gitlab.com/f-nedelec/cytosim/-/blob/master/src/sim/fiber_grid.h
And in figure 10 of "Francois Nedelec and Dietrich Foethke 2007 New J. Phys. 9 427"
to paint a grid with elements based on a maximum range.
"""
struct Painter <: MultiShapeCellList
    grid_start::SVector{3,Float32}

    grid_size::SVector{3,Int32}

    voxel_length::Float32

    inv_voxel_length::Float32

    "Max ranges of each group in grid units"
    max_range::SVector{2, Vector{Float32}}

    "indexed by shape id then (groupid,x,y,z) to get a vector of elements painted"
    grids::SVector{2, Array{Vector{Tuple{Int32, Float32}}, 4}}

    "A place to put stuff that is outside the grid
    indexed by shape id then groupid to get a vector of elements painted outside the grid"
    outsidegrids::SVector{2, Vector{Vector{Tuple{Int32, Float32}}}}

    "in grid units"
    data::Tuple{Vector{Vector{Point}}, Vector{Vector{Line}}}

    "true if element exists"
    exists::Tuple{Vector{Vector{Bool}}, Vector{Vector{Bool}}}
end


function Painter(numpointgroups::Integer, numlinegroups::Integer ;
        grid_start::SVector{3,<:AbstractFloat},
        grid_size::SVector{3,<:Integer},
        voxel_length::AbstractFloat,
        max_range::SVector{2, <:Vector{<:AbstractFloat}},
    )
    @argcheck length(max_range[1]) == numpointgroups
    @argcheck length(max_range[2]) == numlinegroups
    pointgrid = [Tuple{Int32, Float32}[] for i in 1:numpointgroups, j in 1:grid_size[1], k in 1:grid_size[2], l in 1:grid_size[3]]
    linegrid = [Tuple{Int32, Float32}[] for i in 1:numlinegroups, j in 1:grid_size[1], k in 1:grid_size[2], l in 1:grid_size[3]]
    Painter(
        grid_start,
        grid_size,
        voxel_length,
        inv(voxel_length),
        max_range .* inv(voxel_length),
        SA[pointgrid, linegrid],
        SA[[Tuple{Int32, Float32}[] for i in 1:numpointgroups], [Tuple{Int32, Float32}[] for i in 1:numlinegroups]], 
        ([Point[] for i in 1:numpointgroups], [Line[] for i in 1:numlinegroups]),
        ([Bool[] for i in 1:numpointgroups],  [Bool[] for i in 1:numlinegroups]),
    )
end

"""
Return nearest grid point, x is in grid units
Return nothing if outside of the grid
"""
function _nearestGridPoint(grid_size::SVector{3,Int32}, x::SVector{3,Float32})
    x_int = round.(Int32, x)
    if any(x_int .> (grid_size .- Int32(1))) || any(x_int .< Int32(0))
        nothing
    else
        x_int
    end
end

"""
Paint an element onto the grid.
All inputs have grid units.
Max range is extended by `MAX_SAMPLE_SPACING` because of how elements are sampled. 
Any voxel that may contain points within `max_range + MAX_SAMPLE_SPACING` of `element`
will have `i`, and the distance squared, added to that voxel's list.

In grid units, voxel (0,0,0) is a unit cube with center at SA[0.0, 0.0, 0.0]
"""
function _paintElement(m::Painter, groupid, i, element::Simplex{N}, max_range) where {N}
    grid = m.grids[N]
    outsidegrid = m.outsidegrids[N]
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
    bottom_int = round.(Int32, bottom)
    top_int = round.(Int32, top)
    # Go through all voxel center points and check distance to element and paint if in range.
    min_outside_d2 = Inf32
    for voxel_int in Iterators.product(range.(bottom_int,top_int)...)
        # @show voxel_int
        @inline voxel_float = SVector{3,Float32}(voxel_int...)
        @inline d2 = distSqr(SA[voxel_float,], element)
        if d2 ≤ voxelcutoff2
            if any(voxel_int .> (m.grid_size .- Int32(1))) || any(voxel_int .< Int32(0))
                min_outside_d2 = min(min_outside_d2, d2)
            else
                voxelid = voxel_int .+ Int32(1)
                voxel = grid[groupid, voxelid...]
                push!(voxel, (i,d2))
            end
        end
    end
    if isfinite(min_outside_d2)
        push!(outsidegrid[groupid],(i,min_outside_d2))
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
            i = 0
            for stuff in in_group
                i += 1
                group[i] = (stuff .- (m.grid_start,)) .* m.inv_voxel_length
            end
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
add a simplex, return the simplex id.
"""
function addElement(m::Painter, groupid::Integer, element::Simplex{N})::Int32 where {N}
    x = (element .- (m.grid_start,)) .* m.inv_voxel_length
    group = push!(m.data[N][groupid], x)
    push!(m.exists[N][groupid], true)
    i = length(group)
    _paintElement(m, groupid, i, x, m.max_range[N][groupid])
    return i
end

"""
delete an element, the other element ids are stable.
"""
function deleteElement(m::Painter, groupid::Integer, elementid::Integer, elementtype::Type{Simplex{N}})::Nothing where {N} 
    m.exists[N][groupid][elementid] = false
    return
end


"""
This function is to avoid any floating point weirdness
"""
@noinline function _getSampleGridPoint(inv_n::Float32, n::Int32, R::Int32, x::Simplex{2}, grid_size)
    Q::Int32 = n-R
    p1 = Q*x[1] + R*x[2]
    p2 = inv_n*p1
    _nearestGridPoint(grid_size, p2)
end

"""
Get the grid point that is sampled by x, to y.
x and y are in grid units
"""
@noinline function _canonicalGridPoint(x::Simplex{N}, y::Simplex{M}, grid_size, extra)::Union{Nothing,SVector{3, Int32}} where {N,M} 
    if N == 1
        return _nearestGridPoint(grid_size, x[1])
    elseif N == 2
        if isnothing(extra)
            return _nearestGridPoint(grid_size, 1//2*(x[2] + x[1]))
        else
            halfn::Int32, inv_n::Float32 = extra
            d2, smin = distSqr_sMin(x,y)
            n = Int32(2)*halfn
            R::Int32 = clamp(trunc(Int32, smin[1]*halfn),Int32(0),halfn-Int32(1))*Int32(2) + Int32(1)
            return _getSampleGridPoint(inv_n, n, R, x, grid_size)
        end
    else
        error("Triangle not implemented yet")
    end
end

"""
map a function to all elements in `groupid` in range of one simplex

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
        x = (in_x .- (m.grid_start,)) .* m.inv_voxel_length,
        i::Int32=Int32(0),
        filterj=Returns(true),
    ) where {N, M}
    error("Not implemented for those simplexes")
end

#point to other
function mapSimplexElements!(
        f, 
        output, 
        m::Painter, 
        groupid::Integer, 
        in_x::Simplex{1}, 
        elementstype::Type{Simplex{M}}, 
        in_cutoff::Float32
        ;
        x = (in_x .- (m.grid_start,)) .* m.inv_voxel_length,
        i::Int32=Int32(0),
        filterj=Returns(true),
    ) where{M}
    cutoff = in_cutoff * m.inv_voxel_length
    @argcheck cutoff ≤ m.max_range[M][groupid]
    group = m.data[M][groupid]
    exists = m.exists[M][groupid]
    grid = m.grids[M]
    outsidegrid = m.outsidegrids[M]
    voxelcutoff::Float32 = Float32(1//2*√(3)) + cutoff + MAX_SAMPLE_SPACING
    voxelcutoff_sqr = voxelcutoff^2
    cutoff_sqr = cutoff^2
    voxel_length_sqr = m.voxel_length^2
    # sample unique voxels that the element goes in based on sampling points that are at most 1/8 voxel_length away from any point in the element.
    gridpoint = _nearestGridPoint(m.grid_size, x[1])
    voxel = if isnothing(gridpoint)
        outsidegrid[groupid]
    else
        voxelid = gridpoint .+ Int32(1)
        grid[groupid, voxelid...]
    end
    for (j, sqrdist) in voxel
        if filterj(j)
            if sqrdist ≤ voxelcutoff_sqr
                # load other element
                if exists[j]
                    y = group[j]
                    @inline d2 = distSqr(x, y)
                    if d2 ≤ cutoff_sqr
                        in_y = (y .* m.voxel_length) .+ (m.grid_start,)
                        output = f(in_x, in_y, i, j, d2*voxel_length_sqr, output)
                    end
                end
            end
        end
    end
    output
end

#line to other
function mapSimplexElements!(
        f, 
        output, 
        m::Painter, 
        groupid::Integer, 
        in_x::Simplex{2}, 
        elementstype::Type{Simplex{M}}, 
        in_cutoff::Float32
        ;
        x = (in_x .- (m.grid_start,)) .* m.inv_voxel_length,
        i::Int32=Int32(0),
        filterj=Returns(true),
    ) where {M}
    cutoff = in_cutoff * m.inv_voxel_length
    @argcheck cutoff ≤ m.max_range[M][groupid]
    group = m.data[M][groupid]
    exists = m.exists[M][groupid]
    grid = m.grids[M]
    outsidegrid = m.outsidegrids[M]
    voxelcutoff::Float32 = Float32(1//2*√(3)) + cutoff + MAX_SAMPLE_SPACING
    voxelcutoff_sqr = voxelcutoff^2
    cutoff_sqr = cutoff^2
    voxel_length_sqr = m.voxel_length^2
    # sample unique voxels that the element goes in based on sampling points that are at most 1/8 voxel_length away from any point in the element.
    # if the line segment is short enough just us the mid point
    r = x[2] - x[1]
    len2 = r ⋅ r
    midgridpoint = _nearestGridPoint(m.grid_size, 1//2*(x[2] + x[1]))
    usemidpoint = false
    if len2 < (MAX_SAMPLE_SPACING*2)^2
        usemidpoint = true
    else
        # if both ends, moved in by MAX_SAMPLE_SPACING are in the same voxel, just sample that voxel.
        len = √len2
        small_s = MAX_SAMPLE_SPACING/len
        start_pt = (1-small_s)*x[1] + (small_s)*x[2]
        end_pt = (1-small_s)*x[2] + (small_s)*x[1]
        gridpoint1 = _nearestGridPoint(m.grid_size, start_pt)
        gridpoint2 = _nearestGridPoint(m.grid_size, end_pt)
        if gridpoint1 == gridpoint2 && !isnothing(gridpoint1)
            usemidpoint = true
        end
    end
    if usemidpoint
        voxel = if isnothing(midgridpoint)
            outsidegrid[groupid]
        else
            voxelid = midgridpoint .+ Int32(1)
            grid[groupid, voxelid...]
        end
        for (j, sqrdist) in voxel
            if filterj(j)
                if sqrdist ≤ voxelcutoff_sqr
                    # load other element
                    if exists[j]
                        y = group[j]
                        @inline d2 = distSqr(x, y)
                        if d2 ≤ cutoff_sqr
                            in_y = (y .* m.voxel_length) .+ (m.grid_start,)
                            output = f(in_x, in_y, i, j, d2*voxel_length_sqr, output)
                        end
                    end
                end
            end
        end
    else
        # Now more than one voxel needs to be sampled, yikes.
        # I need to ensure that _canonicalGridPoint can find a gridpoint that is actually sampled
        # This is to avoid double counting
        num_sub_div = ceil(inv(2*MAX_SAMPLE_SPACING)*len)
        inv_n::Float32 = inv(2num_sub_div)
        halfn = Int32(num_sub_div)
        extra = (halfn, inv_n)
        n = Int32(2)*halfn
        prevgridpoint = SA[Int32(-1),Int32(-1),Int32(-1)]
        checked_outside = false
        for samplei in Int32(1):Int32(2):n
            gridpoint = _getSampleGridPoint(inv_n, n, samplei, x, m.grid_size)
            if gridpoint != prevgridpoint && !(checked_outside && isnothing(gridpoint))
                voxel = if isnothing(gridpoint)
                    checked_outside = true
                    outsidegrid[groupid]
                else
                    voxelid = gridpoint .+ Int32(1)
                    grid[groupid, voxelid...]
                end
                for (j, sqrdist) in voxel
                    if filterj(j)
                        if sqrdist ≤ voxelcutoff_sqr
                            # load other element
                            if exists[j]
                                y = group[j]
                                @inline d2 = distSqr(x, y)
                                if d2 ≤ cutoff_sqr 
                                    c_gridpoint = _canonicalGridPoint(x, y, m.grid_size, extra)
                                    if gridpoint == c_gridpoint
                                        in_y = (y .* m.voxel_length) .+ (m.grid_start,)
                                        @inline output = f(in_x, in_y, i, j, d2*voxel_length_sqr, output)
                                    end
                                end
                            end
                        end
                    end
                end
                prevgridpoint = gridpoint
            end
        end
    end
    return output
end

"""
map a function to all pairs of elements in the same group in range of each other.

The function f should have the same form as used in CellListMap.jl
Except here `x` and `y` are `SVector{N, SVector{3, Float32}}`, `SVector{N, SVector{3, Float32}}`

    function f(x,y,i,j,d2,output)
        # update output
        return output
    end
"""
function mapPairElements!(f, output, m::Painter, groupid::Integer, elementstype::Type{Simplex{N}}, cutoff::Float32) where {N}
    # loop through all element in groupid
    # and call mapSimplexElements! with extra optional parameters
    # to only run on elements with id > i, to avoid double counting
    group = m.data[N][groupid]
    exists = m.exists[N][groupid]
    n = length(group)
    for i in Int32(1):Int32(n-1)
        if exists[i]
            x = group[i]
            in_x = (x .* m.voxel_length) .+ (m.grid_start,)
            output = mapSimplexElements!(f, output, m, groupid, in_x, elementstype, cutoff;
                i,
                filterj = >(i),
                x,
            )
        end
    end
    return output
end


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