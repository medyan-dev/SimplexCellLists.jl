const MAX_SAMPLE_SPACING = Float32(1//8)

"""$PUBLIC
A grid is painted with element ids based on a max range. Based on the ideas in "cytosim".

    Painter(numpointgroups::Integer, numlinegroups::Integer, numtrianglegroups::Integer;
        grid_start::SVector{3,<:AbstractFloat},
        grid_size::SVector{3,<:Integer},
        voxel_length::AbstractFloat,
        max_range::Union{SVector{3, <:Vector{<:AbstractFloat}}, AbstractFloat},
    )

Where `numpointgroups` is the number of groups of points, `numlinegroups` is the number of groups of lines,
and `numtrianglegroups` is the number of groups of triangles.

The grid is in a box with one corner at `grid_start .- voxel_length/2` 
and another corner at `(grid_start .- voxel_length/2) + grid_size*voxel_length`

`grid_size` is the number of cubic voxels in each dimension and `voxel_length` is the side length of each voxel.

Elements can safely be outside of the grid, but a large number of elements outside the grid will degrade performance.

`max_range` is the max cutoff that can be used when using the mapping functions. It is specified per group, or for all groups.

When `mapElementsElements` is used, the cutoff just has to be less than the `max_range` of the `y` group.

This implementation is based on the idea in https://gitlab.com/f-nedelec/cytosim/-/blob/master/src/sim/fiber_grid.h
And in figure 10 of "Francois Nedelec and Dietrich Foethke 2007 New J. Phys. 9 427"
to paint a grid with elements based on a maximum range.
"""
struct Painter <: SimplexCellList
    grid_start::SVector{3,Float32}

    grid_size::SVector{3,Int32}

    voxel_length::Float32

    inv_voxel_length::Float32

    "Max ranges of each group in grid units"
    max_range::SVector{3, Vector{Float32}}

    "indexed by shape id then (groupid,x,y,z) to get a vector of elements painted"
    grids::SVector{3, Array{Vector{Tuple{Int32, Float32}}, 4}}

    "A place to put stuff that is outside the grid
    indexed by shape id then groupid to get a vector of elements painted outside the grid"
    outsidegrids::SVector{3, Vector{Vector{Tuple{Int32, Float32}}}}

    "in grid units"
    data::Tuple{Vector{Vector{Point}}, Vector{Vector{Line}}, Vector{Vector{Triangle}}}

    "true if element exists"
    exists::Tuple{Vector{Vector{Bool}}, Vector{Vector{Bool}}, Vector{Vector{Bool}}}
end


function Painter(numpointgroups::Integer, numlinegroups::Integer, numtrianglegroups::Integer;
        grid_start::SVector{3,<:AbstractFloat},
        grid_size::SVector{3,<:Integer},
        voxel_length::AbstractFloat,
        max_range::Union{SVector{3, <:Vector{<:AbstractFloat}}, AbstractFloat},
    )
    @argcheck voxel_length > 0
    _max_range = if max_range isa AbstractFloat
        [fill(max_range,n) for n in (numpointgroups, numlinegroups, numtrianglegroups)]
    else
        @argcheck length(max_range[1]) == numpointgroups
        @argcheck length(max_range[2]) == numlinegroups
        @argcheck length(max_range[3]) == numtrianglegroups
        max_range
    end
    pointgrid = [Tuple{Int32, Float32}[] for i in 1:numpointgroups, j in 1:grid_size[1], k in 1:grid_size[2], l in 1:grid_size[3]]
    linegrid = [Tuple{Int32, Float32}[] for i in 1:numlinegroups, j in 1:grid_size[1], k in 1:grid_size[2], l in 1:grid_size[3]]
    trianglegrid = [Tuple{Int32, Float32}[] for i in 1:numtrianglegroups, j in 1:grid_size[1], k in 1:grid_size[2], l in 1:grid_size[3]]
    Painter(
        grid_start,
        grid_size,
        voxel_length,
        inv(voxel_length),
        _max_range .* inv(voxel_length),
        SA[pointgrid, linegrid, trianglegrid],
        SA[[Tuple{Int32, Float32}[] for i in 1:numpointgroups], [Tuple{Int32, Float32}[] for i in 1:numlinegroups], [Tuple{Int32, Float32}[] for i in 1:numtrianglegroups]], 
        ([Point[] for i in 1:numpointgroups], [Line[] for i in 1:numlinegroups], [Triangle[] for i in 1:numtrianglegroups]),
        ([Bool[] for i in 1:numpointgroups], [Bool[] for i in 1:numlinegroups], [Bool[] for i in 1:numtrianglegroups]),
    )
end


"""$PUBLIC
Alternate constructor.
    Painter(numpointgroups::Integer, numlinegroups::Integer, numtrianglegroups::Integer;
        min_point::SVector{3,<:AbstractFloat},
        max_point::SVector{3,<:AbstractFloat},
        voxel_length::AbstractFloat,
        max_range::SVector{3, <:Vector{<:AbstractFloat}},
    )
"""
function Painter(numpointgroups::Integer, numlinegroups::Integer, numtrianglegroups::Integer;
        min_point::SVector{3,<:AbstractFloat},
        max_point::SVector{3,<:AbstractFloat},
        voxel_length::AbstractFloat,
        max_range::SVector{3, <:Vector{<:AbstractFloat}},
    )
    @argcheck all(max_point .≥ min_point)
    grid_size = max.(ceil.(Int32, (max_point-min_point) ./ voxel_length), Int32(1))
    real_grid_size = grid_size .* voxel_length
    grid_center = (1//2)*(max_point + min_point)
    grid_start = grid_center - (1//2)*real_grid_size .+ (1//2)*voxel_length
    Painter(numpointgroups, numlinegroups, numtrianglegroups;
        grid_start,
        grid_size,
        voxel_length,
        max_range,
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
function _paintElement(s::Painter, groupid, i, element::Simplex{N}, max_range) where {N}
    grid = s.grids[N]
    outsidegrid = s.outsidegrids[N]
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
            if any(voxel_int .> (s.grid_size .- Int32(1))) || any(voxel_int .< Int32(0))
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


function setElements!(s::Painter, points, lines, triangles)::Nothing
    @argcheck length(points) == length(s.data[1])
    @argcheck length(lines) == length(s.data[2])
    @argcheck length(triangles) == length(s.data[3])
    in_data_tuple = (points, lines, triangles)
    Base.Cartesian.@nexprs 3 n -> begin
        in_data = in_data_tuple[n]
        for (groupid, in_group) in enumerate(in_data)
            group = s.data[n][groupid]
            groupexists = s.exists[n][groupid]
            resize!(group, length(in_group))
            resize!(groupexists, length(in_group))
            i = 0
            for stuff in in_group
                i += 1
                group[i] = (stuff .- (s.grid_start,)) .* s.inv_voxel_length
            end
            groupexists .= true
        end
    end
    Base.Cartesian.@nexprs 3 n -> begin
        for (groupid, group) in enumerate(s.data[n])
            max_range = s.max_range[n][groupid]
            for (i, element) in enumerate(group)
                _paintElement(s, groupid, i, element, max_range)
            end
        end
    end
end


function addElement!(s::Painter, group_idx::Integer, element::Simplex{N})::Int32 where {N}
    x = (element .- (s.grid_start,)) .* s.inv_voxel_length
    group = push!(s.data[N][group_idx], x)
    push!(s.exists[N][group_idx], true)
    i = length(group)
    _paintElement(s, group_idx, i, x, s.max_range[N][group_idx])
    return i
end


function deactivate!(s::Painter, group_idx::Integer, element_idx::Integer, elementtype::Type{Simplex{N}})::Nothing where {N}
    s.exists[N][group_idx][element_idx] = false
    return
end


function activate!(s::Painter, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Nothing where {N}
    s.exists[N][group_idx][element_idx] = true
    return
end

 
function isActive(s::Painter, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Bool where {N}
    s.exists[N][group_idx][element_idx]
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
    elseif N == 3
        if isnothing(extra)
            return _nearestGridPoint(grid_size, 1//3*(x[3] + x[2] + x[1]))
        else
            d2, smin = distSqr_sMin(x,y)
            pt = (1-(smin[1]+smin[2])).*x[1] + smin[1].*x[2] + smin[2].*x[3]
            return _nearestGridPoint(grid_size, pt)
        end
    else
        error("simplexes with $N points not supported")
    end
end


function mapSimplexElements(
        f,
        output,
        s::Painter,
        in_x::Simplex{N},
        group_idx::Integer,
        elements_type::Type{Simplex{M}},
        in_cutoff::Float32,
        ;
        x = (in_x .- (s.grid_start,)) .* s.inv_voxel_length,
        i::Int32=Int32(0),
        filterj=Returns(true),
    ) where {N, M}
    error("Not implemented for those simplexes")
end

#point to other
function mapSimplexElements(
        f,
        output,
        s::Painter,
        in_x::Simplex{1},
        group_idx::Integer,
        elements_type::Type{Simplex{M}},
        in_cutoff::Float32,
        ;
        x = (in_x .- (s.grid_start,)) .* s.inv_voxel_length,
        i::Int32=Int32(0),
        filterj=Returns(true),
    ) where{M}
    cutoff = in_cutoff * s.inv_voxel_length
    @argcheck cutoff ≤ s.max_range[M][group_idx]
    group = s.data[M][group_idx]
    exists = s.exists[M][group_idx]
    grid = s.grids[M]
    outsidegrid = s.outsidegrids[M]
    voxelcutoff::Float32 = Float32(1//2*√(3)) + cutoff + MAX_SAMPLE_SPACING
    voxelcutoff_sqr = voxelcutoff^2
    cutoff_sqr = cutoff^2
    voxel_length_sqr = s.voxel_length^2
    # sample unique voxels that the element goes in based on sampling points that are at most 1/8 voxel_length away from any point in the element.
    gridpoint = _nearestGridPoint(s.grid_size, x[1])
    voxel = if isnothing(gridpoint)
        outsidegrid[group_idx]
    else
        voxelid = gridpoint .+ Int32(1)
        grid[group_idx, voxelid...]
    end
    for (j, sqrdist) in voxel
        if filterj(j)
            if sqrdist ≤ voxelcutoff_sqr
                # load other element
                if exists[j]
                    y = group[j]
                    @inline d2 = distSqr(x, y)
                    if d2 ≤ cutoff_sqr
                        in_y = (y .* s.voxel_length) .+ (s.grid_start,)
                        output = f(in_x, in_y, i, j, d2*voxel_length_sqr, output)
                    end
                end
            end
        end
    end
    output
end

#line to other
function mapSimplexElements(
        f,
        output,
        s::Painter,
        in_x::Simplex{2},
        group_idx::Integer,
        elements_type::Type{Simplex{M}},
        in_cutoff::Float32,
        ;
        x = (in_x .- (s.grid_start,)) .* s.inv_voxel_length,
        i::Int32=Int32(0),
        filterj=Returns(true),
    ) where {M}
    cutoff = in_cutoff * s.inv_voxel_length
    @argcheck cutoff ≤ s.max_range[M][group_idx]
    group = s.data[M][group_idx]
    exists = s.exists[M][group_idx]
    grid = s.grids[M]
    outsidegrid = s.outsidegrids[M]
    voxelcutoff::Float32 = Float32(1//2*√(3)) + cutoff + MAX_SAMPLE_SPACING
    voxelcutoff_sqr = voxelcutoff^2
    cutoff_sqr = cutoff^2
    voxel_length_sqr = s.voxel_length^2
    # sample unique voxels that the element goes in based on sampling points that are at most 1/8 voxel_length away from any point in the element.
    # if the line segment is short enough just us the mid point
    r = x[2] - x[1]
    len2 = r ⋅ r
    midgridpoint = _nearestGridPoint(s.grid_size, 1//2*(x[2] + x[1]))
    usemidpoint = false
    if len2 < (MAX_SAMPLE_SPACING*2)^2
        usemidpoint = true
    else
        # if both ends, moved in by MAX_SAMPLE_SPACING are in the same voxel, just sample that voxel.
        len = √len2
        small_s = MAX_SAMPLE_SPACING/len
        start_pt = (1-small_s)*x[1] + (small_s)*x[2]
        end_pt = (1-small_s)*x[2] + (small_s)*x[1]
        gridpoint1 = _nearestGridPoint(s.grid_size, start_pt)
        gridpoint2 = _nearestGridPoint(s.grid_size, end_pt)
        if gridpoint1 == gridpoint2 && !isnothing(gridpoint1)
            usemidpoint = true
        end
    end
    if usemidpoint
        voxel = if isnothing(midgridpoint)
            outsidegrid[group_idx]
        else
            voxelid = midgridpoint .+ Int32(1)
            grid[group_idx, voxelid...]
        end
        for (j, sqrdist) in voxel
            if filterj(j)
                if sqrdist ≤ voxelcutoff_sqr
                    # load other element
                    if exists[j]
                        y = group[j]
                        @inline d2 = distSqr(x, y)
                        if d2 ≤ cutoff_sqr
                            in_y = (y .* s.voxel_length) .+ (s.grid_start,)
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
            gridpoint = _getSampleGridPoint(inv_n, n, samplei, x, s.grid_size)
            if gridpoint != prevgridpoint && !(checked_outside && isnothing(gridpoint))
                voxel = if isnothing(gridpoint)
                    checked_outside = true
                    outsidegrid[group_idx]
                else
                    voxelid = gridpoint .+ Int32(1)
                    grid[group_idx, voxelid...]
                end
                for (j, sqrdist) in voxel
                    if filterj(j)
                        if sqrdist ≤ voxelcutoff_sqr
                            # load other element
                            if exists[j]
                                y = group[j]
                                @inline d2 = distSqr(x, y)
                                if d2 ≤ cutoff_sqr 
                                    c_gridpoint = _canonicalGridPoint(x, y, s.grid_size, extra)
                                    if gridpoint == c_gridpoint
                                        in_y = (y .* s.voxel_length) .+ (s.grid_start,)
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

#triangle to other
function mapSimplexElements(
        f,
        output,
        s::Painter,
        in_x::Simplex{3},
        group_idx::Integer,
        elements_type::Type{Simplex{M}},
        in_cutoff::Float32,
        ;
        x = (in_x .- (s.grid_start,)) .* s.inv_voxel_length,
        i::Int32=Int32(0),
        filterj=Returns(true),
    ) where {M}
    cutoff = in_cutoff * s.inv_voxel_length
    @argcheck cutoff ≤ s.max_range[M][group_idx]
    group = s.data[M][group_idx]
    exists = s.exists[M][group_idx]
    grid = s.grids[M]
    outsidegrid = s.outsidegrids[M]
    voxelcutoff::Float32 = Float32(1//2*√(3)) + cutoff + MAX_SAMPLE_SPACING
    voxelcutoff_sqr = voxelcutoff^2
    cutoff_sqr = cutoff^2
    voxel_length_sqr = s.voxel_length^2
    # sample unique voxels that the element goes in based on sampling points that are at most 1/8 voxel_length away from any point in the element.
    # if the line segment is short enough just us the mid point
    mid_pt = 1//3*(x[3] + x[2] + x[1])
    # midgridpoint is nothing if it is outside the grid
    midgridpoint = _nearestGridPoint(s.grid_size, 1//3*(x[3] + x[2] + x[1]))
    usemidpoint = false
    rs = x .- Ref(mid_pt)
    ds = norm.(rs)
    if all(<(MAX_SAMPLE_SPACING), ds)
        usemidpoint = true
    elseif !isnothing(midgridpoint)
        voxel_rs = x .- Ref(midgridpoint)
        voxel_ds = norm.(voxel_rs)
        # if all ends, moved toward the center of the mid point voxel by MAX_SAMPLE_SPACING are in the same voxel, just sample that voxel.
        # TODO move the points into a single voxel more efficently.
        move_ins = min.(MAX_SAMPLE_SPACING ./ voxel_ds, 1.0f0)
        moved_in_rs =  (1 .- move_ins) .* voxel_rs
        moved_in_x = Ref(midgridpoint) .+ moved_in_rs
        gridpoints = _nearestGridPoint.(Ref(s.grid_size), moved_in_x)
        if allequal((gridpoints..., midgridpoint))
            usemidpoint = true
        end
    end
    # @show usemidpoint
    if usemidpoint
        voxel = if isnothing(midgridpoint)
            outsidegrid[group_idx]
        else
            voxelid = midgridpoint .+ Int32(1)
            grid[group_idx, voxelid...]
        end
        for (j, sqrdist) in voxel
            if filterj(j)
                if sqrdist ≤ voxelcutoff_sqr
                    # load other element
                    if exists[j]
                        y = group[j]
                        @inline d2 = distSqr(x, y)
                        if d2 ≤ cutoff_sqr
                            in_y = (y .* s.voxel_length) .+ (s.grid_start,)
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
        # start by getting a axis aligned bounding box around x
        # we need to sample any voxels that x touches or could touch 
        # given floating point rounding.
        bottom = min.(x...)
        top = max.(x...)
        maxabs = max.(abs.(top),abs.(bottom))
        error_est = 32*eps(maximum(maxabs))
        top_extend_int = round.(Int32, top .+ error_est)
        bottom_extend_int = round.(Int32, bottom .- error_est)
        check_outside = false
        if any(top_extend_int .> (s.grid_size .- Int32(1))) || any(bottom_extend_int .< Int32(0))
            check_outside = true
            top_extend_int = clamp.(top_extend_int, Int32(0), (s.grid_size .- Int32(1)))
            bottom_extend_int = clamp.(bottom_extend_int, Int32(0), (s.grid_size .- Int32(1)))
        end
        samplevoxelcutoff = Float32(1//2*√(3)) + error_est + MAX_SAMPLE_SPACING
        samplevoxelcutoff2 = samplevoxelcutoff^2
        for gridpoint_tuple in Iterators.flatten(((nothing,), Iterators.product(range.(bottom_extend_int,top_extend_int)...)))
            gridpoint = isnothing(gridpoint_tuple) ? nothing : SVector(gridpoint_tuple)
            # if gridpoint is outside but we don't need to check outside, skip
            if isnothing(gridpoint) && !check_outside
                continue
            end
            # continue if the triangle could not possibly be in the voxel
            if !isnothing(gridpoint)
                @inline voxel_float = SVector{3,Float32}(gridpoint...)
                @inline sample_voxel_d2 = distSqr(SA[voxel_float,], x)
                if sample_voxel_d2 > samplevoxelcutoff2
                    continue
                end
            end
            # at this point we want to sample gridpoint
            voxel = if isnothing(gridpoint)
                outsidegrid[group_idx]
            else
                voxelid = gridpoint .+ Int32(1)
                grid[group_idx, voxelid...]
            end
            for (j, sqrdist) in voxel
                if filterj(j)
                    if sqrdist ≤ voxelcutoff_sqr
                        # load other element
                        if exists[j]
                            y = group[j]
                            @inline d2 = distSqr(x, y)
                            if d2 ≤ cutoff_sqr 
                                c_gridpoint = _canonicalGridPoint(x, y, s.grid_size, true)
                                if gridpoint == c_gridpoint
                                    in_y = (y .* s.voxel_length) .+ (s.grid_start,)
                                    @inline output = f(in_x, in_y, i, j, d2*voxel_length_sqr, output)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return output
end


function mapPairElements(
        f,
        output,
        s::Painter,
        group_idx::Integer,
        elements_type::Type{Simplex{N}},
        cutoff::Float32,
    ) where {N}
    # loop through all element in groupid
    # and call mapSimplexElements with extra optional parameters
    # to only run on elements with id > i, to avoid double counting
    group = s.data[N][group_idx]
    exists = s.exists[N][group_idx]
    n = length(group)
    for i in Int32(1):Int32(n-1)
        if exists[i]
            x = group[i]
            in_x = (x .* s.voxel_length) .+ (s.grid_start,)
            output = mapSimplexElements(f, output, s, in_x, group_idx, elements_type, cutoff;
                i,
                filterj = >(i),
                x,
            )
        end
    end
    return output
end


function mapElementsElements(
        f,
        output,
        s::Painter,
        x_group_idx::Integer,
        x_type::Type{Simplex{N}},
        y_group_idx::Integer,
        y_type::Type{Simplex{M}},
        cutoff::Float32,
    ) where {N, M}
    # loop through all element in groupid
    # and call mapSimplexElements with extra optional parameters
    x_group = s.data[N][x_group_idx]
    x_exists = s.exists[N][x_group_idx]
    x_n = length(x_group)
    for i in Int32(1):Int32(x_n)
        if x_exists[i]
            x = x_group[i]
            in_x = (x .* s.voxel_length) .+ (s.grid_start,)
            output = mapSimplexElements(f, output, s, in_x, y_group_idx, y_type, cutoff;
                i,
                x,
            )
        end
    end
    return output
end