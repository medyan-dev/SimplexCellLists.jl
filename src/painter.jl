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

    "Max ranges of each group"
    maxrange::SVector{2, Vector{Float32}}

    "indexed by shape id then (groupid,x,y,z) then element id"
    grid::SVector{2, Array{Vector{Tuple{Int32, Float32}}, 4}}

    data::Tuple{Vector{Vector{Point}}, Vector{Vector{Line}}}

    "true if element exists"
    exists::Tuple{Vector{Vector{Bool}}, Vector{Vector{Bool}}}
end


function Painter(numpointgroups::Integer, numlinegroups::Integer ;
        grid_start::SVector{3,AbstractFloat},
        grid_size::SVector{3,Integer},
        voxel_length::AbstractFloat,
        maxrange::SVector{2, Vector{AbstractFloat}},
    )
    pointgrid = [Tuple{Int32, Float32}[] for i in 1:numpointgroups, j in 1:grid_size[1], k in 1:grid_size[2], l in 1:grid_size[3]]
    linegrid = [Tuple{Int32, Float32}[] for i in 1:numlinegroups, j in 1:grid_size[1], k in 1:grid_size[2], l in 1:grid_size[3]]
    Painter(
        grid_start,
        grid_size,
        voxel_length,
        inv(voxel_length),
        maxrange,
        SA[pointgrid, linegrid],
        ([Point[] for i in 1:numpointgroups], [Line[] for i in 1:numlinegroups]),
        ([Bool[] for i in 1:numpointgroups],  [Bool[] for i in 1:numlinegroups]),
    )
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
            group .= in_group
            groupexists .= true
        end
    end
end

"""
add a point, return the point id.
"""
function addElement(m::Painter, groupid::Integer, element::Element{N}) where {N} ::Int64
    group = push!(m.data[N][groupid],element)
    push!(m.exists[N][groupid], true)
    return length(group)
end

"""
delete an element, the other element ids are stable.
"""
function deleteElement(m::Painter, groupid::Integer, elementid::Integer, elementtype::Type{Element{N}}) where {N} ::Nothing
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
function mapElementElements!(f, output, m::Painter, groupid::Integer, element::Element{N}, elementstype::Type{Element{M}}, cutoff_sqr::Float32) where {N, M}
    x = element
    # just loop through all element in groupid
    group = m.data[M][groupid]
    exists = m.exists[M][groupid]
    for (j, y) in enumerate(group)
        if exists[j]
            @inline d2 = distSqr(x, y)
            if d2 ≤ cutoff_sqr
                @inline output = f(x, y, 0, j, d2, output)
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
function mapPairElements!(f, output, m::Painter, groupid::Integer, elementstype::Type{Element{N}}, cutoff_sqr::Float32) where {N}
    # just double loop through all element in groupid
    group = m.data[N][groupid]
    exists = m.exists[N][groupid]
    n = length(group)
    for i in 1:(n-1)
        if exists[i]
            x = group[i]
            for j in (i+1):n
                if exists[j]
                    y = group[j]
                    @inline d2 = distSqr(x, y)
                    if d2 ≤ cutoff_sqr
                        @inline output = f(x, y, i, j, d2, output)
                    end
                end
            end
        end
    end
    return output
end


"""
map a function to all pairs of elements in different groups in range of each other.

The function f should have the same form as used in CellListMap.jl
Except here `x` and `y` are `SVector{N, SVector{3, Float32}}`, `SVector{M, SVector{3, Float32}}`

    function f(x,y,i,j,d2,output)
        # update output
        return output
    end
"""
function mapPairElementsElements!(
        f, 
        output, 
        m::Painter, 
        x_groupid::Integer, 
        x_elementstype::Type{Element{N}}, 
        y_groupid::Integer, 
        y_elementstype::Type{Element{M}}, 
        cutoff_sqr::Float32,
    ) where {N, M}
    # just double loop through all element in groupid
    x_group = m.data[N][x_groupid]
    x_exists = m.exists[N][x_groupid]
    y_group = m.data[M][y_groupid]
    y_exists = m.exists[M][y_groupid]
    xn = length(x_group)
    yn = length(y_group)
    for i in 1:xn
        if x_exists[i]
            x = x_group[i]
            for j in 1:yn
                if y_exists[j]
                    y = y_group[j]
                    @inline d2 = distSqr(x, y)
                    if d2 ≤ cutoff_sqr
                        @inline output = f(x, y, i, j, d2, output)
                    end
                end
            end
        end
    end
    return output
end