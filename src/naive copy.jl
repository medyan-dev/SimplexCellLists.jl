# Naive implementation using double for loops

using StaticArrays

const Element{N} = SVector{N, SVector{3, Float32}}
const Point = Element{1}
const Line = Element{2}


distSqr(x::Point, y::Point) = dist2PointPoint(x, y)
distSqr(x::Point, y::Line)  = dist2PointLine( x, y)
distSqr(x::Line,  y::Point) = dist2PointLine( y, x)
distSqr(x::Line,  y::Line)  = dist2LineLine(  x, y)

struct Naive
    data::Tuple{Vector{Vector{Point}}, Vector{Vector{Line}}}

    "true if element exists"
    exists::Tuple{Vector{Vector{Bool}}, Vector{Vector{Bool}}}
end


function Naive(numpointgroups::Integer, numlinegroups::Integer ;kw...)
    Naive(
        ([Point[] for i in 1:numpointgroups], [Line[] for i in 1:numlinegroups]),
        ([Bool[] for i in 1:numpointgroups],  [Bool[] for i in 1:numlinegroups]),
    )
end


"""
Set the elements stored in the cell list.
"""
function setElements(m::Naive, points, lines)
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
function addElement(m::Naive, groupid::Integer, element::Element{N}) where {N} ::Int64
    group = push!(m.data[N][groupid],element)
    push!(m.exists[N][groupid], true)
    return length(group)
end

"""
delete an element, the other element ids are stable.
"""
function deleteElement(m::Naive, groupid::Integer, elementid::Integer, elementtype::Type{Element{N}}) where {N} ::Nothing
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
function mapElementElements!(f, output, m::Naive, groupid::Integer, element::Element{N}, elementstype::Type{Element{M}}, cutoff_sqr::Float32) where {N, M}
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
function mapPairElements!(f, output, m::Naive, groupid::Integer, elementstype::Type{Element{N}}, cutoff_sqr::Float32) where {N}
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
        m::Naive, 
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