# Naive implementation using double for loops

using StaticArrays

struct Naive
    "points indexed by group id then order added"
    points::Vector{Vector{SVector{1, SVector{3, Float32}}}}

    "true if point exists"
    pointexists::Vector{Vector{Bool}}

    "lines indexed by group id then order added"
    lines::Vector{Vector{SVector{2, SVector{3, Float32}}}}

    "true if line exists"
    lineexists::Vector{Vector{Bool}}
end


function Naive(numpointgroups::Integer, numlinegroups::Integer ;kw...)
    Naive(
        [[] for i in 1:numpointgroups],
        [[] for i in 1:numpointgroups],
        [[] for i in 1:numlinegroups],
        [[] for i in 1:numlinegroups],
    )
end


"""
Set the elements stored in the cell list.
"""
function setElements(m::Naive, points, lines)
    @argcheck length(points) == length(m.points)
    @argcheck length(lines) == length(m.lines)
    for groupid in 1:length(points)
        in_group = points[groupid]
        group = m.points[groupid]
        resize!(group, length(in_group))
        resize!(m.pointexists[groupid], length(in_group))
        group .= in_group
        m.pointexists[groupid] .= true
    end

    for groupid in 1:length(lines)
        in_group = lines[groupid]
        group = m.lines[groupid]
        resize!(group, length(in_group))
        resize!(m.lineexists[groupid], length(in_group))
        group .= in_group
        m.lineexists[groupid] .= true
    end
end


"""
add a point, return the point id.
"""
function addPoint(m::Naive, groupid::Integer, point)::Int64
    @argcheck length(point) == 1
    @argcheck length(point[begin]) == 3
    points = push!(m.points[groupid],point)
    push!(m.pointexists[groupid], true)
    return length(points)
end

"""
add a line, return the line id.
"""
function addLine(m::Naive, groupid::Integer, line)::Int64
    @argcheck length(line) == 2
    @argcheck length(line[begin]) == 3
    @argcheck length(line[begin+1]) == 3
    lines = push!(m.lines[groupid],line)
    push!(m.lineexists[groupid], true)
    return length(lines)
end


"""
delete a point, the other element ids are stable.
"""
function deletePoint(m::Naive, groupid::Integer, pointid::Integer)::Nothing
    m.pointexists[groupid] = false
    return
end

"""
delete a line, the other element ids are stable.
"""
function deleteLine(m::Naive, groupid::Integer, lineid::Integer)::Nothing
    m.lineexists[groupid] = false
    return
end

"""
map a function to all points in `groupid` in range of point

The function f should have the same form as used in CellListMap.jl
Except here `x` and `y` are `SVector{1, SVector{3, Float32}}`

`x` is always `point` converted to 32 bit float and i is always 0.

    function f(x,y,i,j,d2,output)
        # update output
        return output
    end
"""
function mapPointPoints!(f, output, m::Naive, groupid::Integer, point, cutoff_sqr::Float32)
    # convert point to 32bit float form
    x = convert(SVector{1, SVector{3, Float32}}, point)
    # just loop through all points in groupid
    points = m.points[groupid]
    pointexists = m.pointexists[groupid]
    for (j, y) in enumerate(points)
        if pointexists[j]
            @inline d2 = dist2PointPoint(x, y)
            if d2 ≤ cutoff_sqr
                @inline output = f(x, y, 0, j, d2, output)
            end
        end
    end
end



"""
map a function to all line segments in `groupid` in range of point

The function f should have the same form as used in CellListMap.jl
Except here `x::SVector{1, SVector{3, Float32}}` and `y::SVector{2, SVector{3, Float32}}`

`x` is always `point` converted to 32 bit float and i is always 0.

    function f(x,y,i,j,d2,output)
        # update output
        return output
    end
"""
function mapPointLines!(f, output, m::Naive, groupid::Integer, point, cutoff_sqr::Float32)
    # convert point to 32bit float form
    x = convert(SVector{1, SVector{3, Float32}}, point)
    # just loop through all lines in groupid
    lines = m.lines[groupid]
    lineexists = m.lineexists[groupid]
    for (j, y) in enumerate(lines)
        if lineexists[j]
            @inline d2 = dist2PointLine(x, y)
            if d2 ≤ cutoff_sqr
                @inline output = f(x, y, 0, j, d2, output)
            end
        end
    end
end