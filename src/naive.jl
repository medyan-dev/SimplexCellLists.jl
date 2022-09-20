# Naive implementation using double for loops



"""
A simple `SimplexCellList` with no extra options that uses double for loops.

    Naive(numpointgroups::Integer, numlinegroups::Integer, numtrianglegroups::Integer)

Where `numpointgroups` is the number of groups of points, `numlinegroups` is the number of groups of lines,
and `numtrianglegroups` is the number of groups of triangles.
"""
struct Naive <: SimplexCellList
    data::Tuple{Vector{Vector{Point}}, Vector{Vector{Line}}, Vector{Vector{Triangle}}}

    "true if element exists"
    exists::Tuple{Vector{Vector{Bool}}, Vector{Vector{Bool}}, Vector{Vector{Bool}}}
end


function Naive(numpointgroups::Integer, numlinegroups::Integer, numtrianglegroups::Integer)
    Naive(
        (
            [Point[] for i in 1:numpointgroups],
            [Line[] for i in 1:numlinegroups],
            [Triangle[] for i in 1:numtrianglegroups],
        ),
        (
            [Bool[] for i in 1:numpointgroups],
            [Bool[] for i in 1:numlinegroups],
            [Bool[] for i in 1:numtrianglegroups],
        ),
    )
end


"""
Reset the elements stored in `s` in batch:

```julia
setElements!(s::Naive, points, lines, triangles)::Nothing
```

Where `points`, `lines` and `triangles` are collections of collections of objects convertible to 
`Point`, `Line`, and `Triangle` respectively.

For example, each collection in `points` is a group of points that can be mapped over independently from, or together with, other groups.

Added elements will have a group index and element index based on the order of the inputs.
The first group in each type has group index 1, and the first element in each group has element index 1.
"""
function setElements!(s::Naive, points, lines, triangles)::Nothing
    @argcheck length(points) == length(s.data[1])
    @argcheck length(lines) == length(s.data[2])
    @argcheck length(triangles) == length(s.data[3])
    in_data_tuple = (points, lines, triangles)
    Base.Cartesian.@nexprs 3 n -> begin
        in_data = in_data_tuple[n]
        for (group_idx, in_group) in enumerate(in_data)
            group = s.data[n][group_idx]
            group_exists = s.exists[n][group_idx]
            resize!(group, length(in_group))
            resize!(group_exists, length(in_group))
            group .= in_group
            group_exists .= true
        end
    end
    return
end

"""
Add a new element to `s`, and return its element index:

```julia
addElement!(s::Naive, group_idx::Integer, element::Simplex{N})::Int32
```
The new element will be pushed to the end of the specified group.
"""
function addElement!(s::Naive, group_idx::Integer, element::Simplex{N})::Int32 where {N}
    group = push!(s.data[N][group_idx],element)
    push!(s.exists[N][group_idx], true)
    return length(group)
end

"""
Deactivate an existing element in `s`

```julia
deactivate!(s::Naive, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Nothing
```
Inactive elements are not mapped over. Elements are active by default.
"""
function deactivate!(s::Naive, group_idx::Integer, element_idx::Integer, elementtype::Type{Simplex{N}})::Nothing where {N}
    s.exists[N][group_idx][element_idx] = false
    return
end

"""
Re-activate an existing element in `s`

```julia
activate!(s::Naive, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Nothing
```
Inactive elements are not mapped over. Elements are active by default.
"""
function activate!(s::Naive, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Nothing where {N}
    s.exists[N][group_idx][element_idx] = true
    return
end

"""
Return if an existing element in `s` is active.

```julia
isActive(s::Naive, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Bool
```
Inactive elements are not mapped over. Elements are active by default.
"""    
function activate!(s::Naive, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Bool where {N}
    s.exists[N][group_idx][element_idx] = true
    return
end

"""
map a function to all elements in `groupid` within distance `cutoff` of one simplex

The function f should have the same form as used in CellListMap.jl
Except here `x` and `y` are `SVector{N, SVector{3, Float32}}`, `SVector{M, SVector{3, Float32}}`

`x` is always `x` and i is always 0.

    function f(x,y,i,j,d2,output)
        # update output
        return output
    end
"""
function mapSimplexElements(f, output, m::Naive, groupid::Integer, x::Simplex{N}, elementstype::Type{Simplex{M}}, cutoff::Float32) where {N, M}
    # just loop through all element in groupid
    group = m.data[M][groupid]
    exists = m.exists[M][groupid]
    cutoff_sqr = cutoff^2
    for (j, y) in enumerate(group)
        if exists[j]
            @inline d2 = distSqr(x, y)
            if d2 ≤ cutoff_sqr
                @inline output = f(x, y, 0%Int32, j%Int32, d2, output)
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
function mapPairElements(f, output, m::Naive, groupid::Integer, elementstype::Type{Simplex{N}}, cutoff::Float32) where {N}
    # just double loop through all element in groupid
    cutoff_sqr = cutoff^2
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
                        @inline output = f(x, y, i%Int32, j%Int32, d2, output)
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
function mapElementsElements(
        f, 
        output, 
        m::Naive, 
        x_groupid::Integer, 
        x_type::Type{Simplex{N}}, 
        y_groupid::Integer, 
        y_type::Type{Simplex{M}}, 
        cutoff::Float32,
    ) where {N, M}
    cutoff_sqr = cutoff^2
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
                        @inline output = f(x, y, i%Int32, j%Int32, d2, output)
                    end
                end
            end
        end
    end
    return output
end