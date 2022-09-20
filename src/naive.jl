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


function addElement!(s::Naive, group_idx::Integer, element::Simplex{N})::Int32 where {N}
    group = push!(s.data[N][group_idx],element)
    push!(s.exists[N][group_idx], true)
    return length(group)
end


function deactivate!(s::Naive, group_idx::Integer, element_idx::Integer, elementtype::Type{Simplex{N}})::Nothing where {N}
    s.exists[N][group_idx][element_idx] = false
    return
end


function activate!(s::Naive, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Nothing where {N}
    s.exists[N][group_idx][element_idx] = true
    return
end

 
function isActive(s::Naive, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Bool where {N}
    s.exists[N][group_idx][element_idx]
end


function mapSimplexElements(f, output, s::Naive, x::Simplex{N}, group_idx::Integer, elements_type::Type{Simplex{M}}, cutoff::Float32) where {N, M}
    # just loop through all element in groupid
    group = s.data[M][group_idx]
    exists = s.exists[M][group_idx]
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


function mapPairElements(f, output, s::Naive, group_idx::Integer, elements_type::Type{Simplex{N}}, cutoff::Float32) where {N}
    # just double loop through all element in groupid
    cutoff_sqr = cutoff^2
    group = s.data[N][group_idx]
    exists = s.exists[N][group_idx]
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


function mapElementsElements(
        f, 
        output, 
        s::Naive, 
        x_group_idx::Integer, 
        x_type::Type{Simplex{N}}, 
        y_group_idx::Integer, 
        y_type::Type{Simplex{M}}, 
        cutoff::Float32,
    ) where {N, M}
    cutoff_sqr = cutoff^2
    # just double loop through all elements
    x_group = s.data[N][x_group_idx]
    x_exists = s.exists[N][x_group_idx]
    y_group = s.data[M][y_group_idx]
    y_exists = s.exists[M][y_group_idx]
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