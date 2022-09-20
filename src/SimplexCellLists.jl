module SimplexCellLists

include("distances.jl")

using StaticArrays
using ArgCheck

const Simplex{N} = SVector{N, SVector{3, Float32}}
const Point = Simplex{1}
const Line = Simplex{2}
const Triangle = Simplex{3}


distSqr(x::Point,    y::Point)    = dist2PointPoint(      x, y)
distSqr(x::Point,    y::Line)     = dist2PointLine(       x, y)
distSqr(x::Point,    y::Triangle) = dist2PointTriangle(   x, y)
distSqr(x::Line,     y::Point)    = dist2PointLine(       y, x)
distSqr(x::Line,     y::Line)     = dist2LineLine(        x, y)
distSqr(x::Line,     y::Triangle) = dist2LineTriangle(    x, y)
distSqr(x::Triangle, y::Point)    = dist2PointTriangle(   y, x)
distSqr(x::Triangle, y::Line)     = dist2LineTriangle(    y, x)
distSqr(x::Triangle, y::Triangle) = dist2TriangleTriangle(x, y)


distSqr_sMin(x::Point, y::Point)     = dist2PointPoint(x, y), SA_F32[]
distSqr_sMin(x::Point, y::Line)      = dist2PointLine( x, y), SA_F32[]
distSqr_sMin(x::Point, y::Triangle)  = dist2PointTriangle( x, y), SA_F32[]
distSqr_sMin(x::Line, y::Point)  = ((i -> (i[1], SA_F32[i[2]])) ∘ dist2PointLine_t)( y, x)
distSqr_sMin(x::Line, y::Line)   = ((i -> (i[1], SA_F32[i[2]])) ∘ dist2LineLine_s_t)( x, y)
distSqr_sMin(x::Line, y::Triangle)   = ((i -> (i[1], i[2])) ∘ dist2LineTriangle_s_t)( x, y)
distSqr_sMin(x::Triangle, y::Point)    = dist2PointTriangle_t( y, x)
distSqr_sMin(x::Triangle, y::Line)   = ((i -> (i[1], i[3])) ∘ dist2LineTriangle_s_t)( y, x)
distSqr_sMin(x::Triangle, y::Triangle)   = ((i -> (i[1], i[2])) ∘ dist2TriangleTriangle_s_t)( x, y)

"""
Abstract type for all SimplexCellList implementations.

Implementations must have a constructor with signature like

```julia
T(numpointgroups::Integer, numlinegroups::Integer, numtrianglegroups::Integer; kwargs...)::T
```

Where `numpointgroups` is the number of groups of points, `numlinegroups` is the number of groups of lines,
and `numtrianglegroups` is the number of groups of triangles.

`kwargs` are options specific for `T`
"""
abstract type SimplexCellList end

include("naive.jl")
include("painter.jl")


"""
Reset the elements stored in `s` in batch:

```julia
setElements!(s::SimplexCellList, points, lines, triangles)::Nothing
```

Where `points`, `lines` and `triangles` are collections of collections of objects convertible to 
`Point`, `Line`, and `Triangle` respectively.

For example, each collection in `points` is a group of points that can be mapped over independently from, or together with, other groups.

Added elements will have a group index and element index based on the order of the inputs.
The first group in each type has group index 1, and the first element in each group has element index 1.
"""
function setElements! end
export setElements!


"""
Add a new element to `s`, and return its element index:

```julia
addElement!(s::SimplexCellList, group_idx::Integer, element::Simplex{N})::Int32
```
The new element will be pushed to the end of the specified group.
"""
function addElement! end
export addElement!


"""
Deactivate an existing element in `s`

```julia
deactivate!(s::SimplexCellList, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Nothing
```
Inactive elements are not mapped over. Elements are active by default.
"""
function deactivate! end
export deactivate!


"""
Re-activate an existing element in `s`

```julia
activate!(s::SimplexCellList, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Nothing
```
Inactive elements are not mapped over. Elements are active by default.
"""
function activate! end
export activate!


"""
Return if an existing element in `s` is active.

```julia
isActive(s::SimplexCellList, group_idx::Integer, element_idx::Integer, element_type::Type{Simplex{N}})::Bool
```
Inactive elements are not mapped over. Elements are active by default.
"""
function isActive end
export isActive


const mapped_function_docs = """
### Mapped function `f`

The function f should have the same form as used in CellListMap.jl.

`i` is the element index of simplex `x`, `j` is the element index of simplex `y`. 

`d2` is an approximate `Float32` squared distance between `x` and `y`.

Except here `x` and `y` are `Simplex{N}`, `Simplex{M}`

```julia
    function f(x,y,i,j,d2,output)
        # update output
        return output
    end
```
"""


"""
Map `f` to all simplexes in a group close to a single simplex.

```julia
mapSimplexElements(
        f,
        output,
        s::SimplexCellList,
        x::Simplex{N},
        group_idx::Integer,
        elements_type::Type{Simplex{M}},
        cutoff::Float32,
    ) where {N, M}
```

Apply function `f` to all elements in group `group_idx` within the cutoff range of the simplex `x`, and
return the output of the final `f` call.

`x` is always `x` and `i` is always 0, in calls to `f`.

$mapped_function_docs
"""
function mapSimplexElements end
export mapSimplexElements


"""
Map `f` to all pairs of nearby simplexes in a single group.

```julia
mapPairElements(
        f,
        output,
        s::SimplexCellList,
        group_idx::Integer,
        elements_type::Type{Simplex{N}},
        cutoff::Float32,
    ) where {N}
```
Apply function `f` to all unordered pairs of elements in group `group_idx` within cutoff range, and return the output of the final `f` call.

`f` is never called more than once per unordered pair. Which element is `x` and `y` in calls to `f` is implementation dependent.

$mapped_function_docs
"""
function mapPairElements end
export mapPairElements


"""
Map `f` to all pairs of nearby simplexes between two different groups.

```julia
mapElementsElements(
        f, 
        output, 
        s::SimplexCellList, 
        x_group_idx::Integer, 
        x_type::Type{Simplex{N}}, 
        y_group_idx::Integer, 
        y_type::Type{Simplex{M}}, 
        cutoff::Float32,
    ) where {N, M}
```
Apply function `f` to each pair of elements from two different groups that are within cutoff range of each other, and return the output of the final `f` call.

The first element has is an `x_type` in group `x_group_idx` and the second element is a `y_type` in group `y_group_idx`.

`f` is never called more than once per pair.

$mapped_function_docs
"""
function mapElementsElements end
export mapElementsElements

end