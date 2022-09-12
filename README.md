# SimplexCellLists WIP

[![Build Status](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml?query=branch%3Amain)

This Julia package accelerates computations on all pairs of 3D points, line segments, and triangles within a cutoff distance.

This package is largely inspired by [CellListMap.jl](https://github.com/m3g/CellListMap.jl).

However, there is no support for periodic boundary conditions, 2D systems, or types other than Float32.

See [CellListMap.jl](https://github.com/m3g/CellListMap.jl) if you want these features, or higher performance on systems without triangles and line segments.


## Setup

- `Simplex{N}` is `SVector{N, SVector{3, Float32}}`
- `Point` is `Simplex{1}`
- `Line` is `Simplex{2}`
- `Triangle` is `Simplex{3}`

There are multiple algorithms that implement the `SimplexCellList` abstract type interface.
Currently:
- `Naive`: A simple `SimplexCellList` with no extra options that uses double for loops.
- `Painter`: A grid is painted with element ids based on a max range. Based on the ideas in [cytosim](https://gitlab.com/f-nedelec/cytosim/-/blob/af739d2ff768628e4737d3a75457676e1a7f4287/src/sim/fiber_grid.h).

Let `T` be a concrete subtype of `SimplexCellList`

### Constructor
Construct `T`:

```julia
T(numpointgroups::Integer, numlinegroups::Integer, numtrianglegroups::Integer; kwargs...)::T
```

Where `numpointgroups` is the number of groups of points, `numlinegroups` is the number of groups of lines,
and `numtrianglegroups` is the number of groups of triangles.

`kwargs` are options specific for `T`

### `setElements!`
Reset the elements stored in `T` in batch:

```julia
setElements!(m::T, points, lines, triangles)::Nothing
```

Where `points`, `lines` and `triangles` are collections of collections of objects convertible to 
`Point`, `Line`, and `Triangle` respectively.

Added elements will have a group ID and element ID based on the order of the inputs.
The first group in each type has group ID 1, and the first element in each group has element ID 1.

### `addElement!`
Add a new element to `T`, and return its element ID:

```julia
addElement!(m::T, groupid::Integer, element::Simplex{N})::Int32
```
The new element will be pushed to the end of the specified group.

### `deactivate!`
Deactivate an existing element in `T`

```julia
deactivate!(m::T, groupid::Integer, elementid::Integer, elementtype::Type{Simplex{N}})::Nothing
```
Inactive elements are not mapped over. Elements are active by default.

### `activate!`
Re-activate an existing element in `T`

```julia
activate!(m::T, groupid::Integer, elementid::Integer, elementtype::Type{Simplex{N}})::Nothing
```
Inactive elements are not mapped over. Elements are active by default.

### `isActive`
Return if an existing element in `T` is active.

```julia
isActive(m::T, groupid::Integer, elementid::Integer, elementtype::Type{Simplex{N}})::Bool
```
Inactive elements are not mapped over. Elements are active by default.

## Mapping

The following functions allow mapping a custom function over pairs of simplexes within some cutoff.

### Mapped function `f`

The function f should have the same form as used in CellListMap.jl
Except here `x` and `y` are `SVector{N, SVector{3, Float32}}`, `SVector{M, SVector{3, Float32}}`

```julia
    function f(x,y,i,j,d2,output)
        # update output
        return output
    end
```

The order in which pairs of elements in range are mapped is implementation dependent.

The elements passed to `f` may be slightly different from the elements added to `T` due to implementation dependent floating point rounding errors.

If a pair distance is very near the cutoff, it is implementation dependent whether the pair gets mapped or not due to floating point rounding errors.

Therefore, if more precision is needed, add some extra distance to the cutoff, store the elements externally in 64 bit precision, and in `f` use `i` and `j` to get the precise elements and again check distances.

### `mapSimplexElements`

Apply function `f` to all elements in group `groupid` within the cutoff range of `x`, and
return the output of the final `f` call.

```julia
mapSimplexElements(f, output, m::T, groupid::Integer, x::Simplex{N}, elementstype::Type{Simplex{M}}, cutoff::Float32) where {N, M}
```

`x` is always `x` and `i` is always 0, in calls to `f`.

### `mapPairElements`

Apply function `f` to all unordered pairs of elements in group `groupid` within cutoff range, and return the output of the final `f` call.

```julia
mapPairElements(f, output, m::T, groupid::Integer, elementstype::Type{Simplex{N}}, cutoff::Float32) where {N}
```

`f` is never called more than once per unordered pair. Which element is `x` and `y` in calls to `f` is implementation dependent.


### `mapElementsElements`

Apply function `f` to each pair of elements from two different groups that are within cutoff range of each other, and return the output of the final `f` call.


```julia
mapElementsElements(
        f, 
        output, 
        m::T, 
        x_groupid::Integer, 
        x_type::Type{Simplex{N}}, 
        y_groupid::Integer, 
        y_type::Type{Simplex{M}}, 
        cutoff::Float32,
    ) where {N, M}
```

The first element has is an `x_type` in group `x_groupid` and the second element is a `y_type` in group `y_groupid`.

`f` is never called more than once per pair.