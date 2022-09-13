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

For example, each collection in `points` is a group of points that can be mapped over independently from, or together with, other groups.

Added elements will have a group index and element index based on the order of the inputs.
The first group in each type has group index 1, and the first element in each group has element index 1.

### `addElement!`
Add a new element to `T`, and return its element index:

```julia
addElement!(m::T, groupidx::Integer, element::Simplex{N})::Int32
```
The new element will be pushed to the end of the specified group.

### `deactivate!`
Deactivate an existing element in `T`

```julia
deactivate!(m::T, groupidx::Integer, elementidx::Integer, elementtype::Type{Simplex{N}})::Nothing
```
Inactive elements are not mapped over. Elements are active by default.

### `activate!`
Re-activate an existing element in `T`

```julia
activate!(m::T, groupidx::Integer, elementidx::Integer, elementtype::Type{Simplex{N}})::Nothing
```
Inactive elements are not mapped over. Elements are active by default.

### `isActive`
Return if an existing element in `T` is active.

```julia
isActive(m::T, groupidx::Integer, elementidx::Integer, elementtype::Type{Simplex{N}})::Bool
```
Inactive elements are not mapped over. Elements are active by default.

## Mapping

The following functions allow mapping a custom function over pairs of simplexes within some cutoff.

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

The order in which pairs of elements in range are mapped is implementation dependent.

The elements passed to `f` may be slightly different from the elements added to `T` due to implementation dependent floating point rounding errors.

If a pair distance is very near the cutoff, it is implementation dependent whether the pair gets mapped or not due to floating point rounding errors.

Therefore, if more precision is needed, add some extra distance to the cutoff, store the elements externally in 64 bit precision, and in `f` use `i` and `j` to get the precise elements and again check distances.

### `mapSimplexElements`

Map `f` to all simplexes in a group close to a single simplex.

```julia
mapSimplexElements(f, output, m::T, groupidx::Integer, x::Simplex{N}, elementstype::Type{Simplex{M}}, cutoff::Float32) where {N, M}
```

Apply function `f` to all elements in group `groupidx` within the cutoff range of the simplex `x`, and
return the output of the final `f` call.

`x` is always `x` and `i` is always 0, in calls to `f`.

### `mapPairElements`

Map `f` to all pairs of nearby simplexes in a single group.

```julia
mapPairElements(f, output, m::T, groupidx::Integer, elementstype::Type{Simplex{N}}, cutoff::Float32) where {N}
```
Apply function `f` to all unordered pairs of elements in group `groupidx` within cutoff range, and return the output of the final `f` call.

`f` is never called more than once per unordered pair. Which element is `x` and `y` in calls to `f` is implementation dependent.


### `mapElementsElements`

Map `f` to all pairs of nearby simplexes between two different groups.

```julia
mapElementsElements(
        f, 
        output, 
        m::T, 
        x_groupidx::Integer, 
        x_type::Type{Simplex{N}}, 
        y_groupidx::Integer, 
        y_type::Type{Simplex{M}}, 
        cutoff::Float32,
    ) where {N, M}
```
Apply function `f` to each pair of elements from two different groups that are within cutoff range of each other, and return the output of the final `f` call.

The first element has is an `x_type` in group `x_groupidx` and the second element is a `y_type` in group `y_groupidx`.

`f` is never called more than once per pair.
