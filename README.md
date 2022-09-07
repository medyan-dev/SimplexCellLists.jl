# SimplexCellLists WIP

[![Build Status](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml?query=branch%3Amain)

This Julia package accelerates computations on all pairs of 3D points, line segments, and triangles within a cutoff distance.

This package is largely inspired by [CellListMap.jl](https://github.com/m3g/CellListMap.jl).

However, there is no support for periodic boundary conditions, 2D systems, or types other than Float32.

See [CellListMap.jl](https://github.com/m3g/CellListMap.jl) if you want these features, or higher performance on systems without triangles and line segments.


## API

- `Simplex{N}` is `SVector{N, SVector{3, Float32}}`
- `Point` is `SVector{1, SVector{3, Float32}}`
- `Line` is `SVector{2, SVector{3, Float32}}`
- `Triangle` is `SVector{3, SVector{3, Float32}}`

There are multiple algorithms that implement the `SimplexCellList` abstract type interface.
Currently:
- `Naive`: A simple `SimplexCellList` with no extra options that uses double for loops.
- `Painter`: A grid is painted with element ids based on a max range. Based on the ideas in [cytosim](https://gitlab.com/f-nedelec/cytosim/-/blob/af739d2ff768628e4737d3a75457676e1a7f4287/src/sim/fiber_grid.h).

Let `T` be a concrete subtype of `SimplexCellList`

### Constructor
`T` can be constructed with:

```julia
T(numpointgroups::Integer, numlinegroups::Integer, numtrianglegroups::Integer; kwargs...)::T
```

Where `numpointgroups` is the number of groups of points, `numlinegroups` is the number of groups of lines,
and `numtrianglegroups` is the number of groups of triangles.

`kwargs` are options specific for `T`

### `setElements!`
The elements stored in `T` can be reset in batch with:

```julia
setElements!(m::T, points, lines, triangles)::Nothing
```

Where `points`, `lines` and `triangles` are iterators with `length` of the number of groups of each type. Each iterates to get containers of objects convertible to 
`Point`, `Line`, and `Triangle` respectively.
Added elements will have a group id and element id based on the order of the inputs.
Ids are one based.

### `addElement!`
A new element can be added to `T` with:

```julia
addElement!(m::T, groupid::Integer, element::Simplex{N})::Int32
```
The new element id is returned.


### `deleteElement!`
An existing element can be deleted from `T` with:

```julia
SimplexCellLists.deleteElement!(m::T, groupid::Integer, elementid::Integer, elementtype::Type{Simplex{N}})::Int32
```
The other element ids are stable when elements are deleted. To avoid performance issues, with holes, use `setElements!` if you are removing a significant fraction of the elements, to make the ids contiguous again.


### `mapSimplexElements`
