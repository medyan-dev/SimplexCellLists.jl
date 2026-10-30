# SimplexCellLists WIP

[![Build Status](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml?query=branch%3Amain)

This Julia package contains data structures and algorithms for doing computations on pairs of 3D points, line segments, and triangles within a cutoff distance.

It provides:

- Cell lists for points ([`PointCellList`](src/pointcelllist.jl)) and line segments ([`LineSegCellList`](src/linesegcelllist.jl)) with fast nearby-neighbor mapping.
- Squared-distance functions ([`dist_sqr`, `dist_sqr_r_s_t`](src/distances.jl)) for all pairs of points, line segments, and triangles.
- Neighbor list construction with two broad phases: a multi-level grid ([`setup_neighbors_grid!`](src/neighbor-lists-grid.jl)) and sort and sweep ([`setup_neighbors_sort_sweep!`](src/neighbor-lists.jl)).
- Soft-sphere collide forces and energy ([`collide_forces!`](src/collide-forces.jl)) computed from the neighbor lists.

This package is largely inspired by [CellListMap.jl](https://github.com/m3g/CellListMap.jl).

However, there is no support for periodic boundary conditions or 2D systems.

See [CellListMap.jl](https://github.com/m3g/CellListMap.jl) if you want these features.

## Installation

SimplexCellLists requires Julia 1.12 or later.

```julia
using Pkg
Pkg.add("SimplexCellLists")
```

## Quickstart

```julia
using SimplexCellLists, StaticArrays

# 10×10×10 grid of cells, each 1.0 wide, centered on the origin
pcl = PointCellList{Int64,Float32}((10, 10, 10), 1.0f0)
cell_point_add!(pcl, SA[0.1f0, 0.2f0, 0.3f0], 1)
cell_point_add!(pcl, SA[0.3f0, 0.2f0, 0.3f0], 2)

# Count the points within 0.5 of the origin
n = map_nearby_points(pcl, SA[0.0f0, 0.0f0, 0.0f0], 0.5f0, 0) do entry, sep, count
    (count + 1, true)
end
```

`LineSegCellList`, `cell_line_seg_add!`, and `map_nearby_line_segs` are the line segment analogs, and `cell_points_clear!`/`cell_line_segs_clear!` reset a list for reuse.

## Distances

`dist_sqr(x, y)` returns the minimum squared distance between two simplexes,
each an `SVector` of 1, 2, or 3 vertices (a point, line segment, or
triangle). `dist_sqr_r_s_t(x, y)` also returns the separation vector and the
barycentric coordinates of the closest points.

```julia
using SimplexCellLists, StaticArrays

seg = SA[SA[0.0, 0.0, 1.0], SA[1.0, 0.0, 1.0]]
tri = SA[SA[0.0, 0.0, 0.0], SA[1.0, 0.0, 0.0], SA[0.0, 1.0, 0.0]]
dist_sqr(seg, tri) # 1.0
```

The results are accurate in both `Float32` and `Float64`: with `ℓ` the
largest absolute coordinate, `√d²` is within `64*eps(T)*ℓ + floatmin(T)^(1/4)`
of the exact minimum distance, for coordinates below `floatmax(T)^(1/4)/4`
(1e9 for `Float32`). See the `dist_sqr_r_s_t` docstring for the full
contract. This bound was proven in Lean 4 by AI in
[line-segment-triangle-distance-AI-proofs](https://github.com/nhz2/line-segment-triangle-distance-AI-proofs).

## Collide forces

Neighbor lists and collide forces are configured with a `CollidePolicy`.
The policy decides which objects and pairs participate in the neighbor lists
(`filter_object`, `filter_pair`), how per-object parameters combine into
per-pair parameters (`mix_params`), and which force law applies to each edge
(`nl_edge_force`).

`DefaultCollidePolicy` stores a stiffness and collision layer masks per
object, and applies a soft repulsive potential `E = k/2 * (L - d)²` when two
objects overlap, where `L` is the sum of the two radii and `d` is the closest
distance between them. The pair stiffness `k` mixes the two per-object
stiffnesses like springs in series, `k = k₁*k₂/(k₁ + k₂)`, so the two
`10.0f0` stiffnesses in the example below mix to `k = 5`:

```julia
using SimplexCellLists, StaticArrays

policy = DefaultCollidePolicy()

# Two overlapping spheres: radius 0.5, centers 0.5 apart
pos = [SA[0.0, 0.0, 0.0], SA[0.5, 0.0, 0.0]]
params = DefaultObjectParams(10.0f0, UInt32(1), UInt32(0)) # stiffness, layers, no collide mask
inputs = NeighborListInputs(policy;
    points = [PointIdxPart(1), PointIdxPart(2)],
    p_radius = [0.5f0, 0.5f0],
    p_params = [params, params],
)
nl = NeighborLists(policy)
setup_neighbors_sort_sweep!(nl, pos, inputs)

force_energy = ForceEnergyFloat64(length(pos))
collide_forces!(force_energy, pos, nl, Float64)
get_energy(force_energy) # 0.625
get_force(force_energy)  # [-2.5, 0.0, 0.0], [2.5, 0.0, 0.0]
```

`NeighborListInputs` also accepts line segments (`clines`/`lines`) and
triangles, `no_collide_pairs` exclusions, and a `skin` distance so lists can be
reused across steps.

There are three ways to build the neighbor lists, which all produce the same
set of edges:

- `setup_neighbors_grid!` uses a multi-level grid. Larger objects are stored
  at coarser levels and smaller objects at finer levels, so a mix of sizes
  works well. The grid is sized to the number of objects, so memory use
  scales with the number of objects, not the size of the system.
- `setup_neighbors_sort_sweep!` sorts the bounding boxes along the widest
  axis of the system and sweeps over them.
- `setup_neighbors_naive!` checks every pair. It is a reference
  implementation, useful for testing.

The last argument of `collide_forces!` is the floating-point type the forces
are computed in. Positions can be stored in `Float64` and the forces computed
in `Float32`: relative positions are subtracted in the position type before
they are converted. For multithreading, pass the `chunk` and `nthreads`
keywords and merge the per-thread results with `combine_force_energy!`, see
the `collide_forces!` docstring.

A custom policy subtypes `CollidePolicy{ObjectParams, PairParams}` with its
own parameter types and methods for `filter_object`, `filter_pair`, and
`mix_params`, and can define its own `nl_edge_force` methods to change the
force law.
