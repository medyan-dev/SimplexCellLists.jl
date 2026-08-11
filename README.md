# SimplexCellLists WIP

[![Build Status](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml?query=branch%3Amain)

This Julia package contains data structures and algorithms for doing computations on pairs of 3D points, line segments, and triangles within a cutoff distance.

It provides cell lists for points ([`PointCellList`](src/pointcelllist.jl)) and line segments ([`LineSegCellList`](src/linesegcelllist.jl)) with fast nearby-neighbor mapping, squared-distance functions (`dist_sqr`) for all pairs of points, line segments, and triangles, neighbor list construction ([`neighbor-lists.jl`](src/neighbor-lists.jl)) with a sort and sweep broad phase, and soft-sphere collide forces and energy ([`collide-forces.jl`](src/collide-forces.jl)) computed from the neighbor lists.

This package is largely inspired by [CellListMap.jl](https://github.com/m3g/CellListMap.jl).

However, there is no support for periodic boundary conditions or 2D systems.

See [CellListMap.jl](https://github.com/m3g/CellListMap.jl) if you want these features.

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

## Collide forces

Neighbor lists and collide forces are configured with a `CollisionPolicy`.
The policy decides which objects and pairs participate in the neighbor lists
(`filter_object`, `filter_pair`), how per-object parameters combine into
per-pair parameters (`mix_params`), and which force law applies to each edge
(`nl_edge_forces!`).

`DefaultCollisionPolicy` stores a stiffness and collision layer masks per
object, and applies a soft repulsive potential `E = k/2 * (L - d)²` when two
objects overlap, where `L` is the sum of the two radii and `d` is the closest
distance between them:

```julia
using SimplexCellLists, StaticArrays

policy = DefaultCollisionPolicy()

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
reused across steps. `setup_neighbors_naive!` is a reference implementation of
`setup_neighbors_sort_sweep!` useful for testing.

A custom policy subtypes `CollisionPolicy{ObjectParams, PairParams}` with its
own parameter types and methods for `filter_object`, `filter_pair`, and
`mix_params`, and can define its own `nl_edge_forces!` methods to change the
force law.
