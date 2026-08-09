# SimplexCellLists WIP

[![Build Status](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/medyan-dev/SimplexCellLists.jl/actions/workflows/CI.yml?query=branch%3Amain)

This Julia package contains data structures and algorithms for doing computations on pairs of 3D points, line segments, and triangles within a cutoff distance.

It provides cell lists for points ([`PointCellList`](src/pointcelllist.jl)) and line segments ([`LineSegCellList`](src/linesegcelllist.jl)) with fast nearby-neighbor mapping, and squared-distance functions (`dist_sqr`) for all pairs of points, line segments, and triangles.

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
