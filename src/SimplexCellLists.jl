module SimplexCellLists

using LinearAlgebra
using StaticArrays
using ArgCheck

const Vec3 = SVector{3}
const Simplex{N,T} = SVector{N, Vec3{T}}
const Point = Simplex{1}
const LineSeg = Simplex{2}
const Triangle = Simplex{3}
public Vec3, Simplex, Point, LineSeg, Triangle

include("util.jl")
export ball_intersects_unit_cube
export ball_intersects_biunit_cube

include("distances.jl")
export dist_sqr
export dist_sqr_r_s_t

include("pointcelllist.jl")
export PointCellList
export CellPointEntry
export cell_point_add!
export map_nearby_points
export cell_points_clear!

include("linesegcelllist.jl")
export LineSegCellList
export CellLineSegEntry
export cell_line_seg_add!
export cell_line_seg_is_end
export cell_line_seg_tmax
export map_nearby_line_segs
export cell_line_segs_clear!

end