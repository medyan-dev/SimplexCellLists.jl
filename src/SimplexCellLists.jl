module SimplexCellLists

using LinearAlgebra
using StaticArrays
using ArgCheck

const Vec3 = SVector{3}
const Simplex{N,T} = SVector{N, Vec3{T}}
const Point = Simplex{1}
const Line = Simplex{2}
const Triangle = Simplex{3}

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
export cell_point_clear!

end