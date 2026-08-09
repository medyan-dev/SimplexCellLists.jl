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

# The mapping functions are not precompiled because they must be specialized
# on the caller's callback function anyway.
let T = Int64, F = Float32
    precompile(PointCellList{T,F}, (SVector{3,Int}, F))
    precompile(PointCellList{T,F}, (NTuple{3,Int}, F))
    precompile(cell_points_clear!, (PointCellList{T,F},))
    precompile(cell_point_add!, (PointCellList{T,F}, SVector{3,F}, T))
    precompile(LineSegCellList{T,F}, (SVector{3,Int}, F))
    precompile(LineSegCellList{T,F}, (NTuple{3,Int}, F))
    precompile(cell_line_segs_clear!, (LineSegCellList{T,F},))
    precompile(cell_line_seg_add!, (LineSegCellList{T,F}, SVector{3,F}, SVector{3,F}, T))
end

end