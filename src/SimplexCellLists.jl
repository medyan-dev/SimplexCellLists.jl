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

include("force-energy.jl")
export ForceEnergy
export ForceEnergyT
export ForceEnergyFloat64
export ForceEnergyFloat32
export ForceEnergyFixedPoint
export ForceNoEnergyFixedPoint
export DebugForceEnergy
export NullForceEnergy
export add_bead_force!
export zero_bead_force!
export add_energy!
export zero_energy!
export get_nbeads
export zero_force_energy!
export get_force!
export get_bead_force
export get_energy
export get_force
export combine_force_energy!

# Other fixed point fractional bit choices must be specialized on the
# caller's parameters.
for (FE, T) in (
        (ForceEnergyFloat64, Float64),
        (ForceEnergyFixedPoint{30, 30}, Float64),
        (ForceNoEnergyFixedPoint{30}, Float64),
    )
    precompile(FE, (Int,))
    precompile(zero_force_energy!, (FE,))
    precompile(add_bead_force!, (FE, Int, SVector{3,T}))
    precompile(zero_bead_force!, (FE, Int))
    precompile(add_energy!, (FE, T))
    precompile(zero_energy!, (FE,))
    precompile(get_force!, (FE, Vector{SVector{3,T}}))
    precompile(get_bead_force, (FE, Int))
    precompile(get_force, (FE,))
    precompile(combine_force_energy!, (FE, FE))
end
# ForceNoEnergyFixedPoint has no get_energy method.
precompile(get_energy, (ForceEnergyFloat64,))
precompile(get_energy, (ForceEnergyFixedPoint{30, 30},))

end
