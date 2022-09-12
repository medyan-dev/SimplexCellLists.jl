module SimplexCellLists

include("distances.jl")

using StaticArrays
using ArgCheck

const Simplex{N} = SVector{N, SVector{3, Float32}}
const Point = Simplex{1}
const Line = Simplex{2}


distSqr(x::Point, y::Point) = dist2PointPoint(x, y)
distSqr(x::Point, y::Line)  = dist2PointLine( x, y)
distSqr(x::Line,  y::Point) = dist2PointLine( y, x)
distSqr(x::Line,  y::Line)  = dist2LineLine(  x, y)


distSqr_sMin(x::Point, y::Point) = dist2PointPoint(x, y), SA_F32[]
distSqr_sMin(x::Point, y::Line)  = dist2PointLine( x, y), SA_F32[]
distSqr_sMin(x::Line, y::Point)  = ((i -> (i[1], SA_F32[i[2]])) ∘ dist2PointLine_t)( y, x)
distSqr_sMin(x::Line, y::Line)   = ((i -> (i[1], SA_F32[i[2]])) ∘ dist2LineLine_s_t)( x, y)

"""
Abstract type for all SimplexCellList implementations.

Implementations must have a constructor with signature like

`numpointgroups::Integer, numlinegroups::Integer ;kw...`

Fields of implementations are not part of the the interface.


"""
abstract type SimplexCellList end

include("naive.jl")
include("painter.jl")

end
