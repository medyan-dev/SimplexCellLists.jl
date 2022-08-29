module MultiShapeCellLists

include("distances.jl")

using StaticArrays

const Simplex{N} = SVector{N, SVector{3, Float32}}
const Point = Simplex{1}
const Line = Simplex{2}


distSqr(x::Point, y::Point) = dist2PointPoint(x, y)
distSqr(x::Point, y::Line)  = dist2PointLine( x, y)
distSqr(x::Line,  y::Point) = dist2PointLine( y, x)
distSqr(x::Line,  y::Line)  = dist2LineLine(  x, y)

"""
Abstract type for all MultiShapeCellList implementations.

Implementations must have a constructor with signature like

`numpointgroups::Integer, numlinegroups::Integer ;kw...`

Fields of implementations are not part of the the interface.


"""
abstract type MultiShapeCellList end


end
