"""
    CellTriangleEntry{T,F}

A single triangle's entry in one voxel cell. Created by [`cell_triangle_add!`](@ref)
and read by [`map_nearby_triangles`](@ref).

# Fields
- `id::T`: triangle ID.
- `tri::SVector{3, SVector{3, F}}`: triangle vertices.
- `box_ends::UInt8`: where this cell is in the range of cells touched by the
  triangle's bounding box. Bit `a - 1` is set if it is the first cell of the
  range along axis `a`, and bit `a + 2` if it is the last.
"""
struct CellTriangleEntry{T,F}
    id::T
    tri::SVector{3, SVector{3, F}}
    box_ends::UInt8
end

struct CellTriangle{T,F}
    list::Memory{CellTriangleEntry{T,F}}
    len::Int
end
CellTriangle{T,F}() where {T,F} = CellTriangle{T,F}(Memory{CellTriangleEntry{T,F}}(), 0)

"""
    TriangleCellList{T,F}(size::SVector{3, Int}, grid_spacing::F)
    TriangleCellList{T,F}(size::NTuple{3, Int}, grid_spacing::F)

An origin-centered 3D voxel grid for spatial queries on triangles. The grid spans
`[-size * grid_spacing/2, size * grid_spacing/2)` in each dimension.

Insert triangles with [`cell_triangle_add!`](@ref) and query with
[`map_nearby_triangles`](@ref). Use [`cell_triangles_clear!`](@ref) to reset.

Each triangle is stored in every cell its bounding box touches, so triangles
much larger than `grid_spacing` take many entries, growing with the volume of
their bounding box.

Positions are not required to be within the grid bounds: parts of a triangle
outside the grid are stored in the nearest boundary cells, and queries account
for this. Triangles far outside the grid accumulate in boundary cells, which
can slow down queries near those boundaries.

# Arguments
- `size`: number of cells along each axis (x, y, z).
- `grid_spacing::F`: side length of each cubic cell.
"""
struct TriangleCellList{T,F}
    cells::Memory{CellTriangle{T,F}}
    size::SVector{3, Int}
    inv_half_grid_spacing::F
end
function TriangleCellList{T,F}(size::SVector{3, Int}, grid_spacing::F) where {T,F}
    @argcheck all(≥(0), size)
    @argcheck grid_spacing > 0
    inv_half_grid_spacing = 2*inv(grid_spacing)
    @argcheck isfinite(inv_half_grid_spacing)
    n = prod(size)
    cells = Memory{CellTriangle{T,F}}(undef, n)
    fill!(cells, CellTriangle{T,F}())
    TriangleCellList{T,F}(cells, size, inv_half_grid_spacing)
end
TriangleCellList{T,F}(size::NTuple{3, Int}, grid_spacing::F) where {T,F} = TriangleCellList{T,F}(SVector{3,Int}(size...), grid_spacing)

"""
    cell_triangles_clear!(tcl::TriangleCellList)::Nothing

Remove all triangles.
"""
function cell_triangles_clear!(tcl::TriangleCellList{T,F})::Nothing where {T,F}
    for i in eachindex(tcl.cells)
        tcl.cells[i] = CellTriangle{T,F}(tcl.cells[i].list, 0)
    end
end

# Return the `box_ends` bits of axis `a` (1 to 3) for the cell `i` of the
# range `lo:hi`, see `CellTriangleEntry`.
@inline _box_end_bits(i::Int, lo::Int, hi::Int, a::Int)::UInt8 =
    (UInt8(i == lo) << (a - 1)) | (UInt8(i == hi) << (a + 2))

"""
    cell_triangle_add!(tcl::TriangleCellList{T,F}, p0::SVector{3, F}, p1::SVector{3, F}, p2::SVector{3, F}, id::T)::Nothing

Insert the triangle with vertices `p0`, `p1`, and `p2` into the cell list, adding a
[`CellTriangleEntry`](@ref) to every cell its bounding box touches. Each entry
records where its cell is in that range of cells, so
[`map_nearby_triangles`](@ref) reports the triangle exactly once.

Positions are not required to be within the grid bounds: parts of the triangle
outside the grid are stored in the nearest boundary cells, and queries account
for this. Triangles far outside the grid accumulate in boundary cells, which can
slow down queries near those boundaries.

Throws an `ArgumentError` if a coordinate is NaN or Inf.

# Arguments
- `tcl`: the triangle cell list to insert into.
- `p0`, `p1`, `p2`: triangle vertices.
- `id`: ID for this triangle.
"""
function cell_triangle_add!(
        tcl::TriangleCellList{T,F},
        p0::SVector{3, F},
        p1::SVector{3, F},
        p2::SVector{3, F},
        id::T,
    )::Nothing where {T,F}
    h2_inv = tcl.inv_half_grid_spacing
    sz = tcl.size
    @argcheck all(>(0), sz)
    @argcheck all(isfinite, p0) && all(isfinite, p1) && all(isfinite, p2)
    strides = SVector(1, sz[1], sz[1]*sz[2])
    tri = SA[p0, p1, p2]
    # The cells touched by the bounding box, in doubled coordinates as in
    # cell_point_add!. min and max are exact, and scaling, clamping, and
    # round2grid are monotone, so every vertex's cell is in the range.
    lo = round2grid(clamp.(min.(p0, p1, p2) * h2_inv, -F.(sz), prevfloat.(F.(sz))), sz)
    hi = round2grid(clamp.(max.(p0, p1, p2) * h2_inv, -F.(sz), prevfloat.(F.(sz))), sz)
    for iz in lo[3]:hi[3]
        liz = iz*strides[3]
        bz = _box_end_bits(iz, lo[3], hi[3], 3)
        for iy in lo[2]:hi[2]
            liy = iy*strides[2]
            byz = bz | _box_end_bits(iy, lo[2], hi[2], 2)
            for ix in lo[1]:hi[1]
                li = ix + liy + liz
                entry = CellTriangleEntry{T,F}(id, tri, byz | _box_end_bits(ix, lo[1], hi[1], 1))
                # li is inbounds because lo and hi are clamped to the grid,
                # and cells has length prod(sz).
                @inbounds tcl.cells[begin + li] = _cell_push_entry!(tcl.cells[begin + li], entry)
            end
        end
    end
    nothing
end

# Return the `box_ends` bits of axis `a` an entry in the cell `i` must have
# to be reported to a query in the cell `p`: the first cell bit if `i` is
# past `p`, and the last cell bit if `i` is before `p`.
@inline _box_end_mask(i::Int, p::Int, a::Int)::UInt8 =
    (UInt8(i > p) << (a - 1)) | (UInt8(i < p) << (a + 2))

"""
    map_nearby_triangles(f, tcl::TriangleCellList{T,F}, pos::SVector{3, F}, cutoff::F, out) -> out

Call `f` for every triangle in `tcl` whose closest point to `pos` is within `cutoff` distance.
Each triangle is reported exactly once (deduplicated across cells).

The distance is the `d²` of `dist_sqr_r_s_t(SA[pos], entry.tri)`, see
[`dist_sqr_r_s_t`](@ref), which also gives the barycentric coordinates of the
closest point, bit for bit.

# Arguments
- `f(entry::CellTriangleEntry{T,F}, sep::SVector{3,F}, out) -> (out, cont::Bool)`: callback called
  for each nearby triangle. `sep` is the separation vector from `pos` to the closest point on the
  triangle.
  Return `(new_out, true)` to continue, or `(new_out, false)` to stop early.
- `tcl`: the triangle cell list to query.
- `pos`: query point.
- `cutoff`: maximum distance. If zero, returns `out` immediately.
- `out`: initial accumulator value, threaded through `f` and returned.
"""
function map_nearby_triangles(
        f,
        tcl::TriangleCellList{T,F},
        pos::SVector{3, F},
        cutoff::F,
        out,
    ) where {T,F}
    cutoff ≤ 0 && return out
    sz = tcl.size
    all(>(0), sz) || return out
    cutoff2 = cutoff * cutoff
    h2_inv = tcl.inv_half_grid_spacing
    strides = SVector(1, sz[1], sz[1]*sz[2])
    # Clamping the query into the grid box makes the cube test below stay
    # conservative given triangles can extend outside the grid.
    # Clamping in F before converting to Int also keeps far away
    # queries and huge cutoffs from overflowing Int.
    pos_norm = clamp.(pos * h2_inv, -F.(sz), prevfloat.(F.(sz)))
    p = round2grid(pos_norm, sz)
    # Deduplication, like the grid neighbor lists: a triangle within the
    # cutoff is only reported from the cell of the point of its bounding box
    # closest to `pos`. Since round2grid is monotone, along each axis that is
    # the cell of `pos` clamped to the range of cells the box touches, so
    # it's exactly one of the cells the triangle is stored in. Along an axis
    # where a cell is past the cell of `pos`, it must be the first cell of
    # the range, and before, the last, which is a mask of `box_ends` bits.
    #
    # That point is no farther from `pos` than the triangle, so its cell
    # passes the cube test. The test is padded so rounding in the scaling,
    # the cube test, and the distance can't make it miss the cell, for
    # positions within the grid. The distance error is within
    # `64*eps(F)*ℓ`, see dist_sqr_r_s_t.
    r_cells = cutoff * h2_inv
    r_cells += 128 * eps(F) * (2 * r_cells + F(maximum(sz)))
    lo = round2grid(max.(pos_norm .- r_cells, -F.(sz)), sz)
    hi = round2grid(min.(pos_norm .+ r_cells, prevfloat.(F.(sz))), sz)
    r2_cells = r_cells * r_cells
    # The per-axis terms of the ball-to-cube distance test and the mask are
    # hoisted to their loop levels so the z and y contributions are not
    # recomputed for every cell.
    for iz in lo[3]:hi[3]
        liz = iz*strides[3]
        cell_center_z = 2*iz + 1 - sz[3]
        dz2 = relu(abs(pos_norm[3] - cell_center_z) - one(F))^2
        mz = _box_end_mask(iz, p[3], 3)
        for iy in lo[2]:hi[2]
            liy = iy*strides[2]
            cell_center_y = 2*iy + 1 - sz[2]
            dyz2 = dz2 + relu(abs(pos_norm[2] - cell_center_y) - one(F))^2
            myz = mz | _box_end_mask(iy, p[2], 2)
            for ix in lo[1]:hi[1]
                cell_center_x = 2*ix + 1 - sz[1]
                dxyz2 = dyz2 + relu(abs(pos_norm[1] - cell_center_x) - one(F))^2
                dxyz2 ≤ r2_cells || continue
                mask = myz | _box_end_mask(ix, p[1], 1)
                li = ix + liy + liz
                # li is inbounds because lo and hi are clamped to the grid,
                # and cells has length prod(sz). cs.len ≤ length(cs.list) is
                # a _cell_push_entry! invariant.
                cs = @inbounds tcl.cells[begin + li]
                for j in 1:cs.len
                    entry = @inbounds cs.list[j]
                    # ownership check, before the distance
                    entry.box_ends & mask == mask || continue
                    d2, sep, _, _ = dist_sqr_r_s_t(SA[pos], entry.tri)
                    d2 ≤ cutoff2 || continue
                    out, cont = f(entry, sep, out)
                    cont || return out
                end
            end
        end
    end
    out
end
