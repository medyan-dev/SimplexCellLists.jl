"""
    CellPointEntry{T,F}

A single point's entry in one voxel cell. Created by [`cell_point_add!`](@ref)
and read by [`map_nearby_points`](@ref).

# Fields
- `pos::SVector{3,F}`: point position.
- `id::T`: point ID.
"""
struct CellPointEntry{T,F}
    pos::SVector{3,F}
    id::T
end

struct CellPoint{T,F}
    list::Memory{CellPointEntry{T,F}}
    len::Int
end
CellPoint{T,F}() where{T,F} = CellPoint{T,F}(Memory{CellPointEntry{T,F}}(), 0)

function _cell_push_entry!(old::CELLTYPE, entry::ENTRYTYPE) where {CELLTYPE,ENTRYTYPE}
    old_list = old.list
    cap = length(old_list)
    len = old.len
    list = if cap ≤ len
        new_cap = max(cap<<1, 1)
        @assert new_cap > len
        new_list = Memory{ENTRYTYPE}(undef, new_cap)
        unsafe_copyto!(new_list, 1, old_list, 1, cap)
        new_list
    else
        old_list
    end
    # len+1 ≤ length(list) because the branch above grows the list past len.
    @inbounds list[len+1] = entry
    CELLTYPE(list, len+1)
end

"""
    PointCellList{T,F}(size::SVector{3, Int}, grid_spacing::F)
    PointCellList{T,F}(size::NTuple{3, Int}, grid_spacing::F)

An origin-centered 3D voxel grid for spatial queries on points. The grid spans
`[-size * grid_spacing/2, size * grid_spacing/2)` in each dimension.

Insert points with [`cell_point_add!`](@ref) and query with
[`map_nearby_points`](@ref). Use [`cell_points_clear!`](@ref) to reset.

Positions are not required to be within the grid bounds: points outside the
grid are stored in the nearest boundary cell, and queries account for this.
Points far outside the grid accumulate in boundary cells, which can slow down
queries near those boundaries.

# Arguments
- `size`: number of cells along each axis (x, y, z).
- `grid_spacing::F`: side length of each cubic cell.
"""
struct PointCellList{T,F}
    cells::Memory{CellPoint{T,F}}
    size::SVector{3, Int}
    inv_half_grid_spacing::F
end
function PointCellList{T,F}(size::SVector{3, Int}, grid_spacing::F) where {T,F}
    @argcheck all(≥(0), size)
    @argcheck grid_spacing > 0
    inv_half_grid_spacing = 2*inv(grid_spacing)
    @argcheck isfinite(inv_half_grid_spacing)
    n = prod(size)
    cells = Memory{CellPoint{T,F}}(undef, n)
    fill!(cells, CellPoint{T,F}())
    PointCellList{T,F}(cells, size, inv_half_grid_spacing)
end
PointCellList{T,F}(size::NTuple{3, Int}, grid_spacing::F) where {T,F} = PointCellList{T,F}(SVector{3,Int}(size...), grid_spacing)

"""
    cell_points_clear!(pcl::PointCellList)::Nothing

Remove all points.
"""
function cell_points_clear!(pcl::PointCellList{T,F})::Nothing where {T,F}
    for i in eachindex(pcl.cells)
        pcl.cells[i] = CellPoint{T,F}(pcl.cells[i].list, 0)
    end
end

"""
    cell_point_add!(pcl::PointCellList{T,F}, pos::SVector{3, F}, id::T)::Nothing

Insert a point into the cell list, adding a [`CellPointEntry`](@ref) to the cell
containing the point. If the point is outside the grid bounds, it is added to
the nearest boundary cell.

# Arguments
- `pcl`: the point cell list to insert into.
- `pos`: point position.
- `id`: ID for this point.
"""
function cell_point_add!(pcl::PointCellList{T,F}, pos::SVector{3, F}, id::T)::Nothing where {T,F}
    h2_inv = pcl.inv_half_grid_spacing
    sz = pcl.size
    @argcheck all(>(0), sz)
    strides = SVector(1, sz[1], sz[1]*sz[2])
    # Work in doubled coordinates 2*pos/h so the origin-centering offset is
    # added exactly in Int after the floor, keeping cell assignment accurate
    # for any grid size. Clamp in F before converting to Int so far away
    # positions don't overflow Int; the prevfloat upper bound keeps the
    # resulting index ≤ size - 1.
    pos_norm = clamp.(pos*h2_inv, -F.(sz), prevfloat.(F.(sz)))
    i = round2grid(pos_norm, sz)
    li = i ⋅ strides
    entry = CellPointEntry{T,F}(pos, id)
    # li is inbounds because pos_norm is clamped to the grid, and cells has
    # length prod(sz).
    @inbounds pcl.cells[begin + li] = _cell_push_entry!(pcl.cells[begin + li], entry)
    nothing
end

"""
    map_nearby_points(f, pcl::PointCellList{T,F}, pos::SVector{3, F}, cutoff::F, out) -> out

Call `f` for every point in `pcl` within `cutoff` distance of `pos`.

# Arguments
- `f(entry::CellPointEntry{T,F}, out) -> (out, cont::Bool)`: callback called for each nearby point.
  Return `(new_out, true)` to continue, or `(new_out, false)` to stop early.
- `pcl`: the point cell list to query.
- `pos`: query point.
- `cutoff`: maximum distance. If zero, returns `out` immediately.
- `out`: initial accumulator value, threaded through `f` and returned.
"""
function map_nearby_points(
        f,
        pcl::PointCellList{T,F},
        pos::SVector{3, F},
        cutoff::F,
        out,
    ) where {T,F}
    cutoff ≤ 0 && return out
    sz = pcl.size
    all(>(0), sz) || return out
    cutoff2 = cutoff * cutoff
    h2_inv = pcl.inv_half_grid_spacing
    strides = SVector(1, sz[1], sz[1]*sz[2])
    # Clamping the query into the grid box makes the cube test below treat
    # boundary cells as extending outward to infinity, so points clamped into
    # boundary cells by cell_point_add! are never pruned away. For any stored
    # point q in cell c, dist(clamped pos_norm, cube of c) ≤ |pos - q| * 2/h,
    # so the test stays conservative. Clamping in F before converting to Int
    # also keeps far away queries and huge cutoffs from overflowing Int.
    pos_norm = clamp.(pos * h2_inv, -F.(sz), prevfloat.(F.(sz)))
    r_cells = cutoff * h2_inv
    lo = round2grid(max.(pos_norm .- r_cells, -F.(sz)), sz)
    hi = round2grid(min.(pos_norm .+ r_cells, prevfloat.(F.(sz))), sz)
    r2_cells = r_cells * r_cells
    # The per-axis terms of the ball-to-cube distance test are hoisted to
    # their loop levels so the z and y contributions are not recomputed for
    # every cell.
    for iz in lo[3]:hi[3]
        liz = iz*strides[3]
        cell_center_z = 2*iz + 1 - sz[3]
        dz2 = relu(abs(pos_norm[3] - cell_center_z) - one(F))^2
        for iy in lo[2]:hi[2]
            liy = iy*strides[2]
            cell_center_y = 2*iy + 1 - sz[2]
            dyz2 = dz2 + relu(abs(pos_norm[2] - cell_center_y) - one(F))^2
            for ix in lo[1]:hi[1]
                cell_center_x = 2*ix + 1 - sz[1]
                dxyz2 = dyz2 + relu(abs(pos_norm[1] - cell_center_x) - one(F))^2
                dxyz2 ≤ r2_cells || continue
                li = ix + liy + liz
                # li is inbounds because lo and hi are clamped to the grid,
                # and cells has length prod(sz). cs.len ≤ length(cs.list) is
                # a _cell_push_entry! invariant.
                cs = @inbounds pcl.cells[begin + li]
                for j in 1:cs.len
                    entry = @inbounds cs.list[j]
                    diff = entry.pos - pos
                    d2 = diff ⋅ diff
                    d2 ≤ cutoff2 || continue
                    out, cont = f(entry, out)
                    cont || return out
                end
            end
        end
    end
    out
end
