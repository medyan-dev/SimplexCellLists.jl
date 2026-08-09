
"""
    CellLineSegEntry{T,F}

A single line segment's contribution to one voxel cell. Created by [`cell_line_seg_add!`](@ref)
and read by [`map_nearby_line_segs`](@ref).

# Fields
- `id::T`: line segment ID.
- `p0::SVector{3, F}`: line segment start point.
- `d_hat::SVector{3, F}`: unit direction vector along the line segment. Zero for zero-length line segments.
- `tmin::F`: start of the owned parametric range along the line segment (length units from `p0`).
- `packed_tmax::F`: end of the owned parametric range (length units from `p0`). The sign bit
  encodes the `is_end` flag; use [`cell_line_seg_tmax`](@ref) for the absolute value and
  [`cell_line_seg_is_end`](@ref) for the flag.
"""
struct CellLineSegEntry{T,F}
    id::T
    p0::SVector{3, F}
    d_hat::SVector{3, F}
    tmin::F
    packed_tmax::F
end
function CellLineSegEntry{T,F}(id::T, is_end::Bool,
                      p0::SVector{3, F}, d_hat::SVector{3, F},
                      tmin::F, tmax::F) where {T,F}
    # Pack is_end flag into sign bit of tmax (-0.0 encodes is_end when tmax==0)
    packed_tmax = ifelse(is_end, -abs(tmax), abs(tmax))
    CellLineSegEntry{T,F}(id, p0, d_hat, tmin, packed_tmax)
end

"""
    cell_line_seg_is_end(e::CellLineSegEntry)::Bool

Return true if `e` is the entry containing the end point of its line segment.
"""
cell_line_seg_is_end(e::CellLineSegEntry) = signbit(e.packed_tmax)

"""
    cell_line_seg_tmax(e::CellLineSegEntry{T,F})::F

Return the end of the owned parametric range of `e` (length units from `e.p0`).
"""
cell_line_seg_tmax(e::CellLineSegEntry) = abs(e.packed_tmax)

struct CellLineSeg{T,F}
    list::Memory{CellLineSegEntry{T,F}}
    len::Int
end
CellLineSeg{T,F}() where {T,F} = CellLineSeg{T,F}(Memory{CellLineSegEntry{T,F}}(), 0)

"""
    LineSegCellList{T,F}(size::SVector{3, Int}, grid_spacing::F)
    LineSegCellList{T,F}(size::NTuple{3, Int}, grid_spacing::F)

An origin-centered 3D voxel grid for spatial queries on line segments. The grid spans
`[-size * grid_spacing/2, size * grid_spacing/2)` in each dimension.

Insert line segments with [`cell_line_seg_add!`](@ref) and query with
[`map_nearby_line_segs`](@ref). Use [`cell_line_segs_clear!`](@ref) to reset.

Positions are not required to be within the grid bounds: parts of a line segment
outside the grid are stored in the nearest boundary cells, and queries account
for this. Line segments far outside the grid accumulate in boundary cells, which
can slow down queries near those boundaries.

# Arguments
- `size`: number of cells along each axis (x, y, z).
- `grid_spacing::F`: side length of each cubic cell.
"""
struct LineSegCellList{T,F}
    cells::Memory{CellLineSeg{T,F}}
    size::SVector{3, Int}
    inv_half_grid_spacing::F
end
function LineSegCellList{T,F}(size::SVector{3, Int}, grid_spacing::F) where {T,F}
    @argcheck all(≥(0), size)
    @argcheck grid_spacing > 0
    inv_half_grid_spacing = 2*inv(grid_spacing)
    @argcheck isfinite(inv_half_grid_spacing)
    n = prod(size)
    cells = Memory{CellLineSeg{T,F}}(undef, n)
    fill!(cells, CellLineSeg{T,F}())
    LineSegCellList{T,F}(cells, size, inv_half_grid_spacing)
end
LineSegCellList{T,F}(size::NTuple{3, Int}, grid_spacing::F) where {T,F} = LineSegCellList{T,F}(SVector{3,Int}(size...), grid_spacing)

"""
    cell_line_segs_clear!(scl::LineSegCellList)::Nothing

Remove all line segments.
"""
function cell_line_segs_clear!(scl::LineSegCellList{T,F})::Nothing where {T,F}
    for i in eachindex(scl.cells)
        scl.cells[i] = CellLineSeg{T,F}(scl.cells[i].list, 0)
    end
end

"""
    cell_line_seg_add!(scl::LineSegCellList{T,F}, p0::SVector{3, F}, p1::SVector{3, F}, id::T)::Nothing

Insert a line segment into the cell list, adding a [`CellLineSegEntry`](@ref) to every cell the
line segment passes through. Each entry records the parametric range of the line segment owned
by its cell; the ranges partition the line segment so [`map_nearby_line_segs`](@ref) reports it
exactly once.

Uses the Amanatides & Woo "Fast Voxel Traversal" algorithm (1987) for 3D DDA traversal.

Positions are not required to be within the grid bounds: parts of the line segment outside the
grid are owned by the nearest boundary cells, and queries account for this. Line segments far
outside the grid accumulate in boundary cells, which can slow down queries near those boundaries.

Throws an `ArgumentError` if the length of the line segment is not finite: for example if a
coordinate is NaN or Inf, or if the squared length overflows `F`.

# Arguments
- `scl`: the line segment cell list to insert into.
- `p0`: line segment start point.
- `p1`: line segment end point.
- `id`: ID for this line segment.
"""
function cell_line_seg_add!(scl::LineSegCellList{T,F}, p0::SVector{3, F}, p1::SVector{3, F}, id::T)::Nothing where {T,F}
    h2_inv = scl.inv_half_grid_spacing
    sz = scl.size
    @argcheck all(>(0), sz)
    strides = SVector(1, sz[1], sz[1]*sz[2])
    d = p1 - p0
    seg_len = norm_fast(d)
    @argcheck isfinite(seg_len)
    invseglen = inv(seg_len)
    d_hat = ifelse(isfinite(invseglen), d * invseglen, zero(d))
    # Work in doubled coordinates 2*pos/h so the origin-centering offset is
    # added exactly in Int after the floor, keeping cell assignment accurate
    # for any grid size. Clamp in F before converting to Int so far away
    # positions don't overflow Int; the prevfloat upper bound keeps the
    # resulting index ≤ size - 1.
    q0 = p0 * h2_inv
    i = round2grid(clamp.(q0, -F.(sz), prevfloat.(F.(sz))), sz)
    step_signbit = signbit.(d)
    step = ifelse.(step_signbit, -1, +1)
    # tDelta: how much t to cross one cell (2 doubled units).
    tDelta = 2 .* inv.(abs.(d_hat * h2_inv)) # h = |d̂ᵢ| * tΔᵢ, tΔᵢ = h/|d̂ᵢ|, tΔᵢ = 2/(|d̂ᵢ*2*inv(h)|)
    # tMax: t at which the next cell boundary is crossed per axis. Boundary
    # cells own everything outside the grid (matching the clamping above), so
    # an axis already at the wall in its step direction never crosses again.
    at_wall = ifelse.(step_signbit, i .== 0, i .== sz .- 1)
    next_boundary = 2 .* (i .+ .!step_signbit) .- sz
    tMax_raw = abs.(next_boundary .- q0) .* tDelta .* F(0.5)
    # Replace NaN with Inf (from 0*Inf when the segment is parallel to an axis
    # exactly at a cell boundary)
    tMax_1, tMax_2, tMax_3 = ifelse.(at_wall .| isnan.(tMax_raw), F(Inf), tMax_raw)
    # DDA loop
    t_prev = zero(F) # parametric t at entry to current cell
    li = i ⋅ strides
    X, Y, Z = i
    while true
        t_exit = min(tMax_1, tMax_2, tMax_3, seg_len)
        is_last = !(t_exit < seg_len)
        entry = CellLineSegEntry{T,F}(id, is_last, p0, d_hat, t_prev, t_exit)
        # li is inbounds because the start cell is clamped into the grid and
        # the wall checks above stop stepping at the boundaries, so X, Y, Z
        # stay in [0, sz-1].
        @inbounds scl.cells[begin + li] = _cell_push_entry!(scl.cells[begin + li], entry)
        is_last && break
        # Step to next cell: advance the axis with smallest tMax.
        t_prev = t_exit
        if tMax_1 ≤ tMax_2 && tMax_1 ≤ tMax_3
            X += step[1]
            li += step[1] * strides[1]
            tMax_1 = ifelse(ifelse(step_signbit[1], X == 0, X == sz[1]-1), F(Inf), tMax_1 + tDelta[1])
        elseif tMax_2 ≤ tMax_3
            Y += step[2]
            li += step[2] * strides[2]
            tMax_2 = ifelse(ifelse(step_signbit[2], Y == 0, Y == sz[2]-1), F(Inf), tMax_2 + tDelta[2])
        else
            Z += step[3]
            li += step[3] * strides[3]
            tMax_3 = ifelse(ifelse(step_signbit[3], Z == 0, Z == sz[3]-1), F(Inf), tMax_3 + tDelta[3])
        end
    end
    nothing
end


"""
    map_nearby_line_segs(f, scl::LineSegCellList{T,F}, pos::SVector{3, F}, cutoff::F, out) -> out

Call `f` for every line segment in `scl` whose closest point to `pos` is within `cutoff` distance.
Each line segment is reported exactly once (deduplicated across cells).

# Arguments
- `f(entry::CellLineSegEntry{T,F}, out) -> (out, cont::Bool)`: callback called for each nearby line segment.
  Return `(new_out, true)` to continue, or `(new_out, false)` to stop early.
- `scl`: the line segment cell list to query.
- `pos`: query point.
- `cutoff`: maximum distance. If zero, returns `out` immediately.
- `out`: initial accumulator value, threaded through `f` and returned.
"""
function map_nearby_line_segs(
        f,
        scl::LineSegCellList{T,F},
        pos::SVector{3, F},
        cutoff::F,
        out,
    ) where {T,F}
    cutoff ≤ 0 && return out
    sz = scl.size
    all(>(0), sz) || return out
    cutoff2 = cutoff * cutoff
    h2_inv = scl.inv_half_grid_spacing
    strides = SVector(1, sz[1], sz[1]*sz[2])
    # Clamping the query into the grid box makes the cube test below stay
    # conservative given line segments can extend outside the grid.
    # Clamping in F before converting to Int also keeps far away
    # queries and huge cutoffs from overflowing Int.
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
                cs = @inbounds scl.cells[begin + li]
                for j in 1:cs.len
                    entry = @inbounds cs.list[j]
                    # closest point parameter on the whole line segment
                    tmax = cell_line_seg_tmax(entry)
                    pos_m_p0 = pos - entry.p0
                    t = clamp(pos_m_p0 ⋅ entry.d_hat, zero(F), tmax)
                    # distance check
                    diff = t * entry.d_hat - pos_m_p0
                    d2 = diff ⋅ diff
                    d2 ≤ cutoff2 || continue
                    # ownership check: tmin ≤ t < tmax, or is_end
                    t ≥ entry.tmin || continue
                    (t < tmax || cell_line_seg_is_end(entry)) || continue
                    out, cont = f(entry, out)
                    cont || return out
                end
            end
        end
    end
    out
end
