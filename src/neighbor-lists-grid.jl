# Cells are at least this many times the smallest `r + skin` wide along each axis.
const GRID_MIN_CELL_WIDTH_FACTOR = 4
# The grid has at most this many cells per object that passes `filter_object`.
const GRID_MAX_CELLS_PER_OBJECT = 1
# The grid has at most this many cells along each axis, so a cell coordinate
# fits in a `UInt8`.
const GRID_MAX_AXIS_CELLS = 256
# Each cell is split into `GRID_CELL_UNITS = 2^GRID_CELL_BITS` unit
# intervals along each axis, so a unit interval index is a `UInt16` with the
# cell index in the high byte and the unit interval within the cell in the
# low byte.
const GRID_CELL_BITS = 0x08
const GRID_CELL_UNITS = 1 << GRID_CELL_BITS

# A padded bounding box of one object quantized to the grid, 128 bits.
# Each axis is split into `GRID_CELL_UNITS * n_cells` unit intervals,
# numbered from zero, and the box is assigned the unit intervals it touches,
# from `lo` to `hi`. At level `k` it touches the cells `_grid_cells(lo, k)`
# to `_grid_cells(hi, k)`.
# Thanks to conservative rounding, two boxes overlap on an axis iff
#     a.lo ≤ b.hi && b.lo ≤ a.hi
struct GridAABB
    lo::SVector{3, UInt16}
    hi::SVector{3, UInt16}
    # Index into the object array of this object's type
    index::UInt32
end

# Return the first and last unit intervals touched by the bounds `lo` and
# `hi` on one axis, padded by `pad`. The axis starts at `base` and has
# `scale` unit intervals per length, numbered from zero to `last_unit`.
# Rounding down is conservative since it is monotone, clamping to the axis.
# NaN geometry errors.
function _quantize_grid_bounds(lo, hi, pad, base, scale, last_unit)::NTuple{2, UInt16}
    # Check the sign before scaling, so zero times the infinite scale of a
    # zero width axis isn't NaN.
    unit(x) = x ≤ 0 ? zero(UInt16) : trunc(UInt16, min(x * scale, last_unit))
    (unit((lo - base) - pad), unit((hi - base) + pad))
end

# Return the bounding box of `obj` padded by `pad` and quantized to the grid.
Base.@propagate_inbounds function _grid_aabb(
        pos,
        obj,
        pad,
        base::SVector{3},
        scale::SVector{3},
        last_unit::SVector{3},
        index::UInt32
    )::GridAABB
    local bounds = ntuple(Val(3)) do axis
        local lo, hi = _load_axis_bounds(pos, obj, axis)
        _quantize_grid_bounds(lo, hi, pad, base[axis], scale[axis], last_unit[axis])
    end
    GridAABB(SVector(map(first, bounds)), SVector(map(last, bounds)), index)
end

# Return the cell coordinates at level `k` of the unit intervals `units`.
@inline function _grid_cells(units::SVector{3, UInt16}, k::UInt8)::SVector{3, Int}
    Int.(units .>> (GRID_CELL_BITS + k))
end

# A `GridAABB` stored in a cell as 8 `UInt16` lanes,
#     (hi..., ~lo..., low half of index, high half of index)
# so checking a box against it is a few operations on all lanes at once,
# see `_grid_candidate`.
const GridEntry = SIMD.Vec{8, UInt16}

@inline function _grid_entry(a::GridAABB)::GridEntry
    GridEntry((a.hi[1], a.hi[2], a.hi[3], ~a.lo[1], ~a.lo[2], ~a.lo[3], a.index % UInt16, (a.index >> 16) % UInt16))
end

@inline function _grid_entry_index(e::GridEntry)::UInt32
    UInt32(e[7]) | (UInt32(e[8]) << 16)
end

# Return `(lo..., ~hi..., 0, 0)` of a box. Compared with an entry `b`, the
# first lanes are `lo ≤ b.hi`, the next are `~hi ≤ ~b.lo` so `b.lo ≤ hi`,
# and the last are always true, so `all(query ≤ b)` is the overlap.
@inline function _grid_query(a::GridAABB)::GridEntry
    GridEntry((a.lo[1], a.lo[2], a.lo[3], ~a.hi[1], ~a.hi[2], ~a.hi[3], 0x0000, 0x0000))
end
@inline function _grid_query(e::GridEntry)::GridEntry
    SIMD.shufflevector(~e, zero(GridEntry), Val((3, 4, 5, 0, 1, 2, 8, 8)))
end

# Shift right by `s`, which is 16 at the coarsest level of an axis with 256
# cells. That gives zero like Julia's `>>`, but isn't defined in LLVM.
@inline function _grid_entry_shift(x::GridEntry, s::UInt8)::GridEntry
    s ≥ 0x10 ? zero(x) : x >> s
end

# Lanes 3 to 5 hold the cell coordinates in the dedup check, the others are
# always true.
const GRID_DEDUP_LANES = SIMD.Vec{8, Bool}((true, true, true, false, false, false, true, true))

# Return the lanes where the lower cell of a box at the shift `s` is the cell
# `c`, given the lower bound of the box in lanes 3 to 5 of `lo_lanes` and the
# cell coordinates in lanes 3 to 5 of `c`.
@inline function _grid_lower_in_cell(lo_lanes::GridEntry, s::UInt8, c::GridEntry)::SIMD.Vec{8, Bool}
    (_grid_entry_shift(lo_lanes, s) == c) | GRID_DEDUP_LANES
end

# Return whether the entry `b` overlaps the box with `query`, and the cell
# `c` contains the lower corner of their intersection. A pair can touch
# several common cells, so it is only checked in that one. Both boxes touch
# `c`, so on each axis the larger of their lower cells is `c` iff either is.
@inline function _grid_candidate(query::GridEntry, a_in_cell::SIMD.Vec{8, Bool}, b::GridEntry, s::UInt8, c::GridEntry)::Bool
    all((query <= b) & (_grid_lower_in_cell(~b, s, c) | a_in_cell))
end
# TODO we could try to sneak in the results of _grid_lower_in_cell in 3 bits of idx

# One pair type, like in `setup_neighbors_naive!`: the neighbor list `nl`
# its edges go to, the exclusions `excl`, whether both objects are the same
# type `self`, and the `a` and `b` objects as
# `(object type, objects, radii, params)`.
struct GridPairType{NL, E, A, B}
    nl::NL
    excl::E
    self::Bool
    a::A
    b::B
end

# Initial number of candidate pairs buffered before the narrow phase runs,
# small enough that the buffer stays in cache.
const GRID_CANDIDATE_CHUNK = 1024
# Pairs checked at once by the narrow phase distance check.
const GRID_SIMD_WIDTH = 8

# The first `n` of `pairs` are the candidate pairs of one pair type found by
# the broad phase, waiting for the narrow phase, as `i => j` with `i` the
# index of the `a` object and `j` of the `b` object, in the order of
# `NeighborListEdge`, which for pairs of the same type is the smaller index
# first. The last `GRID_SIMD_WIDTH - 1` of `pairs` are padding, so the
# narrow phase can always check whole groups of pairs, see
# `_grid_distance_pass!`. `d2` is scratch space for the narrow phase, as
# long as `pairs`, with the element type `T` of the positions.
mutable struct GridCandidates{T}
    pairs::Memory{Pair{UInt32, UInt32}}
    d2::Memory{T}
    n::Int
end
GridCandidates{T}() where T = GridCandidates{T}(GRID_CANDIDATE_CHUNK)
# Room for `capacity` candidates
function GridCandidates{T}(capacity::Int) where T
    local len = capacity + GRID_SIMD_WIDTH - 1
    GridCandidates{T}(Memory{Pair{UInt32, UInt32}}(undef, len), Memory{T}(undef, len), 0)
end

# Return the number of candidates `cands` has room for, without the padding.
_grid_capacity(cands::GridCandidates)::Int = length(cands.pairs) - (GRID_SIMD_WIDTH - 1)

# Make room for `extra` more candidates of `pair_type`, by running the
# narrow phase if they don't fit, and growing the buffer if they don't fit
# even when it is empty. The broad phase calls this before each run of
# entries it checks, so it never needs to check for room inside the run.
# The overlap test is unpredictable, so the broad phase appends without
# branching on it: it writes every pair it checks at `n + 1` and only
# advances `n` for overlapping pairs, which needs room for the whole run.
@inline function _grid_make_room!(cands::GridCandidates, extra::Int, pos, policy, extra_cutoff, pair_type::GridPairType)
    if cands.n + extra > _grid_capacity(cands)
        _grid_narrow_phase!(cands, pos, policy, extra_cutoff, pair_type)
        @assert iszero(cands.n)
        # Over allocate, so later runs from crowded cells fit too. More than
        # doubles the buffer, since `extra` is more than its capacity.
        if extra > _grid_capacity(cands)
            local grown = typeof(cands)(2 * extra)
            cands.pairs = grown.pairs
            cands.d2 = grown.d2
        end
    end
    nothing
end

# Return `pairs[k:k + W - 1]` as a tuple.
@inline function _grid_load_pairs(pairs::Memory{Pair{UInt32, UInt32}}, k::Int, ::Val{W}) where W
    ntuple(l -> @inbounds(pairs[k + l - 1]), Val(W))
end

# Return the simplices of the objects `objs[idx]` as one simplex whose
# coordinates have a lane per object. The indexes were checked when the
# objects were quantized.
@inline function _grid_lanes(pos, objs, idx::Tuple{UInt32, Vararg{UInt32, M}}) where M
    local simplices = ntuple(l -> @inbounds(load_positions(pos, objs[idx[l]])), Val(M + 1))
    SVector(ntuple(Val(length(first(simplices)))) do v
        SVector(ntuple(Val(3)) do axis
            SIMD.Vec(ntuple(l -> simplices[l][v][axis], Val(M + 1)))
        end)
    end)
end

# Keep the buffered candidates within the cutoff distance, compacting them
# to the first `m` for the returned `m`. The distance check is the one of
# `setup_neighbors_naive!`, done `GRID_SIMD_WIDTH` pairs at a time, which
# gives identical results since the distance functions do the same
# operations on every lane. The last group is padded with copies of the
# last pair.
function _grid_distance_pass!(cands::GridCandidates, pos, extra_cutoff, a_objs, a_r, b_objs, b_r)::Int
    local pairs = cands.pairs
    local d2s = cands.d2
    local n = cands.n
    iszero(n) && return 0
    local last_pair = @inbounds pairs[n]
    for k in (n + 1):(n + GRID_SIMD_WIDTH - 1)
        @inbounds pairs[k] = last_pair
    end
    # The distances go to their own loop and buffer, so the groups are
    # independent. Checking the cutoff and compacting in this loop makes
    # each group wait for the previous one, and makes LLVM load the
    # indexes of a group as a vector and extract them one by one before
    # loading the positions, both about twice as slow.
    for k in 1:GRID_SIMD_WIDTH:n
        local group = _grid_load_pairs(pairs, k, Val(GRID_SIMD_WIDTH))
        local d2 = dist_sqr(_grid_lanes(pos, a_objs, map(first, group)), _grid_lanes(pos, b_objs, map(last, group)))
        @inbounds SIMD.vstore(d2, d2s, k)
    end
    # Branchless compaction, the writes never pass the reads. The cutoff
    # check is per pair with scalar loads of the radii, since LLVM turns
    # the vector version into slow gather instructions.
    local m = 0
    for k in 1:n
        local (i, j) = @inbounds pairs[k]
        local cutoff = @inbounds a_r[i] + b_r[j] + extra_cutoff
        @inbounds pairs[m + 1] = i => j
        m += @inbounds(d2s[k]) < cutoff^2
    end
    m
end

# Push the edges of the buffered candidates that pass the distance check,
# `filter_pair`, and the exclusions, like `setup_neighbors_naive!`, and
# empty the buffer. The distance check comes first, which `filter_pair`
# allows.
function _grid_narrow_phase!(cands::GridCandidates, pos, policy, extra_cutoff, pair_type::GridPairType)
    local (; nl, excl) = pair_type
    local (a_type, a_objs, a_r, a_params) = pair_type.a
    local (b_type, b_objs, b_r, b_params) = pair_type.b
    local m = _grid_distance_pass!(cands, pos, extra_cutoff, a_objs, a_r, b_objs, b_r)
    for k in 1:m
        local (i, j) = @inbounds cands.pairs[k]
        local r_i = a_r[i]
        local r_j = b_r[j]
        local params_i = a_params[i]
        local params_j = b_params[j]
        filter_pair(policy, params_i, params_j, a_type, b_type, i, j, r_i, r_j) || continue
        (i => j) in excl && continue
        local pair_params = mix_params(
            policy, params_i, params_j, a_type, b_type, i, j, r_i, r_j)
        push!(nl, NeighborListEdge(a_objs[i], b_objs[j], r_i + r_j, pair_params))
    end
    cands.n = 0
    nothing
end

# Store the quantized padded bounding boxes of the objects of one type that
# pass `filter_object` in `aabbs[range]`, in index order.
function _grid_quantize_type!(aabbs, policy, pos, (obj_type, a_objs, a_r, a_params), range,
        skin_with_eps, base, scale, last_unit,
    )
    local k = first(range)
    for index in UInt32(1):UInt32(length(a_objs))
        local r = a_r[index]
        filter_object(policy, a_params[index], obj_type, index, r) || continue
        aabbs[k] = _grid_aabb(pos, a_objs[index], r + skin_with_eps, base, scale, last_unit, index)
        k += 1
    end
    # filter_object must be pure
    @assert k == last(range) + 1
    nothing
end

# The grid has levels, level `k` has cells `2^k` times wider than level 0
# along each axis, so its cell of the unit interval `x` is
# `x >> (GRID_CELL_BITS + k)`.
# Each object is stored in at most 8 adjacent cells.
struct GridLevel
    n_cells::SVector{3, Int}
    # The entries of cell `c` start after `values[cell_offsets[begin + c]]`.
    cell_offsets::Vector{Int}
    # Entries of type `t` in cell `c` are from `type_offsets[t, begin + c]+1`
    # to `type_offsets[t + 1, begin + c]` within the cell's entries. Local to
    # the cell, so a `UInt32` is enough, and a cell's column is 32 bytes.
    type_offsets::Matrix{UInt32}
    values::Vector{GridEntry}
end

Base.@propagate_inbounds function _grid_type_range(level::GridLevel, t::Int, cell::Int)::UnitRange{Int}
    local start = level.cell_offsets[begin + cell]
    (start + level.type_offsets[t, begin + cell] + 1):(start + level.type_offsets[t + 1, begin + cell])
end

# Return the level a box `aabb` is added to: the finest level whose cells are
# wider than the box on every axis, capped at the coarsest level, which has
# one cell. So the box touches at most two cells along each axis wherever it
# is. A box aligned with the cells can touch at most two cells at a finer
# level too, e.g. a box 300 units wide from unit 0 touches two cells at level
# 0, but is added to level 1.
@inline function _grid_level(aabb::GridAABB, n_levels::UInt8)::UInt8
    # Less than 256 cells, so fits in a `UInt8`
    local size = (maximum(aabb.hi .- aabb.lo) >> GRID_CELL_BITS) % UInt8
    min(Base.top_set_bit(size) % UInt8, n_levels - 0x01)
end

# Return the levels of the grid with each object binned into the cells it
# touches at its level `level_of[k]`, type by type.
function _grid_levels(aabbs, level_of, type_ranges, n_cells, n_levels::UInt8)::Vector{GridLevel}
    # Each object is in a cell at most once, so the `UInt32` local type
    # offsets can't overflow
    @argcheck length(aabbs) ≤ typemax(UInt32)
    local levels = map(0x00:(n_levels - 0x01)) do k
        # The last cell at level 0 is `n_cells - 1`, so at level `k` it is
        # `(n_cells - 1) >> k`
        local n_k = ((n_cells .- 1) .>> k) .+ 1
        GridLevel(n_k,
            zeros(Int, prod(n_k)),
            zeros(UInt32, N_COLLIDE_OBJECT_TYPES + 1, prod(n_k)),
            GridEntry[],
        )
    end
    # Count the objects in each cell they touch, which is one or two cells
    # along each axis at their level. Which one is hard to predict, so this
    # is branchless: always do the 8 corner increments, with x fastest,
    # adding zero for corners the box doesn't touch, clamped to a cell it
    # does touch. Cells are in bounds since `_grid_cells(hi, l) < n_cells`.
    # TODO try cells in Morton order, (might help locality of multi level searches)
    for k in eachindex(aabbs)
        local aabb = aabbs[k]
        local l = level_of[k]
        local level = levels[begin + l]
        local n_cells_l = level.n_cells
        local lo = _grid_cells(aabb.lo, l)
        local span = _grid_cells(aabb.hi, l) .- lo
        @assert all(span .≤ 1)
        for dz in 0:1
            local plane = n_cells_l[2] * (lo[3] + min(dz, span[3]))
            for dy in 0:1
                local row = n_cells_l[1] * (lo[2] + min(dy, span[2]) + plane)
                for dx in 0:1
                    local cell = row + lo[1] + min(dx, span[1])
                    local touched = (dx ≤ span[1]) & (dy ≤ span[2]) & (dz ≤ span[3])
                    @inbounds level.cell_offsets[begin + cell] += touched
                end
            end
        end
    end
    # Each cell starts where the previous one ends
    for level in levels
        local stop = 0
        for c in 0:prod(level.n_cells)-1
            local count = level.cell_offsets[begin + c]
            level.cell_offsets[begin + c] = stop
            stop += count
        end
        resize!(level.values, stop)
    end
    # Fill type by type. Row `t + 1` starts as a copy of row `t` and is the
    # fill position of type `t`, so it ends as the start of type `t + 1`.
    for t in 1:N_COLLIDE_OBJECT_TYPES
        for level in levels
            @views level.type_offsets[t + 1, :] .= level.type_offsets[t, :]
        end
        for k in type_ranges[t]
            local aabb = aabbs[k]
            local entry = _grid_entry(aabb)
            local l = level_of[k]
            local level = levels[begin + l]
            local n_cells_l = level.n_cells
            local lo = _grid_cells(aabb.lo, l)
            local hi = _grid_cells(aabb.hi, l)
            # Not branchless like counting, since the extra stores of
            # untouched corners cost more than the mispredictions here
            for cz in lo[3]:hi[3]
                local plane = n_cells_l[2] * cz
                for cy in lo[2]:hi[2]
                    local row = n_cells_l[1] * (cy + plane)
                    for cx in lo[1]:hi[1]
                        local cell = row + cx
                        local e = level.type_offsets[t + 1, begin + cell] += 0x00000001
                        level.values[level.cell_offsets[begin + cell] + e] = entry
                    end
                end
            end
        end
    end
    levels
end

# Buffer the candidate pairs of one pair type between objects at the same
# level, by checking the pairs of entries in each cell. For pairs of the
# same type, each pair of entries is checked once, and the pair is from the
# smaller index, like `i < j` in `setup_neighbors_naive!`.
function _grid_same_level_pairs!(levels, counts, cands, pos, policy, extra_cutoff, pair_type::GridPairType)
    local self = pair_type.self
    local A = Int(first(pair_type.a))
    local B = Int(first(pair_type.b))
    for k in 0x00:UInt8(length(levels) - 1)
        (iszero(counts[A, begin + k]) || iszero(counts[B, begin + k])) && continue
        local level = levels[begin + k]
        local n_cells = level.n_cells
        local s = GRID_CELL_BITS + k
        for cz in 0:n_cells[3]-1
            local plane = n_cells[2] * cz
            for cy in 0:n_cells[2]-1
                local row = n_cells[1] * (cy + plane)
                for cx in 0:n_cells[1]-1
                    local cell = row + cx
                    local a_range = _grid_type_range(level, A, cell)
                    local b_range = _grid_type_range(level, B, cell)
                    isempty(b_range) && continue
                    # Cell coordinates in the lower bound lanes
                    local c = GridEntry((0x0000, 0x0000, 0x0000, cx % UInt16, cy % UInt16, cz % UInt16, 0x0000, 0x0000))
                    for p in a_range
                        local a_entry = level.values[p]
                        local a_query = _grid_query(a_entry)
                        local a_in_cell = _grid_lower_in_cell(~a_entry, s, c)
                        local i = _grid_entry_index(a_entry)
                        local q_range = self ? ((p + 1):last(a_range)) : b_range
                        _grid_make_room!(cands, length(q_range), pos, policy, extra_cutoff, pair_type)
                        local n = cands.n
                        for q in q_range
                            local b_entry = level.values[q]
                            local hit = _grid_candidate(a_query, a_in_cell, b_entry, s, c)
                            local j = _grid_entry_index(b_entry)
                            local (i_a, i_b) = self ? minmax(i, j) : (i, j)
                            @inbounds cands.pairs[n + 1] = i_a => i_b
                            n += hit
                        end
                        cands.n = n
                    end
                end
            end
        end
    end
    nothing
end

# Buffer the candidate pairs of `pair_type` between `x_aabb` at a coarser
# level and the entries of type `Y` in the cells it touches at the finer
# level `k`. The pair is from `x` if `x_is_a`, or from the smaller index if
# both are the same type.
function _grid_check_finer!(cands, pos, policy, extra_cutoff, pair_type::GridPairType,
        level::GridLevel, k::UInt8, x_aabb, Y::Int, x_is_a::Bool,
    )
    local self = pair_type.self
    local i = x_aabb.index
    local x_query = _grid_query(x_aabb)
    # The lower bound of `x` in the lower bound lanes
    local x_lo = GridEntry((0x0000, 0x0000, 0x0000, x_aabb.lo[1], x_aabb.lo[2], x_aabb.lo[3], 0x0000, 0x0000))
    local s = GRID_CELL_BITS + k
    local lo = _grid_cells(x_aabb.lo, k)
    local hi = _grid_cells(x_aabb.hi, k)
    local n_cells = level.n_cells
    for cz in lo[3]:hi[3]
        local plane = n_cells[2] * cz
        for cy in lo[2]:hi[2]
            local row = n_cells[1] * (cy + plane)
            for cx in lo[1]:hi[1]
                local cell = row + cx
                # Cell coordinates in the lower bound lanes
                local c = GridEntry((0x0000, 0x0000, 0x0000, cx % UInt16, cy % UInt16, cz % UInt16, 0x0000, 0x0000))
                local x_in_cell = _grid_lower_in_cell(x_lo, s, c)
                local q_range = _grid_type_range(level, Y, cell)
                _grid_make_room!(cands, length(q_range), pos, policy, extra_cutoff, pair_type)
                local n = cands.n
                for q in q_range
                    local y_entry = level.values[q]
                    local hit = _grid_candidate(x_query, x_in_cell, y_entry, s, c)
                    local j = _grid_entry_index(y_entry)
                    local x_first = self ? (i < j) : x_is_a
                    @inbounds cands.pairs[n + 1] = x_first ? (i => j) : (j => i)
                    n += hit
                end
                cands.n = n
            end
        end
    end
    nothing
end

# Buffer the candidate pairs of `pair_type` between objects of its `a` type
# if `x_is_a` or its `b` type otherwise, and the entries of the other type
# at finer levels.
function _grid_check_finer_levels!(levels, counts, level_of, aabbs, type_ranges, cands, pos, policy, extra_cutoff,
        pair_type::GridPairType, x_is_a::Bool,
    )
    local a_type = first(pair_type.a)
    local b_type = first(pair_type.b)
    local X = Int(x_is_a ? a_type : b_type)
    local Y = Int(x_is_a ? b_type : a_type)
    # if X doesn't have any boxes at a non finest level then skip
    all(iszero, @view(counts[X, begin+1:end])) && return
    for k in type_ranges[X]
        # TODO ideally we go through objects in the order they have been sorted into by the grid.
        local m = level_of[k]
        iszero(m) && continue
        local x_aabb = aabbs[k]
        for l in 0x00:(m - 0x01)
            iszero(counts[Y, begin + l]) && continue
            _grid_check_finer!(cands, pos, policy, extra_cutoff, pair_type,
                levels[begin + l], l, x_aabb, Y, x_is_a,
            )
        end
    end
    nothing
end

# Buffer the candidate pairs of one pair type between objects at different
# levels, found from the coarser object.
function _grid_cross_level_pairs!(levels, counts, level_of, aabbs, type_ranges, cands, pos, policy, extra_cutoff, pair_type::GridPairType)
    _grid_check_finer_levels!(levels, counts, level_of, aabbs, type_ranges, cands, pos, policy, extra_cutoff,
        pair_type, true,
    )
    if !pair_type.self
        _grid_check_finer_levels!(levels, counts, level_of, aabbs, type_ranges, cands, pos, policy, extra_cutoff,
            pair_type, false,
        )
    end
    nothing
end

# Push the edges of one pair type. The broad phase buffers the pairs of
# objects whose padded bounding boxes overlap: pairs at the same level are
# found within each cell, and pairs at different levels from the coarser
# object, so each pair is found once. The narrow phase pushes the edges of
# the buffered pairs, whenever the buffer is full and once at the end.
function _grid_pairs!(levels, counts, level_of, aabbs, type_ranges, cands, pos, policy, extra_cutoff, pair_type::GridPairType)
    _grid_same_level_pairs!(levels, counts, cands, pos, policy, extra_cutoff, pair_type)
    _grid_cross_level_pairs!(levels, counts, level_of, aabbs, type_ranges, cands, pos, policy, extra_cutoff, pair_type)
    _grid_narrow_phase!(cands, pos, policy, extra_cutoff, pair_type)
    nothing
end

# Return the smallest radius and the number of objects that pass `filter_object`.
function _filtered_min_radius(policy, obj_type::CollideObjectTypes, radii, params)::Tuple{Float32, Int}
    local min_r = Inf32
    local n = 0
    for index in UInt32(1):UInt32(length(radii))
        local r = radii[index]
        filter_object(policy, params[index], obj_type, index, r) || continue
        min_r = min(min_r, r)
        n += 1
    end
    (min_r, n)
end

# Return the number of cells `n` along each axis of a grid spanning `width`,
# such that `1 ≤ n[i] ≤ GRID_MAX_AXIS_CELLS`,
# `n[i] ≤ max(1, floor(width[i] / min_cell_width))`, and
# `prod(n) ≤ max_cells`. Cells don't need to be cubic, since objects are binned
# by their bounding boxes. The minimum is a 1x1x1 grid.
#
# Starting from `n_max`, conceptually increase a cell width h, merging cells
# on every axis where they are narrower than h, until `prod(n) ≤ max_cells`.
# Axes with equal cell widths merge together, so this can leave some of
# `max_cells` unused.
function _grid_size(width::SVector{3}, min_cell_width, max_cells::Int)::SVector{3, Int}
    @argcheck max_cells ≥ 1
    # This also rejects NaN
    @argcheck 0 < min_cell_width < Inf
    # Infinite widths, from overflow in the bounds, get a 1x1x1 grid.
    all(isfinite, width) || return SA[1, 1, 1]
    # The largest n with cells at least h wide, within the axis limits.
    # w / h is never NaN for h ≥ min_cell_width since w ≥ 0 is never NaN and
    # h is finite and positive.
    merged(h) = map(w -> Int(clamp(floor(w / h), 1, GRID_MAX_AXIS_CELLS)), width)
    n_max = merged(min_cell_width)
    # No merging needed, this is the common case for dense systems.
    prod(n_max) ≤ max_cells && return n_max
    # Bisect for the smallest h where the merged grid fits.
    # h_lo gives n_max, which doesn't fit, and h_hi gives a 1x1x1 grid, which
    # does.
    h_lo, h_hi = promote(min_cell_width, maximum(width))
    n_lo = n_max
    n = SA[1, 1, 1]
    # Stop once n_lo and n differ by at most one cell on each axis, and pick
    # the smaller n. This avoids bisecting down to adjacent floats when
    # several axes merge at the same h. The h_mid check guards against
    # getting stuck on rounding.
    while any(n_lo .- n .> 1)
        h_mid = sqrt(h_lo) * sqrt(h_hi)
        h_lo < h_mid < h_hi || break
        n_mid = merged(h_mid)
        if prod(n_mid) ≤ max_cells
            h_hi, n = h_mid, n_mid
        else
            h_lo, n_lo = h_mid, n_mid
        end
    end
    n
end

"""
    setup_neighbors_grid!(s::NeighborLists, pos, inputs::NeighborListInputs)::Nothing

Build the neighbor lists `s` using a hierarchy of uniform grids. Each object
is stored at the finest level whose cells are wider than its padded bounding
box, so it touches at most 8 cells there. Produces the
same set of edges as [`setup_neighbors_naive!`](@ref) and
[`setup_neighbors_sort_sweep!`](@ref).

Throws an `ArgumentError` if `inputs.skin` is negative or NaN. May also throw
an `ArgumentError` if the radius of any object that passes `filter_object` is
negative or NaN.
"""
function setup_neighbors_grid!(s::NeighborLists, pos, inputs::NeighborListInputs)::Nothing
    _argcheck_neighbor_list_input(inputs)
    @argcheck inputs.skin ≥ 0 "skin must be non-negative"
    (; policy,
       points, p_radius, p_params,
       clines, c_radius, c_params,
       lines, l_radius, l_params,
       line_points, lp_radius, lp_params,
       triangles, t_radius, t_params,
       triangle_lines, tl_radius, tl_params,
       triangle_points, tp_radius, tp_params,
       no_collide_pairs, skin, nthreads) = inputs

    _prepare_neighbor_lists!(s, inputs)
    if length(pos) < 2
        # Only one or zero points
        return
    end

    PTYPE = eltype(eltype(pos))

    # Get simulation bounds
    # Using the bounds of `pos` instead of the padded bounding boxes is safe
    # because cell ranges are clamped into the grid, and clamping preserves
    # interval overlap.
    min_p::SVector{3, PTYPE} = SVector{3, PTYPE}(Inf, Inf, Inf)
    max_p::SVector{3, PTYPE} = SVector{3, PTYPE}(-Inf, -Inf, -Inf)
    for p in pos
        all(isfinite, p) || error("positions must be finite")
        min_p = min.(min_p, p)
        max_p = max.(max_p, p)
    end
    box_width::SVector{3, PTYPE} = max_p - min_p
    all(isfinite, box_width) || error("positions span too large a range")

    # Objects of each type
    P = (CollidePoint, points, p_radius, p_params)
    C = (CollideCLine, clines, c_radius, c_params)
    L = (CollideLine, lines, l_radius, l_params)
    LP = (CollideLinePoint, line_points, lp_radius, lp_params)
    T = (CollideTriangle, triangles, t_radius, t_params)
    TL = (CollideTriangleLine, triangle_lines, tl_radius, tl_params)
    TP = (CollideTrianglePoint, triangle_points, tp_radius, tp_params)

    # Smallest radius and number of objects that pass filter_object, in
    # CollideObjectTypes order. The radii and params have the same types for
    # every object type.
    filtered_stats = map((
        (CollidePoint, p_radius, p_params),
        (CollideCLine, c_radius, c_params),
        (CollideLine, l_radius, l_params),
        (CollideLinePoint, lp_radius, lp_params),
        (CollideTriangle, t_radius, t_params),
        (CollideTriangleLine, tl_radius, tl_params),
        (CollideTrianglePoint, tp_radius, tp_params),
    )) do (obj_type, a_r, a_params)
        _filtered_min_radius(policy, obj_type, a_r, a_params)
    end
    min_r = minimum(first, filtered_stats)
    # min propagates NaN, so this also rejects NaN radii
    @argcheck min_r ≥ 0 "radii of objects that pass filter_object must be non-negative"
    n_filtered = sum(last, filtered_stats)
    if n_filtered < 2
        # Only one or zero objects, so no pairs
        return
    end

    # Clamp the minimum cell width to the finite positive Float32 range, so
    # zero or infinite `min_r + skin` is handled.
    min_cell_width = clamp(
        GRID_MIN_CELL_WIDTH_FACTOR * (min_r + skin),
        nextfloat(0.0f0),
        floatmax(Float32),
    )
    n_cells = _grid_size(
        box_width,
        min_cell_width,
        GRID_MAX_CELLS_PER_OBJECT * n_filtered,
    )

    # Quantize the padded bounding boxes of the filtered objects
    # Include extra epsilon padding to account for floating-point rounding errors
    skin_with_eps = skin + 8 * eps(maximum(box_width))
    n_units = GRID_CELL_UNITS .* n_cells
    # For a zero width axis this is Inf, which is safe: every bounding box
    # then starts at 0, so all pairs overlap on that axis.
    scale = n_units ./ box_width
    last_unit = PTYPE.(n_units .- 1)
    # Filtered objects of type `t` are `aabbs[type_ranges[t]]`
    type_counts = map(last, filtered_stats)
    type_ranges = map((stop, n) -> (stop - n + 1):stop, cumsum(type_counts), type_counts)
    aabbs = Vector{GridAABB}(undef, n_filtered)
    _grid_quantize_type!(aabbs, policy, pos, P, type_ranges[Int(CollidePoint)], skin_with_eps, min_p, scale, last_unit)
    _grid_quantize_type!(aabbs, policy, pos, C, type_ranges[Int(CollideCLine)], skin_with_eps, min_p, scale, last_unit)
    _grid_quantize_type!(aabbs, policy, pos, L, type_ranges[Int(CollideLine)], skin_with_eps, min_p, scale, last_unit)
    _grid_quantize_type!(aabbs, policy, pos, LP, type_ranges[Int(CollideLinePoint)], skin_with_eps, min_p, scale, last_unit)
    _grid_quantize_type!(aabbs, policy, pos, T, type_ranges[Int(CollideTriangle)], skin_with_eps, min_p, scale, last_unit)
    _grid_quantize_type!(aabbs, policy, pos, TL, type_ranges[Int(CollideTriangleLine)], skin_with_eps, min_p, scale, last_unit)
    _grid_quantize_type!(aabbs, policy, pos, TP, type_ranges[Int(CollideTrianglePoint)], skin_with_eps, min_p, scale, last_unit)

    # Assign each object to the finest level whose cells are wider than its
    # box, see `_grid_level`. The coarsest possible level has one cell.
    n_levels = UInt8(1 + maximum(n -> Base.top_set_bit(n - 1), n_cells))
    level_of = UInt8[_grid_level(aabb, n_levels) for aabb in aabbs]
    # Number of objects of type `t` at level `k` is `counts[t, begin + k]`
    counts = zeros(Int, N_COLLIDE_OBJECT_TYPES, n_levels)
    for t in 1:N_COLLIDE_OBJECT_TYPES, k in type_ranges[t]
        counts[t, begin + level_of[k]] += 1
    end
    levels = _grid_levels(aabbs, level_of, type_ranges, n_cells, n_levels)

    # Each pair type, like `setup_neighbors_naive!`
    extra_cutoff = 2skin
    cands = GridCandidates{PTYPE}()
    args = (levels, counts, level_of, aabbs, type_ranges, cands, pos, policy, extra_cutoff)
    ncp = no_collide_pairs
    _grid_pairs!(args..., GridPairType(s.PPNL, ncp[Int(CollidePoint_Point)], true, P, P))
    _grid_pairs!(args..., GridPairType(s.PCNL, ncp[Int(CollidePoint_CLine)], false, P, C))
    _grid_pairs!(args..., GridPairType(s.PLNL, ncp[Int(CollidePoint_Line)], false, P, L))
    _grid_pairs!(args..., GridPairType(s.PTNL, ncp[Int(CollidePoint_Triangle)], false, P, T))
    _grid_pairs!(args..., GridPairType(s.PTNL, ncp[Int(CollideLinePoint_Triangle)], false, LP, T))
    _grid_pairs!(args..., GridPairType(s.PTNL, ncp[Int(CollideTrianglePoint_Triangle)], false, TP, T))
    _grid_pairs!(args..., GridPairType(s.CCNL, ncp[Int(CollideCLine_CLine)], true, C, C))
    _grid_pairs!(args..., GridPairType(s.CLNL, ncp[Int(CollideCLine_Line)], false, C, L))
    _grid_pairs!(args..., GridPairType(s.CLNL, ncp[Int(CollideCLine_TriangleLine)], false, C, TL))
    _grid_pairs!(args..., GridPairType(s.LLNL, ncp[Int(CollideLine_Line)], true, L, L))
    _grid_pairs!(args..., GridPairType(s.LLNL, ncp[Int(CollideLine_TriangleLine)], false, L, TL))
    _grid_pairs!(args..., GridPairType(s.LLNL, ncp[Int(CollideTriangleLine_TriangleLine)], true, TL, TL))
    nothing
end
