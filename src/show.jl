# `text/plain` `show` methods for the container types, printing a summary
# instead of the default dump of every field.

# Print `n` followed by the singular or plural form of a noun.
function _print_count(io::IO, n::Integer, one::AbstractString, many::AbstractString)
    print(io, n, " ", isone(n) ? one : many)
end

function Base.show(io::IO, ::MIME"text/plain", pcl::PointCellList)
    lens = Iterators.map(c -> c.len, pcl.cells)
    summary(io, pcl)
    print(io, " with ")
    _print_count(io, sum(lens; init=0), "point", "points")
    print(io, ":\n  size: ", join(pcl.size, "×"), " cells")
    print(io, "\n  grid spacing: ")
    show(IOContext(io, :compact => true), 2*inv(pcl.inv_half_grid_spacing))
    print(io, "\n  most points in one cell: ", maximum(lens; init=0))
end

function Base.show(io::IO, ::MIME"text/plain", scl::LineSegCellList)
    lens = Iterators.map(c -> c.len, scl.cells)
    # Each line segment has exactly one entry with the `is_end` flag.
    nsegs = sum(scl.cells; init=0) do c
        count(cell_line_seg_is_end, view(c.list, 1:c.len))
    end
    summary(io, scl)
    print(io, " with ")
    _print_count(io, nsegs, "line segment", "line segments")
    print(io, ":\n  size: ", join(scl.size, "×"), " cells")
    print(io, "\n  grid spacing: ")
    show(IOContext(io, :compact => true), 2*inv(scl.inv_half_grid_spacing))
    print(io, "\n  cell entries: ", sum(lens; init=0))
    print(io, "\n  most entries in one cell: ", maximum(lens; init=0))
end

function Base.show(io::IO, ::MIME"text/plain", tcl::TriangleCellList)
    lens = Iterators.map(c -> c.len, tcl.cells)
    # Each triangle has exactly one entry in the first cell of its range
    # along every axis.
    ntris = sum(tcl.cells; init=0) do c
        count(e -> e.box_ends & 0x07 == 0x07, view(c.list, 1:c.len))
    end
    summary(io, tcl)
    print(io, " with ")
    _print_count(io, ntris, "triangle", "triangles")
    print(io, ":\n  size: ", join(tcl.size, "×"), " cells")
    print(io, "\n  grid spacing: ")
    show(IOContext(io, :compact => true), 2*inv(tcl.inv_half_grid_spacing))
    print(io, "\n  cell entries: ", sum(lens; init=0))
    print(io, "\n  most entries in one cell: ", maximum(lens; init=0))
end

# Print the forces of all beads as type `F`, one bead per line, for the
# `text/plain` `show` of accumulators that store forces. Like printing a
# `Vector`, a limited `io` shows only the first and last beads so the output
# fits the display, including the `header_lines` already printed.
function _show_bead_forces(io::IO, fe::ForceEnergy, ::Type{F}, header_lines::Int) where F
    n = get_nbeads(fe)
    iszero(n) && return
    print(io, "\n  forces:")
    max_lines = if get(io, :limit, false)
        # Leave 5 lines: the "forces:" line, the "⋮" line, and 3 free lines
        # for the prompt, the same as a long `Vector` leaves free.
        max(displaysize(io)[1] - header_lines - 5, 2)
    else
        n
    end
    beads = if n ≤ max_lines
        (1:n,)
    else
        (1:cld(max_lines, 2), n-fld(max_lines, 2)+1:n)
    end
    width = ndigits(n)
    # The `:typeinfo` stops a `Rational{BigInt}[...]` style prefix.
    compact_io = IOContext(io, :compact => true, :typeinfo => SVector{3, F})
    for (k, range) in enumerate(beads)
        isone(k) || print(io, "\n    ", lpad("⋮", width))
        for i in range
            print(io, "\n    ", lpad(i, width), ": ")
            show(compact_io, get_bead_force(fe, i, F))
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", fe::ForceEnergyT{T}) where T
    summary(io, fe)
    print(io, " with ")
    _print_count(io, get_nbeads(fe), "bead", "beads")
    print(io, ":\n  energy: ")
    show(IOContext(io, :compact => true), fe.energy)
    _show_bead_forces(io, fe, T, 2)
end

function Base.show(io::IO, ::MIME"text/plain", fe::ForceEnergyFixedPoint)
    summary(io, fe)
    print(io, " with ")
    _print_count(io, get_nbeads(fe), "bead", "beads")
    print(io, ":\n  energy: ")
    show(IOContext(io, :compact => true), get_energy(fe))
    _show_bead_forces(io, fe, Float64, 2)
end

function Base.show(io::IO, ::MIME"text/plain", fe::ForceNoEnergyFixedPoint)
    summary(io, fe)
    print(io, " with ")
    _print_count(io, get_nbeads(fe), "bead", "beads")
    print(io, " and no energy:")
    _show_bead_forces(io, fe, Float64, 1)
end

function Base.show(io::IO, ::MIME"text/plain", fe::DebugForceEnergy{T}) where T
    summary(io, fe)
    print(io, " with ")
    _print_count(io, get_nbeads(fe), "bead", "beads")
    print(io, ":\n  ")
    _print_count(io, length(fe.added_energy), "energy contribution", "energy contributions")
    print(io, ", ")
    _print_count(io, length(fe.added_bead_force), "force contribution", "force contributions")
    print(io, "\n  energy: ")
    show(IOContext(io, :compact => true), get_energy(fe, T))
    _show_bead_forces(io, fe, T, 3)
end

function Base.show(io::IO, ::MIME"text/plain", fe::NullForceEnergy)
    summary(io, fe)
    print(io, " with ")
    _print_count(io, get_nbeads(fe), "bead", "beads")
    print(io, ", discarding all forces and energy")
end

function Base.show(io::IO, ::MIME"text/plain", inputs::NeighborListInputs)
    compact_io = IOContext(io, :compact => true)
    summary(io, inputs)
    print(io, ":\n  policy: ")
    show(compact_io, inputs.policy)
    # Name and lengths of the object, radius, and params vectors of each
    # object type, which must match for the inputs to be valid.
    objects = (
        ("points", "p", length(inputs.points), length(inputs.p_radius), length(inputs.p_params)),
        ("clines", "c", length(inputs.clines), length(inputs.c_radius), length(inputs.c_params)),
        ("lines", "l", length(inputs.lines), length(inputs.l_radius), length(inputs.l_params)),
        ("line_points", "lp", length(inputs.line_points), length(inputs.lp_radius), length(inputs.lp_params)),
        ("triangles", "t", length(inputs.triangles), length(inputs.t_radius), length(inputs.t_params)),
        ("triangle_lines", "tl", length(inputs.triangle_lines), length(inputs.tl_radius), length(inputs.tl_params)),
        ("triangle_points", "tp", length(inputs.triangle_points), length(inputs.tp_radius), length(inputs.tp_params)),
    )
    name_width = maximum(o -> length(o[1]), objects) + 1
    count_width = maximum(o -> ndigits(o[3]), objects)
    for (name, prefix, n, n_radius, n_params) in objects
        print(io, "\n  ", rpad(name * ":", name_width), " ", lpad(n, count_width))
        if n_radius != n || n_params != n
            print(io, ", but ", n_radius, " ", prefix, "_radius and ", n_params, " ", prefix, "_params")
        end
    end
    print(io, "\n  no_collide_pairs: ", sum(length, inputs.no_collide_pairs))
    for pair_type in instances(CollidePairs)
        n = length(inputs.no_collide_pairs[Int(pair_type)])
        iszero(n) || print(io, "\n    ", pair_type, ": ", n)
    end
    print(io, "\n  skin: ")
    show(compact_io, inputs.skin)
    print(io, "\n  nthreads: ", inputs.nthreads)
end

function Base.show(io::IO, ::MIME"text/plain", s::NeighborLists)
    lists = (
        ("PPNL", "point–point", length(s.PPNL)),
        ("PCNL", "point–cline", length(s.PCNL)),
        ("PLNL", "point–line", length(s.PLNL)),
        ("PTNL", "point–triangle", length(s.PTNL)),
        ("CCNL", "cline–cline", length(s.CCNL)),
        ("CLNL", "cline–line", length(s.CLNL)),
        ("LLNL", "line–line", length(s.LLNL)),
    )
    summary(io, s)
    print(io, " with ")
    _print_count(io, sum(l -> l[3], lists), "edge", "edges")
    print(io, ":\n  policy: ")
    show(IOContext(io, :compact => true), s.policy)
    desc_width = maximum(l -> textwidth(l[2]), lists) + 1
    count_width = maximum(l -> ndigits(l[3]), lists)
    for (name, desc, n) in lists
        # `rpad` pads by `textwidth`, so the en dashes align.
        print(io, "\n  ", name, " ", rpad(desc * ":", desc_width), " ", lpad(n, count_width))
    end
end
