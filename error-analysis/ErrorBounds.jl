"""
First order rounding error bounds for the distance kernels, valid for all
inputs, found by running the kernels on [`EB`](@ref) numbers.

An `EB` holds upper bounds as polynomials in size symbols: `mag` on the
absolute value of the exact result, and `err` on the rounding error, in
units of `u`. Each operation applies a rule that holds for all inputs, so
the bounds are functions of the sizes, not of a particular input. Products
of two errors, `O(u²)`, are dropped.

A sum or difference of two exact values, such as `y[1] - x[1]`, becomes a
new size symbol, for example `|Pᵢ|`, rounded once. So the bounds are in the
sizes of the differences, not of the coordinates.

In `:backward` mode, the rounded differences of input vertices are instead
treated as exact differences of slightly moved vertices. The kernel is then
analysed on the moved vertices, and [`vertex_moves`](@ref) gives how far
each vertex moved, which changes the exact distance by at most the largest
move on each simplex. This is how the kernels use the differences: the same
rounded `E0`, `E1` everywhere, so the result is the accurate distance of a
nearby triangle, even when it is a sliver.

The kernels are branch free, so tracing evaluates both arms of every
select. A select joins the bounds of its arms, and adds a jump term when
the computed condition can disagree with the exact one:

- A running minimum, `cand[1] < best[1]`, has no jump, since `min` is
  1-Lipschitz. `r`, `s` and `t` are compared with the exact values of the
  candidate that was chosen.
- An `inside` test has a jump of `Σ err(v) * L(v)` over its distinct
  variables `v`, where `L(v)` comes from [`JUMPS`](@ref): a hand derived
  bound on how far the boundary is, in distance, per unit of `v`.
"""
module ErrorBounds

using SimplexCellLists: SimplexCellLists, SIMD
using StaticArrays: SVector

export EB, SPoly, Sym, input_simplex, normbound, vecnormbound, maxpoly, boundvalue,
    reset!, vertex_moves, substitute_abscross, NAMES, CROSS_DEFS, ABSCROSS_DEFS, JUMPS, ASSUMPTIONS, NOTES

const Q = Rational{BigInt}

# A size symbol: `|base_i|` for component `i` in 1:3, or `‖base‖` when `i == 0`.
struct Sym
    base::Symbol
    i::Int
end
Base.isless(a::Sym, b::Sym) = isless((a.base, a.i), (b.base, b.i))

# Sorted `Sym => exponent` pairs, with nonzero exponents.
const Mono = Vector{Pair{Sym, Int}}

function monomul(a::Mono, b::Mono)
    d = Dict{Sym, Int}()
    for (s, e) in Iterators.flatten((a, b))
        d[s] = get(d, s, 0) + e
    end
    sort!([s => e for (s, e) in d if e != 0]; by=first)
end

"""
    SPoly

A polynomial with nonnegative coefficients in size symbols, used as an
upper bound. Exponents may be negative for norms, as in `‖R‖⁻²`.
"""
struct SPoly
    t::Dict{Mono, Q}
end
SPoly() = SPoly(Dict{Mono, Q}())
SPoly(c::Real) = iszero(c) ? SPoly() : SPoly(Dict(Mono() => Q(c)))
SPoly(s::Sym, e::Int=1) = SPoly(Dict([s => e] => one(Q)))

Base.iszero(p::SPoly) = isempty(p.t)
function Base.:+(a::SPoly, b::SPoly)
    t = copy(a.t)
    for (m, c) in b.t
        t[m] = get(t, m, zero(Q)) + c
    end
    SPoly(t)
end
function Base.:*(a::SPoly, b::SPoly)
    t = Dict{Mono, Q}()
    for (ma, ca) in a.t, (mb, cb) in b.t
        m = monomul(ma, mb)
        t[m] = get(t, m, zero(Q)) + ca*cb
    end
    SPoly(t)
end
Base.:*(c::Real, p::SPoly) = iszero(c) ? SPoly() : SPoly(Dict(m => Q(c)*v for (m, v) in p.t))

"""
    maxpoly(a::SPoly, b::SPoly)

An upper bound on `max(a, b)`: the largest coefficient of each monomial.
"""
function maxpoly(a::SPoly, b::SPoly)
    t = copy(a.t)
    for (m, c) in b.t
        t[m] = max(get(t, m, zero(Q)), c)
    end
    SPoly(t)
end

# The square root of a monomial with even exponents and coefficient one.
function _sqrtmag(p::SPoly)
    length(p.t) == 1 || return nothing
    (m, c), = p.t
    isone(c) && all(iseven ∘ last, m) || return nothing
    SPoly(Dict(Mono([s => e ÷ 2 for (s, e) in m]) => one(Q)))
end

# The component indices that appear in `p`.
indices(p::SPoly) = Set(s.i for m in keys(p.t) for (s, _) in m if s.i > 0)

const SUPERSCRIPTS = Dict(zip("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
function Base.show(io::IO, p::SPoly)
    iszero(p) && return print(io, "0")
    terms = sort!(collect(p.t); by=x -> string(x[1]))
    for (k, (m, c)) in enumerate(terms)
        k > 1 && print(io, " + ")
        (isone(c) && !isempty(m)) || print(io, isinteger(c) ? numerator(c) : c)
        for (s, e) in m
            print(io, s.i == 0 ? "‖$(s.base)‖" : "|$(s.base)$(('₀':'₉')[s.i + 1])|")
            e == 1 || print(io, map(ch -> SUPERSCRIPTS[ch], string(e)))
        end
    end
end

"""
    boundvalue(p::SPoly, values)

Evaluate `p` with `values[s]` for each symbol `s`.
"""
boundvalue(p::SPoly, values) = sum((c*prod((values[s]^e for (s, e) in m); init=big(1)) for (m, c) in p.t); init=big(0))

"""
    normbound(p::SPoly)

Bound a scalar that sums over components by norms. For each group of terms
that differ only in their component index `i`, with the largest coefficient
`c`, `Σᵢ c|Xᵢ||Yᵢ|… ≤ c‖X‖‖Y‖…` (Cauchy-Schwarz, then `|Zᵢ| ≤ ‖Z‖`) and
`Σᵢ c|Xᵢ| ≤ c√n‖X‖` for `n` indices. Other terms use `|Xⱼ| ≤ ‖X‖` alone.
"""
function normbound(p::SPoly)
    out = SPoly()
    groups = Dict{Tuple{Vector{Pair{Symbol, Int}}, Mono}, Tuple{Q, Set{Int}}}()
    for (m, c) in p.t
        comps = [s => e for (s, e) in m if s.i > 0]
        glob = Mono([s => e for (s, e) in m if s.i == 0])
        @assert all(>(0), last.(comps)) "negative exponent on a component"
        idx = Set(s.i for (s, _) in comps)
        if isempty(comps)
            out += SPoly(Dict(m => c))
        elseif length(idx) == 1
            key = (sort!([s.base => e for (s, e) in comps]; by=first), glob)
            cmax, seen = get(groups, key, (zero(Q), Set{Int}()))
            groups[key] = (max(cmax, c), union(seen, idx))
        else
            out += SPoly(Dict(monomul(glob, Mono(sort!([Sym(s.base, 0) => e for (s, e) in comps]; by=first))) => c))
        end
    end
    for ((bases, glob), (cmax, seen)) in groups
        norms = Mono(sort!([Sym(b, 0) => e for (b, e) in bases]; by=first))
        degree = sum(last, bases)
        # √n for a sum of single components, rounded up to a rational.
        factor = degree ≥ 2 ? one(Q) : (length(seen) == 1 ? one(Q) : length(seen) == 2 ? Q(3//2) : Q(7//4))
        out += SPoly(Dict(monomul(glob, norms) => cmax*factor))
    end
    out
end

"""
    vecnormbound(ps)

Bound the Euclidean norm of a vector whose component `i` is at most `ps[i]`,
by Minkowski: group the terms by the part with index `i`, take the largest
coefficient, and bound the norm of `(|Xᵢ||Yᵢ|…)ᵢ` by `‖X‖‖Y‖…`. A term
without index `i` is a constant vector, with norm `√3` times its value.
"""
function vecnormbound(ps)
    groups = Dict{Tuple{Vector{Pair{Symbol, Int}}, Mono}, Q}()
    for (i, p) in enumerate(ps)
        rests = Dict{Vector{Pair{Symbol, Int}}, SPoly}()
        for (m, c) in p.t
            own = sort!([s.base => e for (s, e) in m if s.i == i]; by=first)
            rest = SPoly(Dict(Mono([s => e for (s, e) in m if s.i != i]) => c))
            rests[own] = get(rests, own, SPoly()) + rest
        end
        for (own, rest) in rests
            # The rest no longer depends on `i`, so reduce it to norms.
            for (m, c) in normbound(rest).t
                key = (own, m)
                groups[key] = max(get(groups, key, zero(Q)), c)
            end
        end
    end
    out = SPoly()
    for ((own, m), c) in groups
        norms = Mono(sort!([Sym(b, 0) => e for (b, e) in own]; by=first))
        out += SPoly(Dict(monomul(norms, m) => isempty(own) ? Q(7//4)*c : c))
    end
    out
end

"""
    EB

A traced number: `mag` bounds the absolute value of the exact result, and
`err` bounds the rounding error in units of `u`. If `rel` is set, the error
is also at most `rel*u` times the exact value. `exact` means `mag` is the
absolute value of the exact result, not just a bound, as division needs.

`root`, if set, is `(abs, rel)` for a nonnegative result `ρ²`, such as a
squared distance: `|√computed - ρ| ≤ u*(abs + rel*ρ)`.

`sqscale`, if set, records that this is the product `r*w` with `w` a
nonnegative exact square, so that a later `*r` is recognized as a squared
distance, as in `(hN*invΔ)*hN`.
"""
struct EB <: Real
    mag::SPoly
    err::SPoly
    rel::Union{Nothing, Q}
    nonneg::Bool
    exact::Bool
    input::Union{Nothing, Sym}
    con::Union{Nothing, Q}
    root::Union{Nothing, Tuple{SPoly, Q}}
    sqscale::Union{Nothing, Tuple{EB, EB}}
end

EB(mag, err, rel, nonneg, exact, input, con, root) =
    EB(mag, err, rel, nonneg, exact, input, con, root, nothing)
EB(c::Real) = EB(SPoly(abs(c)), SPoly(), zero(Q), c ≥ 0, true, nothing, Q(c), nothing)
withrel(mag, rel, nonneg, exact) = EB(mag, rel*mag, rel, nonneg, exact, nothing, nothing, nothing)
withabs(mag, err, nonneg, exact) = EB(mag, err, nothing, nonneg, exact, nothing, nothing, nothing)
_withroot(x::EB, root) = EB(x.mag, x.err, x.rel, x.nonneg, x.exact, x.input, x.con, root)

Base.zero(::Type{EB}) = EB(0)
Base.one(::Type{EB}) = EB(1)
Base.zero(::EB) = EB(0)
Base.one(::EB) = EB(1)
Base.convert(::Type{EB}, c::Real) = EB(c)
Base.convert(::Type{EB}, x::EB) = x
Base.promote_rule(::Type{EB}, ::Type{<:Real}) = EB
function Base.show(io::IO, x::EB)
    print(io, "EB(mag = ", x.mag, ", err = (", x.err, ")·u")
    x.rel === nothing || print(io, ", rel = ", x.rel, "·u")
    x.root === nothing || print(io, ", root err = (", x.root[1], ")·u + ", x.root[2], "·u·ρ")
    print(io, ")")
end

"""
    NAMES

Names for sums and differences of exact values, keyed by
`(op, base_a, base_b)`, for example `(:-, :a, :x) => :P` for
`P = a - x`. Missing names are made up.
"""
const NAMES = Dict{Tuple{Symbol, Symbol, Symbol}, Symbol}()
const DEFS = Dict{Symbol, Tuple{Symbol, Symbol, Symbol}}()
"Names for the results of `cross3_accurate`, used in order, and their definitions."
const CROSS_NAMES = Symbol[]
const CROSS_DEFS = Dict{Symbol, Tuple{Symbol, Symbol}}()
"""
For a plain `cross3(a, b)`, `X => (a, b)`: `|Xᵢ| = |aⱼ||bₖ| + |aₖ||bⱼ|` bounds
each component, and `‖X‖ ≤ √2‖a‖‖b‖`.
"""
const ABSCROSS_DEFS = Dict{Symbol, Tuple{Symbol, Symbol}}()
"For each `inside` select, in order, the distance per unit of each distinct condition variable."
const JUMPS = Vector{Vector{SPoly}}()
const MOVES = Dict{Symbol, Tuple{Symbol, Symbol}}()
const MODE = Ref(:backward)
const ASSUMPTIONS = String[]
const NOTES = String[]

"""
    reset!(; names, cross=Symbol[], jumps=Vector{SPoly}[], mode=:backward)

Clear the trace state and set the names and jump constants for one kernel.
"""
function reset!(; names=Dict{Tuple{Symbol, Symbol, Symbol}, Symbol}(), cross=Symbol[], jumps=Vector{SPoly}[], mode=:backward)
    empty!(NAMES); merge!(NAMES, names)
    empty!(DEFS); empty!(CROSS_DEFS); empty!(ABSCROSS_DEFS); empty!(MOVES)
    empty!(CROSS_NAMES); append!(CROSS_NAMES, cross)
    empty!(JUMPS); append!(JUMPS, jumps)
    empty!(ASSUMPTIONS); empty!(NOTES)
    MODE[] = mode
    nothing
end

"""
    input_simplex(names...)

Return a simplex of exact input vertices named `names`, with each
coordinate's size the symbol `|nameᵢ|`.
"""
function input_simplex(names::Symbol...)
    SVector(map(names) do n
        SVector(ntuple(i -> EB(SPoly(Sym(n, i)), SPoly(), zero(Q), false, true, Sym(n, i), nothing, nothing), 3))
    end)
end

"""
    vertex_moves(v::Symbol)

In `:backward` mode, a bound on how far input vertex `v` moved, as one
`SPoly` per component in units of `u`, following the rounded differences
back to a vertex that did not move.
"""
function vertex_moves(v::Symbol)
    haskey(MOVES, v) || return [SPoly() for _ in 1:3]
    anchor, name = MOVES[v]
    [m + SPoly(Sym(name, i)) for (i, m) in enumerate(vertex_moves(anchor))]
end

# The symbol of an exact value that is a single size symbol, or `nothing`.
function _single(x::EB)
    iszero(x.err) && x.exact && length(x.mag.t) == 1 || return nothing
    (m, c), = x.mag.t
    isone(c) && length(m) == 1 && m[1][2] == 1 && m[1][1].i > 0 ? m[1][1] : nothing
end

# One rounding of the exact result `x`.
function _round(x::EB)
    x.con === nothing || return x
    x.rel === nothing ? withabs(x.mag, x.err + x.mag, x.nonneg, x.exact) : withrel(x.mag, x.rel + 1, x.nonneg, x.exact)
end

function _neg(a::EB)
    EB(a.mag, a.err, a.rel, iszero(a.mag), a.exact, nothing, a.con === nothing ? nothing : -a.con, nothing)
end

function _add_exact(a::EB, b::EB)
    a.con == 0 && return b
    b.con == 0 && return a
    mag = a.mag + b.mag
    nonneg = a.nonneg && b.nonneg
    exact = nonneg && a.exact && b.exact
    # Each nonnegative term has relative error at most its `rel`, so their
    # sum has relative error at most the larger one.
    if nonneg && a.rel !== nothing && b.rel !== nothing
        withrel(mag, max(a.rel, b.rel), true, exact)
    else
        withabs(mag, a.err + b.err, nonneg, exact)
    end
end

# A sum or difference of two exact values with single symbols of the same
# component is a new symbol, rounded once, so its error is relative to it.
function _addsub(op::Symbol, a::EB, b::EB)
    sa, sb = _single(a), _single(b)
    if sa !== nothing && sb !== nothing && sa.i == sb.i
        name = get!(NAMES, (op, sa.base, sb.base), Symbol(sa.base, op == :- ? :_m_ : :_p_, sb.base))
        DEFS[name] = (op, sa.base, sb.base)
        if op == :- && a.input !== nothing && b.input !== nothing && MODE[] == :backward
            # The rounded difference is the exact difference of `a` moved.
            get!(MOVES, sa.base, (sb.base, name)) == (sb.base, name) || error("vertex $(sa.base) moved twice")
            return EB(SPoly(Sym(name, sa.i)), SPoly(), zero(Q), false, true, nothing, nothing, nothing)
        end
        return withrel(SPoly(Sym(name, sa.i)), one(Q), false, true)
    end
    _round(_add_exact(a, op == :- ? _neg(b) : b))
end

# Reduce a scalar that sums over components to norms before it scales a
# single component.
function _reduce(a::EB, b::EB)
    length(indices(a.mag) ∪ indices(a.err)) > 1 && length(indices(b.mag) ∪ indices(b.err)) == 1 || return a
    EB(normbound(a.mag), normbound(a.err), nothing, a.nonneg, false, nothing, nothing, nothing)
end

function _mul_exact(a::EB, b::EB)
    (a.con == 0 || b.con == 0) && return EB(0)
    a.con == 1 && return b
    b.con == 1 && return a
    nonneg = a === b || (a.nonneg && b.nonneg)
    a, b = _reduce(a, b), _reduce(b, a)
    mag = a.mag * b.mag
    exact = a.exact && b.exact
    if a.rel !== nothing && b.rel !== nothing
        withrel(mag, a.rel + b.rel, nonneg, exact)
    else
        withabs(mag, a.err*b.mag + a.mag*b.err, nonneg, exact)
    end
end

# `1/|B|` for the exact value of a divisor: a single monomial, or `Σᵢ Xᵢ²`
# for the squared norm `‖X‖²`.
function _inv_mag(p::SPoly)
    if length(p.t) == 1
        (m, c), = p.t
        return SPoly(Dict(Mono([s => -e for (s, e) in m]) => inv(c)))
    end
    bases = unique(s.base for m in keys(p.t) for (s, _) in m)
    sq(b) = Dict([Sym(b, i) => 2] => one(Q) for i in 1:3)
    length(bases) == 1 && p.t == sq(only(bases)) && return SPoly(Sym(only(bases), 0), -2)
    error("cannot divide by an exact value of $p")
end

function _div_exact(a::EB, b::EB)
    b.nonneg && b.exact && b.rel !== nothing || error("divisor must be nonnegative with a known exact size and relative error")
    invb = _inv_mag(b.mag)
    mag = a.mag * invb
    if a.rel !== nothing
        withrel(mag, a.rel + b.rel, a.nonneg, a.exact)
    else
        # First order: a/(B(1 + β)) = a/B - (a/B)β.
        withabs(mag, a.err*invb + b.rel*a.mag*invb, a.nonneg, a.exact)
    end
end

function Base.:*(a::EB, b::EB)
    x = _round(_mul_exact(a, b))
    if a === b
        # √(x̂²(1 + δ)) = |x̂|(1 + δ/2).
        return _withroot(x, a.rel === nothing ? (a.err, Q(1//2)) : (SPoly(), a.rel + Q(1//2)))
    end
    # A squared distance times a nonnegative exact square, such as `1/‖N‖²`.
    for (s, w) in ((a, b), (b, a))
        s.root === nothing && continue
        w.nonneg && w.exact && w.rel !== nothing || continue
        sq = _sqrtmag(w.mag)
        sq === nothing && continue
        return _withroot(x, (s.root[1]*sq, s.root[2] + w.rel/2 + Q(1//2)))
    end
    # `(r*w)*r` with the recorded exact square `w`: the same two roundings as
    # `(r*r)*w`, so the same distance bound as the two rules above combined.
    for (s, o) in ((a, b), (b, a))
        s.sqscale === nothing && continue
        r, w = s.sqscale
        o === r || continue
        sq = _sqrtmag(w.mag)
        rt = r.rel === nothing ? (r.err, Q(1//2)) : (SPoly(), r.rel + Q(1//2))
        return _withroot(x, (rt[1]*sq, rt[2] + w.rel/2 + Q(1//2)))
    end
    # Remember `r*w` with `w` a nonnegative exact square for the rule above.
    for (r, w) in ((a, b), (b, a))
        w.nonneg && w.exact && w.rel !== nothing || continue
        _sqrtmag(w.mag) === nothing && continue
        return EB(x.mag, x.err, x.rel, x.nonneg, x.exact, x.input, x.con, x.root, (r, w))
    end
    x
end
Base.:-(a::EB, b::EB) = _addsub(:-, a, b)
Base.:+(a::EB, b::EB) = _addsub(:+, a, b)
Base.:-(a::EB) = _neg(a)
Base.:/(a::EB, b::EB) = _round(_div_exact(a, b))
Base.inv(b::EB) = EB(1)/b
Base.fma(a::EB, b::EB, c::EB) = _round(_add_exact(_mul_exact(a, b), c))

function SimplexCellLists.dot3(a::SVector{3, EB}, b::SVector{3, EB})
    x = fma(a[3], b[3], fma(a[2], b[2], a[1]*b[1]))
    a === b || return x
    # The sum of squares of the computed `a` has relative error γ₃, and
    # `‖â‖` is within `‖â - a‖` of `‖a‖`.
    _withroot(x, x.rel === nothing ? (vecnormbound([ai.err for ai in a]), Q(3//2)) : (SPoly(), x.rel/2))
end

function SimplexCellLists.cross3_accurate(a::SVector{3, EB}, b::SVector{3, EB})
    sa, sb = _single.(a), _single.(b)
    if all(!isnothing, sa) && all(!isnothing, sb) && allequal(s.base for s in sa) && allequal(s.base for s in sb)
        # Each component is Kahan's 2×2 determinant of exact values, within
        # 2u of the exact component (Jeannerod, Louvet, Muller 2013).
        name = isempty(CROSS_NAMES) ? Symbol(sa[1].base, :_x_, sb[1].base) : popfirst!(CROSS_NAMES)
        CROSS_DEFS[name] = (sa[1].base, sb[1].base)
        return SVector(ntuple(i -> withrel(SPoly(Sym(name, i)), Q(2), false, true), 3))
    end
    det(p, q, u, v) = (x = _add_exact(_mul_exact(p, q), _neg(_mul_exact(u, v))); withabs(x.mag, x.err + 2*x.mag, false, false))
    SVector(det(a[2], b[3], a[3], b[2]), det(a[3], b[1], a[1], b[3]), det(a[1], b[2], a[2], b[1]))
end

# The symbol and relative error of a value that is a single size symbol
# times `1 + rel*u`, or `nothing`.
function _symrel(x::EB)
    x.exact && x.rel !== nothing && length(x.mag.t) == 1 || return nothing
    (m, c), = x.mag.t
    isone(c) && length(m) == 1 && m[1][2] == 1 && m[1][1].i > 0 ? (m[1][1], x.rel) : nothing
end

function SimplexCellLists.cross3(a::SVector{3, EB}, b::SVector{3, EB})
    sa, sb = _symrel.(a), _symrel.(b)
    if all(!isnothing, sa) && all(!isnothing, sb) && allequal(first(s).base for s in sa) && allequal(first(s).base for s in sb)
        # Component i is fma(aⱼ, bₖ, -fl(aₖbⱼ)). With |Xᵢ| = |aⱼ||bₖ| + |aₖ||bⱼ|,
        # the input errors give (rel_a + rel_b)|Xᵢ|, the product rounding at
        # most |Xᵢ|, and the fma rounding at most |cᵢ| ≤ |Xᵢ|.
        name = Symbol(:X_, first(sa[1]).base, :_, first(sb[1]).base)
        ABSCROSS_DEFS[name] = (first(sa[1]).base, first(sb[1]).base)
        rel = maximum(last, sa) + maximum(last, sb) + 2
        return SVector(ntuple(i -> withabs(SPoly(Sym(name, i)), rel*SPoly(Sym(name, i)), false, false), 3))
    end
    SVector(fma(a[2], b[3], -(a[3]*b[2])), fma(a[3], b[1], -(a[1]*b[3])), fma(a[1], b[2], -(a[2]*b[1])))
end

"""
    substitute_abscross(p::SPoly)

Replace each `‖X‖` of a plain `cross3(a, b)` by its bound `√2‖a‖‖b‖`,
rounded up to `3/2‖a‖‖b‖`. Call after [`normbound`](@ref).
"""
function substitute_abscross(p::SPoly)
    out = SPoly()
    for (m, c) in p.t
        term = SPoly(c)
        for (s, e) in m
            if s.i == 0 && haskey(ABSCROSS_DEFS, s.base)
                e > 0 || error("‖$(s.base)‖ has a negative exponent")
                a, b = ABSCROSS_DEFS[s.base]
                term = term * SPoly(Q(3//2)^e) * SPoly(Sym(a, 0), e) * SPoly(Sym(b, 0), e)
            else
                term = term * SPoly(s, e)
            end
        end
        out += term
    end
    out
end

abstract type AbstractCond end
struct Cond <: AbstractCond
    op::Symbol
    lhs::EB
    rhs::EB
end
struct CondAnd <: AbstractCond
    parts::Vector{Cond}
end
_parts(c::Cond) = [c]
_parts(c::CondAnd) = c.parts
Base.:<(a::EB, b::EB) = Cond(:lt, a, b)
Base.:≤(a::EB, b::EB) = Cond(:le, a, b)
Base.:&(a::AbstractCond, b::AbstractCond) = CondAnd([_parts(a); _parts(b)])

struct NaNCond
    x::EB
end
Base.isnan(x::EB) = NaNCond(x)
function SIMD.vifelse(c::NaNCond, a, b)
    push!(ASSUMPTIONS, "a clamped value is never NaN (no zero length divisor)")
    b
end
SIMD.vifelse(c::AbstractCond, a, b) = error("select on $c outside `_select_d2_r_s_t` is not traced")

# Clamping is 1-Lipschitz, so it does not increase the error.
function Base.FastMath.max_fast(x::EB, z::EB)
    z.con == 0 || error("max_fast is only traced against zero")
    x.nonneg ? x : withabs(x.mag, x.err, true, false)
end
function Base.FastMath.min_fast(x::EB, o::EB)
    o.con == 1 && x.nonneg && return withabs(SPoly(1), x.err, true, false)
    # |min(x̂, ô) - min(x, o)| ≤ max(|x̂ - x|, |ô - o|).
    withabs(maxpoly(x.mag, o.mag), maxpoly(x.err, o.err), x.nonneg && o.nonneg, false)
end

# Join two arms: each bound is the larger of the two.
function _join(a::EB, b::EB)
    root = a.root !== nothing && b.root !== nothing ? (maxpoly(a.root[1], b.root[1]), max(a.root[2], b.root[2])) : nothing
    rel = a.rel !== nothing && b.rel !== nothing ? max(a.rel, b.rel) : nothing
    EB(maxpoly(a.mag, b.mag), maxpoly(a.err, b.err), rel, a.nonneg && b.nonneg, false, nothing, nothing, root)
end

function SimplexCellLists._select_d2_r_s_t(c::AbstractCond, a::Tuple, b::Tuple)
    d2 = _join(a[1], b[1])
    if c isa Cond && c.op == :lt && c.lhs === a[1] && c.rhs === b[1]
        # A running minimum: |min(â, b̂) - min(a, b)| ≤ max(|â - a|, |b̂ - b|),
        # and either arm that can be chosen has a distance within the error
        # of the minimum, so `rel*ρ` stays in terms of the minimum.
        push!(NOTES, "running minimum: no jump in the distance; r, s, t are those of the chosen candidate")
    else
        vars = EB[]
        for p in _parts(c)
            v = p.lhs.con === nothing ? p.lhs : p.rhs
            (p.lhs.con === nothing) ⊻ (p.rhs.con === nothing) || error("condition must compare a variable with a constant")
            any(w -> w === v, vars) || push!(vars, v)
        end
        isempty(JUMPS) && error("no jump constants for this select")
        Ls = popfirst!(JUMPS)
        length(Ls) == length(vars) || error("$(length(vars)) condition variables, $(length(Ls)) jump constants")
        jump = sum(v.err * L for (v, L) in zip(vars, Ls))
        d2.root === nothing && error("select arms have no distance bound")
        d2 = _withroot(d2, (d2.root[1] + jump, d2.root[2]))
        push!(NOTES, "inside select: jump $(jump) in the distance; r, s, t are those of the chosen arm")
    end
    r = SVector(ntuple(i -> _join(a[2][i], b[2][i]), 3))
    (d2, r, _join(a[3], b[3]), _join(a[4], b[4]))
end

end
