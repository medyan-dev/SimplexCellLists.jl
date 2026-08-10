"""
    abstract type ForceEnergy

Accumulator interface for summing per-bead forces and total energy during a
force calculation. Concrete strategies trade off precision, thread
combinability, and whether energy is tracked at all — see
[`ForceEnergyT`](@ref), [`ForceEnergyFixedPoint`](@ref),
[`ForceNoEnergyFixedPoint`](@ref), [`DebugForceEnergy`](@ref), and
[`NullForceEnergy`](@ref).

# Implementing a custom accumulator

A concrete `FE <: ForceEnergy` provides a constructor `FE(nbeads)` and the
write/read methods below.

Required, called in the hot force loop:
  - [`add_bead_force!`](@ref)`(fe, i, force)` — add an `SVector{3}` to bead `i`.
    Annotate `Base.@propagate_inbounds`.
  - [`zero_bead_force!`](@ref)`(fe, i)` — set force of bead `i` to zero. Useful for pinning.
    Annotate `Base.@propagate_inbounds`.
  - [`add_energy!`](@ref)`(fe, energy)` — add a scalar to the total energy.
  - [`zero_energy!`](@ref)`(fe)` — set the total energy to zero.

Required, for setup and reading results (accumulators that discard data, like
[`NullForceEnergy`](@ref) and [`ForceNoEnergyFixedPoint`](@ref), leave
the readers for the discarded data unimplemented so reading errors instead of
returning made-up values):
  - [`get_nbeads`](@ref)`(fe)` — number of beads.
  - [`zero_force_energy!`](@ref)`(fe)` — reset all forces and energy to zero.
  - [`get_force!`](@ref)`(fe, out)` — write per-bead `SVector{3}` forces into `out`.
  - [`get_bead_force`](@ref)`(fe, i, F=Float64)` — one bead's force.
  - [`get_energy`](@ref)`(fe, F=Float64)` — the accumulated energy.

Provided generically from the above, and may be specialized for speed or
bit-exactness: [`get_force`](@ref),
[`combine_force_energy!`](@ref) (merge one accumulator into another of matching
size, e.g. to reduce per-thread partials).
"""
abstract type ForceEnergy end

"""
    add_bead_force!(force_energy::ForceEnergy, i::Integer, force)::Nothing

Add `force`, a 3-element vector, to the accumulated force of bead `i`.
"""
function add_bead_force! end

"""
    zero_bead_force!(force_energy::ForceEnergy, i::Integer)::Nothing

Set the accumulated force of bead `i` to zero.
"""
function zero_bead_force! end

"""
    add_energy!(force_energy::ForceEnergy, energy)::Nothing

Add the scalar `energy` to the accumulated total energy.
"""
function add_energy! end

"""
    zero_energy!(force_energy::ForceEnergy)::Nothing

Set the accumulated total energy to zero, leaving forces unchanged.
"""
function zero_energy! end

"""
    get_nbeads(force_energy::ForceEnergy)::Int

Return the number of beads.
"""
function get_nbeads end

"""
    zero_force_energy!(force_energy::ForceEnergy)::Nothing

Set all accumulated forces and the total energy to zero.
"""
function zero_force_energy! end

"""
    get_force!(force_energy::ForceEnergy, force_out)::Nothing

Write the accumulated per-bead forces into `force_out`, a one-based vector of
3-element vectors of length [`get_nbeads`](@ref)`(force_energy)`.
"""
function get_force! end

"""
    get_bead_force(force_energy::ForceEnergy, i::Integer, F::Type=Float64)::SVector{3, F}

Return the accumulated force of bead `i`, converted to element type `F`.
"""
function get_bead_force end

"""
    get_energy(force_energy::ForceEnergy, F::Type=Float64)::F

Return the accumulated total energy, converted to type `F`.
"""
function get_energy end

"""
    get_force(force_energy::ForceEnergy, F::Type=Float64)::Vector{SVector{3, F}}

Return a new vector of the accumulated per-bead forces.
"""
function get_force(force_energy::ForceEnergy, ::Type{F}=Float64) where F
    force_out = zeros(SVector{3, F}, get_nbeads(force_energy))
    get_force!(force_energy, force_out)
    force_out
end

"""
    combine_force_energy!(fe_out::ForceEnergy, fe_in::ForceEnergy) -> fe_out

Merge `fe_in` into `fe_out`, which must have the same number of beads, by
adding its accumulated forces and energy. Useful to reduce per-thread partial
accumulators. The generic fallback round-trips through `Float64`; matching
fixed-point types merge bit-exactly.
"""
function combine_force_energy!(fe_out::ForceEnergy, fe_in::ForceEnergy)
    @argcheck get_nbeads(fe_in) == get_nbeads(fe_out)
    add_energy!(fe_out, get_energy(fe_in))
    for i in 1:get_nbeads(fe_in)
        @inbounds add_bead_force!(fe_out, i, get_bead_force(fe_in, i))
    end
    fe_out
end

"""
    ForceEnergyT{T}(nbeads::Integer)

A [`ForceEnergy`](@ref) accumulator storing per-bead forces and the total
energy as numbers of type `T`, summed in place. If `T` is an inexact type like
`Float64`, results can depend on the order contributions are added in, so
parallel reductions are not reproducible bit-for-bit.
"""
mutable struct ForceEnergyT{T} <: ForceEnergy
    const forces::Memory{SVector{3,T}}
    energy::T
end
"""
    ForceEnergyFloat64

Alias for [`ForceEnergyT`](@ref)`{Float64}`.
"""
const ForceEnergyFloat64 = ForceEnergyT{Float64}
"""
    ForceEnergyFloat32

Alias for [`ForceEnergyT`](@ref)`{Float32}`.
"""
const ForceEnergyFloat32 = ForceEnergyT{Float32}
function ForceEnergyT{T}(nbeads::Integer) where T
    forces = Memory{SVector{3,T}}(undef, nbeads)
    fill!(forces, zero(SVector{3,T}))
    ForceEnergyT{T}(forces, zero(T))
end
function zero_force_energy!(force_energy::ForceEnergyT{T}) where T
    fill!(force_energy.forces, zero(SVector{3, T}))
    force_energy.energy = zero(T)
    nothing
end
function get_nbeads(force_energy::ForceEnergyT)
    length(force_energy.forces)
end
function add_energy!(force_energy::ForceEnergyT, energy)
    force_energy.energy += energy
    nothing
end
function zero_energy!(force_energy::ForceEnergyT{T}) where {T}
    force_energy.energy = zero(T)
    nothing
end
Base.@propagate_inbounds function add_bead_force!(force_energy::ForceEnergyT, i::Integer, force)
    force_energy.forces[i] += force
    nothing
end
Base.@propagate_inbounds function zero_bead_force!(force_energy::ForceEnergyT{T}, i::Integer) where {T}
    force_energy.forces[i] = zero(SVector{3,T})
    nothing
end
function get_force!(force_energy::ForceEnergyT, force_out)
    Base.require_one_based_indexing(force_out)
    @argcheck length(force_out) == length(force_energy.forces)
    @inbounds for i in eachindex(force_out)
        force_out[i] = force_energy.forces[i]
    end
    nothing
end
Base.@propagate_inbounds function get_bead_force(force_energy::ForceEnergyT{T}, i::Integer, ::Type{F}=Float64)::SVector{3, F} where {T, F}
    SVector{3, F}(force_energy.forces[i])
end
function get_energy(force_energy::ForceEnergyT, ::Type{F}=Float64) where F
    F(force_energy.energy)
end


@inline convert_to_fixed_point(x, fbits) = unsafe_trunc(Int64, x*(oftype(x, 2)^fbits))
@inline convert_to_floating_point(T, x, fbits) = T(x)*(convert(T, 2)^-fbits)

"""
    ForceEnergyFixedPoint{F_FBITS, E_FBITS}(nbeads::Integer)

A [`ForceEnergy`](@ref) accumulator storing per-bead forces and the total
energy as `Int64` fixed-point numbers with `F_FBITS` and `E_FBITS` fractional
bits respectively. Accumulation is exact once inputs are quantized, so results
are independent of the order contributions are added in, and accumulators with
matching parameters merge bit-exactly with [`combine_force_energy!`](@ref)
(e.g. to reduce per-thread partials reproducibly).

Added values are quantized by truncation toward zero, losing up to
`2^-F_FBITS` (or `2^-E_FBITS`) of magnitude per contribution. Negating an
input exactly negates its quantized value, so equal-and-opposite contributions
cancel exactly. Inputs must be finite with magnitude less than
`2^(63 - FBITS)`; non-finite or out-of-range values silently corrupt the
accumulator (`unsafe_trunc`), as does overflow of the running `Int64` sums.
"""
mutable struct ForceEnergyFixedPoint{F_FBITS, E_FBITS} <: ForceEnergy
    const forces::Memory{Int64}
    energy::Int64
end
function ForceEnergyFixedPoint{F_FBITS, E_FBITS}(nbeads::Integer) where {F_FBITS, E_FBITS}
    forces = Memory{Int64}(undef, 3*nbeads)
    fill!(forces, Int64(0))
    ForceEnergyFixedPoint{F_FBITS, E_FBITS}(forces, Int64(0))
end
function zero_force_energy!(force_energy::ForceEnergyFixedPoint)
    fill!(force_energy.forces, Int64(0))
    force_energy.energy = 0
    nothing
end
function get_nbeads(force_energy::ForceEnergyFixedPoint)
    length(force_energy.forces)÷3
end
function add_energy!(force_energy::ForceEnergyFixedPoint{F_FBITS, E_FBITS}, energy) where {F_FBITS, E_FBITS}
    force_energy.energy += convert_to_fixed_point(energy, E_FBITS)
    nothing
end
function zero_energy!(force_energy::ForceEnergyFixedPoint)
    force_energy.energy = 0
    nothing
end
Base.@propagate_inbounds function add_bead_force!(force_energy::ForceEnergyFixedPoint{F_FBITS, E_FBITS}, i::Integer, force) where {F_FBITS, E_FBITS}
    for j in 1:3
        force_energy.forces[3(i-1) + j] += convert_to_fixed_point(force[j], F_FBITS)
    end
    nothing
end
Base.@propagate_inbounds function zero_bead_force!(force_energy::ForceEnergyFixedPoint, i::Integer)
    for j in 1:3
        force_energy.forces[3(i-1) + j] = 0
    end
    nothing
end
function get_force!(force_energy::ForceEnergyFixedPoint{F_FBITS, E_FBITS}, force_out) where {F_FBITS, E_FBITS}
    Base.require_one_based_indexing(force_out)
    @argcheck 3*length(force_out) == length(force_energy.forces)
    @inbounds for i in eachindex(force_out)
        force_out[i] = @SVector [
            convert_to_floating_point(eltype(eltype(force_out)), force_energy.forces[3(i-1) + j], F_FBITS)
            for j in 1:3
        ]
    end
    nothing
end
Base.@propagate_inbounds function get_bead_force(force_energy::ForceEnergyFixedPoint{F_FBITS, E_FBITS}, i::Integer, ::Type{F}=Float64)::SVector{3, F} where {F_FBITS, E_FBITS, F}
    base = 3*(Int(i) - 1)
    @SVector [
        convert_to_floating_point(F, force_energy.forces[base + j], F_FBITS)
        for j in 1:3
    ]
end
function get_energy(force_energy::ForceEnergyFixedPoint{F_FBITS, E_FBITS}, ::Type{F}=Float64) where {F_FBITS, E_FBITS, F}
    convert_to_floating_point(F, force_energy.energy, E_FBITS)
end
function combine_force_energy!(
        fe_out::ForceEnergyFixedPoint{F, E},
        fe_in::ForceEnergyFixedPoint{F, E},
    ) where {F, E}
    @argcheck length(fe_out.forces) == length(fe_in.forces)
    fe_out.energy += fe_in.energy
    @inbounds for i in eachindex(fe_in.forces)
        fe_out.forces[i] += fe_in.forces[i]
    end
    fe_out
end


"""
    ForceNoEnergyFixedPoint{F_FBITS}(nbeads::Integer)

Identical fixed-point force storage to [`ForceEnergyFixedPoint`](@ref), but
with no energy accumulator: [`add_energy!`](@ref) is a no-op and
[`get_energy`](@ref) is intentionally not implemented, so reading the energy
errors instead of returning a made-up value.
"""
struct ForceNoEnergyFixedPoint{F_FBITS} <: ForceEnergy
    forces::Memory{Int64}
end
function ForceNoEnergyFixedPoint{F_FBITS}(nbeads::Integer) where {F_FBITS}
    forces = Memory{Int64}(undef, 3*nbeads)
    fill!(forces, Int64(0))
    ForceNoEnergyFixedPoint{F_FBITS}(forces)
end
function zero_force_energy!(force_energy::ForceNoEnergyFixedPoint)
    fill!(force_energy.forces, Int64(0))
    nothing
end
function get_nbeads(force_energy::ForceNoEnergyFixedPoint)
    length(force_energy.forces)÷3
end
function add_energy!(force_energy::ForceNoEnergyFixedPoint, energy)
    nothing
end
function zero_energy!(force_energy::ForceNoEnergyFixedPoint)
    nothing
end
Base.@propagate_inbounds function add_bead_force!(force_energy::ForceNoEnergyFixedPoint{F_FBITS}, i::Integer, force) where {F_FBITS}
    for j in 1:3
        force_energy.forces[3(i-1) + j] += convert_to_fixed_point(force[j], F_FBITS)
    end
    nothing
end
Base.@propagate_inbounds function zero_bead_force!(force_energy::ForceNoEnergyFixedPoint, i::Integer)
    for j in 1:3
        force_energy.forces[3(i-1) + j] = 0
    end
    nothing
end
function get_force!(force_energy::ForceNoEnergyFixedPoint{F_FBITS}, force_out) where {F_FBITS}
    Base.require_one_based_indexing(force_out)
    @argcheck 3*length(force_out) == length(force_energy.forces)
    @inbounds for i in eachindex(force_out)
        force_out[i] = @SVector [
            convert_to_floating_point(eltype(eltype(force_out)), force_energy.forces[3(i-1) + j], F_FBITS)
            for j in 1:3
        ]
    end
    nothing
end
Base.@propagate_inbounds function get_bead_force(force_energy::ForceNoEnergyFixedPoint{F_FBITS}, i::Integer, ::Type{F}=Float64)::SVector{3, F} where {F_FBITS, F}
    base = 3*(Int(i) - 1)
    @SVector [
        convert_to_floating_point(F, force_energy.forces[base + j], F_FBITS)
        for j in 1:3
    ]
end
function combine_force_energy!(
        fe_out::ForceNoEnergyFixedPoint{F},
        fe_in::ForceNoEnergyFixedPoint{F},
    ) where {F}
    @argcheck length(fe_out.forces) == length(fe_in.forces)
    @inbounds for i in eachindex(fe_in.forces)
        fe_out.forces[i] += fe_in.forces[i]
    end
    fe_out
end


"""
    DebugForceEnergy{T}(nbeads::Integer)

A [`ForceEnergy`](@ref) accumulator for testing that records every
[`add_energy!`](@ref) and [`add_bead_force!`](@ref) contribution individually
as type `T` instead of summing in place. Readers sum the recorded
contributions on demand, and [`zero_bead_force!`](@ref) removes the recorded
contributions to that bead. All bead indices are bounds-checked.
"""
struct DebugForceEnergy{T} <: ForceEnergy
    nbeads::Int64
    added_energy::Vector{T}
    added_bead_force::Vector{Tuple{Int64, SVector{3, T}}}
end
function DebugForceEnergy{T}(nbeads::Integer) where T
    DebugForceEnergy{T}(nbeads, [], [])
end
function zero_force_energy!(force_energy::DebugForceEnergy)
    empty!(force_energy.added_energy)
    empty!(force_energy.added_bead_force)
    nothing
end
function get_nbeads(force_energy::DebugForceEnergy)
    force_energy.nbeads
end
function add_energy!(force_energy::DebugForceEnergy, energy)
    push!(force_energy.added_energy, energy)
    nothing
end
function zero_energy!(force_energy::DebugForceEnergy)
    empty!(force_energy.added_energy)
    nothing
end
function add_bead_force!(force_energy::DebugForceEnergy, i::Integer, force)
    checkbounds(1:force_energy.nbeads, i)
    push!(force_energy.added_bead_force, (i, force))
    nothing
end
function zero_bead_force!(force_energy::DebugForceEnergy, i::Integer)
    checkbounds(1:force_energy.nbeads, i)
    filter!(((k, f),) -> k != i, force_energy.added_bead_force)
    nothing
end
function get_force!(force_energy::DebugForceEnergy, force_out)
    Base.require_one_based_indexing(force_out)
    @argcheck length(force_out) == force_energy.nbeads
    fill!(force_out, zero(eltype(force_out)))
    for (i, f) in force_energy.added_bead_force
        force_out[i] += f
    end
    nothing
end
function get_bead_force(force_energy::DebugForceEnergy{T}, i::Integer, ::Type{F}=Float64)::SVector{3, F} where {T, F}
    checkbounds(1:force_energy.nbeads, i)
    acc = zero(SVector{3, F})
    for (k, f) in force_energy.added_bead_force
        if k == i
            acc += SVector{3, F}(f)
        end
    end
    acc
end
function get_energy(force_energy::DebugForceEnergy, ::Type{F}=Float64) where F
    F(sum(force_energy.added_energy))
end


"""
    NullForceEnergy(nbeads::Integer)

A fake [`ForceEnergy`](@ref) if no forces or energy need to be calculated. All
writes are no-ops. Nothing is accumulated, so the getters
([`get_force!`](@ref), [`get_bead_force`](@ref), [`get_energy`](@ref), ...)
are intentionally not implemented and reading results errors instead of
returning made-up values.
"""
struct NullForceEnergy <: ForceEnergy
    nbeads::Int64
end
function zero_force_energy!(force_energy::NullForceEnergy)
    nothing
end
function get_nbeads(force_energy::NullForceEnergy)
    force_energy.nbeads
end
function add_energy!(force_energy::NullForceEnergy, energy)
    nothing
end
function zero_energy!(force_energy::NullForceEnergy)
    nothing
end
function add_bead_force!(force_energy::NullForceEnergy, i::Integer, force)
    nothing
end
function zero_bead_force!(force_energy::NullForceEnergy, i::Integer)
    nothing
end
