using Chairmarks
using Random
using SimplexCellLists
using SimplexCellLists: nl_edge_forces!

module Ring
include("ring-system.jl")
end
module Isosphere
include("isosphere-system.jl")
end

# Count the edges whose force is nonzero
function count_interacting(pos, nl_list, policy)
    count(nl_list) do edge
        fe = DebugForceEnergy{Float64}(length(pos))
        nl_edge_forces!(fe, pos, edge, policy, Float64)
        !isempty(fe.added_energy)
    end
end

systems = [
    ("ring", Ring, Ring.pos),
    ("isosphere", Isosphere, Isosphere.pos),
    # The isosphere is so regular that no edges are in range, so also move its vertices randomly.
    ("isosphere jittered", Isosphere, Isosphere.pos .+ 3f0 .* randn(Xoshiro(42), eltype(Isosphere.pos), length(Isosphere.pos))),
]

for (name, sys, pos) in systems
    local nl = NeighborLists(sys.policy)
    setup_neighbors_grid!(nl, pos, sys.inputs)
    println("$name: $(length(pos)) beads")
    for f in (:PTNL, :CCNL, :LLNL)
        local list = getfield(nl, f)
        isempty(list) && continue
        println("  $f: $(length(list)) edges, $(count_interacting(pos, list, sys.policy)) interacting")
    end
    for (FE, calc_type) in (
            (ForceEnergyFixedPoint{30, 30}, Float64),
            (ForceEnergyFixedPoint{30, 30}, Float32),
            (ForceEnergyFloat64, Float64),
        )
        local fe = FE(length(pos))
        # Positions in the calculation type.
        local pos_calc = map(p -> SimplexCellLists.Vec3{calc_type}(p), pos)
        print("  collide_forces! $FE $calc_type:")
        display(@b (zero_force_energy!($fe); collide_forces!($fe, $pos_calc, $nl, $calc_type)) seconds=2)
    end
end
