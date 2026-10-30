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

# Each system has copies of its positions, which the benchmark runs in turn.
# When the same positions are run over and over, the branch predictor learns
# which groups of edges have one in range, if there are fewer than about 6000
# such groups, and the time is too low. So the ring is jittered into copies,
# as its positions would change in a simulation. The isosphere systems have
# enough groups already.
systems = [
    ("ring", Ring, [Ring.pos .+ 0.3f0 .* randn(Xoshiro(s), eltype(Ring.pos), length(Ring.pos)) for s in 1:32]),
    ("isosphere", Isosphere, [Isosphere.pos]),
    # The isosphere is so regular that no edges are in range, so also move its vertices randomly.
    ("isosphere jittered", Isosphere, [Isosphere.pos .+ 3f0 .* randn(Xoshiro(42), eltype(Isosphere.pos), length(Isosphere.pos))]),
]

for (name, sys, copies) in systems
    local nls = map(copies) do pos
        nl = NeighborLists(sys.policy)
        setup_neighbors_grid!(nl, pos, sys.inputs)
        nl
    end
    local K = length(copies)
    println("$name: $(length(first(copies))) beads", K > 1 ? ", $K copies, edge counts are means" : "")
    for f in (:PTNL, :CCNL, :LLNL)
        isempty(getfield(first(nls), f)) && continue
        local nedges = sum(nl -> length(getfield(nl, f)), nls)
        local ninteracting = sum(((pos, nl),) -> count_interacting(pos, getfield(nl, f), sys.policy), zip(copies, nls))
        println("  $f: $(round(Int, nedges/K)) edges, $(round(Int, ninteracting/K)) interacting")
    end
    for (FE, calc_type) in (
            (ForceEnergyFixedPoint{30, 30}, Float64),
            (ForceEnergyFixedPoint{30, 30}, Float32),
            (ForceEnergyFloat64, Float64),
        )
        local fe = FE(length(first(copies)))
        # Positions in the calculation type.
        local pos_calc = map(pos -> map(p -> SimplexCellLists.Vec3{calc_type}(p), pos), copies)
        # Each sample is one call on the next copy.
        local k = Ref(0)
        print("  collide_forces! $FE $calc_type:")
        display(@b (k[] = mod1(k[] + 1, $K)) (zero_force_energy!($fe); collide_forces!($fe, $pos_calc[_], $nls[_], $calc_type)) evals=1 seconds=2)
    end
end
