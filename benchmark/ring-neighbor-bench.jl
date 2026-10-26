using Chairmarks
using SimplexCellLists

include("ring-system.jl")

nl = NeighborLists(policy)
setup_neighbors_sort_sweep!(nl, pos, inputs)
nl_grid = NeighborLists(policy)
setup_neighbors_grid!(nl_grid, pos, inputs)
@assert is_neighbor_list_subset(nl, nl_grid) && is_neighbor_list_subset(nl_grid, nl)
println("$(length(clines)) clines on $(length(pos)) points, radius $radius nm, skin $skin nm, $(length(nl.CCNL)) neighbor pairs")
print("  setup_neighbors_sort_sweep!:")
display(@b setup_neighbors_sort_sweep!($nl, $pos, $inputs) evals=1 seconds=2)
print("  setup_neighbors_grid!:")
display(@b setup_neighbors_grid!($nl_grid, $pos, $inputs) evals=1 seconds=2)
print("  setup_neighbors_naive!:")
display(@b setup_neighbors_naive!($nl, $pos, $inputs) evals=1)
