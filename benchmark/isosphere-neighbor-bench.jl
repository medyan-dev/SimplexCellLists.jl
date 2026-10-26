using Chairmarks
using SimplexCellLists

include("isosphere-system.jl")

nl = NeighborLists(policy)
setup_neighbors_sort_sweep!(nl, pos, inputs)
nl_grid = NeighborLists(policy)
setup_neighbors_grid!(nl_grid, pos, inputs)
@assert is_neighbor_list_subset(nl, nl_grid) && is_neighbor_list_subset(nl_grid, nl)
println("isosphere radius $sphere_radius nm, $(length(triangles)) triangles, $(length(triangle_lines)) edges, $(length(pos)) vertices, radius $radius nm, skin $skin nm")
println("  $(length(nl.PTNL)) point-triangle pairs, $(length(nl.LLNL)) line-line pairs")
print("  setup_neighbors_sort_sweep!:")
display(@b setup_neighbors_sort_sweep!($nl, $pos, $inputs) evals=1 seconds=2)
print("  setup_neighbors_grid!:")
display(@b setup_neighbors_grid!($nl_grid, $pos, $inputs) evals=1 seconds=2)
