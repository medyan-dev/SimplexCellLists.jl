using Chairmarks
using Random
using SimplexCellLists
using StaticArrays

include("isosphere-system.jl")

# Benchmark building and query time of a TriangleCellList for the triangles
# of a sphere surface 1000 nm in radius, with edges about 26 nm long.
# Query time is measured for random points in the bounding box of the
# sphere, and for random points within 20 nm of the surface.

tris = [SimplexCellLists.Triangle{Float32}(pos[f[1]], pos[f[2]], pos[f[3]]) for f in faces]

function build!(tcl::TriangleCellList{Int64,Float32}, tris)
    cell_triangles_clear!(tcl)
    for (i, tri) in enumerate(tris)
        cell_triangle_add!(tcl, tri[1], tri[2], tri[3], i)
    end
end

function count_nearby(tcl::TriangleCellList{Int64,Float32}, points, cutoff::Float32)
    total = 0
    for p in points
        total = map_nearby_triangles((entry, sep, out) -> (out + 1, true), tcl, p, cutoff, total)
    end
    total
end

const grid_dims = SA[2100.0f0, 2100.0f0, 2100.0f0]

Random.seed!(1234)
n_points = 10^5
random_points = [(rand(SimplexCellLists.Vec3{Float32}) .- 0.5f0) .* grid_dims for _ in 1:n_points]
surface_points = map(1:n_points) do _
    v = normalize(randn(SimplexCellLists.Vec3{Float32}))
    v * (1000.0f0 + 40.0f0 * (rand(Float32) - 0.5f0))
end

for grid_spacing in Float32[25, 50, 100]
    grid_size = ceil.(Int, grid_dims ./ grid_spacing)
    tcl = TriangleCellList{Int64, Float32}(Tuple(grid_size), grid_spacing)
    println("$(length(tris)) triangles, grid spacing $grid_spacing nm, grid size $(Tuple(grid_size))")
    print("  build:")
    display(@b build!($tcl, $tris) evals=1)
    build!(tcl, tris)
    println("  cell entries: ", sum(c.len for c in tcl.cells))
    for cutoff in Float32[10, 25, 50]
        n = count_nearby(tcl, random_points, cutoff)
        print("  cutoff $cutoff nm, $n_points random points ($n hits):")
        display(@b count_nearby($tcl, $random_points, $cutoff) evals=1)
        n = count_nearby(tcl, surface_points, cutoff)
        print("  cutoff $cutoff nm, $n_points surface points ($n hits):")
        display(@b count_nearby($tcl, $surface_points, $cutoff) evals=1)
    end
    println()
end
