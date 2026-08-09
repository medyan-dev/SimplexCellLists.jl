using Chairmarks
using Random
using SimplexCellLists
using StaticArrays

function load_cylinders()
    data = read(joinpath(@__DIR__, "ring-system-line-segs-f32.bin"))
    cylinders = collect(reinterpret(SimplexCellLists.LineSeg{Float32}, data))
    # center
    cylinders .+= (SA[SA[-2000,-2000,-200], SA[-2000,-2000,-200]],)
    cylinders
end

# Grid dimensions are 4000 nm by 4000 nm by 400 nm
# Benchmark building and query time with 50 nm, 100 nm, and 200 nm grid spacing and 10 nm, 50 nm, 100 nm, and 200 nm cutoff.
# Query time is measured for random points on the grid, and for cylinder end points

function build!(scl::LineSegCellList{Int64,Float32}, cylinders)
    cell_line_segs_clear!(scl)
    for (i, cyl) in enumerate(cylinders)
        cell_line_seg_add!(scl, cyl[1], cyl[2], i)
    end
end

function count_nearby(scl::LineSegCellList{Int64,Float32}, points, cutoff::Float32)
    total = 0
    for pos in points
        total = map_nearby_line_segs((entry, out) -> (out + 1, true), scl, pos, cutoff, total)
    end
    total
end

function build!(pcl::PointCellList{Int64,Float32}, points)
    cell_points_clear!(pcl)
    for (i, pos) in enumerate(points)
        cell_point_add!(pcl, pos, i)
    end
end

function count_nearby(pcl::PointCellList{Int64,Float32}, points, cutoff::Float32)
    total = 0
    for pos in points
        total = map_nearby_points((entry, out) -> (out + 1, true), pcl, pos, cutoff, total)
    end
    total
end

const grid_dims = SA[4000.0f0, 4000.0f0, 400.0f0]

cylinders = load_cylinders()
end_points = [cyl[2] for cyl in cylinders]
Random.seed!(1234)
random_points = [(rand(SimplexCellLists.Vec3{Float32}) .- 0.5f0) .* grid_dims for _ in eachindex(cylinders)]

for grid_spacing in Float32[100, 200]
    grid_size = round.(Int, grid_dims ./ grid_spacing)
    scl = LineSegCellList{Int64, Float32}(Tuple(grid_size), grid_spacing)
    println("$(length(cylinders)) line segments, grid spacing $grid_spacing nm, grid size $(Tuple(grid_size))")
    print("  build:")
    display(@b build!($scl, $cylinders) evals=1)
    build!(scl, cylinders)
    for cutoff in Float32[10, 50, 100]
        n = count_nearby(scl, random_points, cutoff)
        print("  cutoff $cutoff nm, $(length(random_points)) random points ($n hits):")
        display(@b count_nearby($scl, $random_points, $cutoff) evals=1)
        n = count_nearby(scl, end_points, cutoff)
        print("  cutoff $cutoff nm, $(length(end_points)) end points ($n hits):")
        display(@b count_nearby($scl, $end_points, $cutoff) evals=1)
    end
    println()
end

for grid_spacing in Float32[100, 200]
    grid_size = round.(Int, grid_dims ./ grid_spacing)
    pcl = PointCellList{Int64, Float32}(Tuple(grid_size), grid_spacing)
    println("$(length(end_points)) points, grid spacing $grid_spacing nm, grid size $(Tuple(grid_size))")
    print("  build:")
    display(@b build!($pcl, $end_points) evals=1)
    build!(pcl, end_points)
    for cutoff in Float32[10, 50, 100]
        n = count_nearby(pcl, random_points, cutoff)
        print("  cutoff $cutoff nm, $(length(random_points)) random points ($n hits):")
        display(@b count_nearby($pcl, $random_points, $cutoff) evals=1)
        n = count_nearby(pcl, end_points, cutoff)
        print("  cutoff $cutoff nm, $(length(end_points)) end points ($n hits):")
        display(@b count_nearby($pcl, $end_points, $cutoff) evals=1)
    end
    println()
end
