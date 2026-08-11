using Chairmarks
using SimplexCellLists
using StaticArrays

function load_cylinders()
    data = read(joinpath(@__DIR__, "ring-system-line-segs-f32.bin"))
    cylinders = collect(reinterpret(SimplexCellLists.LineSeg{Float32}, data))
    # center
    cylinders .+= (SA[SA[-2000,-2000,-200], SA[-2000,-2000,-200]],)
    cylinders
end

# Convert the line segments into positions and CLines,
# merging neighboring end points that match.
function load_cline_system()
    cylinders = load_cylinders()
    pos = SimplexCellLists.Vec3{Float32}[]
    clines = CLineIdxPart[]
    for cyl in cylinders
        if isempty(pos) || last(pos) != cyl[1]
            push!(pos, cyl[1])
        end
        push!(pos, cyl[2])
        push!(clines, CLineIdxPart(length(pos) - 1))
    end
    pos, clines
end

pos, clines = load_cline_system()

radius = 3.0f0
skin = 7.0f0
stiffness = 500.0f0

# Exclude bonded pairs: consecutive clines that share a merged point
no_collide_pairs = empty_no_collide_pairs()
for k in 1:length(clines)-1
    if clines[k].i + UInt32(1) == clines[k+1].i
        push!(no_collide_pairs[Int(SimplexCellLists.CollideCLine_CLine)], UInt32(k) => UInt32(k+1))
    end
end

policy = DefaultCollisionPolicy()
inputs = NeighborListInputs(policy;
    clines,
    c_radius = fill(radius, length(clines)),
    c_params = fill(DefaultObjectParams(stiffness, UInt32(1), UInt32(0)), length(clines)),
    no_collide_pairs,
    skin,
)

nl = NeighborLists(policy)
setup_neighbors_sort_sweep!(nl, pos, inputs)
println("$(length(clines)) clines on $(length(pos)) points, radius $radius nm, skin $skin nm, $(length(nl.CCNL)) neighbor pairs")
print("  setup_neighbors_sort_sweep!:")
display(@b setup_neighbors_sort_sweep!($nl, $pos, $inputs) evals=1 seconds=2)
print("  setup_neighbors_naive!:")
display(@b setup_neighbors_naive!($nl, $pos, $inputs) evals=1)
