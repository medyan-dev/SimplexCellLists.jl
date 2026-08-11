using Aqua: Aqua
using SimplexCellLists: SimplexCellLists
Aqua.test_all(SimplexCellLists)

include("test-mindistance.jl")

include("test-pointcelllist.jl")

include("test-linesegcelllist.jl")

include("test-force-energy.jl")

include("test-neighbor-lists.jl")

include("test-collide-forces.jl")

