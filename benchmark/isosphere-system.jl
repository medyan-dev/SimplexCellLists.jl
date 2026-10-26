# Defines `pos` and `inputs` for a triangulated sphere membrane.
using LinearAlgebra
using SimplexCellLists
using StaticArrays

# `icosphere_mesh` and `sphere_slerp` are copied from MEDYAN.jl `src/membrane/make_memb_sphere.jl`

"""
    icosphere_mesh(divisions::Int) -> (verts::Vector{SVector{3,Float64}}, faces::Vector{SVector{3,Int}})

Return a geodesic icosphere mesh on the unit sphere, with each icosahedron edge divided into `divisions` segments.
Triangles are oriented counterclockwise when viewed from outside.
"""
function icosphere_mesh(divisions::Int)
    divisions ≥ 1 || throw(ArgumentError("divisions must be at least 1"))
    # Golden ratio
    ϕ = (1 + sqrt(5.0)) / 2
    a, b = normalize(SA[1, ϕ])
    # Create icosahedron vertices
    verts = [
        SA[-a,  b,  0],
        SA[ a,  b,  0],
        SA[-a, -b,  0],
        SA[ a, -b,  0],
        SA[ 0, -a,  b],
        SA[ 0,  a,  b],
        SA[ 0, -a, -b],
        SA[ 0,  a, -b],
        SA[ b,  0, -a],
        SA[ b,  0,  a],
        SA[-b,  0, -a],
        SA[-b,  0,  a],
    ]
    # Faces of icosahedron (1-based)
    faces = SVector{3,Int}[
        SA[1,12,6],
        SA[1,6,2],
        SA[1,2,8],
        SA[1,8,11],
        SA[1,11,12],

        SA[2,6,10],
        SA[6,12,5],
        SA[12,11,3],
        SA[11,8,7],
        SA[8,2,9],

        SA[4,10,5],
        SA[4,5,3],
        SA[4,3,7],
        SA[4,7,9],
        SA[4,9,10],

        SA[5,10,6],
        SA[3,5,12],
        SA[7,3,11],
        SA[9,7,8],
        SA[10,9,2],
    ]

    n = divisions
    # Each icosahedron face is split into a triangular grid of n² triangles.
    # Grid point (i, j) of face (A, B, C) is i steps from edge AC towards B and j steps from edge AB towards C.
    # Vertices on icosahedron edges are placed by spherical interpolation, so that
    # neighboring faces agree on them. Each vertex gets a unique index from its edge or face,
    # so no lookup table of previously created vertices is needed.
    edges = Tuple{Int,Int}[]
    for f in faces, (u, v) in ((f[1], f[2]), (f[2], f[3]), (f[3], f[1]))
        u < v && push!(edges, (u, v))
    end
    sizehint!(verts, 10n^2 + 2)
    for (u, v) in edges, t in 1:n-1
        push!(verts, sphere_slerp(verts[u], verts[v], t/n))
    end
    # index of the vertex t steps along the icosahedron edge from u to v
    function edge_vert(u::Int, v::Int, t::Int)::Int
        t == 0 && return u
        t == n && return v
        e = findfirst(==(minmax(u, v)), edges)
        12 + (e-1)*(n-1) + (u < v ? t : n-t)
    end
    newfaces = SVector{3,Int}[]
    sizehint!(newfaces, 20n^2)
    for f in faces
        A, B, C = f
        # interior vertices of this face, in order of i then j
        # For corner A, take the grid line through (i, j) that is parallel to the opposite edge BC.
        # Its ends are on edges AB and AC, and the point is placed along the great-circle arc
        # between those ends, j steps out of i+j from the AB end. Doing the same from corners
        # B and C gives two more points. On the sphere the three points don't quite coincide,
        # so they are averaged to keep the grid symmetric.
        start = length(verts)
        pA, pB, pC = verts[A], verts[B], verts[C]
        for i in 1:n-2, j in 1:n-1-i
            k = n - i - j
            from_A = sphere_slerp(sphere_slerp(pA, pB, (i+j)/n), sphere_slerp(pA, pC, (i+j)/n), j/(i+j))
            from_B = sphere_slerp(sphere_slerp(pB, pC, (j+k)/n), sphere_slerp(pB, pA, (j+k)/n), k/(j+k))
            from_C = sphere_slerp(sphere_slerp(pC, pA, (k+i)/n), sphere_slerp(pC, pB, (k+i)/n), i/(k+i))
            push!(verts, normalize(from_A + from_B + from_C))
        end
        # index of the vertex at grid point (i, j)
        function vid(i::Int, j::Int)::Int
            if j == 0
                edge_vert(A, B, i)
            elseif i == 0
                edge_vert(A, C, j)
            elseif i + j == n
                edge_vert(B, C, j)
            else
                start + (i-1)*(n-2) - (i-1)*(i-2)÷2 + j
            end
        end
        for i in 0:n-1, j in 0:n-1-i
            push!(newfaces, SA[vid(i, j), vid(i+1, j), vid(i, j+1)])
            if i + j ≤ n - 2
                push!(newfaces, SA[vid(i+1, j), vid(i+1, j+1), vid(i, j+1)])
            end
        end
    end
    verts, newfaces
end

"""
Spherical linear interpolation between unit vectors `p` and `q`.

Angles are measured from the midpoint direction `(p + q)/|p + q|`, towards the
perpendicular direction `(q - p)/|q - p|`. `p` and `q` are at angles `∓θ/2`, where
`|p + q| = 2cos(θ/2)` and `|q - p| = 2sin(θ/2)`.
"""
function sphere_slerp(p, q, t)
    a = norm(p + q)
    b = norm(q - p)
    iszero(b) && return p
    s, c = sincos((2t - 1)*atan(b, a))
    c*(p + q)/a + s*(q - p)/b
end

# A membrane is added to the neighbor lists the same way as in MEDYAN.jl:
# its triangles, its edges as triangle lines, and its vertices as triangle points.
sphere_radius = 1000.0 # nm
divisions = 40 # edges about 26 nm long
verts, faces = icosphere_mesh(divisions)
pos = [SimplexCellLists.Vec3{Float32}(sphere_radius*v) for v in verts]
triangles = [TriangleIdxPart(f[1], f[2], f[3]) for f in faces]
edges = sort!(unique!([minmax(f[k], f[mod1(k+1, 3)]) for f in faces for k in 1:3]))
triangle_lines = [LineIdxPart(u, v) for (u, v) in edges]
triangle_points = [PointIdxPart(i) for i in eachindex(pos)]

radius = 10.0f0 # MEDYAN default membrane `radius_collide`
skin = 7.0f0
stiffness = 500.0f0

# Exclude vertices from their own triangles, and edges that share a vertex.
no_collide_pairs = empty_no_collide_pairs()
for (t, f) in enumerate(faces), v in f
    push!(no_collide_pairs[Int(SimplexCellLists.CollideTrianglePoint_Triangle)], UInt32(v) => UInt32(t))
end
vert_edges = [UInt32[] for _ in eachindex(pos)]
for (e, (u, v)) in enumerate(edges)
    push!(vert_edges[u], e)
    push!(vert_edges[v], e)
end
for es in vert_edges, i in eachindex(es), j in i+1:length(es)
    push!(no_collide_pairs[Int(SimplexCellLists.CollideTriangleLine_TriangleLine)], Pair(minmax(es[i], es[j])...))
end

policy = DefaultCollidePolicy()
params = DefaultObjectParams(stiffness, UInt32(1), UInt32(0))
inputs = NeighborListInputs(policy;
    triangles,
    t_radius = fill(radius, length(triangles)),
    t_params = fill(params, length(triangles)),
    triangle_lines,
    tl_radius = fill(radius, length(triangle_lines)),
    tl_params = fill(params, length(triangle_lines)),
    triangle_points,
    tp_radius = fill(radius, length(triangle_points)),
    tp_params = fill(params, length(triangle_points)),
    no_collide_pairs,
    skin,
)
