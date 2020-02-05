import halfedge_mesh

# .off are supported
mesh = halfedge_mesh.HalfedgeMesh("strange.off")

# Returns a list of Vertex type (in order of file)--similarly for halfedges,
# and facets
mesh.vertices

# The number of facets in the mesh
print("number of facets : " + str(len(mesh.facets)))

# Get the 10th halfedge
#mesh.halfedges[10]

# Get the halfedge that starts at vertex 25 and ends at vertex 50
#mesh.get_halfedge(25, 50)

mesh.geodesique(mesh.vertices[2])
