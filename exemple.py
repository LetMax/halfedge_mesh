import halfedge_mesh

# .off are supported
titre = "strange"
mesh = halfedge_mesh.HalfedgeMesh("figures_normales/"+titre+".off")

# Returns a list of Vertex type (in order of file)--similarly for halfedges,
# and facets
# mesh.vertices

# The number of facets in the mesh
# print("number of facets : " + str(len(mesh.facets)))

# Get the 10th halfedge
#mesh.halfedges[10]

# Get the halfedge that starts at vertex 25 and ends at vertex 50
#mesh.get_halfedge(25, 50)
mesh.color_geodesique(titre)
# mesh.set_composantes_connexes()
#mesh.color_composante(titre)
# mesh.geodesique(mesh.vertices[2])
#for f in mesh.facets:
#    print(f.vertex)
