from halfedge_mesh import *
from halfedge_mesh.halfedge_mesh import *

# .off are supported
titre = "tritore"
mesh = HalfedgeMesh("figures_normales/"+titre+".off")

# Returns a list of Vertex type (in order of file)--similarly for halfedges,
# and facets
# mesh.vertices

# The number of facets in the mesh
# print("number of facets : " + str(len(mesh.facets)))

# Get the 10th halfedge
#mesh.halfedges[10]

# Get the halfedge that starts at vertex 25 and ends at vertex 50
#mesh.get_halfedge(25, 50)
# print(mesh.set_composantes_connexes())

mesh.color_composante(titre)
mesh.color_geodesique(mesh.vertices[0], titre)
mesh.color_genre(titre)

# mesh.set_composantes_connexes()
# for v in mesh.vertices:
#     if v.composante == -1:
#         print(v.composante)
# calul_time(mesh.genre_by_composante, [])
# mesh.geodesique(mesh.vertices[2])

# for f in mesh.facets:
#     f.adjacent_faces
#     f.adjacent_vertices
#     f.adjacent_halfedges
# for f in mesh.halfedges:
#     f.adjacent_faces
#     f.adjacent_vertices
#     f.adjacent_halfedges
# for f in mesh.vertices:
#     f.adjacent_faces
#     f.adjacent_vertices
#     f.adjacent_halfedges
