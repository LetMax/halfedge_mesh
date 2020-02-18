from halfedge_mesh import *
from halfedge_mesh.halfedge_mesh import *

# .off are supported


titres = ["strange", "CC_genre"]
titres = ["cube_bord1"]
for titre in titres:
    mesh = HalfedgeMesh("figures_normales/"+titre+".off")


    mesh.color_geodesique(mesh.vertices[2], titre)
    mesh.color_genre(titre)
    mesh.color_composante(titre)

    # calul_time(mesh.geodesique, [mesh.vertices[0]])
    # calul_time(mesh.set_composantes_connexes, [])
    # calul_time(mesh.genre_by_composante, [])

for f in mesh.facets:
    print(f.adjacent_vertices())
