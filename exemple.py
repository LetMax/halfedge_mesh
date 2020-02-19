from halfedge_mesh import *
from halfedge_mesh.halfedge_mesh import *

# .off are supported

titres = ["strange", "CC_genre"]
titres = ["strange"]
for titre in titres:
    mesh = HalfedgeMesh("figures_normales/"+titre+".off")
    print(len(mesh.vertices))

    mesh.color_geodesique(mesh.vertices[2], titre)
    mesh.color_genre(titre)
    mesh.color_composante(titre)

    calcul_time(mesh.geodesique, [mesh.vertices[0]])
    calcul_time(mesh.set_composantes_connexes, [])
    calcul_time(mesh.genre_by_composante, [])
