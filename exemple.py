from halfedge_mesh import *
from halfedge_mesh.halfedge_mesh import *

# .off are supported

titres = ["CC_genre", "adapter", "bonhomme", "strange"]
# titres = ["adapter"]
for titre in titres:
    mesh = HalfedgeMesh("figures_normales/"+titre+".off")

    mesh.color_geodesique(mesh.vertices[2], titre)
    mesh.color_genre(titre)
    mesh.color_composante(titre)
    mesh.color_classe(2, [calcul_perimetre], [1], titre)
    mesh.color_composante_with_class(2, [calcul_perimetre,], [1], titre)
