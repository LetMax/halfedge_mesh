from halfedge_mesh import *
from halfedge_mesh.halfedge_mesh import *

# .off are supported

titres = ["CC_genre", "adapter", "bonhomme", "strange"]
# titres = ["bonhomme"]
for titre in titres:
    mesh = HalfedgeMesh("figures_normales/"+titre+".off")
    # print(len(mesh.vertices))
    #
    mesh.color_geodesique(mesh.vertices[2], titre)
    mesh.color_genre(titre)
    mesh.color_composante(titre)
    mesh.color_classe(2, [calcul_perimetre], [1], titre)
    mesh.color_composante_with_class(10, [calcul_perimetre], [1], titre)

    # mesh.color_composante(titre)
    # mesh.write_mesh(True, titre)
# for f in mesh.facets:
#     f.calcul_aire()

# tab = [0] * 50
# for f in mesh.facets :
#     tab[f.classe] += 1
# print(len(mesh.facets))
# print(len(mesh.vertices))
# print(len(mesh.halfedges))
# print(tab)
