from halfedge_mesh import *
from halfedge_mesh.halfedge_mesh import *

# .off are supported

titres = ["adapter", "bonhomme"]
# titres = ["adapter"]
for titre in titres:
    mesh = HalfedgeMesh("figures_normales/"+titre+".off")
    # print(len(mesh.vertices))
    #
    # mesh.color_geodesique(mesh.vertices[2], titre)
    # mesh.color_genre(titre)
    # mesh.color_composante(titre)
    #
    # calcul_time(mesh.geodesique, [mesh.vertices[0]])
    # calcul_time(mesh.set_composantes_connexes, [])
    # calcul_time(mesh.genre_by_composante, [])
    # mesh.classification(2)
    mesh.color_classe(3, titre)

tab = [0] * 50
for f in mesh.facets :
    tab[f.classe] += 1
print(len(mesh.facets))
print(len(mesh.vertices))
print(len(mesh.halfedges))
print(tab)
