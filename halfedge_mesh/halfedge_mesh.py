import sys
from . import config
import math
import random
import functools
import time
from .facet import Facet
from .vertex import Vertex
from .halfedge import Halfedge
from .functions import *

# python3 compatibility
try:
    xrange
except NameError:
    xrange = range
try:
    dict.iteritems
except AttributeError:
    # Python 3
    def itervalues(d):
        return iter(d.values())
    def iteritems(d):
        return iter(d.items())
else:
    # Python 2
    def itervalues(d):
        return d.itervalues()
    def iteritems(d):
        return d.iteritems()


# TODO: Reorder functions

class HalfedgeMesh:

    def __init__(self, filename=None, vertices=[], halfedges=[], facets=[]):
        """Make an empty halfedge mesh.

           filename   - a string that holds the directory location and name of
               the mesh
            vertices  - a list of Vertex types
            halfedges - a list of HalfEdge types
            facets    - a list of Facet types
        """

        self.vertices = vertices
        self.halfedges = halfedges
        self.facets = facets
        self.filename = filename
        # dictionary of all the edges given indexes
        # TODO: Figure out if I need halfedges or if I should just use edges
        # Which is faster?
        self.edges = None

        if filename:
            self.vertices, self.halfedges, self.facets, self.edges = \
                    self.read_file(filename)


    def classification(self, nb_classe):
        color_tab = [0] * nb_classe
        for i in range(nb_classe):
            color_tab[i] = random_color()
        min = float("inf")
        max = 0

        for face in self.facets :
            tmp = calcul_perimetre(face)
            if tmp < min :
                min = tmp
            if tmp > max :
                max = tmp

        tab_classe = [0] * nb_classe
        ecart = (max-min)/nb_classe
        for i in range(len(tab_classe)) :
            tab_classe[i] = min + (i) * ecart

        for face in self.facets :
            for i,j in enumerate(tab_classe):
                if face.perimetre >= j :
                    face.classe = i
                    face.color = color_tab[i]

    def calcul_compar(self, compar_functions):
        tab_min = [float("inf")] * len(compar_functions)
        tab_max = [0] * len(compar_functions)

        for face in self.facets :
            face.classe = -1
            face.compar = []
            for i, function in enumerate(compar_functions):
                tmp = function(face)
                face.compar.append(tmp)
                if tmp < tab_min[i] :
                    tab_min[i] = tmp
                if tmp > tab_max[i] :
                    tab_max[i] = tmp
        return tab_min, tab_max

    def recalcul_classe(self, nb_classe, nb_functions):
        sum = [0] * nb_classe
        tmps = [0] * nb_functions

        for i in range(nb_classe):
            sum[i] = tmps[:]
        tab_classe = sum[:]
        count = [0] * nb_classe
        for face in self.facets :
            for i in range(nb_functions):
                sum[face.classe][i] += face.compar[i]
            count[face.classe] += 1
        for i in range(nb_classe):
            for j in range(nb_functions):
                if count[i] != 0:
                    tab_classe[i][j] = sum[i][j]/count[i]
                else:
                    tab_classe[i][j] = 0
        return tab_classe

    def k_moyenne(self, nb_classe, compar_functions, poids_functions):
        color_tab = [0] * nb_classe
        for i in range(nb_classe):
            color_tab[i] = random_color()

        tab_min, tab_max = self.calcul_compar(compar_functions)

        tab_classe = init_classe(nb_classe, len(compar_functions), tab_min, tab_max)

        change = True
        while change:
            change = False
            for face in self.facets:
                for i, j in enumerate(tab_classe):
                    ecart = face.calcul_ecart(j, poids_functions)
                    if face.class_assignation(ecart, i):
                        change = True

            tab_classe = self.recalcul_classe(nb_classe, len(compar_functions))

            for face in self.facets:
                face.ecart = face.calcul_ecart(tab_classe[face.classe], poids_functions)

        for face in self.facets:
            face.color = color_tab[face.classe]

    def color_classe(self, nb_classe, compar_functions, poids_functions, titre) :
        self.k_moyenne(nb_classe, compar_functions, poids_functions)
        self.write_mesh(True, "figures_classe/" + titre + "ClasseColor.off")

    def color_composante_with_class(self, nb_classe, compar_functions, poids_functions, titre):
        self.k_moyenne(nb_classe, compar_functions, poids_functions)
        self.set_composantes_connexes_with_classes()
        self.write_mesh(True, "figures_classe_compo/" + titre + "ClasseCompoColor.off")

    def new_composante(self, composante, s):
        composante += 1
        tab_color = random_color()
        s.composante = composante
        return tab_color, composante

    def set_composantes_connexes_with_classes(self):
        tab_color = random_color()
        composante_tmp = []
        composante = 1

        nb_facet = len(self.facets)
        for f in self.facets :
            f.init_for_composante()

        s = self.facets[0]
        s.set_in_composante(composante, tab_color)
        s.traiter = True
        nb_facet -= 1

        while nb_facet > 0 :
            voisins = s.adjacent_faces()
            for v in voisins:
                if v.vu == False and v.classe == s.classe:
                    composante_tmp.append(v)
                    v.set_in_composante(composante, tab_color)

            if len(composante_tmp) == 0 :
                s = self.retrieve_face()
                tab_color, composante = self.new_composante(composante, s)
            else:
                s = composante_tmp[0]
                del composante_tmp[0]

            if s.traiter == False:
                s.traiter = True
                nb_facet -= 1
        return composante

    def set_composantes_connexes(self) :
        tab_color = random_color()
        composante_tmp = []
        composante = 1

        nb_vert = len(self.vertices)
        for v in self.vertices :
            v.init_for_composante()

        s = self.vertices[0]
        s.set_in_composante(composante, tab_color)
        s.traiter = True
        nb_vert -= 1

        while nb_vert > 0 :
            voisins = s.adjacent_vertices()
            for v in voisins:
                if v.vu == False:
                    composante_tmp.append(v)
                    v.set_in_composante(composante, tab_color)

            if len(composante_tmp) == 0 :
                s = self.retrieve_vert()
                tab_color, composante = self.new_composante(composante, s)
            else:
                s = composante_tmp[0]
                del composante_tmp[0]

            if s.traiter == False:
                s.traiter = True
                nb_vert -= 1
        return composante

    def count_composante(self, nb_composante):
        tab_halfedge = [0] * nb_composante
        count = [0] * nb_composante
        tmp = [0] * 3
        tmp2 = []
        for i in range(nb_composante):
            count[i] = tmp[:]
            tab_halfedge[i] = tmp2[:]

        for v in self.vertices:
            count[v.composante-1][0] += 1

        for f in self.facets:
            count[f.halfedge.vertex.composante-1][1] += 1

        for h in self.halfedges:
            if h.opposite == None:
                tab_halfedge[h.vertex.composante-1].append(h)
            count[h.vertex.composante-1][2] += 1
        return count, tab_halfedge

    def genre_by_composante(self):
        nb_composante = self.set_composantes_connexes()
        tab = [0] * nb_composante
        count, tab_halfedge  = self.count_composante(nb_composante)
        for i in range(len(count)):

            euler = count[i][0] + count[i][1] - ((count[i][2]+len(tab_halfedge[i]))/2)
            nb_bord = count_bord(tab_halfedge[i])
            euler += nb_bord
            tab[i] = int((2 - euler)/2)
        return tab

    def set_color_genre(self):
        tab = self.genre_by_composante()
        genre0 = [255, 255, 255]
        genre1 = [127, 127, 255]
        genre2 = [0, 255, 255]
        genre3 = [255, 127, 127]
        genre4 = [255, 0, 255]
        genre5 = [127, 255, 127]
        genre6 = [255, 255, 0]
        genre7 = [0, 0, 0]
        color = [genre0, genre1, genre2, genre3, genre4, genre5, genre6, genre7]
        tab = self.genre_by_composante()
        for i in range(len(tab)):
            tab[i] = color[tab[i]]

        for vert in self.vertices:
            vert.color = tab[vert.composante-1]


    def retrieve_face(self):
        for face in self.facets:
            if face.traiter == False:
                return face

    def retrieve_vert(self):
        for vert in self.vertices:
            if vert.traiter == False:
                return vert

    def geodesique(self, s) :
        inf = float('inf')
        dist_max = 0
        continu = True
        nb_vert = len(self.vertices)
        voisins_calcules = []
        for v in self.vertices :
            v.traiter = False
            v.dist = inf

        s.traiter = True
        nb_vert -= 1
        s.dist = 0

        while nb_vert > 0 and continu:
            voisins = s.adjacent_vertices()
            for v in voisins :

                nouvelle_dist = s.dist + s.distance(v)
                if v.dist > nouvelle_dist :
                    v.dist = nouvelle_dist
                    if v.traiter == False:
                        voisins_calcules.append([v, nouvelle_dist])

            voisins_calcules.sort(key = lambda voisins_calcules : voisins_calcules[1])

            while s.traiter == True:
                if len(voisins_calcules) == 0:
                    continu = False
                    break
                else:
                    s = voisins_calcules[0][0]
                    del voisins_calcules[0]

            s.traiter = True
            if s.dist > dist_max:
                dist_max = s.dist
            nb_vert -= 1
        return dist_max

    def set_color_dist(self, sommet):
        inf = float("inf")
        dist_max = self.geodesique(sommet)
        for vert in self.vertices:
            if vert.dist == 0:
                vert.color = [0, 0, 0]
            elif vert.dist == inf:
                vert.color = [0, 255, 0]
            elif vert.dist == dist_max:
                vert.color = [127, 255, 255]
            else:
                tmp = 255 - (255 * (vert.dist/dist_max))
                vert.color = [255, tmp, tmp]

    def color_composante(self, titre):
        nb_composantes = self.set_composantes_connexes()
        self.write_mesh(True, "figures_compo/" + titre + "CompColor.off")

    def color_genre(self, titre):
        self.set_color_genre()
        self.write_mesh(True, "figures_genre/" + titre + "GenreColor.off")

    def color_geodesique(self, sommet, titre):
        self.set_color_dist(sommet)
        self.write_mesh(True, "figures_geo/" + titre + "GeoColor.off")

    def write_mesh(self, color, titre):
        file = create_file(titre, color)

        multiple_write(file, [len(self.vertices), len(self.facets), len(self.edges)])
        file.write("\n")
        for v in self.vertices:
            v.write_vertex(file)

        for face in self.facets:
            face.write_face(file)

        close_file(file)

    def __eq__(self, other):
        return (isinstance(other, type(self)) and
            (self.vertices, self.halfedges, self.facets) ==
            (other.vertices, other.halfedges, other.facets))

    def __hash__(self):
        return (hash(str(self.vertices)) ^ hash(str(self.halfedges)) ^ hash(str(self.facets)) ^
            hash((str(self.vertices), str(self.halfedges), str(self.facets))))

    def read_file(self, filename):
        """Determine the type of file and use the appropriate parser.

        Returns a HalfedgeMesh
        """
        try:
            with open(filename, 'r') as file:

                first_line = file.readline().strip().upper()

                if first_line != "OFF":
                    raise ValueError("Filetype: " + first_line + " not accepted")

                # TODO: build OBJ, PLY parsers
                parser_dispatcher = {"OFF": self.parse_off}

                return parser_dispatcher[first_line](file)

        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            return
        except ValueError as e:
            print("Value error: {0}:".format(e))
            return

    def read_off_vertices(self, file_object, number_vertices):
        """Read each line of the file_object and return a list of Vertex types.
        The list will be as [V1, V2, ..., Vn] for n vertices

        Return a list of vertices.
        """
        vertices = []

        # Read all the vertices in
        for index in xrange(number_vertices):
            line = file_object.readline().split()

            try:
                # convert strings to floats
                line = list(map(float, line))
            except ValueError as e:
                raise ValueError("vertices " + str(e))

            vertices.append(Vertex(line[0], line[1], line[2], index))

        return vertices

    def parse_build_halfedge_off(self, file_object, number_facets, vertices):
        """Link to the code:
        http://stackoverflow.com/questions/15365471/initializing-half-edge-
        data-structure-from-vertices

        Pseudo code:
        map< pair<unsigned int, unsigned int>, HalfEdge* > Edges;

        for each face F
        {
            for each edge (u,v) of F
            {
                Edges[ pair(u,v) ] = new HalfEdge();
                Edges[ pair(u,v) ]->face = F;
            }
            for each edge (u,v) of F
            {
                set Edges[ pair(u,v) ]->nextHalfEdge to next half-edge in F
                if ( Edges.find( pair(v,u) ) != Edges.end() )
                {
                    Edges[ pair(u,v) ]->oppositeHalfEdge = Edges[ pair(v,u) ];
                    Edges[ pair(v,u) ]->oppositeHalfEdge = Edges[ pair(u,v) ];
            }
        }

        """
        Edges = {}
        facets = []
        halfedge_count = 0
        #TODO Check if vertex index out of bounds

        # For each facet
        for index in xrange(number_facets):
            line = file_object.readline().split()

            # convert strings to ints
            line = list(map(int, line))
            nb_vert = line[0]
            # TODO: make general to support non-triangular meshes
            # Facets vertices are in counter-clockwise order
            vertex = [0] * nb_vert
            for i in range(nb_vert):
                vertex[i] = line[i+1]
            facet = Facet(line[1], line[2], line[3], index, vertex)
            facets.append(facet)

            # create pairing of vertices for example if the vertices are
            # verts = [1,2,3] then zip(verts, verts[1:]) = [(1,2),(2,3)]
            # note: we skip line[0] because it represents the number of vertices
            # in the facet.
            all_facet_edges = list(zip(line[1:], line[2:]))
            all_facet_edges.append((line[nb_vert], line[1]))

            # For every halfedge around the facet
            for i in xrange(nb_vert):
                Edges[all_facet_edges[i]] = Halfedge()
                Edges[all_facet_edges[i]].facet = facet
                Edges[all_facet_edges[i]].vertex = vertices[
                    all_facet_edges[i][1]]
                vertices[all_facet_edges[i][1]].halfedge = Edges[all_facet_edges[i]]
                halfedge_count +=1

            facet.halfedge = Edges[all_facet_edges[0]]

            for i in xrange(nb_vert):
                Edges[all_facet_edges[i]].next = Edges[
                    all_facet_edges[(i + 1) % nb_vert]]
                Edges[all_facet_edges[i]].prev = Edges[
                    all_facet_edges[(i - 1) % nb_vert]]

                # reverse edge ordering of vertex, e.g. (1,2)->(2,1)
                if all_facet_edges[i][2::-1] in Edges:
                    Edges[all_facet_edges[i]].opposite = \
                        Edges[all_facet_edges[i][2::-1]]

                    Edges[all_facet_edges[i][2::-1]].opposite = \
                        Edges[all_facet_edges[i]]

        return facets, Edges

    def parse_off(self, file_object):
        """Parses OFF files

        Returns a HalfedgeMesh
        """
        facets, halfedges, vertices = [], [], []

        # TODO Make ability to discard # lines
        vertices_faces_edges_counts = list(map(int, file_object.readline().split()))

        number_vertices = vertices_faces_edges_counts[0]
        vertices = self.read_off_vertices(file_object, number_vertices)

        number_facets = vertices_faces_edges_counts[1]
        facets, Edges = self.parse_build_halfedge_off(file_object,
                                                      number_facets, vertices)

        i = 0
        for key, value in iteritems(Edges):
            value.index = i
            halfedges.append(value)
            i += 1

        return vertices, halfedges, facets, Edges

    def get_halfedge(self, u, v):
        """Retrieve halfedge with starting vertex u and target vertex v

        u - starting vertex
        v - target vertex

        Returns a halfedge
        """
        return self.edges[(u, v)]

    def update_vertices(self, vertices):
        # update vertices
        vlist = []
        i = 0
        for v in vertices:
            vlist.append(Vertex(v[0], v[1], v[2], i))
            i += 1
        self.vertices = vlist

        hlist = []
        # update all the halfedges
        for he in self.halfedges:
            vi = he.vertex.index
            hlist.append(Halfedge(None, None, None, self.vertices[vi], None,
                he.index))

        flist = []
        # update neighboring halfedges
        for f in self.facets:
            hi = f.halfedge.index
            flist.append(Facet(f.a, f.b, f.c, f.index, f.vertex, hlist[hi]))
        self.facets = flist


        i = 0
        for he in self.halfedges:
            nextid = he.next.index
            oppid = he.opposite.index
            previd = he.prev.index

            hlist[i].next = hlist[nextid]
            hlist[i].opposite = hlist[oppid]
            hlist[i].prev = hlist[previd]


            fi = he.facet.index
            hlist[i].facet = flist[fi]
            i += 1

        self.halfedges = hlist
