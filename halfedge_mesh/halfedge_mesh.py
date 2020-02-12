import sys
from . import config
import math
import random
import functools
import time
from .facet import Facet
from .vertex import Vertex
from .halfedge import Halfedge

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

    def set_composantes_connexes(self) :
        composante = 1
        for v in self.vertices :
            v.composante = -1
            v.vu = False

        vertices = self.vertices[:]
        s = vertices[0]
        vertices.remove(s)
        s.vu = True
        s.composante = composante
        composante_tmp = []

        while len(vertices) != 0 :
            voisins = s.adjacent_vertices()
            for v in voisins :
                if composante_tmp.count(v) == 0 and v.vu == False:
                    composante_tmp.append(v)
                    v.composante = composante

            if len(composante_tmp) == 0 :
                s = vertices[0]
                composante += 1
                s.composante = composante
            else:
                s = composante_tmp[0]
                del composante_tmp[0]

            s.vu = True
            vertices.remove(s)
        print(composante)
        return composante


    def geodesique(self, s) :
        debut = time.time()
        print(debut)
        inf = float('inf')
        vertices = []
        for v in self.vertices :
            v.vu = False
            v.dist = inf
        sommets = self.vertices[:]
        sommets.remove(s)
        s.vu = True
        s.dist = 0

        while len(sommets) != 0 :
            voisins = s.adjacent_vertices()
            for v in voisins :
                nouvelle_dist = s.dist + s.distance(v)

                if v.vu == False:
                    vertices.append([v, nouvelle_dist])

                if v.dist > nouvelle_dist :
                    v.dist = nouvelle_dist

            vertices.sort(key = lambda vertices : vertices[1])
            s = vertices[0][0]
            while s.vu == True:
                del vertices[0]
                s = vertices[0][0]

            del vertices[0]
            s.vu = True
            sommets.remove(s)

        fin = time.time()
        print(fin-debut)

    def create_file(self, titre, color):
        f = open(titre, "w")
        if color:
            f.write("COFF\n")
        else:
            f.write("OFF\n")
        return f

    def close_file(self, f):
	   # f.write("</svg>")
       f.close()

    def color_composante(self, titre):

        nb_composantes = self.set_composantes_connexes()
        tab = []
        for i in range(nb_composantes):
            tmp = []
            tmp.append(random.uniform(0, 1))
            tmp.append(random.uniform(0, 1))
            tmp.append(random.uniform(0, 1))
            tab.append(tmp[:])
        print(tab)
        f = self.create_file("figures_compo/" + titre + "CompColor.off", True)

        self.write3(f, len(self.vertices), len(self.facets), len(self.edges))
        f.write("\n")
        for v in self.vertices:
            self.write3(f, v.x, v.y, v.z)
            f.write(" ")
            self.write3(f, 255 * tab[v.composante-1][0], 255 * tab[v.composante-1][1], 255 * tab[v.composante-1][2])
            f.write("\n")

        for face in self.facets:
            f.write("3 ")
            self.write3(f, face.a, face.b, face.c)
            f.write("\n")
        self.close_file(f)


    def color_geodesique(self, titre):
        self.geodesique(self.vertices[2])
        dist_max = 0
        inf = float('inf')
        for v in self.vertices:
            if dist_max < v.dist:
                dist_max = v.dist
        f = self.create_file("figures_geo/" + titre + "GeoColor.off", True)

        self.write3(f, len(self.vertices), len(self.facets), len(self.edges))
        f.write("\n")
        for v in self.vertices:
            self.write3(f, v.x, v.y, v.z)
            f.write(" ")
            if v.dist == 0 or v.dist == inf:
                self.write3(f, 0, 0, 0)
            else:
                tmp = 255 - (255 * (v.dist/dist_max))
                self.write3(f, 255, tmp, tmp)
            f.write("\n")

        for face in self.facets:
            f.write("3 ")
            self.write3(f, face.a, face.b, face.c)
            f.write("\n")
        self.close_file(f)

    def write3(self, file, one, two, three):
        file.write(str(one))
        file.write(" ")
        file.write(str(two))
        file.write(" ")
        file.write(str(three))



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

            # TODO: make general to support non-triangular meshes
            # Facets vertices are in counter-clockwise order
            facet = Facet(line[1], line[2], line[3], index, [line[1], line[2], line[3]])
            facets.append(facet)

            # create pairing of vertices for example if the vertices are
            # verts = [1,2,3] then zip(verts, verts[1:]) = [(1,2),(2,3)]
            # note: we skip line[0] because it represents the number of vertices
            # in the facet.
            all_facet_edges = list(zip(line[1:], line[2:]))
            all_facet_edges.append((line[3], line[1]))

            # For every halfedge around the facet
            for i in xrange(3):
                Edges[all_facet_edges[i]] = Halfedge()
                Edges[all_facet_edges[i]].facet = facet
                Edges[all_facet_edges[i]].vertex = vertices[
                    all_facet_edges[i][1]]
                vertices[all_facet_edges[i][1]].halfedge = Edges[all_facet_edges[i]]
                halfedge_count +=1

            facet.halfedge = Edges[all_facet_edges[0]]

            for i in xrange(3):
                Edges[all_facet_edges[i]].next = Edges[
                    all_facet_edges[(i + 1) % 3]]
                Edges[all_facet_edges[i]].prev = Edges[
                    all_facet_edges[(i - 1) % 3]]

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

def allclose(v1, v2):
    """Compare if v1 and v2 are close

    v1, v2 - any numerical type or list/tuple of numerical types

    Return bool if vectors are close, up to some epsilon specified in config.py
    """

    v1 = make_iterable(v1)
    v2 = make_iterable(v2)

    elementwise_compare = list(map(
        (lambda x, y: abs(x - y) < config.EPSILON), v1, v2))
    return functools.reduce((lambda x, y: x and y), elementwise_compare)


def make_iterable(obj):
    """Check if obj is iterable, if not return an iterable with obj inside it.
    Otherwise just return obj.

    obj - any type

    Return an iterable
    """
    try:
        iter(obj)
    except:
        return [obj]
    else:
        return obj


def dot(v1, v2):
    """Dot product(inner product) of v1 and v2

    v1, v2 - python list

    Return v1 dot v2
    """
    elementwise_multiply = list(map((lambda x, y: x * y), v1, v2))
    return functools.reduce((lambda x, y: x + y), elementwise_multiply)


def norm(vec):
    """ Return the Euclidean norm of a 3d vector.

    vec - a 3d vector expressed as a list of 3 floats.
    """
    return math.sqrt(functools.reduce((lambda x, y: x + y * y), vec, 0.0))


def normalize(vec):
    """Normalize a vector

    vec - python list

    Return normalized vector
    """
    if norm(vec) < 1e-6:
        return [0 for i in xrange(len(vec))]
    return list(map(lambda x: x / norm(vec), vec))


def cross_product(v1, v2):
    """ Return the cross product of v1, v2.

    v1, v2 - 3d vector expressed as a list of 3 floats.
    """
    x3 = v1[1] * v2[2] - v2[1] * v1[2]
    y3 = -(v1[0] * v2[2] - v2[0] * v1[2])
    z3 = v1[0] * v2[1] - v2[0] * v1[1]
    return [x3, y3, z3]

def create_vector(p1, p2):
    """Contruct a vector going from p1 to p2.

    p1, p2 - python list wth coordinates [x,y,z].

    Return a list [x,y,z] for the coordinates of vector
    """
    return list(map((lambda x,y: x-y), p2, p1))
