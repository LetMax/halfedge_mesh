import sys
from . import config
import math
import random
import functools
import time

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

def create_file(titre, color):
    file = open(titre, "w")
    if color:
        file.write("COFF\n")
    else:
        file.write("OFF\n")
    return file

def close_file(file):
    file.close()

def multiple_write(file, param):
    for p in param:
        file.write(str(p))
        file.write(" ")

def count_bord(tab_halfedge):
    nb_bord = 0

    if tab_halfedge != []:
        print(" Pour la composante il y a au moins un bord")
        while tab_halfedge != []:
            first = tab_halfedge[0]
            del tab_halfedge[0]
            next = first.next_in_bord()
            while next != first:
                tab_halfedge.remove(next)
                next = next.next_in_bord()
            nb_bord += 1
    return nb_bord

def random_color():
    tmp = []
    tmp.append(random.uniform(0, 1) * 255)
    tmp.append(random.uniform(0, 1) * 255)
    tmp.append(random.uniform(0, 1) * 255)
    return tmp

def calcul_perimetre(face):
    perimetre = 0
    for v in face.adjacent_halfedges():
        perimetre += v.calcul_distance()
    face.perimetre = perimetre
    return perimetre

def return_z(face):
    return face.halfedge.vertex.z

def init_classe(nb_classe, nb_functions, tab_min, tab_max):
    tab_classe = [0] * nb_classe
    tmps = [0] * nb_functions
    ecart_tab = [0] * nb_functions
    for j in range(len(ecart_tab)):
        ecart_tab[j] = (tab_max[j]-tab_min[j])/nb_classe
    for i in range(nb_classe) :
        tab_classe[i] = tmps[:]
        for j in range(nb_functions):
            tab_classe[i][j] = tab_min[j] + i * ecart_tab[j]
    return tab_classe

def calcul_aire(face):
    aire = 0
    first = face.adjacent_vertices()[0]
    list1 = face.adjacent_vertices()
    list2 = list1[:]
    del list1[len(list1)-1]
    del list1[0]
    del list2[0]
    del list2[0]
    for i,j in zip(list1, list2):
        t = abs(first.determinant(i,j))
        aire += (t/2)
    face.aire = aire
    return aire

def calcul_time(fonction, params):
    debut = time.time()
    if params == []:
        fonction()
    else:
        fonction(params[0])
    fin = time.time()
    print("\n La fonction\n", fonction, "\n met", fin - debut, "secondes pour s'executer")
    return fin - debut
