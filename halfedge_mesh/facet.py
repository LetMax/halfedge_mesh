import sys
from . import config
import math
import functools
from .vertex import Vertex
from .functions import *

class Facet:
    def __init__(self, a=-1, b=-1, c=-1, index=None, vertex = [], halfedge=None, color=[]):
        """Create a facet with the given index with three vertices.

        a, b, c - indices for the vertices in the facet, counter clockwise.
        index - index of facet in the mesh
        halfedge - a Halfedge that belongs to the facet
        """
        self.vertex = vertex
        self.a = a
        self.b = b
        self.c = c
        self.index = index
        # halfedge going ccw around this facet.
        self.halfedge = halfedge
        self.vu = False
        self.traiter = False
        self.perimetre = -1
        self.composante = -1
        self.air = -1
        self.compar = []
        self.classe = -1
        self.ecart = float("inf")
        self.color = color

    def init_for_composante(self):
        self.composante = -1
        self.traiter = False
        self.vu = False

    def set_in_composante(self, composante, tab_color):
        self.vu = True
        self.composante = composante
        self.color = tab_color

    def calcul_ecart(self, classe, poids_functions):
        ecart = 0
        for i, j, poids in zip(self.compar, classe, poids_functions):
            ecart += poids*abs(i-j)
        return ecart

    def class_assignation(self, ecart, i):
        if ecart < self.ecart :
            self.ecart = ecart
            if self.classe != i:
                self.classe = i
                return True

    def adjacent_vertices(self):
        halfedges = self.adjacent_halfedges()
        tab = []
        for halfedge in halfedges:
            tab.append(halfedge.vertex)
        return tab

    def adjacent_halfedges(self):
        first = self.halfedge
        tab = [first]
        tmp = self.halfedge
        while tmp.next != first:
            tmp = tmp.next
            tab.append(tmp)

        return tab

    def adjacent_faces(self):
        facets = []
        adj = self.adjacent_halfedges()
        for i in adj :
            facets.append(i.opposite.facet)
        return facets

    def __eq__(self, other):
        return self.vertex == other.vertex and self.index == other.index and self.halfedge == other.halfedge

    def __hash__(self):
        res = hash(self.halfedge)
        for vert in self.vertex:
            res = res ^ hash(vert)
        res = res ^ hash(self.index) ^ hash(self.halfedges, self.vertex, self.index)
        return res
        # return hash(self.halfedge) ^ hash(self.a) ^ hash(self.b) ^ \
        #     hash(self.c) ^ hash(self.index) ^ \
        #     hash((self.halfedges, self.a, self.b, self.c, self.index))

    def get_normal(self):
        """Calculate the normal of facet

        Return a python list that contains the normal
        """
        vertex_a = [self.halfedge.vertex.x, self.halfedge.vertex.y,
                    self.halfedge.vertex.z]

        vertex_b = [self.halfedge.next.vertex.x, self.halfedge.next.vertex.y,
                    self.halfedge.next.vertex.z]

        vertex_c = [self.halfedge.prev.vertex.x, self.halfedge.prev.vertex.y,
                    self.halfedge.prev.vertex.z]

        # create edge 1 with vector difference
        edge1 = [u - v for u, v in zip(vertex_b, vertex_a)]
        edge1 = normalize(edge1)
        # create edge 2 ...
        edge2 = [u - v for u, v in zip(vertex_c, vertex_b)]
        edge2 = normalize(edge2)

        # cross product
        normal = cross_product(edge1, edge2)

        normal = normalize(normal)

        return normal

    def write_face(self, file):
        file.write(str(len(self.vertex)))
        file.write(" ")
        for v in self.vertex:
            file.write(str(v))
            file.write(" ")

        multiple_write(file,self.color)
        file.write("\n")
