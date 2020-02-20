import sys
from . import config
import math
import functools
from .facet import Facet
from .vertex import Vertex
from .functions import *

class Halfedge:

    def __init__(self, next=None, opposite=None, prev=None, vertex=None,
                 facet=None, index=None):
        """Create a halfedge with given index.
        """
        self.opposite = opposite
        self.next = next
        self.prev = prev
        self.vertex = vertex
        self.facet = facet
        self.index = index

    def calcul_distance(self) :
        S1 = self.vertex
        S2 = self.opposite.vertex
        deltaX = abs(S1.x - S2.x)
        deltaY = abs(S1.y - S2.y)
        deltaZ = abs(S1.z - S2.z)
        distcarre = pow(deltaX,2) + pow(deltaY,2) + pow(deltaZ, 2)
        return math.sqrt(distcarre)

    def adjacent_faces(self):
        return self.vertex

    def adjacent_vertices(self):
        return self.opposite

    def adjacent_halfedges(self):
        return self.facet

    def next_in_bord(self):
        next = self.next
        while next.opposite != None:
            next = next.opposite.next
        return next

    def __eq__(self, other):
        # TODO Test more
        return  (type(self) == type(other)) and \
                (self.vertex == other.vertex) and \
                (self.prev.vertex == other.prev.vertex) and \
                (self.index == other.index)



    def __hash__(self):
        return hash(self.opposite) ^ hash(self.next) ^ hash(self.prev) ^ \
                hash(self.vertex) ^ hash(self.facet) ^ hash(self.index) ^ \
                hash((self.opposite, self.next, self.prev, self.vertex,
                    self.facet, self.index))

    def get_angle_normal(self):
        """Calculate the angle between the normals that neighbor the edge.

        Return an angle in radians
        """
        a = self.facet.get_normal()
        b = self.opposite.facet.get_normal()

        dir = [self.vertex.x - self.prev.vertex.x,
               self.vertex.y - self.prev.vertex.y,
               self.vertex.z - self.prev.vertex.z]
        dir = normalize(dir)

        ab = dot(a, b)

        args = ab / (norm(a) * norm(b))

        if allclose(args, 1):
            args = 1
        elif allclose(args, -1):
            args = -1

        assert (args <= 1.0 and args >= -1.0)

        angle = math.acos(args)

        if not (angle % math.pi == 0):
            e = cross_product(a, b)
            e = normalize(e)

            vec = dir
            vec = normalize(vec)

            if (allclose(vec, e)):
                return angle
            else:
                return -angle
        else:
            return 0
