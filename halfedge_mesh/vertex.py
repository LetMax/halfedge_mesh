import sys
from . import config
import math
import functools
from .functions import *

class Vertex:

    def __init__(self, x=0, y=0, z=0, index=None, halfedge=None, color=[]):
        """Create a vertex with given index at given point.

        x        - x-coordinate of the point
        y        - y-coordinate of the point
        z        - z-coordinate of the point
        index    - integer index of this vertex
        halfedge - a halfedge that points to the vertex
        """

        self.x = x
        self.y = y
        self.z = z
        self.vu = False
        self.traiter = False
        self.dist = 0
        self.composante = -1
        self.color = color

        self.index = index

        self.halfedge = halfedge

    def distance(self, S2) :
        S1 = self
        deltaX = abs(S1.x - S2.x)
        deltaY = abs(S1.y - S2.y)
        deltaZ = abs(S1.z - S2.z)
        distcarre = pow(deltaX,2) + pow(deltaY,2) + pow(deltaZ, 2)
        return math.sqrt(distcarre)

    def adjacent_vertices(self):
        adj = self.adjacent_halfedges()
        tab = []
        for i in adj:
            if i .opposite != None:
                tab.append(i.opposite.vertex)
        return tab

    def adjacent_faces(self):
        adj = self.adjacent_halfedges()
        tab = []
        for i in adj:
            tab.append(i.facet)
        return tab

    def adjacent_halfedges(self):

        continu = True
        first = self.halfedge
        tab = []
        next =  first.next.opposite

        if next != None:
            tab.append(next)
            while next != first and next != None:
                next = next.next.opposite
                if next == None:
                    break
                tab.append(next)

        if first == next:
            continu = False

        if continu:
            tab.append(first)
            if first.opposite != None:
                prev =  first.opposite.prev
                if prev != None:
                    tab.append(prev)
                    while prev != first and prev != None:

                        if prev.opposite == None:
                            break
                        prev = prev.opposite.prev
                        tab.append(prev)
        return tab

    def __eq__(x, y):
        return x.__key() == y.__key() and type(x) == type(y)

    def __key(self):
        return (self.x, self.y, self.z, self.index)

    def __hash__(self):
        return hash(self.__key())

    def get_vertex(self):
        return [self.x, self.y, self.z]

    def write_vertex(self, file):
        multiple_write(file, [self.x, self.y, self.z])
        file.write(" ")
        multiple_write(file,self.color)
        file.write("\n")
