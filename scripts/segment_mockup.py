#!/usr/bin/env python

# This is a Python mockup of the refactored C++ flexible linker geometry code

# the actual residue data object
class Residue(object):
    def __init__(self, label, chain):
        self.label = label
        self.chain = chain


# utility class for the graph
class Graph(object):
    def __init__(self):
        self.adjacency_map = {} # edges
        self.flexibility_map = {} # whether edges are flexible

    def add_edge(self, i, j, flexible):
        if i not in self.adjacency_map:
            self.adjacency_map[i] = set()
        self.adjacency_map[i].add(j)

        if j not in self.adjacency_map:
            self.adjacency_map[j] = set()
        self.adjacency_map[j].add(i)

        self.flexibility_map[(i, j)] = flexible
        self.flexibility_map[(j, i)] = flexible

    def vertices(self):
        return self.adjacency_map.keys()

    def edges(self):
        return [(i, j) for (i, j) in self.flexibility_map.keys() if i < j]

    def neighbours(self, i):
        return self.adjacency_map[i]

    def flexible(self, i, j):
        return self.flexibility_map[(i, j)]

    # TODO also add method for returning residues for flex move
    # TODO also add method for returning residues for crankshaft move


# classes for the pre-generated lists of geometric features
class Bond(object):
    def __init__(self, i, j):
        self.i = i
        self.j = j


class Angle(object):
    def __init__(self, i, j, k):
        self.i = i
        self.j = j
        self.k = k


class Torsion(object):
    def __init__(self, i, j, k, l):
        self.i = i
        self.j = j
        self.k = k
        self.l = l

residues = []
segments = [(72, 73, 74, 75), (72 + 76, 73 + 76, 74 + 76, 75 + 76)] # indices (starting from zero)

# Create diubiquitin

for i in range(76):
    residues.append(Residue(i + 1, "A"))

for i in range(76):
    residues.append(Residue(i + 77, "B"))

# create graph

graph = Graph()

vertices = list(range(len(residues)))

for i, j in zip(vertices[:-1], vertices[1:])
    if residues[i].chain == residues[j].chain:
        graph.add_edge(i, j, false) # molecule is not all flexible

for s in segments:
    for i, j in zip(segments[:-1], segments[1:]):
        graph.add_edge(i, j, true) # these are flexible

# create sets of geometric features for potential calculation

bonds = set()
angles = set()
torsions = set()

for (i, j) in graph.edges():
    bonds.add(Bond(i, j))

for j in graph.vertices():
    neighbours = graph.neighbours(j)
    if len(neighbours) > 1: # j is between at least two other residues
        for (i, k) in zip(neighbours[:-1], neighbours[1:]):
            if graph.flexible(i, j) or graph.flexible(j, k):
                angles.add(Angle(i, j, k))

for (j, k) in graph.edges():
    neighbours_j = graph.neighbours(j)
    neighbours_k = graph.neighbours(k)

    if len(neighbours_j) > 1 and len(neighbours_k) > 1: # jk is between at least two other edges
        for i in neighbours_j - set(k):
            for l in neighbours_k - set(j):
                if graph.flexible(i, j) or graph.flexible(j, k) or graph.flexible(k, l):
                    torsions.add(Torsion(i, j, k, l))

# create sets of residues available for MC moves

mc_local_residues = set()
mc_flex_residues = set()

for i in graph.vertices():
    neighbours = graph.neighbours(i)

    if all(graph.flexible(i, j) for j in neighbours):
        mc_local_residues.add(i) # local translation or crankshaft is allowed (should residues be excluded from cranksaft if they have > 2 neighbours? Maybe?)

    if any(graph.flexible(i, j) for j in neighbours):
        mc_flex_residues.add(i) # flex move about this residue is allowed
