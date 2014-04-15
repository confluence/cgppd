#!/usr/bin/env python

import random
import sys
from collections import defaultdict

# This is a Python mockup of the refactored C++ flexible linker geometry code

# the actual residue data object
class Residue(object):
    def __init__(self, label, chain):
        self.label = label
        self.chain = chain

    def __str__(self):
        return "Residue(%d, %s)" % (self.label, self.chain)
    
    
# classes for the pre-generated lists of geometric features
class Bond(object):
    def __init__(self, i, j):
        self.i = i
        self.j = j

    def __str__(self):
        return "Bond(%d, %d)" % (self.i, self.j)


class Angle(object):
    def __init__(self, i, j, k):
        self.i = i
        self.j = j
        self.k = k

    def __str__(self):
        return "Angle(%d, %d, %d)" % (self.i, self.j, self.k)


class Torsion(object):
    def __init__(self, i, j, k, l):
        self.i = i
        self.j = j
        self.k = k
        self.l = l

    def __str__(self):
        return "Torsion(%d, %d, %d, %d)" % (self.i, self.j, self.k, self.l)


def pretty_print(lst):
    print ", ".join([str(i) for i in lst])


# utility class for the graph
class Graph(object):
    def __init__(self):
        self.adjacency_map = {} # edges
        self.flexibility_map = {} # whether edges are flexible
        
        self.bonds = set()
        self.angles = set()
        self.torsions = set()

        self.mc_local_residues = set()
        self.mc_crankshaft_residues = set()
        self.mc_flex_residues = set()
        
        self.bonds_for_residue = defaultdict(set)
        self.angles_for_residue = defaultdict(set)
        self.torsions_for_residue = defaultdict(set)
        
        self.rigid_domains = []
        self.segment_bonds = []
        self.indirect_neighbours = set()

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
    
    def parse(self, residues, all_flexible, segments):
        vertices = list(range(len(residues)))

        for i, j in zip(vertices[:-1], vertices[1:]):
            if residues[i].chain == residues[j].chain:
                self.add_edge(i, j, False) # molecule is not all flexible

        for s in segments:
            for i, j in zip(s[:-1], s[1:]):
                self.add_edge(i, j, True) # these are flexible

        # create sets of geometric features for potential calculation

        for (i, j) in self.edges():
            if self.flexible(i, j):
                self.bonds.add(Bond(i, j))

        for j in self.vertices():
            neighbours = list(self.neighbours(j))
            if len(neighbours) > 1: # j is between at least two other residues
                for i in neighbours:
                    for k in neighbours:
                        if i < k and (self.flexible(i, j) or self.flexible(j, k)):
                            self.angles.add(Angle(i, j, k))

        for (j, k) in self.edges():
            neighbours_j = self.neighbours(j)
            neighbours_k = self.neighbours(k)

            if len(neighbours_j) > 1 and len(neighbours_k) > 1: # jk is between at least two other edges
                for i in (neighbours_j - set([k])):
                    for l in (neighbours_k - set([j])):
                        if self.flexible(i, j) or self.flexible(j, k) or self.flexible(k, l):
                            self.torsions.add(Torsion(i, j, k, l))

        # create sets of residues available for MC moves

        for i in self.vertices():
            neighbours = self.neighbours(i)

            if all(self.flexible(i, j) for j in neighbours):
                self.mc_local_residues.add(residues[i]) # local translation is allowed

                if len(neighbours) > 1:
                    self.mc_crankshaft_residues.add(residues[i]) # crankshaft is allowed

            if any(self.flexible(i, j) for j in neighbours):
                self.mc_flex_residues.add(residues[i]) # flex move about this residue is allowed
                
                
        # create set of bonds, angles and torsions for each residue
        
        for b in self.bonds:
            self.bonds_for_residue[b.i].add(b)
            self.bonds_for_residue[b.j].add(b)

        for a in self.angles:
            self.angles_for_residue[a.i].add(a)
            self.angles_for_residue[a.j].add(a)
            self.angles_for_residue[a.k].add(a)

        for t in self.torsions:
            self.torsions_for_residue[t.i].add(t)
            self.torsions_for_residue[t.j].add(t)
            self.torsions_for_residue[t.k].add(t)
            self.torsions_for_residue[t.l].add(t)
                
        # create set of rigid domains
        
        vertices_to_process = list(self.vertices())
        
        while vertices_to_process:
            i = vertices_to_process[0]
            domain = self.rigid_domain_around(i)
            
            for v in domain:
                vertices_to_process.remove(v)
                
            if len(domain) > 1:
                self.rigid_domains.append(domain)
                
        # find set of segment bonds
        
        for b in self.bonds:
            if residues[b.i].chain != residues[b.j].chain or b.j - b.i > 1 : # not neighbours on the backbone
                self.segment_bonds.append(b)
                
        # find set of indirect neighbour pairs to subtract from unbonded total
        
        for b in self.segment_bonds:
            torsions_around_bond = self.torsions_for_residue[b.i].intersection(self.torsions_for_residue[b.j]) # all torsions containing both i and j
            
            def add_indirect_neighbour(n):
                i, j = n
                if (residues[i].chain != residues[j].chain or abs(i - j) >= 4): # check that residues are not also direct neighbours
                    self.indirect_neighbours.add(n)

            for t in torsions_around_bond:
                add_indirect_neighbour((t.i, t.l))
                if (t.i, t.j) == (b.i, b.j) or (t.j, t.i) == (b.i, b.j):
                    add_indirect_neighbour((t.i, t.k))
                if (t.k, t.l) == (b.i, b.j) or (t.l, t.k) == (b.i, b.j):
                    add_indirect_neighbour((t.j, t.l))

    def rigid_domain_around(self, vertex, visited_edges=None):
        if visited_edges is None:
            visited_edges = set()
            
        i = vertex
        domain = set([i])
        
        out_vertices = set(j for j in self.neighbours(i) if not self.flexible(i, j) and (j, i) not in visited_edges)
        
        for j in out_vertices:
            domain |= self.rigid_domain_around(j, visited_edges | set([(i, j)]))
        
        return domain
    
    def branch(self, edge, visited_edges=None):
        if visited_edges is None:
            visited_edges = set()

        i, j = edge
        branch = set(edge)
        
        out_edges = set((j, k) for k in self.neighbours(j) if k != i and (j, k) not in visited_edges)
        
        for out_edge in out_edges:
            branch |= self.branch(out_edge, visited_edges | set([edge]))

        return branch
                
    def random_branch(self, i):
        j = random.choice(list(self.neighbours[i]))
        return self.branch((i, j))
    
    def random_neighbour_pair(self, i):
        return set(random.sample(list(self.neighbours(i)), 2))
    
if __name__ == "__main__":
    residues = []
    segments = [(72, 73, 74, 75, 123), (72 + 76, 73 + 76, 74 + 76, 75 + 76)] # indices (starting from zero)

    # Create diubiquitin

    for i in range(76):
        residues.append(Residue(i, "A"))

    for i in range(76):
        residues.append(Residue(i + 76, "B"))

    # create graph

    graph = Graph()
    graph.parse(residues, False, segments)

    # see if we did it right

    pretty_print(graph.bonds)
    pretty_print(graph.angles)
    pretty_print(graph.torsions)

    pretty_print(graph.mc_local_residues)
    pretty_print(graph.mc_crankshaft_residues)
    pretty_print(graph.mc_flex_residues)

    print "branch starting from (123, 122)"
    pretty_print(graph.branch((123, 122)))

    print "branch starting from (75, 123)"
    pretty_print(graph.branch((75, 123)))

    print "branch starting from (148, 149)"
    pretty_print(graph.branch((148, 149)))

    print "random pairs of neighbours for 123 (not that we can crankshaft 123)"
    pretty_print(graph.random_neighbour_pair(123))
    pretty_print(graph.random_neighbour_pair(123))
    pretty_print(graph.random_neighbour_pair(123))
    pretty_print(graph.random_neighbour_pair(123))

    print "random pair of neighbours for 74 (should always be the same)"
    pretty_print(graph.random_neighbour_pair(74))

    print "number of rigid domains"
    print len(graph.rigid_domains)
    
    for i, d in enumerate(graph.rigid_domains):
        print "domain %d" % i
        pretty_print(d)
        
    print "segment bonds"
    pretty_print(graph.segment_bonds)
    
    print "indirect neighbour pairs"
    pretty_print(graph.indirect_neighbours)
