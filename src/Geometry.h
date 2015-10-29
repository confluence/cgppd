#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "structures.h"
#include "Residue.h"
#include <vector>
#include <map>
#include <utility>
#include <set>
#include <tuple>
#include <algorithm>

#include "easylogging++.h"

// Simple structures describing links between residues in a molecule
// i, j, etc. are indices in the residue array, which should be the same as indices into the original residue vector

struct Pair
{
    // residue indices
    int i;
    int j;

    Pair(): i(0), j(0) {}
    Pair(int i, int j): i(i), j(j) {}

    bool operator<(const Pair &p) const
    {
        return make_tuple(i, j) < make_tuple(p.i, p.j);
    }

    bool operator==(const Pair &p) const
    {
        return make_tuple(i, j) == make_tuple(p.i, p.j);
    }

    operator std::string() const {
        char str_char[256];
        sprintf(str_char, "Pair(%d, %d)", i, j);
        return string(str_char);
    }
};

struct Bond
{
    // residue indices
    int i;
    int j;

    Bond(): i(0), j(0) {}
    Bond(int i, int j): i(i), j(j) {}

    bool operator<(const Bond &b) const
    {
        return make_tuple(i, j) < make_tuple(b.i, b.j);
    }

    bool operator==(const Bond &b) const
    {
        return make_tuple(i, j) == make_tuple(b.i, b.j);
    }

    operator std::string() const {
        char str_char[256];
        sprintf(str_char, "Bond(%d, %d)", i, j);
        return string(str_char);
    }
};

struct Angle
{
    // residue indices
    int i;
    int j;
    int k;

    Angle(): i(0), j(0), k(0) {}
    Angle(int i, int j, int k): i(i), j(j), k(k) {}

    bool operator<(const Angle &a) const
    {
        return make_tuple(i, j, k) < make_tuple(a.i, a.j, a.k);
    }

    bool operator==(const Angle &a) const
    {
        return make_tuple(i, j, k) == make_tuple(a.i, a.j, a.k);
    }

    operator std::string() const {
        char str_char[256];
        sprintf(str_char, "Angle(%d, %d, %d)", i, j, k);
        return string(str_char);
    }
};

struct Torsion
{
    // residue indices
    int i;
    int j;
    int k;
    int l;

    Torsion(): i(0), j(0), k(0), l(0) {}
    Torsion(int i, int j, int k, int l): i(i), j(j), k(k), l(l) {}

    bool operator<(const Torsion &t) const
    {
        return make_tuple(i, j, k, l) < make_tuple(t.i, t.j, t.k, t.l);
    }

    bool operator==(const Torsion &t) const
    {
        return make_tuple(i, j, k, l) == make_tuple(t.i, t.j, t.k, t.l);
    }

    operator std::string() const {
        char str_char[256];
        sprintf(str_char, "Torsion(%d, %d, %d, %d)", i, j, k, l);
        return string(str_char);
    }

    bool contains(const int &ri) const {
        return ri == i || ri == j || ri == k || ri == l;
    }

    bool contains(const Pair &p) const { // TODO replace with a Pair
        return contains(p.i) && contains(p.j);
    }
};

class Graph
{
    set<int> vertices;
    set<pair<int, int> > edges;
    map<int, set<int> > adjacency_map;
    map<pair<int, int>, bool> flexibility_map;

    void add_edge(int i, int j, bool flexible);
    bool is_flexible(int i, int j);
public:
    Graph();
    ~Graph();
    // INPUT
    void init(vector<Residue> residues, bool all_flexible, vector<segdata> segments, int num_chains);

    // OUTPUT -- calculate these once during init

    int num_chains;
    
    set<Bond> bonds;
    set<Angle> angles;
    set<Torsion> torsions;

    // These have to be vectors because we index them when we pick them randomly

    vector<int> MC_local_residues; // residue array indices
    vector<int> MC_crankshaft_residues; // residue array indices
    vector<int> MC_flex_residues; // residue array indices

    // These have to be vectors because their indices are mapped to UIDs
    vector<set<int> > rigid_domains;
    vector<Pair> segment_bonds; // bond indexes; these bonds connect residues which are not adjacent on the backbone

    set<Pair> indirect_neighbours; // pairs of residues to be excluded (or subtracted) from potential sum because they are close neighbours on opposite sides of a segment bond

    vector<int> neighbours(int i); // this has to be a vector because we pick a random neighbour
    set<int> branch(int i, int j, set<pair<int, int> > visited_edges=set<pair<int, int> >());
    set<int> rigid_domain_around(int i, set<pair<int, int> > visited_edges=set<pair<int, int> >());

    // These are replica-wide unique identifiers for chains, domains and segment bonds
    map<int, int> chain_uid; // chains aren't stored in the graph, but we can put this map here for consistency
    map<int, int> domain_uid;
    map<int, int> bond_uid;
    
    void assign_uids(Residue * residues, int & chain_offset, int & domain_offset, int & bond_offset);
};

#endif
