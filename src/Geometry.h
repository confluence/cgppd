#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "structures.h"
#include "Residue.h"
#include <vector>
#include <map>
#include <utility>
#include <set>
#include <tr1/tuple>

// Simple structures describing links between residues in a molecule
// i, j, etc. are indices in the residue array, which should be the same as indices into the original residue vector

struct Bond
{
    // residue indices
    int i;
    int j;

    // cache for tests
    double length;
    double potential;

    Bond(): i(0), j(0), length(0.0f), potential(0.0f) {}
    Bond(int i, int j, double length=0.0f, double potential=0.0f)
    {
        this->i = i;
        this->j = j;
        this->length = length;
        this->potential = potential;
    }

    bool operator<(const Bond &b) const
    {
        return tr1::make_tuple(i, j) < tr1::make_tuple(b.i, b.j);
    }

    bool operator==(const Bond &b) const
    {
        return tr1::make_tuple(i, j) == tr1::make_tuple(b.i, b.j);
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

    // cache for tests
    double theta;
    double potential;

    Angle(): i(0), j(0), k(0), theta(0.0f), potential(0.0f) {}
    Angle(int i, int j, int k, double theta=0.0f, double potential=0.0f)
    {
        this->i = i;
        this->j = j;
        this->k = k;
        this->theta = theta;
        this->potential = potential;
    }

    bool operator<(const Angle &a) const
    {
        return tr1::make_tuple(i, j, k) < tr1::make_tuple(a.i, a.j, a.k);
    }

    bool operator==(const Angle &a) const
    {
        return tr1::make_tuple(i, j, k) == tr1::make_tuple(a.i, a.j, a.k);
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

    // cache for tests
    double phi;
    double potential;

    Torsion(): i(0), j(0), k(0), l(0), phi(0.0f), potential(0.0f) {}
    Torsion(int i, int j, int k, int l, double phi=0.0f, double potential=0.0f)
    {
        this->i = i;
        this->j = j;
        this->k = k;
        this->l = l;
        this->phi = phi;
        this->potential = potential;
    }

    bool operator<(const Torsion &t) const
    {
        return tr1::make_tuple(i, j, k, l) < tr1::make_tuple(t.i, t.j, t.k, t.l);
    }

    bool operator==(const Torsion &t) const
    {
        return tr1::make_tuple(i, j, k, l) == tr1::make_tuple(t.i, t.j, t.k, t.l);
    }

    operator std::string() const {
        char str_char[256];
        sprintf(str_char, "Torsion(%d, %d, %d, %d)", i, j, k, l);
        return string(str_char);
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

    // INPUT
    void init(vector<Residue> residues, bool all_flexible, vector<segdata> segments);

    // OUTPUT -- calculate these once during init
    set<Bond> bonds;
    set<Angle> angles;
    set<Torsion> torsions;

    set<int> MC_local_residues;
    set<int> MC_crankshaft_residues;
    set<int> MC_flex_residues;

    set<int> branch(int i, int j, set<pair<int, int> > visited_edges=set<pair<int, int> >());
};

#endif
