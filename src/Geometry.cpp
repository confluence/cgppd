#include "Geometry.h"

void Graph::add_edge(int i, int j, bool flexible)
{
    vertices.insert(i);
    vertices.insert(j);

    if (i < j) {
        edges.insert(make_pair(i, j));
    } else {
        edges.insert(make_pair(j, i));
    }

    adjacency_map[i].insert(j);
    adjacency_map[j].insert(i);

    flexibility_map[make_pair(i, j)] = flexible;
    flexibility_map[make_pair(j, i)] = flexible;
}

bool Graph::is_flexible(int i, int j)
{
    return flexibility_map[make_pair(i, j)];
}

void Graph::init(vector<Residue> residues, bool all_flexible, vector<segdata> segments)
{
    // import residue chains
    for (int i = 0; i < residues.size() - 1; i++) {
        int j = i + 1;

        if (residues[i].chainId == residues[j].chainId) {
            add_edge(i, j, all_flexible); // flexible if the whole molecule is flexible
        }
    }

    // import segment data
    for (vector<segdata>::iterator s = segments.begin(); s != segments.end(); s++) {
        for (int i = 0; i < s->residue_indices.size() - 1; i++) {
            int j = i + 1;

            add_edge(s->residue_indices[i], s->residue_indices[j], true); // all of these are flexible
        }
    }

    // Temporary sets used for construction
    set<Bond> s_bonds;
    set<Angle> s_angles;
    set<Torsion> s_torsions;

    set<int> s_MC_local_residues;
    set<int> s_MC_crankshaft_residues;
    set<int> s_MC_flex_residues;

    for (set<pair<int, int> >::iterator e = edges.begin(); e != edges.end(); e++) {
        const int & j = e->first;
        const int & k = e->second;

        if (flexibility_map[*e]) {
            s_bonds.insert(Bond(j, k));
        }

        set<int> & neighbours_j = adjacency_map[j];
        set<int> & neighbours_k = adjacency_map[k];

        // TODO: is the outer check actually necessary?
        if (neighbours_j.size() > 1 && neighbours_k.size() > 1) { // jk is between at least two other edges
            for (set<int>::iterator i = neighbours_j.begin(); i != neighbours_j.end(); i++) {
                if (*i != k) {
                    for (set<int>::iterator l = neighbours_k.begin(); l != neighbours_k.end(); l++) {
                        if (*l != j) {
                            if (is_flexible(*i, j) || is_flexible(j, k) || is_flexible(k, *l)) {
                                s_torsions.insert(Torsion(*i, j, k, *l));
                            }
                        }
                    }
                }
            }
        }
    }

    for (set<int>::iterator v = vertices.begin(); v != vertices.end(); v++) {
        const int & j = *v;
        set<int> & neighbours = adjacency_map[j];

        bool has_rigid_neighbour = false;
        bool has_flexible_neighbour = false; 


        for (set<int>::iterator i = neighbours.begin(); i != neighbours.end(); i++) {
            if (is_flexible(*i, j)) {
                has_flexible_neighbour = true;
            } else {
                has_rigid_neighbour = true;
            }

            for (set<int>::iterator k = neighbours.begin(); k != neighbours.end(); k++) {
                if (*i < *k) {
                    if (is_flexible(*i, j) || is_flexible(j, *k)) {
                        s_angles.insert(Angle(*i, j, *k));
                    }
                }
            }
        }

        if (!has_rigid_neighbour) {
            s_MC_local_residues.insert(j);

            if (neighbours.size() > 1) {
                s_MC_crankshaft_residues.insert(j);
            }
        }

        if (has_flexible_neighbour) {
                s_MC_flex_residues.insert(j);
        }
    }

    // Copy the sets to vectors
    std::copy(s_bonds.begin(), s_bonds.end(), std::back_inserter(bonds));
    std::copy(s_angles.begin(), s_angles.end(), std::back_inserter(angles));
    std::copy(s_torsions.begin(), s_torsions.end(), std::back_inserter(torsions));

    std::copy(s_MC_local_residues.begin(), s_MC_local_residues.end(), std::back_inserter(MC_local_residues));
    std::copy(s_MC_crankshaft_residues.begin(), s_MC_crankshaft_residues.end(), std::back_inserter(MC_crankshaft_residues));
    std::copy(s_MC_flex_residues.begin(), s_MC_flex_residues.end(), std::back_inserter(MC_flex_residues));

    // Cache bonds, etc. for each residue

    for (int bi = 0; bi < bonds.size(); bi++) {
        Bond & b = bonds[bi];
        bonds_for_residue[b.i].insert(bi);
        bonds_for_residue[b.j].insert(bi);
    }

    for (int ai = 0; ai < angles.size(); ai++) {
        Angle & a = angles[ai];
        angles_for_residue[a.i].insert(ai);
        angles_for_residue[a.j].insert(ai);
        angles_for_residue[a.k].insert(ai);
    }

    for (int ti = 0; ti < torsions.size(); ti++) {
        Torsion & t = torsions[ti];
        torsions_for_residue[t.i].insert(ti);
        torsions_for_residue[t.j].insert(ti);
        torsions_for_residue[t.k].insert(ti);
        torsions_for_residue[t.l].insert(ti);
    }


}

void Graph::copy(Graph g) {
    vertices = g.vertices;
    edges = g.edges;
    adjacency_map = g.adjacency_map;
    flexibility_map = g.flexibility_map;

    bonds = g.bonds;
    angles = g.angles;
    torsions = g.torsions;

    MC_local_residues = g.MC_local_residues;
    MC_crankshaft_residues = g.MC_crankshaft_residues;
    MC_flex_residues = g.MC_flex_residues;

    bonds_for_residue = g.bonds_for_residue; 
    angles_for_residue = g.angles_for_residue;  
    torsions_for_residue = g.torsions_for_residue;  
}

set<int> Graph::branch(int i, int j, set<pair<int, int> > visited_edges)
{
    set<pair<int, int> > out_edges;
    set<int> & neighbours = adjacency_map[j];

    for (set<int>::iterator k = neighbours.begin(); k != neighbours.end(); k++) {
        pair<int, int> jk = make_pair(j, *k);
        if (*k != i && !visited_edges.count(jk)) {
            out_edges.insert(jk);
        }
    }

    set<int> this_branch;
    this_branch.insert(i);
    this_branch.insert(j);

    set<pair<int, int> > new_visited_edges = visited_edges;
    new_visited_edges.insert(make_pair(i, j));

    for (set<pair<int, int> >::iterator oe = out_edges.begin(); oe != out_edges.end(); oe++) {
        const set<int> & out_branch = branch(oe->first, oe->second, new_visited_edges);
        this_branch.insert(out_branch.begin(), out_branch.end());
    }

    return this_branch;
}

vector<int> Graph::neighbours(int i)
{
    set<int> & s_neighbours = adjacency_map[i];
    vector<int> v_neighbours;
    
    std::copy(s_neighbours.begin(), s_neighbours.end(), std::back_inserter(v_neighbours));
    return v_neighbours;
}
