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

    for (set<pair<int, int> >::iterator e = edges.begin(); e != edges.end(); e++) {
        const int & j = e->first;
        const int & k = e->second;

        if (flexibility_map[*e]) {
            bonds.insert(Bond(j, k));
            int & b_i = bonds.size() - 1;
            bonds_for_residue[j].insert(b_i);
            bonds_for_residue[k].insert(b_i);
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
                                torsions.insert(Torsion(*i, j, k, *l));
                                int & t_i = torsions.size() - 1;
                                torsions_for_residue[*i].insert(t_i);
                                torsions_for_residue[j].insert(t_i);
                                torsions_for_residue[k].insert(t_i);
                                torsions_for_residue[*l].insert(t_i);
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
                        angles.insert(Angle(*i, j, *k));
                        int & a_i = angles.size() - 1;
                        angles_for_residue[*i].insert(a_i);
                        angles_for_residue[j].insert(a_i);
                        angles_for_residue[*k].insert(a_i);
                    }
                }
            }
        }

        if (!has_rigid_neighbour) {
            MC_local_residues.insert(j);

            if (neighbours.size() > 1) {
                MC_crankshaft_residues.insert(j);
            }
        }

        if (has_flexible_neighbour) {
                MC_flex_residues.insert(j);
        }
    }
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
