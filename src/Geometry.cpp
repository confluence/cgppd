#include "Geometry.h"

Graph::Graph() : num_chains(0)
{
}

Graph::~Graph()
{
}

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

void Graph::init(vector<Residue> residues, bool all_flexible, vector<segdata> segments, int num_chains)
{
    this->num_chains = num_chains;

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

    set<int> s_MC_local_residues;
    set<int> s_MC_crankshaft_residues;
    set<int> s_MC_flex_residues;

    for (set<pair<int, int> >::iterator e = edges.begin(); e != edges.end(); e++) {
        const int & j = e->first;
        const int & k = e->second;

        if (flexibility_map[*e]) {
            bonds.insert(Bond(j, k));
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

    std::copy(s_MC_local_residues.begin(), s_MC_local_residues.end(), std::back_inserter(MC_local_residues));
    std::copy(s_MC_crankshaft_residues.begin(), s_MC_crankshaft_residues.end(), std::back_inserter(MC_crankshaft_residues));
    std::copy(s_MC_flex_residues.begin(), s_MC_flex_residues.end(), std::back_inserter(MC_flex_residues));
    
    // Detect rigid domains
    
    // TODO TODO TODO: make sure this doesn't edit the vertices
    set<int> vertices_to_process(vertices);
    
    while (vertices_to_process.size()) {
        set<int>::iterator i = vertices_to_process.begin();
        set<int> domain = rigid_domain_around(*i);

        for (set<int>::iterator v = domain.begin(); v != domain.end(); v++) {
            vertices_to_process.erase(*v); // should erase by key value 
        }

        if (domain.size() > 1) {
            rigid_domains.push_back(domain);
        }
    }
    
    // Find segment bonds
    
    for (set<Bond>::iterator b = bonds.begin(); b != bonds.end(); b++) {
        if (residues[b->i].chainId != residues[b->j].chainId || abs(b->i - b->j) >= 4) {
            segment_bonds.push_back(Pair(b->i, b->j));
        }
    }
    
    // Find indirect neighbours
    
    for (vector<Pair>::iterator p = segment_bonds.begin(); p != segment_bonds.end(); p++) {
        set<Torsion> torsions_around_bond;

        for (set<Torsion>::iterator t = torsions.begin(); t != torsions.end(); t++) {
            if (t->contains(*p)) {
                torsions_around_bond.insert(*t);
            }
        }

        set<Pair> proposed;

        for (set<Torsion>::iterator t = torsions_around_bond.begin(); t != torsions_around_bond.end(); t++) {
            proposed.insert(Pair(t->i, t->l));

            if ((t->i, t->j) == (p->i, p->j) || (t->j, t->i) == (p->i, p->j)) {
                proposed.insert(Pair(t->i, t->k));
            }
            if ((t->k, t->l) == (p->i, p->j) || (t->l, t->k) == (p->i, p->j)) {
                proposed.insert(Pair(t->j, t->l));
            }
        }

        for (set<Pair>::iterator pr = proposed.begin(); pr != proposed.end(); pr++) {
            if (residues[pr->i].chainId != residues[pr->j].chainId || abs(pr->i - pr->j) >= 4) {
                if (pr->i < pr->j) {
                    indirect_neighbours.insert(*pr);
                } else {
                    indirect_neighbours.insert(Pair(pr->j, pr->i));
                }
            }
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

vector<int> Graph::neighbours(int i)
{
    set<int> & s_neighbours = adjacency_map[i];
    vector<int> v_neighbours;
    
    std::copy(s_neighbours.begin(), s_neighbours.end(), std::back_inserter(v_neighbours));
    return v_neighbours;
}

set<int> Graph::rigid_domain_around(int i, set<pair<int, int> > visited_edges)
{
    set<int> out_vertices;
    set<int> & neighbours = adjacency_map[i];

    for (set<int>::iterator j = neighbours.begin(); j != neighbours.end(); j++) {
        pair<int, int> ji = make_pair(*j, i);
        if (!is_flexible(i, *j) && !visited_edges.count(ji)) {
            out_vertices.insert(*j);
        }
    }
        
    set<int> domain;
    domain.insert(i);

    for (set<int>::iterator j = out_vertices.begin(); j != out_vertices.end(); j++) {
        set<pair<int, int> > new_visited_edges = visited_edges;
        new_visited_edges.insert(make_pair(i, *j));
        
        const set<int> & out_domain = rigid_domain_around(*j, new_visited_edges);
        domain.insert(out_domain.begin(), out_domain.end());
    }
    
    return domain;
}

void Graph::assign_uids(Residue * residues, int & chain_offset, int & domain_offset, int & bond_offset)
{
    // Assign the replica-wide uids

    for (int c = 0; c < num_chains; c++) {
        chain_uid[c] = chain_offset;
        chain_offset++;
    }
    
    for (int d = 0; d < rigid_domains.size(); d++) {
        domain_uid[d] = domain_offset;
        domain_offset++;
    }

    for (int b = 0; b < segment_bonds.size(); b++) {
        bond_uid[b] = bond_offset;
        bond_offset++;
    }
    
    // Assign domain UIDs to residues
    
    for (int d = 0; d < rigid_domains.size(); d++) {
        set<int> & domain = rigid_domains[d];
        for (set<int>::iterator i = domain.begin(); i != domain.end(); i++) {
            Residue & r = residues[*i];
            r.rigid_domain_UID = domain_uid[d];
        }
    }
        
    // Assign bond UIDs to residues
    
    for (int b = 0; b < segment_bonds.size(); b++) {
        Pair & sb = segment_bonds[b];
        residues[sb.i].segment_bond_UID = bond_uid[b];
        residues[sb.j].segment_bond_UID = bond_uid[b];
    }

    // Assign chain UIDs to residues, and pack float values to be transferred to the GPU, if necessary
    
    for (int i = 0; i < vertices.size(); i++) {
        residues[i].chain_UID = chain_uid[residues[i].chain_int_id];
        
#if USING_CUDA
        residues[i].pack_GPU_floats(bond_offset, chain_offset);
#endif
    }
}

