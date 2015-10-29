#include "Molecule.h"

using namespace std;

// TODO: we probably need a default constructor for dynamic arrays. :/
Molecule::Molecule() : residueCount(0), chainCount(0), length(0.0f), contiguous(false), bounding_value(0), moleculeRoleIdentifier(0.0f), index(-2), rotation(Quaternion(1.0f, 0, 0, 0)), only_molecule(false),
#if FLEXIBLE_LINKS
is_flexible(false),
#endif // FLEXIBLE_LINKS
volume(0.0f)
{
    memset(filename, 0, 256);
    memset(name, 0, 256);
    memset(last_MC_move, 0, 256);
}

Molecule::~Molecule()
{
    if (!contiguous && residueCount > 0)
    {
        delete [] Residues;
    }
};

void Molecule::init(const moldata mol, AminoAcids &a, int index, const float bounding_value)
{
    AminoAcidsData = a;

#if FLEXIBLE_LINKS
    torsion_data.loadData(TORSIONALPAIRDATA, a);
#endif // FLEXIBLE_LINKS

    this->index = index;
    this->bounding_value = bounding_value;

    vector<Residue> vResidues = initFromPDB(mol.pdbfilename);

    calculateVolume();
    calculate_length();
    
    if (mol.translate) {
        if (mol.px || mol.py || mol.pz) {
            translate(Vector3f(mol.px, mol.py, mol.pz));
        }
    } else {
        setPosition(Vector3f(mol.px, mol.py, mol.pz));
    }

    if (mol.ra) {
        Vector3double v = Vector3double(mol.rx,mol.ry,mol.rz);
        v.normalizeInPlace();
        rotate(v,mol.ra);
    }

    if (mol.crowder) {
        setMoleculeRoleIdentifier(CROWDER_IDENTIFIER);
    }

    graph.init(vResidues, mol.all_flexible, mol.segments, chainCount);

    if (mol.all_flexible || mol.segments.size()) {
        is_flexible = true;
    }
}

void Molecule::copy(const Molecule& m, Residue * contiguous_residue_offset)
{
    contiguous = true;
    memcpy(contiguous_residue_offset, m.Residues, m.residueCount*sizeof(Residue));
    Residues = contiguous_residue_offset;
    residueCount = m.residueCount;
    graph = m.graph; // TODO TODO TODO make sure this is OK. It should be fine for all molecule copies to use the same graph; it's read-only once constructed.
    bounding_value = m.bounding_value;
    AminoAcidsData = m.AminoAcidsData;
    center = m.center;
    rotation = m.rotation;
    length = m.length;
    moleculeRoleIdentifier = m.moleculeRoleIdentifier;
    volume = m.volume;
    index = m.index;
    strcpy(filename, m.filename);
    strcpy(name, m.name);

#if FLEXIBLE_LINKS
    is_flexible = m.is_flexible;
    torsion_data = m.torsion_data;
#endif // FLEXIBLE_LINKS
}

void Molecule::log_info(int index)
{
    VLOG(1) << "\t\t" << index << " " << filename << "; " << residueCount << " residues; centre: " << center << "; volume: " << volume << " A^3";
}

vector<string> Molecule::residue_sequence()
{
    vector<string> sequence;
    
    for (int i = 0; i < residueCount; i++) {
        sequence.push_back(string(AminoAcidsData.get(Residues[i].aminoAcidIndex).getSNAME()));
    }

    return sequence;
}

void Molecule::MC_backup_restore(const Molecule* m)
{
    memcpy(Residues, m->Residues, sizeof(Residue) * m->residueCount);

    residueCount = m->residueCount; // yes, we DO need to copy this, so that this function can be used in both directions.

    center = m->center;
    rotation = m->rotation;
    length = m->length;
    volume = m->volume;

    strcpy(last_MC_move, m->last_MC_move);
}

void Molecule::init_saved_molecule(int max_residue_count)
{
    Residues = new Residue[max_residue_count];
}

void Molecule::setMoleculeRoleIdentifier(float moleculeRoleIdentifier)
{
    this->moleculeRoleIdentifier = moleculeRoleIdentifier;
    if(moleculeRoleIdentifier == CROWDER_IDENTIFIER)
    {
#if REPULSIVE_CROWDING
        for (int i =0; i<residueCount; i++)
        {
            Residues[i].aminoAcidIndex = int(CROWDER_IDENTIFIER);
        }
#endif // REPULSIVE_CROWDING
    }
    else if (moleculeRoleIdentifier == PADDER_IDENTIFIER)
    {
        for (int i =0; i<residueCount; i++)
        {
            Residues[i].aminoAcidIndex = int(PADDER_IDENTIFIER);
        }
    }
}

void Molecule::translate(const Vector3f v)
{
    Vector3f new_centre = center + v;
    
#if BOUNDING_METHOD == BOUNDING_SPHERE
    if (centre_outside_boundary(new_centre)) {
        return;
    }
#elif BOUNDING_METHOD == PERIODIC_BOUNDARY
    new_centre = wrapped_centre(new_centre);
#endif // BOUNDING_METHOD

    setPosition(new_centre);
}

void Molecule::recalculate_relative_positions()
{
    // TODO: ugh, no, don't do this.
    for (size_t i=0; i<residueCount; i++)
    {
        Residues[i].relativePosition = Residues[i].position - center;
    }
}

Vector3f Molecule::recalculate_center()
{
    // TODO: this does not appear to be used. Why?
    Vector3f accumulator(0, 0, 0);

    for (size_t i=0; i<residueCount; i++)
    {
        accumulator += Residues[i].position;
    }

    return accumulator / residueCount;
}

void Molecule::setPosition(Vector3f v)
{
    center = v;

    for (size_t i = 0; i < residueCount; i++)
    {
        Residues[i].position = Residues[i].relativePosition + center;
    }
}

void Molecule::rotate(const Vector3double Raxis, const double angle)
{
    // No boundary checking is necessary because rotation does not change the centre.
    Quaternion q(angle, Raxis);
    setRotation(q);
}

void Molecule::setRotation(Quaternion q)
{
    q.normalize();
    rotation = q * rotation; // TODO: get rid of this; keeping a cumulative rotation value is obsolete because the molecule changes shape.
    for (size_t i=0; i<residueCount; i++)
    {
        Residues[i].set_rotation_about_center(q, center);
    }
}

Vector3f Molecule::normalised_random_vector_f(gsl_rng * rng)
{
    // TODO: do we actually need double brackets?
    return (Vector3f(gsl_rng_uniform(rng) - 0.5, 0.5 - gsl_rng_uniform(rng), gsl_rng_uniform(rng) - 0.5)).normalize();
}

Vector3double Molecule::normalised_random_vector_d(gsl_rng * rng)
{
    Vector3double x(gsl_rng_uniform(rng) - 0.5, 0.5 - gsl_rng_uniform(rng), gsl_rng_uniform(rng) - 0.5);
    x.normalizeInPlace();
    return x;
}

bool Molecule::centre_outside_boundary(Vector3f c)
{
    bool outside = c.sumSquares() < bounding_value * bounding_value;

    if (outside) {
        VLOG(1) << "Move rejected because it would move molecule beyond boundary sphere.";
    }

    return outside;
}

Vector3f Molecule::wrapped_centre(Vector3f c)
{
    Vector3f w_c;
    w_c.x = fmod(c.x, bounding_value);
    w_c.y = fmod(c.y, bounding_value);
    w_c.z = fmod(c.z, bounding_value);
    return w_c;
}

void Molecule::wrap_if_necessary()
{
    // Check if wrapping is required
    Vector3f wrapped_new_centre = wrapped_centre(center);

    // If wrapping is required, we need to update all the positions
    // Relative positions are not affected by wrapping; it's just a whole-molecule translation
    if ((wrapped_new_centre - center).magnitude()) {
        setPosition(wrapped_new_centre);
    }
}

void Molecule::rotate(gsl_rng * rng, const double rotate_step)
{
    // no boundary conditions required, because the centre doesn't change
    Vector3double v = normalised_random_vector_d(rng);
    rotate(v, rotate_step);
}

void Molecule::translate(gsl_rng * rng, const double translate_step)
{
    Vector3f v = translate_step * normalised_random_vector_f(rng);
    translate(v);
}

#if FLEXIBLE_LINKS
Vector3f Molecule::recalculate_center(Vector3f difference)
{
    return ((center * residueCount) + difference) / residueCount;
}

void Molecule::translate(Vector3f v, const int ri)
{
    Vector3f v_fraction = v / residueCount;
    Vector3f new_centre = center + v_fraction;

#if BOUNDING_METHOD == BOUNDING_SPHERE
    if (centre_outside_boundary(new_centre)) {
        return;
    }
#endif // BOUNDING_METHOD == BOUNDING_SPHERE

    // One position changes
    Residues[ri].position += v;
    Residues[ri].relativePosition += v;

    // Centre and all relative positions change
    center = new_centre;

    for (int i = 0; i < residueCount; i++) {
        Residues[i].relativePosition -= v_fraction;
    }

    calculate_length();

#if BOUNDING_METHOD == PERIODIC_BOUNDARY
    wrap_if_necessary();
#endif // BOUNDING_METHOD == PERIODIC_BOUNDARY
}

void Molecule::crankshaft(double angle, const bool flip_angle, const int ri)
{
    // calculate axis from neighbouring residues
    Vector3double raxis(Residues[ri + 1].position - Residues[ri - 1].position);

    // normalise axis
    raxis.normalizeInPlace();

    // TODO: clean this up -- do it in make_local_moves and just pass in the angle.
    if (flip_angle) {
        angle = -angle;
    }

    // calculate quaternion from angle and axis
    Quaternion q(angle, raxis);
    // TODO: is this necessary, or is it already normalised?
    q.normalize();

    // apply rotation to residue
    Vector3f relative_position = Residues[ri].position - Residues[ri - 1].position;
    relative_position = q.rotateVector(relative_position);
    
    Vector3f v = relative_position + Residues[ri - 1].position - Residues[ri].position;

    // Now this is reduced to a local translation of residue ri by vector v
    translate(v, ri);
}

void Molecule::flex(const Vector3double raxis, const double angle, const int ri, const int neighbour)
{
    // calculate quaternion from angle and axis
    Quaternion q(angle, raxis);

    // get the neighbours again
    const vector<int> & neighbours = graph.neighbours(ri);

    // select residues to rotate: the branch starting from the edge ri -> neighbour
    const set<int> & main_branch = graph.branch(ri, neighbour);
    set<int> branch;
    branch.insert(main_branch.begin(), main_branch.end());

    // possibly backtrack select other residues from the same rigid domain

    const set<int> & ri_domain = graph.rigid_domain_around(ri);

    // if the edge we're following is in a rigid domain
    if (ri_domain.find(neighbour) != ri_domain.end()) {
        for (vector<int>::const_iterator i = neighbours.begin(); i != neighbours.end(); i++) {
            // for every other out edge in the same domain
            if (*i != neighbour && ri_domain.find(*i) != ri_domain.end()) {
                // add that branch to the residues to be rotated
                const set<int> & additional_branch = graph.branch(ri, *i);
                branch.insert(additional_branch.begin(), additional_branch.end());
            }
        }
    }
    
    map<int, Vector3f> difference_vectors;
    Vector3f accumulated_difference(0,0,0);

    for (set<int>::iterator i = branch.begin(); i != branch.end(); i++) {
        Vector3f relative_position = Residues[*i].position - Residues[ri].position;
        relative_position = q.rotateVector(relative_position);
        difference_vectors[*i] = relative_position + Residues[ri].position - Residues[*i].position;
        accumulated_difference += difference_vectors[*i];
    }
    
    Vector3f v_fraction = accumulated_difference / residueCount;
    Vector3f new_centre = center + v_fraction;
    
    #if BOUNDING_METHOD == BOUNDING_SPHERE
    if (centre_outside_boundary(new_centre)) {
        return;
    }
    #endif // BOUNDING_METHOD == BOUNDING_SPHERE
    
    // apply the move

    ostringstream debug_log;
    debug_log << "Rotating neighbours ";
    
    for (set<int>::iterator i = branch.begin(); i != branch.end(); i++) {
        debug_log << *i << " ";
        Residues[*i].position += difference_vectors[*i];
        Residues[*i].relativePosition += difference_vectors[*i];
    }

    VLOG(1) << debug_log.str();
    
    // Centre and all relative positions change
    
    center = new_centre;

    for (int i = 0; i < residueCount; i++) {
        Residues[i].relativePosition -= v_fraction;
    }
    
    calculate_length();
    
    #if BOUNDING_METHOD == PERIODIC_BOUNDARY
    wrap_if_necessary();
    #endif // BOUNDING_METHOD == PERIODIC_BOUNDARY
}

void Molecule::flex(gsl_rng * rng, const double rotate_step)
{
    Vector3double raxis = normalised_random_vector_d(rng);

    vector<int> & index_set = graph.MC_flex_residues;
    int ri = index_set[gsl_rng_uniform_int(rng, index_set.size())];

    VLOG(1) << "Selected residue " << ri << " for flex move.";

    // pick a random neighbour
    const vector<int> & neighbours = graph.neighbours(ri);
    int neighbour = neighbours[gsl_rng_uniform_int(rng, neighbours.size())];

    VLOG(1) << "Selected neighbour " << neighbour;

    flex(raxis, rotate_step, ri, neighbour);
}

void Molecule::make_local_moves(gsl_rng * rng, const double rotate_step, const double translate_step)
{
    for (size_t i = 0; i < NUM_LOCAL_MOVES; i++) {
        //TODO: if linker too short for crankshaft, only return translate?
        //TODO: if crankshaft disabled, only return translate?
        int move = gsl_ran_bernoulli(rng, LOCAL_TRANSLATE_BIAS);

        switch (move)
        {
            case MC_LOCAL_TRANSLATE:
            {
                strcat(last_MC_move, "T");
                vector<int> & index_set = graph.MC_local_residues;
                int ri = index_set[gsl_rng_uniform_int(rng, index_set.size())];

                VLOG(1) << "Selected residue " << ri << " for local translation.";

                Vector3f v = LOCAL_TRANSLATE_STEP_SCALING_FACTOR * translate_step * normalised_random_vector_f(rng);
                translate(v, ri);
                break;
            }
            case MC_LOCAL_CRANKSHAFT: // TODO: remove if crankshaft disabled
            {
                strcat(last_MC_move, "C");
                vector<int> & index_set = graph.MC_crankshaft_residues;
                int ri = index_set[gsl_rng_uniform_int(rng, index_set.size())];

                VLOG(1) << "Selected residue " << ri << " for crankshaft.";

                bool flip = (bool) gsl_ran_bernoulli(rng, 0.5);
                crankshaft(rotate_step, flip, ri);
                break;
            }
            default:
                break;
        }
    }
}
#endif // FLEXIBLE_LINKS

int Molecule::get_MC_mutation_type(gsl_rng * rng)
{
#if FLEXIBLE_LINKS
    if (is_flexible) {
        if (!only_molecule) {
            return (int) gsl_ran_discrete(rng, MC_discrete_table);
        } else { // this is the only molecule; no translations or rotations
            return (int) gsl_ran_discrete(rng, MC_discrete_table) + 2; // shift by two integers to get values of flex and local mutations
        }
    } else {
        return (int) gsl_ran_bernoulli(rng, translate_rotate_bernoulli_bias);
    }
#else // if not FLEXIBLE_LINKS
    return (int) gsl_ran_bernoulli(rng, translate_rotate_bernoulli_bias);
#endif // FLEXIBLE_LINKS
}

void Molecule::make_MC_move(gsl_rng * rng, const double rotate_step, const double translate_step)
{
    int mutationType = get_MC_mutation_type(rng);
    memset(last_MC_move, 0, 256);

    switch (mutationType)
    {
        case MC_ROTATE:
        {
            strcat(last_MC_move, "Rotate");
            rotate(rng, rotate_step);
            break;
        }
        case MC_TRANSLATE:
        {
            strcat(last_MC_move, "Transl");
            translate(rng, translate_step);
            break;
        }
#if FLEXIBLE_LINKS
        case MC_FLEX:
        {
            strcat(last_MC_move, "Flex");
            flex(rng, rotate_step);
            break;
        }
        case MC_LOCAL:
        {
            strcat(last_MC_move, "Local ");
            make_local_moves(rng, rotate_step, translate_step);
            break;
        }
#endif // FLEXIBLE_LINKS
        default:
            break;
    }
}

vector<Residue> Molecule::initFromPDB(const char* pdbfilename)
{
    vector<Residue> vResidues;

    strcpy(filename, pdbfilename);
    // TODO: we only store this for printing.  Put it at the end?

    ifstream input(pdbfilename);
    if (!input.good())
    {
        LOG(INFO) << "Failed to open file: " << pdbfilename << " ("<<filename<<")";
        exit(0);
    }

    // for eachline in the file
    char line[512] = {' '};
    bool reached_end_marker = false;
    while (!reached_end_marker && !input.eof())
    {
        input.getline(line, 512);

        if (strncmp(line, "ATOM", 4) == 0)
        {
            char name[5] = {' '};

            sscanf(line + 12, "%s", name);

            // push the atom onto the vector of our structure if it is a CA
            if (strcmp(name, "CA") == 0) // it is a CA atom. center of our bead
            {
                char resName[4] = {' '};
                char chainID = ' ';
                int resSeq(0);
                float x(0);
                float y(0);
                float z(0);
                float occupancy(0);

                sscanf(line + 17, "%s", resName);
                sscanf(line + 21, "%c", &chainID);
                sscanf(line + 22, "%d", &resSeq);
                sscanf(line + 30, "%f", &x);
                sscanf(line + 38, "%f", &y);
                sscanf(line + 46, "%f", &z);
                sscanf(line + 54, "%f", &occupancy);

                Residue R;

                R.aminoAcidIndex = AminoAcidsData.getAminoAcidIndex(resName);
                R.electrostaticCharge = float(AminoAcidsData.get(R.aminoAcidIndex).electrostaticCharge); // TODO already a float; remove
                R.vanderWaalRadius = AminoAcidsData.data[R.aminoAcidIndex].vanderWaalRadius; // TODO use get method consistently or remove it

                R.position = Vector3f(x,y,z);
                R.relativePosition = Vector3f(0,0,0);
                R.chainId = chainID;
                R.chain_int_id = chainCount; // TODO: test this
                R.resSeq = resSeq;
                //R.sasa = 1.0f;
                //R.vanderWaalRadius = float(AminoAcidsData.data[R.aminoAcidIndex].vanderWaalRadius);
                R.moleculeId = index;
                vResidues.push_back(R);
            }
        }
        else if (strncmp(line, "TER", 3) == 0) // Optional,Mandatory if ATOM records exist.
        {
            chainCount++;
        }
        else if (strncmp(line, "END", 3) == 0)
        {
            reached_end_marker = true;
        }
    }
    input.close();
    
    if (vResidues.size() && !chainCount) {
        LOG(WARNING) << "Warning: check file " << pdbfilename << " for missing TER. Automatically setting number of chains to 1.";
        chainCount++;
    }

    // copy the residues into an array
    residueCount = vResidues.size();
    Residues = new Residue[residueCount];

    // TODO: use recalculate_center here, and also don't do this inside this function because it should go in init.

    // set the residues and accumulate the positions to determine the center of the molecule
    for(size_t r = 0; r < residueCount; r++)
    {
        Residue R = vResidues[r];
        memcpy (&Residues[r], &R, sizeof(R));
        center.x += R.position.x;
        center.y += R.position.y;
        center.z += R.position.z;

    }
    //calculate the center
    center = center / residueCount;

    // set the relative positions of all residues
    for(size_t r = 0; r < residueCount; r++)
    {
        Residues[r].relativePosition.x = Residues[r].position.x - center.x;
        Residues[r].relativePosition.y = Residues[r].position.y - center.y;
        Residues[r].relativePosition.z = Residues[r].position.z - center.z;
    }

    return vResidues;
}

void Molecule::calculateVolume()
{
    // simple volume, todo: surface enclosed calc which is nontrivial
    volume = 0.0f;

    for (size_t i=0; i<residueCount; i++)
    {
        float r(AminoAcidsData.get(Residues[i].aminoAcidIndex).vanderWaalRadius);
        volume += 4.0/3.0*r*r*r*M_PIl ;  // 4/3 * pi * r^3
    }
}

void Molecule::calculate_length()
{
    // Absolute distance, ignoring periodic boundary
    length = (float) Residues[0].distance(Residues[residueCount - 1]);
}

#if FLEXIBLE_LINKS



Potential Molecule::E(bool include_LJ_and_DH)
{
    Potential potential;

    // We may be using the GPU to calculate the internal LJ and DH, or not calculating it at all
    if (is_flexible) { // TODO TODO TODO are we still doing this?
        if (include_LJ_and_DH) {
            // We're calculating the LJ and DH on the CPU
            // LJ and DH between all residue pairs within molecule

            for (size_t i = 0; i < residueCount; i++) {
                for (size_t j = i + 1; j < residueCount; j++) {
#if !ASSUME_POLYMER_FOLDING_TEST
                    /* ignore pairs which are:
                        - in the same rigid domain
                        - close neighbours on the backbone in the same chain
                        - bonded to each other
                        - close neighbours across a bond */

                    if (!Residues[i].same_rigid_domain_as(Residues[j]) && !Residues[i].chain_neighbour_of(Residues[j]) && !Residues[i].bonded_to(Residues[j]) && !graph.indirect_neighbours.count(Pair(i, j))) {
#else // ASSUME_POLYMER_FOLDING_TEST
                    // More efficient handling of this special case.  Assume we are folding a single polymer.
                    // There's only one chain; we calculate all pairs except close neighbours on the backbone.
                    // Because there's only one chain we can just compare residue indices. And we know that j > i.
                    if ((j - i) >= 4) {
#endif
                        double r(Residues[i].distance(Residues[j], bounding_value) + EPS);
#if !LJ_OFF
                        /* Calculate LJ-type potential for each residue pair; increment molecule total. */
                        potential.increment_LJ(calculate_LJ(Residues[i], Residues[j], r, AminoAcidsData));
#endif
                        /* Calculate electrostatic potential for each residue pair; increment molecule total. */
                        potential.increment_DH(calculate_DH(Residues[i], Residues[j], r));
                    }
                }
            }
        } else {
            // We calculated the LJ and DH on the GPU; here we need to *subtract* the indirect neighbour total, which is too complex to do on the GPU for now
            for (set<Pair>::iterator p = graph.indirect_neighbours.begin(); p != graph.indirect_neighbours.end(); p++) {
                double r(Residues[p->i].distance(Residues[p->j], bounding_value) + EPS);
#if !LJ_OFF
                /* Calculate LJ-type potential for each residue pair; DECREMENT molecule total. */
                potential.increment_LJ(-calculate_LJ(Residues[p->i], Residues[p->j], r, AminoAcidsData));
#endif
                /* Calculate electrostatic potential for each residue pair; DECREMENT molecule total. */
                potential.increment_DH(-calculate_DH(Residues[p->i], Residues[p->j], r));
            }
        }

        for (set<Bond>::iterator b = graph.bonds.begin(); b != graph.bonds.end(); b++) {
            potential.increment_bond(calculate_bond(Residues, *b, bounding_value));
        }

        for (set<Angle>::iterator a = graph.angles.begin(); a != graph.angles.end(); a++) {
            potential.increment_angle(calculate_angle(Residues, *a));
        }

        for (set<Torsion>::iterator t = graph.torsions.begin(); t != graph.torsions.end(); t++) {
            potential.increment_torsion(calculate_torsion(Residues, *t, torsion_data));
        }
    }

    return potential;
}

#endif // FLEXIBLE_LINKS
