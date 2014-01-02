#include "Molecule.h"

using namespace std;

// TODO: we probably need a default constructor for dynamic arrays. :/
Molecule::Molecule() : residueCount(0), length(0.0f), contiguous(false), bounding_value(0), moleculeRoleIdentifier(0.0f), index(-2), rotation(Quaternion(1.0f, 0, 0, 0)),
#if FLEXIBLE_LINKS
linkCount(0), segmentCount(0), linkerCount(0), LJ(0), DH(0), update_LJ_and_DH(true), local_move_successful(false), is_flexible(false),
#endif // FLEXIBLE_LINKS
volume(0.0f)
{
    memset(filename, 0, 256);
    memset(name, 0, 256);
    memset(last_MC_move, 0, 256);
}

Molecule::~Molecule()
{
#if FLEXIBLE_LINKS
    if (linkerCount > 0)
    {
        delete [] Linkers;
    }

    if (segmentCount > 0)
    {
        delete [] Segments;
    }

    if (linkCount > 0)
    {
        delete [] Links;
    }
#endif // FLEXIBLE_LINKS

    if (!contiguous && residueCount > 0)
    {
        delete [] Residues;
    }
};

void Molecule::init(const char* pdbfilename, AminoAcids &a, int index, const float bounding_value)
{
    AminoAcidsData = a;
#if FLEXIBLE_LINKS
    torsions.loadData(TORSIONALPAIRDATA, a);
#endif // FLEXIBLE_LINKS

    initFromPDB(pdbfilename);
    this->index = index;
    this->bounding_value = bounding_value;

    calculate_length();
}

void Molecule::copy(const Molecule& m, Residue * contiguous_residue_offset)
{
    contiguous = true;
    memcpy(contiguous_residue_offset, m.Residues, m.residueCount*sizeof(Residue));
    Residues = contiguous_residue_offset;
    residueCount = m.residueCount;

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
    Links = new Link[m.linkCount];
    memcpy(Links, m.Links, m.linkCount*sizeof(Link));
    linkCount = m.linkCount;

    Segments = new Segment[m.segmentCount];
    memcpy(Segments, m.Segments, m.segmentCount*sizeof(Segment));
    segmentCount = m.segmentCount;

    Linkers = new int[m.linkerCount];
    memcpy(Linkers, m.Linkers, m.linkerCount*sizeof(int));
    linkerCount = m.linkerCount;

    torsions = m.torsions;

    LJ = m.LJ;
    DH = m.DH;
    update_LJ_and_DH = m.update_LJ_and_DH;
#endif // FLEXIBLE_LINKS
}

void Molecule::log_info(int index, int verbosity)
{
    LOG(verbosity > 1, "\t\t%d %s; %d residues; centre: (%0.3f, %0.3f, %0.3f); volume: %.3f A^3\n", index, filename, residueCount, center.x, center.y, center.z, volume);
}

void Molecule::MC_backup_restore(const Molecule* m)
{
    memcpy(Residues, m->Residues, sizeof(Residue) * m->residueCount);
    residueCount = m->residueCount;

    center = m->center;
    rotation = m->rotation;
    length = m->length;

#if FLEXIBLE_LINKS
    // we need to save / restore all the cached potentials
    LJ = m->LJ;
    DH = m->DH;
    update_LJ_and_DH = m->update_LJ_and_DH;

    memcpy(Links, m->Links, sizeof(Link) * m->linkCount);
    linkCount = m->linkCount;

    memcpy(Segments, m->Segments, sizeof(Segment) * m->segmentCount);
    segmentCount = m->segmentCount;
#endif // FLEXIBLE_LINKS
}

#if FLEXIBLE_LINKS
void Molecule::init_saved_molecule(int max_residue_count, int max_segment_count)
{
    Residues = new Residue[max_residue_count];
    Links = new Link[max_residue_count - 1];
    Segments = new Segment[max_segment_count]; // TODO: make this a vector; then we don't have to know the size
}
#else // if not FLEXIBLE_LINKS
void Molecule::init_saved_molecule(int max_residue_count)
{
    Residues = new Residue[max_residue_count];
}
#endif // FLEXIBLE_LINKS

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
    new_center = center + v;

    if (test_boundary_conditions())
    {
        setPosition(new_center_wrapped);
    }
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

    for (size_t i=0; i<residueCount; i++)
    {
        // more efficient to do it this way than with the overloaded +
        // TODO: is that actually true?
        Residues[i].position.x = Residues[i].relativePosition.x + center.x;
        Residues[i].position.y = Residues[i].relativePosition.y + center.y;
        Residues[i].position.z = Residues[i].relativePosition.z + center.z;
    }
}

void Molecule::rotate(const Vector3double Raxis, const double angle)
{
    Quaternion q(angle, Raxis);
    setRotation(q);
}

//TODO: see if component-wise addition is faster
void Molecule::setRotation(Quaternion q)
{
    q.normalize();
    rotation = q * rotation; // TODO: get rid of this
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

bool Molecule::test_boundary_conditions()
{
#if BOUNDING_METHOD == BOUNDING_SPHERE
    if (new_center.sumSquares() < bounding_value * bounding_value)
    {
        // Don't do anything to the new centre; just say that it's OK
        new_center_wrapped = new_center;
        return true;
    }
    else
    {
        LOG(DEBUG_BOUNDARY, "Move rejected because it would move molecule beyond boundary sphere.");
        return false;
    }
#elif BOUNDING_METHOD == PERIODIC_BOUNDARY
    new_center_wrapped.x = fmod(new_center.x, bounding_value);
    new_center_wrapped.y = fmod(new_center.y, bounding_value);
    new_center_wrapped.z = fmod(new_center.z, bounding_value);

    center_wrap_delta = new_center_wrapped - new_center;

    return true;
#endif // BOUNDING_METHOD
}

void Molecule::wrap_if_necessary()
{
    if (center_wrap_delta.magnitude())
    {
        // If the molecule has wrapped, we have to move everything
        LOG(DEBUG_LOCAL, "Molecule has wrapped around the periodic boundary.\t");
        setPosition(new_center_wrapped);
    }
    else
    {
        // otherwise just update the center
        center = new_center_wrapped;
    }
}

void Molecule::rotate(gsl_rng * rng, const double rotate_step)
{
    // no boundary conditions required, because the centre doesn't change
    rotate(normalised_random_vector_d(rng), rotate_step);
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

void Molecule::mark_cached_potentials_for_update(const int ri, const bool crankshaft)
{
    // For any move, update the LJ and DH for the entire molecule and for all segments which contain the residue

    update_LJ_and_DH = true;

    for (int i = 0; i < segmentCount; i++) {
        if (ri >= Segments[i].start && ri <= Segments[i].end) {
            Segments[i].update_LJ_and_DH = true;
        }
    }

    // We don't need to update bonds or angles for crankshaft moves; only torsions.
    // We're not checking whether the marked residues and links fall within flexible segments, because the potential calculation ignores unnecessary segments.

    if (ri > 1) {
        Links[ri - 2].update_e_torsion = true;
    }

    if (ri > 0) {
        if (!crankshaft) {
            Links[ri - 1].update_e_bond = true;
            Residues[ri - 1].update_e_angle = true;
        }

        Links[ri - 1].update_e_torsion = true;
    }

    if (!crankshaft) {
        Residues[ri].update_e_angle = true;
    }

    if (ri + 1 < residueCount) {
        if (!crankshaft) {
            Links[ri].update_e_bond = true;
            Residues[ri + 1].update_e_angle = true;
        }

        Links[ri].update_e_torsion = true;
    }

    if (ri + 2 < residueCount) {
        Links[ri + 1].update_e_torsion = true;
    }
}

void Molecule::translate(Vector3f v, const int ri)
{
    new_center = recalculate_center(v);

    // apply boundary conditions, possibly rejecting the move
    if (test_boundary_conditions())
    {
        wrap_if_necessary();

        // translate one residue
        LOG(DEBUG_LOCAL, "Residue: %d Old: (%f, %f, %f)\t", ri, Residues[ri].position.x, Residues[ri].position.y, Residues[ri].position.z);
        Residues[ri].position += v;
        LOG(DEBUG_LOCAL, "New: (%f, %f, %f)\n", Residues[ri].position.x, Residues[ri].position.y, Residues[ri].position.z);

        mark_cached_potentials_for_update(ri, false);
        local_move_successful = true;
    }
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
    Residues[ri].new_position = relative_position + Residues[ri - 1].position;

    new_center = recalculate_center(Residues[ri].new_position - Residues[ri].position);

    // apply boundary conditions, possibly rejecting the move
    if (test_boundary_conditions())
    {
        wrap_if_necessary();

        // apply the move
        LOG(DEBUG_LOCAL, "Residue: %d Old: (%f, %f, %f)\t", ri, Residues[ri].position.x, Residues[ri].position.y, Residues[ri].position.z);
        Residues[ri].position = Residues[ri].new_position + center_wrap_delta;
        LOG(DEBUG_LOCAL, "New: (%f, %f, %f)\n", Residues[ri].position.x, Residues[ri].position.y, Residues[ri].position.z);

        mark_cached_potentials_for_update(ri, true);
        local_move_successful = true;
    }
}

void Molecule::flex(const Vector3double raxis, const double angle, const int ri, const bool before)
{
    // calculate quaternion from angle and axis
    Quaternion q(angle, raxis);

    // apply rotation to everything before or after residue
    int start, end;

    if (before) {
        start = 0;
        end = ri;
    } else {
        start = ri + 1;
        end = residueCount;
    }

    Vector3f accumulated_difference(0,0,0);

    for (int i = start; i < end; i++) {
        Vector3f relative_position = Residues[i].position - Residues[ri].position;
        relative_position = q.rotateVector(relative_position);
        Residues[i].new_position = relative_position + Residues[ri].position;
        accumulated_difference += Residues[i].new_position - Residues[i].position;
    }

    // apply boundary conditions, possibly rejecting the move
    new_center = recalculate_center(accumulated_difference);

    if (test_boundary_conditions())
    {

        wrap_if_necessary();

        // apply the move
        for (int i = start; i < end; i++) {
            Residues[i].position = Residues[i].new_position + center_wrap_delta;
        }

        mark_cached_potentials_for_update(ri, false);
        recalculate_relative_positions();
        calculate_length();
    }
}

void Molecule::flex(gsl_rng * rng, const double rotate_step)
{
    Vector3double raxis = normalised_random_vector_d(rng);
    uint li = random_linker_index(rng);
    uint ri = random_residue_index(rng, li);
    bool before = (bool) gsl_ran_bernoulli(rng, 0.5);

    flex(raxis, rotate_step, ri, before);
}

void Molecule::make_local_moves(gsl_rng * rng, const double rotate_step, const double translate_step)
{

    local_move_successful = false;
    uint li = random_linker_index(rng);
    for (size_t i = 0; i < NUM_LOCAL_MOVES; i++) {
        //TODO: if linker too short for crankshaft, only return translate?
        //TODO: if crankshaft disabled, only return translate
        uint move = gsl_ran_bernoulli(rng, LOCAL_TRANSLATE_BIAS);

        switch (move)
        {
            case MC_LOCAL_TRANSLATE:
            {
                strcat(last_MC_move, "T");
                uint ri = random_residue_index(rng, li);
                Vector3f v = LOCAL_TRANSLATE_STEP_SCALING_FACTOR * translate_step * normalised_random_vector_f(rng);
                translate(v, ri);
                break;
            }
            case MC_LOCAL_CRANKSHAFT: // TODO: remove if crankshaft disabled
            {
                strcat(last_MC_move, "C");
                uint ri = random_residue_index_middle(rng, li);
                bool flip = (bool) gsl_ran_bernoulli(rng, 0.5);
                crankshaft(rotate_step, flip, ri);
                break;
            }
            default:
                break;
        }
    }

    // recalculate positions of residues relative to centre
    if (local_move_successful)
    {
        recalculate_relative_positions();
        calculate_length();
    }
}
#endif // FLEXIBLE_LINKS

uint Molecule::get_MC_mutation_type(gsl_rng * rng)
{
#if FLEXIBLE_LINKS
    if (linkerCount > 0)
    {
        return gsl_ran_discrete(rng, MC_discrete_table);
    }
    else
    {
        return gsl_ran_bernoulli(rng, translate_rotate_bernoulli_bias);
    }
#else // if not FLEXIBLE_LINKS
    return gsl_ran_bernoulli(rng, translate_rotate_bernoulli_bias);
#endif // FLEXIBLE_LINKS
}

void Molecule::make_MC_move(gsl_rng * rng, const double rotate_step, const double translate_step)
{
    uint mutationType = get_MC_mutation_type(rng);
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
            strcat(last_MC_move, "Translate");
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

bool Molecule::initFromPDB(const char* pdbfilename)
{
    vector<Residue> vResidues;
#if FLEXIBLE_LINKS
    vector<Link> vLinks;
#endif // FLEXIBLE_LINKS

    chainCount = 1; // there must be at least one
    strcpy(filename, pdbfilename);
    // TODO: we only store this for printing.  Put it at the end?

    ifstream input(pdbfilename);
    if (!input.good())
    {
        cout << "Failed to open file: " << pdbfilename << " ("<<filename<<")\n";
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

                LOG(DEBUG_PARSING, "residue: %s; chain: %c; sequence no.: %d; position: (%f, %f, %f); occupancy: %f\n", resName, chainID, resSeq, x, y, z, occupancy);

                Residue R;

                R.aminoAcidIndex = AminoAcidsData.getAminoAcidIndex(resName);
                R.electrostaticCharge = float(AminoAcidsData.get(R.aminoAcidIndex).electrostaticCharge);// * E_charge;
                R.vanderWaalRadius = AminoAcidsData.data[R.aminoAcidIndex].vanderWaalRadius;// * E_charge;
                if (R.vanderWaalRadius != AminoAcidsData.get(R.aminoAcidIndex).vanderWaalRadius)
                {
                    cout << R.aminoAcidIndex << " differing vdw radii" << endl;
                }

                R.position = Vector3f(x,y,z);
                R.relativePosition = Vector3f(0,0,0);
                R.chainId = chainID;
                R.resSeq = resSeq;
                //R.sasa = 1.0f;
                //R.vanderWaalRadius = float(AminoAcidsData.data[R.aminoAcidIndex].vanderWaalRadius);
                R.moleculeID = index;
                vResidues.push_back(R);

#if FLEXIBLE_LINKS
                Link L;
                // occupancy set to 1 on a residue indicates that the following link is flexible.
                L.flexible = bool(occupancy);
                vLinks.push_back(L);
#endif // FLEXIBLE_LINKS
            }
        }
        else if (strncmp(line, "TER", 3) == 0) // Optional,Mandatory if ATOM records exist.
        {
            chainCount++;
#if FLEXIBLE_LINKS
            // The last link isn't a real link
            vLinks.back().dummy = true;
            // Set it to non-flexible regardless of the occupancy of the last residue
            vLinks.back().flexible = false;
#endif // FLEXIBLE_LINKS
        }
        else if (strncmp(line, "END", 3) == 0)
        {
            reached_end_marker = true;
        }
    }
    input.close();

    // copy the resides and links read into the arrays.
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

    calculateVolume();

    // TODO: eliminate all the arrays; we can just keep using vectors.
#if FLEXIBLE_LINKS
    linkCount = vLinks.size() - 1;
    Links = new Link[linkCount];
    vector<Segment> vSegments;
    vector<int> vLinkers;

    for(size_t l = 0; l < linkCount; l++)
    {
        Link L = vLinks[l];
        memcpy (&Links[l], &L, sizeof(L));

        // End last segment
        if (l && L.flexible != vLinks[l-1].flexible || L.dummy)
        {
            vSegments.back().end = l;
        }

        // Start new segment
        if (!l || L.flexible != vLinks[l-1].flexible || vLinks[l-1].dummy)
        {
            Segment S;
            S.start = l;
            S.flexible = L.flexible;
            vSegments.push_back(S);
        }

        // End last segment at end of molecule
        if (l == linkCount - 1)
        {
            vSegments.back().end = linkCount;
            vSegments.back().size = vSegments.back().end - vSegments.back().start + 1;
        }
    }

    segmentCount = vSegments.size();
    Segments = new Segment[segmentCount];

    for (size_t s = 0; s < segmentCount; s++)
    {
        Segment S = vSegments[s];
        memcpy(&Segments[s], &S, sizeof(S));

        if (S.flexible) {
            vLinkers.push_back(s);
            is_flexible = true;
        }
    }

    linkerCount = vLinkers.size();
    Linkers = new int[linkerCount];

    for (size_t l = 0; l < linkerCount; l++)
    {
        Linkers[l] = vLinkers[l];
    }

#endif // FLEXIBLE_LINKS

    return true;
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

uint Molecule::random_residue_index(gsl_rng * rng)
{
    return (int) gsl_rng_uniform_int(rng, residueCount);
}

#if FLEXIBLE_LINKS

// TODO: add comments

Potential Molecule::E(bool include_LJ_and_DH)
{
    Potential potential;

#define iRes Residues[i]
#define jRes Residues[j]
#define iSeg Segments[si]
#define jSeg Segments[sj]

    // We may be using the GPU to calculate the internal LJ and DH, or not calculating it at all
    if (is_flexible && include_LJ_and_DH)
    {
        // LJ and DH between all residue pairs within molecule
        if (update_LJ_and_DH)
        {
            // TODO: we don't need this anymore. Eliminate it.
            potential.reset_LJ_subtotal();
            potential.reset_DH_subtotal();

            for (size_t i = 0; i < residueCount; i++)
            {
                for (size_t j = i + 1; j < residueCount; j++) {
                    double r(iRes.distance(jRes, bounding_value) + EPS);
                    // apply intra-linker condition to eliminate overlap and closeness
                    if (j - i >= 4) // TODO: merge into for
                    {
                        /* Calculate LJ-type potential for each residue pair; increment molecule total. */
                        potential.increment_LJ_subtotal(calculate_LJ(iRes, jRes, r, AminoAcidsData));
                        /* Calculate electrostatic potential for each residue pair; increment molecule total. */
                        potential.increment_DH_subtotal(calculate_DH(iRes, jRes, r));
                    }
                }
            }

            /* Cache new values on the molecule */
            LJ = potential.LJ_subtotal;
            DH = potential.DH_subtotal;

            update_LJ_and_DH = false;
        }

        /* Add molecule totals to potential totals */
        /* TODO: what impact do side calculations like this have on kahan sum accuracy? */

        potential.increment_LJ(LJ);
        potential.increment_DH(DH);
    }

    for (size_t si = 0; si < segmentCount; si++)
    {
        // Flexible linker potentials
        if (iSeg.flexible)
        {
            for (size_t i = iSeg.start; i <= iSeg.end; i++)
            {
                // Pseudo-bond
                if (i < iSeg.end)
                {
                    potential.increment_bond(calculate_bond(iRes, Links[i], Residues[i+1], bounding_value));
                }

                if (i > 0 && !Links[i-1].dummy)
                {

                    // Pseudo-angle
                    if (i < linkCount)
                    {
                        potential.increment_angle(calculate_angle(Residues[i-1], iRes, Residues[i+1]));
                    }

                    // Pseudo-torsion
                    if (i < linkCount - 1 && !Links[i+1].dummy)
                    {
                        potential.increment_torsion(calculate_torsion(Residues[i-1], iRes, Links[i], Residues[i+1], Residues[i+2], torsions));
                    }
                }
            }
        }
    }

    return potential;
}

uint Molecule::random_linker_index(gsl_rng * rng)
{
    return (int) gsl_rng_uniform_int(rng, linkerCount);
}

uint Molecule::random_residue_index(gsl_rng * rng, int li)
{
    return Segments[Linkers[li]].random_residue_index(rng);
}

uint Molecule::random_residue_index_middle(gsl_rng * rng, int li)
{
    return Segments[Linkers[li]].random_residue_index_middle(rng);
}


#endif // FLEXIBLE_LINKS