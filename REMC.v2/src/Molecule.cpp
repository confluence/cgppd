#include "Molecule.h"

using namespace std;

// TODO: we probably need a default constructor for dynamic arrays. :/
Molecule::Molecule()
{
    translationalStep = INITIAL_TRANSLATIONAL_STEP;
    rotationalStep = INITIAL_ROTATIONAL_STEP;
    residueCount = 0;
    moleculeRoleIdentifier = 0.0f;
    index = -2;
    rotation = Quaternion(1.0f,0,0,0);
    volume = 0.0f;
    hasFilename = false;
    amIACrowder = false;
#if FLEXIBLE_LINKS
    linkCount = 0;
    segmentCount = 0;
    LJ = 0.0f;
    DH = 0.0f;
    update_LJ_and_DH = true;
#endif
}

Molecule::Molecule(const char* pdbfilename) // constructor for molecule loaded from file; not final parameter list
{
}

Molecule::Molecule(const int size) // constructor for blank molecule of a particular size (for saving before MC move); not final parameter list
{
}

Molecule::~Molecule()
{
    if(hasFilename)
        delete [] filename;
    // TODO: shouldn't we delete residues and other dynamically allocated arays here?  Use a flag to see if they have been allocated?
    // yes, we should
    if(residueCount > 0)
    {
        delete [] Residues;
    }
};

void Molecule::init_amino_acid_data(AminoAcids &a)
{
    AminoAcidsData = a;
#if FLEXIBLE_LINKS
//     torsions.setAminoAcidLookupData(AminoAcidsData);
    torsions.loadData("data/torsional_pair_potentials", a);
#endif
}

// TODO: are we using this?
Molecule::Molecule(const Molecule& m)
{
    translationalStep = m.translationalStep;
    rotationalStep = m.rotationalStep;

    if(residueCount > 0 && residueCount != m.residueCount)
    {
        delete [] Residues;
        Residues = new Residue[m.residueCount];
    }
    else if (residueCount == 0)
    {
        Residues = new Residue[m.residueCount];
    }

    residueCount = m.residueCount;
    memcpy(Residues,m.Residues,sizeof(Residue)*residueCount);

    AminoAcidsData = m.AminoAcidsData;
    center = m.center;
    rotation = m.rotation;
    index = m.index;
    volume = m.volume;
    moleculeRoleIdentifier = m.moleculeRoleIdentifier;
    amIACrowder = m.amIACrowder;
    hasFilename = false;
#if FLEXIBLE_LINKS
    torsions = m.torsions;
    linkCount = m.linkCount;
    segmentCount = m.segmentCount;
    LJ = m.LJ;
    DH = m.DH;
    update_LJ_and_DH = m.update_LJ_and_DH;
#endif
}

// TODO: are we using this?
void Molecule::copy(const Molecule& m)
{
    translationalStep = m.translationalStep;
    rotationalStep = m.rotationalStep;

    if(residueCount > 0 && residueCount != m.residueCount)
    {
        delete [] Residues;
        Residues = new Residue[m.residueCount];
    }
    else if (residueCount == 0)
    {
        try
        {
            delete [] Residues;
        } catch (char * err) {}
        Residues = new Residue[m.residueCount];
    }

    residueCount = m.residueCount;
    memcpy(Residues,m.Residues,sizeof(Residue)*residueCount);

    AminoAcidsData = m.AminoAcidsData;
    center = m.center;
    rotation = m.rotation;
    volume = m.volume;

    index = m.index;
    moleculeRoleIdentifier = m.moleculeRoleIdentifier;
    amIACrowder = m.amIACrowder;
#if FLEXIBLE_LINKS
    torsions = m.torsions;
    linkCount = m.linkCount;
    segmentCount = m.segmentCount;
    LJ = m.LJ;
    DH = m.DH;
    update_LJ_and_DH = m.update_LJ_and_DH;
#endif
}

void Molecule::saveBeforeStateChange(const Molecule* m)
{
    memcpy(Residues,m->Residues,sizeof(Residue)*m->residueCount);
    residueCount = m->residueCount; // make sure we copy the right number of residues back, and also delete the residues in destructor
    center = m->center;
    rotation = m->rotation;
}

void Molecule::undoStateChange(const Molecule* m)
{
    // TODO: how did this ever work?! If residueCount wasn't set in the save, it would have always been zero, and no bytes would have been copied!
    memcpy(Residues,m->Residues,sizeof(Residue)*residueCount);
    center = m->center;
    rotation = m->rotation;
}

void Molecule::reserveResidueSpace(int size)
{
    /*try
    {
    	delete [] Residues;
    } catch (char * err) { cout << "Molecule already sized, deleted and resized" <<endl; }*/
    Residues = new Residue[size]; // TODO make sure array is deleted in destructor
}

void Molecule::setMoleculeRoleIdentifier(float moleculeRoleIdentifier)
{
    this->moleculeRoleIdentifier = moleculeRoleIdentifier;
    if(moleculeRoleIdentifier == CROWDER_IDENTIFIER)
    {
        amIACrowder = true;

#if REPULSIVE_CROWDING
        for (int i =0; i<residueCount; i++)
        {
            //Residues[i].isCrowder = true;
            Residues[i].aminoAcidIndex = int(CROWDER_IDENTIFIER);
        }
#endif
    }
    else if (moleculeRoleIdentifier == PADDER_IDENTIFIER)
    {
        for (int i =0; i<residueCount; i++)
        {
            //Residues[i].isCrowder = false;
            Residues[i].aminoAcidIndex = int(PADDER_IDENTIFIER);
        }
    }
}

void Molecule::translate(const Vector3f v)
{
    //center = center+v;

    center.x += v.x;
    center.y += v.y;
    center.z += v.z;

    for (size_t i=0; i<residueCount; i++)
    {
        //Residues[i].position = Residues[i].position + v;

        // faster. no oo overhead
        Residues[i].position.x += v.x;
        Residues[i].position.y += v.y;
        Residues[i].position.z += v.z;

        //even faster, use sse?

    }
}

void Molecule::recalculate_relative_positions()
{
    for (size_t i=0; i<residueCount; i++)
    {
        Residues[i].relativePosition = center + Residues[i].position;
    }
}

void Molecule::recalculate_center()
{
    Vector3f accumulator(0,0,0);
    for (size_t i=0; i<residueCount; i++)
    {
        accumulator = accumulator + Residues[i].position;
    }
    center = accumulator / residueCount;
    recalculate_relative_positions();
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
    rotation = q * rotation;
    for (size_t i=0; i<residueCount; i++)
    {
//         Residues[i].relativePosition = q.rotateVector(Residues[i].relativePosition);
//         Residues[i].position.x = Residues[i].relativePosition.x + center.x;
//         Residues[i].position.y = Residues[i].relativePosition.y + center.y;
//         Residues[i].position.z = Residues[i].relativePosition.z + center.z;
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

void Molecule::rotate(gsl_rng * rng)
{
    rotate(normalised_random_vector_d(rng), INITIAL_ROTATIONAL_STEP);
}

// TODO: add boundary conditions to everything?
void Molecule::translate(gsl_rng * rng, const float bounding_value)
{
    Vector3f v = INITIAL_TRANSLATIONAL_STEP * normalised_random_vector_f(rng);

#if BOUNDING_METHOD == BOUNDING_SPHERE
    if ((center + v).sumSquares() < bounding_value * bounding_value)
    {
        translate(v);
    }
    else
    {
    }
#elif BOUNDING_METHOD == PERIODIC_BOUNDARY
    Vector3f new_pos = center + v;
    new_pos.x = fmod(new_pos.x + bounding_value, bounding_value);
    new_pos.y = fmod(new_pos.y + bounding_value, bounding_value);
    new_pos.z = fmod(new_pos.z + bounding_value, bounding_value);
    setPosition(new_pos);
#endif
}

#if FLEXIBLE_LINKS
void Molecule::recalculate_center(Vector3f difference)
{
    // TODO: doing this after each local move is expensive -- can we get away with calling it once from Replica at the end of all local moves?
    // YES, because none of the local moves care about the centre.
    center = ((center * residueCount) + difference) / residueCount;
    recalculate_relative_positions();
}

void Molecule::mark_cached_potentials_for_update(const int ri)
{
    update_LJ_and_DH = true;

    for (int i = 0; i < segmentCount; i++) {
        if (ri >= Segments[i].start && ri <= Segments[i].end) {
            Segments[i].update_LJ_and_DH = true;
        }
    }

    // Not checking segment boundaries here, because they are checked when linker potentials are calculated

    if (ri > 1) {
        Links[ri - 2].update_e_torsion = true;
    }

    if (ri > 0) {
        Links[ri - 1].update_e_bond = true;
        Residues[ri - 1].update_e_angle = true;
        Links[ri - 1].update_e_torsion = true;
    }

    Links[ri].update_e_bond = true;
    Residues[ri].update_e_angle = true;
    Links[ri].update_e_torsion = true;

    if (ri < residueCount - 1) {
        Residues[ri + 1].update_e_angle = true;
        Links[ri + 1].update_e_torsion = true;
    }
}

//TODO TODO TODO boundary conditions -- see how they're done for other MC moves

void Molecule::translate(Vector3f v, const int ri)
{
    // translate one residue
    Residues[ri].position.x += v.x;
    Residues[ri].position.y += v.y;
    Residues[ri].position.z += v.z;

    // recalculate centre (just one residue changed)
    recalculate_center(v);

    // mark potential values for update
    mark_cached_potentials_for_update(ri);
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
    Vector3f old_position = Residues[ri].position;
    Vector3f relative_position = Residues[ri].position - Residues[ri - 1].position;
    relative_position = q.rotateVector(relative_position);
    Residues[ri].position = relative_position + Residues[ri - 1].position;

    // recalculate centre (just one residue changed)
    recalculate_center(Residues[ri].position - old_position);

    // mark potential values for update
    mark_cached_potentials_for_update(ri);
}

void Molecule::rotate_domain(const Vector3double raxis, const double angle, const int ri, const bool before)
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

    for (int i = start; i < end; i++) {
        Vector3f relative_position = Residues[i].position - Residues[ri].position;
        relative_position = q.rotateVector(relative_position);
        Residues[i].position = relative_position + Residues[ri].position;
    }

    // recalculate centre (everything)
    recalculate_center();

    // mark potential values for update
    mark_cached_potentials_for_update(ri);
}

void Molecule::rotate_domain(gsl_rng * rng)
{
    Vector3double raxis = normalised_random_vector_d(rng);
    uint li = random_linker_index(rng);
    uint ri = random_residue_index(rng, li);
    bool before = (bool) gsl_ran_bernoulli(rng, 0.5);

    rotate_domain(raxis, angle, ri, before);
    // TODO TODO TODO recalculate the centre afterwards
}

void Molecule::make_local_moves(gsl_rng * rng)
{
    uint li = random_linker_index(rng);
    for (size_t i = 0; i <= NUM_LOCAL_MOVES; i++) {
        uint ri = random_residue_index(rng, li);
        //TODO: if crankshaft disabled, only return translate
        uint move = gsl_ran_bernoulli(rng, LOCAL_TRANSLATE_BIAS);

        switch (move)
        {
            case MC_LOCAL_TRANSLATE:
            {
                Vector3f v = INITIAL_TRANSLATIONAL_STEP * normalised_random_vector_f(rng);
                translate(v, ri);
                break;
            }
            case MC_LOCAL_CRANKSHAFT: // TODO: remove if crankshaft disabled
            {
                bool flip = (bool) gsl_ran_bernoulli(rng, 0.5);
                // TODO: we actually need to select a residue that is not at the end of the linker!
                crankshaft(INITIAL_ROTATIONAL_STEP, flip, ri);
                break;
            }
            default:
                break;
        }
    }
    // TODO TODO TODO recalculate the centre afterwards
}
#endif // FLEXIBLE_LINKS

bool Molecule::initFromPDB(const char* pdbfilename)
{
    vector<Residue> vResidues;
#if FLEXIBLE_LINKS
    vector<Link> vLinks;
#endif

    chainCount = 1; // there must be at least one
    hasFilename = true;
    filename = new char[256];
    strcpy(filename, pdbfilename);
    // TODO: we only store this for printing.  Put it at the end?

    ifstream input(pdbfilename);
    if (!input.good())
    {
        cout << "Failed to open file: " << pdbfilename << " ("<<filename<<")\n";
        exit(0);
    }

#ifdef VERY_VERBOSE
    //cout << "\tParsing: " << pdbfilename << "\n";
#endif

    // for eachline in the file
    char line[512] = {' '};
    bool reached_end_marker = false;
    while (!reached_end_marker && !input.eof())
    {
        input.getline(line,512);

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
#endif
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
#endif
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

    calculateSASA();
    calculateVolume();

#if FLEXIBLE_LINKS
    linkCount = vLinks.size() - 1;
    Links = new Link[linkCount];
    vector<Segment> vSegments;
    vector<Segment*> vLinkers;

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
            vLinkers.push_back(&Segments[s]);
        }
    }

    linkerCount = vLinkers.size();
    Linkers = new Segment*[linkerCount];

    for (size_t l = 0; l < linkerCount; l++)
    {
        memcpy(&Linkers[l], &vLinkers[l], sizeof(Segment*));
    }

#endif

    return true;
}

float Molecule::calculateSASA()
{
    // todo, not required for now
    return 1.0f;
}

float Molecule::getVolume()
{
    return volume;
}

float Molecule::calculateVolume()
{
    // simple volume, todo: surface enclosed calc which is nontrivial
    volume = 0.0f;
    for (size_t i=0; i<residueCount; i++)
    {
        float r(AminoAcidsData.get(Residues[i].aminoAcidIndex).vanderWaalRadius);
        volume += 4.0/3.0*r*r*r*M_PIl ;  // 4/3 * pi * r^3
    }
    return volume;
}

uint Molecule::random_residue_index(gsl_rng * rng)
{
    return (int) gsl_rng_uniform_int(rng, residueCount);
}

#if FLEXIBLE_LINKS

// TODO: add comments

Potential Molecule::E(const float bounding_value)
{
    Potential potential;

#define iRes Residues[i]
#define jRes Residues[j]
#define iSeg Segments[si]
#define jSeg Segments[sj]

    // LJ and DH between segments within molecule
    if (update_LJ_and_DH)
    {
        potential.reset_LJ_subtotal();
        potential.reset_DH_subtotal();

        for (size_t si = 0; si < segmentCount; si++)
        {
            for (size_t sj = si + 1; sj < segmentCount; sj++)
            {
                for (size_t i = iSeg.start; i <= iSeg.end; i++)
                {
                    for (size_t j = jSeg.start; j <= iSeg.end; j++) {
                        double r(iRes.distance(jRes, bounding_value) + EPS);
                        /* Calculate LJ-type potential for each residue pair; increment molecule total. */
                        potential.increment_LJ_subtotal(calculate_LJ(iRes, jRes, r, AminoAcidsData));
                        /* Calculate electrostatic potential for each residue pair; increment molecule total. */
                        potential.increment_DH_subtotal(calculate_DH(iRes, jRes, r));
                    }
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

    for (size_t si = 0; si < segmentCount; si++)
    {
        // Flexible linker potentials
        if (iSeg.flexible)
        {
            for (size_t i = iSeg.start; i <= iSeg.end; i++)
            {
                // LJ and DH within linker
                if (iSeg.update_LJ_and_DH)
                {
                    potential.reset_LJ_subtotal();
                    potential.reset_DH_subtotal();

                    for (size_t j = i + 1; j <= iSeg.end; j++)
                    {
                        double r(iRes.distance(jRes, bounding_value) + EPS);
                        /* Calculate LJ-type potential for each residue pair; increment segment total. */
                        potential.increment_LJ_subtotal(calculate_LJ(iRes, jRes, r, AminoAcidsData));
                        /* Calculate electrostatic potential if residues separated by more than 3 residues (kim2008 p. 1429); increment segment total. */
                        if (j - i >= 4)
                        {
                            potential.increment_DH_subtotal(calculate_DH(iRes, jRes, r));
                        }
                    }

                    /* Cache new values on the segment */
                    iSeg.LJ = potential.LJ_subtotal;
                    iSeg.DH = potential.DH_subtotal;

                    iSeg.update_LJ_and_DH = false;
                }

                /* Add segment totals to potential totals */
                potential.increment_LJ(iSeg.LJ);
                potential.increment_DH(iSeg.DH);

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
    return Linkers[li]->random_residue_index(rng);
}


#endif // FLEXIBLE_LINKS
