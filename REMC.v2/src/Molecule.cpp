#include "Molecule.h"

using namespace std;

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
#ifdef FLEXIBLE_LINKS
    linkCount = 0;
    segmentCount = 0;
    LJ = 0.0f;
    DH = 0.0f;
    update_LJ_and_DH = true;
#endif
}

Molecule::~Molecule()
{
    if(hasFilename)
        delete [] filename;
};

void Molecule::init_amino_acid_data(AminoAcids &a)
{
    AminoAcidsData = a;
#ifdef FLEXIBLE_LINKS
//     torsions.setAminoAcidLookupData(AminoAcidsData);
    torsions.loadData("data/torsional_pair_potentials", a);
#endif
}

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
    position = m.position;
    rotation = m.rotation;
    index = m.index;
    volume = m.volume;
    moleculeRoleIdentifier = m.moleculeRoleIdentifier;
    amIACrowder = m.amIACrowder;
    hasFilename = false;
#ifdef FLEXIBLE_LINKS
    torsions = m.torsions;
    linkCount = m.linkCount;
    segmentCount = m.segmentCount;
    LJ = m.LJ;
    DH = m.DH;
    update_LJ_and_DH = m.update_LJ_and_DH;
#endif
}

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
    position = m.position;
    rotation = m.rotation;
    volume = m.volume;

    index = m.index;
    moleculeRoleIdentifier = m.moleculeRoleIdentifier;
    amIACrowder = m.amIACrowder;
#ifdef FLEXIBLE_LINKS
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
    center = m->center;
    position = m->position;
    rotation = m->rotation;
}

void Molecule::undoStateChange(const Molecule* m)
{
    memcpy(Residues,m->Residues,sizeof(Residue)*residueCount);
    center = m->center;
    position = m->position;
    rotation = m->rotation;
}

void Molecule::reserveResidueSpace(int size)
{
    /*try
    {
    	delete [] Residues;
    } catch (char * err) { cout << "Molecule already sized, deleted and resized" <<endl; }*/
    Residues = new Residue[size];
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

void Molecule::deleteResidues()
{
    delete [] Residues;
}

bool Molecule::translate(const Vector3f v)
{
    //center = center+v;
    //position = position+v;

    center.x += v.x;
    center.y += v.y;
    center.z += v.z;
    position.x += v.x;
    position.y += v.y;
    position.z += v.z;

    for (size_t i=0; i<residueCount; i++)
    {
        //Residues[i].position = Residues[i].position + v;

        // faster. no oo overhead
        Residues[i].position.x += v.x;
        Residues[i].position.y += v.y;
        Residues[i].position.z += v.z;

        //even faster, use sse?

    }
    return true;
}

Vector3f Molecule::calculateCenter()
{
    Vector3f accumulator(0,0,0);
    for (size_t i=0; i<residueCount; i++)
    {
        accumulator = accumulator + Residues[i].position;
    }
    center = accumulator / residueCount;
    return center;
}

// void Molecule::saveAsPDB(const char* filename)
// {
//     FILE * output;
//     output = fopen (filename,"w");
//     fprintf(output,"REMARK %s \n",filename);
//     fprintf(output,"REMARK Translation relative to input Q(w,x,y,z): %f %f %f %f\n",rotation.w,rotation.x,rotation.y,rotation.z);
//     fprintf(output,"REMARK Rotation relative to input +P(x,y,z): %f %f %f\n",center.x,center.y,center.z);
//
//     size_t i=0;
//     int itemcount = 0;
//     while (i<residueCount)
//     {
//         itemcount++;
//         fprintf(output,"ATOM  %5d %4s%C%3s %C%4d%C  %8.3f%8.3f%8.3f%6.2f%6.2f\n",itemcount,"CA",' ',AminoAcidsData.get(Residues[i].aminoAcidIndex).getSNAME(),Residues[i].chainId,Residues[i].resSeq,' ',Residues[i].position.x,Residues[i].position.y,Residues[i].position.z,1.0f,1.0f);
//         i++;
//     }
//
//     fprintf(output,"END \n");
//     fflush(output);
//     fclose(output);
//
// }

void Molecule::setPosition(Vector3f v)
{
    center = v;
    position = v;
    for (size_t i=0; i<residueCount; i++)
    {
        // more efficient to do it this way than with the overloaded +
        Residues[i].position.x = Residues[i].relativePosition.x + center.x;
        Residues[i].position.y = Residues[i].relativePosition.y + center.y;
        Residues[i].position.z = Residues[i].relativePosition.z + center.z;
    }
}

bool Molecule::rotateQ(const Vector3double Raxis, const double angle)
{
    double sina = sin(angle/2.0);
    double cosa = cos(angle/2.0);

    Quaternion q(cosa,sina*Raxis.x,sina*Raxis.y,sina*Raxis.z);

    // track the global change from initial conditions.
    q.normalize();
    rotation = q * rotation;

    for (size_t i=0; i<residueCount; i++)
    {
        Residues[i].relativePosition = q.rotateVector(Residues[i].relativePosition);
        Residues[i].position.x = Residues[i].relativePosition.x + center.x;
        Residues[i].position.y = Residues[i].relativePosition.y + center.y;
        Residues[i].position.z = Residues[i].relativePosition.z + center.z;
    }

    return true;
}

void Molecule::setRotation(Quaternion q)
{
    q.normalize();
    rotation = q * rotation;
    for (size_t i=0; i<residueCount; i++)
    {
        Residues[i].relativePosition = q.rotateVector(Residues[i].relativePosition);
        Residues[i].position.x = Residues[i].relativePosition.x + center.x;
        Residues[i].position.y = Residues[i].relativePosition.y + center.y;
        Residues[i].position.z = Residues[i].relativePosition.z + center.z;
    }
}

bool Molecule::initFromPDB(const char* pdbfilename)
{
    vector<Residue> vResidues;
#ifdef FLEXIBLE_LINKS
    vector<Link> vLinks;
#endif

    chainCount = 1; // there must be at least one
    hasFilename = true;
    filename = new char[256];
    strcpy(filename, pdbfilename);

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

#ifdef FLEXIBLE_LINKS
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
#ifdef FLEXIBLE_LINKS
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

#ifdef FLEXIBLE_LINKS
    linkCount = vLinks.size() - 1;
    Links = new Link[linkCount];
    vector<Segment> vSegments;

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
        }
    }

    segmentCount = vSegments.size();
    Segments = new Segment[segmentCount];

    for (size_t s = 0; s < segmentCount; s++)
    {
        Segment S = vSegments[s];
        memcpy (&Segments[s], &S, sizeof(S));
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

#ifdef FLEXIBLE_LINKS

// TODO: add comments

Potential Molecule::E()
{
    Potential potential;
    double c_lj(0.0f);
    double c_dh(0.0f);

    // LJ and DH between segments within molecule
    if (update_LJ_and_DH)
    {
        LJ = 0.0f;
        DH = 0.0f;
        for (size_t si = 0; si < segmentCount; si++)
        {
            for (size_t sj = si + 1; sj < segmentCount; sj++)
            {
                for (size_t i = Segments[si].start; i <= Segments[si].end; i++)
                {
                    for (size_t j = Segments[sj].start; j <= Segments[si].end; j++) {
                        // TODO: do we need to wrap around the bounding box for this calculation?
                        double r((Residues[i].position - Residues[j].position).magnitude() + EPS);
                        /* Calculate LJ-type potential for each residue pair; increment molecule total. */
                        potential.increment_LJ(Residues[i], Residues[j], r, AminoAcidsData, LJ, c_lj);
                        /* Calculate electrostatic potential for each residue pair; increment molecule total. */
                        potential.increment_DH(Residues[i], Residues[j], r, DH, c_dh);
                    }
                }
            }
        }

        update_LJ_and_DH = false;
    }

    /* Add molecule totals to potential totals */
    /* TODO: what impact do side calculations like this have on kahan sum accuracy? */
    potential.increment_LJ(LJ);
    potential.increment_DH(DH);

    for (size_t si = 0; si < segmentCount; si++)
    {
        // Flexible linker potentials
        if (Segments[si].flexible)
        {
            for (size_t i = Segments[si].start; i <= Segments[si].end; i++)
            {
                // LJ and DH within linker
                if (Segments[si].update_LJ_and_DH)
                {
                    c_lj = 0.0f;
                    c_dh = 0.0f;
                    Segments[si].LJ = 0.0f;
                    Segments[si].DH = 0.0f;
                    for (size_t j = i + 1; j <= Segments[si].end; j++)
                    {
                        // TODO: do we need to wrap around the bounding box for this calculation?
                        double r((Residues[i].position - Residues[j].position).magnitude() + EPS);
                        /* Calculate LJ-type potential for each residue pair; increment segment total. */
                        potential.increment_LJ(Residues[i], Residues[j], r, AminoAcidsData, Segments[si].LJ, c_lj);
                        /* Calculate electrostatic potential if residues separated by more than 3 residues (kim2008 p. 1429); increment segment total. */
                        if (j - i >= 4)
                        {
                            potential.increment_DH(Residues[i], Residues[j], r, Segments[si].DH, c_dh);
                        }
                        Segments[si].update_LJ_and_DH = false;
                    }
                }

                /* Add segment totals to potential totals */
                potential.increment_LJ(Segments[si].LJ);
                potential.increment_DH(Segments[si].DH);

                // Pseudo-bond
                if (i < Segments[si].end)
                {
                    potential.increment_bond(Residues[i], Links[i], Residues[i+1]);
                }

                if (i > 0 && !Links[i-1].dummy)
                {

                    // Pseudo-angle
                    if (i < linkCount)
                    {
                        potential.increment_angle(Residues[i-1], Residues[i], Residues[i+1]);
                    }

                    // Pseudo-torsion
                    if (i < linkCount - 1 && !Links[i+1].dummy)
                    {
                        potential.increment_torsion(Residues[i-1], Residues[i], Links[i], Residues[i+1], Residues[i+2], torsions);
                    }
                }
            }
        }
    }

    return potential;
}

#endif
