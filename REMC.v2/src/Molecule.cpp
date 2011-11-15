#include "Molecule.h"

using namespace std;

Molecule::Molecule()
{
    translationalStep = INITIAL_TRANSLATIONAL_STEP;
    rotationalStep = INITIAL_ROTATIONAL_STEP;
    residueCount = 0;
    linkCount = 0;
    moleculeRoleIdentifier = 0.0f;
    index = -2;
    rotation = Quaternion(1.0f,0,0,0);
    volume = 0.0f;
    hasFilename = false;
    amIACrowder = false;
}

Molecule::~Molecule()
{
    if(hasFilename)
        delete [] filename;
};

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
    //temperature = m.temperature;
    //label = m.label;
    moleculeRoleIdentifier = m.moleculeRoleIdentifier;
    amIACrowder = m.amIACrowder;
    hasFilename = false;
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

    /*index = m.index;
    temperature = m.temperature;
    label = m.label;*/
    moleculeRoleIdentifier = m.moleculeRoleIdentifier;
    amIACrowder = m.amIACrowder;
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

void Molecule::saveAsPDB(const char* filename)
{
    FILE * output;
    output = fopen (filename,"w");
    fprintf(output,"REMARK %s \n",filename);
    fprintf(output,"REMARK Translation relative to input Q(w,x,y,z): %f %f %f %f\n",rotation.w,rotation.x,rotation.y,rotation.z);
    fprintf(output,"REMARK Rotation relative to input +P(x,y,z): %f %f %f\n",center.x,center.y,center.z);

    size_t i=0;
    int itemcount = 0;
    while (i<residueCount)
    {
        itemcount++;
        fprintf(output,"ATOM  %5d %4s%C%3s %C%4d%C  %8.3f%8.3f%8.3f%6.2f%6.2f\n",itemcount,"CA",' ',AminoAcidsData.get(Residues[i].aminoAcidIndex).getSNAME(),Residues[i].chainId,Residues[i].resSeq,' ',Residues[i].position.x,Residues[i].position.y,Residues[i].position.z,1.0f,1.0f);
        i++;
    }

    fprintf(output,"END \n");
    fflush(output);
    fclose(output);

}

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
    strcpy(filename,pdbfilename);

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

            char field[9] = {' '};

            strncpy(field,line+12,4);
            sscanf (field,"%s",name);

            // push the atom onto the vector of our structure if it is a CA
            if (strcmp(name, "CA") == 0) // it is a CA atom. center of our bead
            {
                char resName[4] = {' '};
                char chainID = ' ';
                int resSeq(0);
                float x(0);
                float y(0);
                float z(0);
                float occupancy;

                strncpy(field,line+17,3);
                sscanf (field,"%s",resName);

                strncpy(field,line+21,1);
                sscanf (field,"%c",&chainID);

                strncpy(field,line+22,4);
                sscanf (field,"%d",&resSeq);

                strncpy (field,line+30,8);
                sscanf (field,"%f",&x);

                strncpy (field,line+38,8);
                sscanf (field,"%f",&y);

                strncpy (field,line+46,8);
                sscanf (field,"%f",&z);

                strncpy (field,line+54,6);
                sscanf (field,"%f",&occupancy);

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
            vLinks.back().terminal = true;
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

    for(size_t l = 0; l < linkCount; l++)
    {
        Link L = vLinks[l];
        memcpy (&Links[l], &L, sizeof(L));
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
// TODO: do all these calculations in one loop?

float Molecule::E(AminoAcids *a)
{
    double epotential = 0.0f;
    double LJAccumulator = 0.0f;
    double DHAccumulator = 0.0f;
    double DH_constant_component =  DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f ;

    // TODO: apply this to residues joined by flexible links only
    // TODO: revisit caching indices of flexible segments somehow
    for (size_t i = 0; i < linkCount; i++)
    {
        for (size_t j = i + 1; j < residueCount; j++)
        {
            // TODO: do we need to wrap around the bounding box for this calculation?
            double r((Residues[i].position - Residues[j].position).magnitude() + EPS);

            /* Calculate electrostatic potential if residues separated by more than 3 residues (kim2008 p. 1429). */
            if (j - i >= 4)
            {
                DHAccumulator += (Residues[i].electrostaticCharge * Residues[j].electrostaticCharge * expf(-r / Xi) / r);
            }

            /* Calculate LJ-type potential for each residue pair. */
            float Eij(lambda * (a->LJpotentials[Residues[i].aminoAcidIndex][Residues[j].aminoAcidIndex] - e0));
            float sigmaij(0.5f * (Residues[i].vanderWaalRadius + Residues[j].vanderWaalRadius));
            double LJtmp(powf(sigmaij / r, 6.0f)); //sigT*sigT*sigT*sigT*sigT*sigT;
            double LJ(-4.0f * Eij * LJtmp * (LJtmp - 1.0f));
            if (Eij > 0.0f && r < sigmaij * 1.122462048309372981433533049679f)  // attractive pairs
            {
                LJ = -LJ + 2.0f * Eij;
            }
            LJAccumulator += LJ;
        }
    }
}

float Molecule::E_bond()
{
    float e_bond = 0.0;
    for (size_t i = 0; i < linkCount; i++)
    {
        /* Bond length calculated for every flexible link. */
        if (Links[i].flexible) {
            if (Links[i].update_e_bond)
            {
                // eqn 9: kim2008
                float rmag = float((Residues[i+1].position - Residues[i].position).magnitude());  // reduces the number of sqrts by 1
                Links[i].pseudo_bond = rmag;
                Links[i].e_bond = (rmag - R0 * Angstrom) * (rmag - R0 * Angstrom);
                Links[i].update_e_bond = false;
            }
            e_bond += Links[i].e_bond;
        }
    }
    return 0.5 * K_spring * e_bond;
}

float Molecule::E_angle()
{
    double e_angle = 1.0;
    Vector3f ab;
    Vector3f cb;
    float theta;

    for (size_t i = 1; i < linkCount; i++)
    {
        /* Angle calculated for residue before every flexible link and
         * first residue after flexible link, unless we are at the beginning or
         * end of a chain or the molecule. */
        if (Links[i].flexible && !Links[i-1].terminal || Links[i-1].flexible && !Links[i].terminal) {
            if (Links[i].update_e_angle)
            {
                ab = Residues[i-1].position - Residues[i].position;
                cb = Residues[i+1].position - Residues[i].position;
                theta = ab.angle(cb);
                Links[i].pseudo_angle = theta;
                // eqn 10: kim2008
                Links[i].e_angle = exp(-GammaAngle * (KAlpha * (theta - ThetaAlpha) * (theta - ThetaAlpha) + EpsilonAlpha)) +
                            exp(-GammaAngle * (KBeta *(theta - ThetaBeta) *(theta - ThetaBeta)));
                Links[i].update_e_angle = false;
            }
            e_angle *= Links[i].e_angle;
        }
    }
    return -GammaAngleReciprocal * log(e_angle);
}

float Molecule::E_torsion()
{
    float e_torsion = 0.0;

    Vector3f b1;
    Vector3f b2;
    Vector3f b3;
    Vector3f b2xb3;
    float phi;

    int r1;
    int r2;

    for (size_t i = 1; i < linkCount - 1; i++)
    {
        /* Torsion calculated for each flexible link except first, last, and
         * first/last of chain. */
        if (Links[i].flexible && !Links[i-1].terminal && !Links[i+1].terminal) {
            if (Links[i].update_e_torsion)
            {
                // TODO: is this the most efficient way?

                b1 = Residues[i].position - Residues[i-1].position;
                b2 = Residues[i+1].position - Residues[i].position;
                b3 = Residues[i+2].position - Residues[i+1].position;
                b2xb3 = b2.cross(b3);
                phi = atan2((b2.magnitude() * b1.dot(b2xb3)), b1.cross(b2).dot(b2xb3));

                Links[i].pseudo_torsion = phi;
                // eqn 11: kim2008
                r1 = Residues[i].aminoAcidIndex;
                r2 = Residues[i+1].aminoAcidIndex;
                Links[i].e_torsion = (1 + cos(phi - torsions.getSigma(r1, r2, 1))) * torsions.getV(r1, r2, 1) +
                                    (1 + cos(2 * phi - torsions.getSigma(r1, r2, 2))) * torsions.getV(r1, r2, 2) +
                                    (1 + cos(3 * phi - torsions.getSigma(r1, r2, 3))) * torsions.getV(r1, r2, 3) +
                                    (1 + cos(4 * phi - torsions.getSigma(r1, r2, 4))) * torsions.getV(r1, r2, 4);
                Links[i].update_e_torsion = false;
            }
            e_torsion += Links[i].e_torsion;
        }
    }
    return e_torsion;
}

#endif
