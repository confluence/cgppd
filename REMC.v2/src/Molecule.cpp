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

float Molecule::getMoleculeRoleIdentifier()
{
    return moleculeRoleIdentifier;
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

Vector3f Molecule::getPosition()
{
    return position;

}

int Molecule::getIndex()
{
    return index;
}

Quaternion Molecule::getRotation()
{
    return rotation;
}

int Molecule::getLength()
{
    return residueCount;
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

                Residue R;
#ifdef FLEXIBLE_LINKS
                Link L;
#endif
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
                vLinks.push_back(L);
#endif
            }
        }
        else if (strncmp(line, "TER", 3) == 0) // Optional,Mandatory if ATOM records exist.
        {
            chainCount++;
#ifdef FLEXIBLE_LINKS
            vLinks[vLinks.size()-1].terminal = true;
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
//         cout << l << endl;
        Link L = vLinks[l];
        //TODO is this the same relative position?!
//         cout << Residues[l + 1].relativePosition.magnitude() << ", " << Residues[l].relativePosition.magnitude() << endl;
        L.Angle = (float) Residues[l + 1].relativePosition.angle(Residues[l].relativePosition);
//         cout << L.Angle << endl;
        L.TorsionAngle = 0; // twist relative to previous bead; TODO
//         cout << L.TorsionAngle << endl;
        L.BondLength = Residues[l + 1].relativePosition.magnitude(); // TODO: is this right?
//         cout << L.BondLength << endl;
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

//TODO: calculate these just for the flexible linker residues

// total bond energy for this molecule
// float Molecule::Ebond()
// {
//     float ebond = 0;
//     for (size_t i=1; i<residueCount; i++)
//     {
//         float rmag = float((Residues[i].position-Residues[i-1].position).magnitude());  // reduces the number of sqrts by 1
//         // eqn 9: kim2008
//         ebond += (rmag - R0*Angstrom)*(rmag - R0*Angstrom);
//     }
//     return 0.5 * K_spring * ebond;   // reduce ops by multiplying at the end by constants.
// }

float Molecule::Ebond()
{
    float e_bond = 0.0;
    // TODO: iterate over flexible links only
    for (size_t i=0; i < linkCount; i++)
    {
        if (Links[i].update_e_bond)
        {
            // eqn 9: kim2008
            float rmag = float((Residues[i+1].position - Residues[i].position).magnitude());  // reduces the number of sqrts by 1
            Links[i].e_bond = (rmag - R0 * Angstrom) * (rmag - R0 * Angstrom);
            Links[i].update_e_bond = false;
        }
        e_bond += Links[i].e_bond;
    }
    return 0.5 * K_spring * e_bond;
}

// total angular bond energy for this molecule
// float Molecule::Eangle()
// {
//     double eangle = 0;
//     Vector3f ab;
//     Vector3f cb;
//     float theta;
//
//     for (size_t i=1; i<residueCount-1; i++)
//     {
//         ab = Residues[i-1].position - Residues[i].position;
//         cb = Residues[i+1].position - Residues[i].position;
//         theta = ab.angle(cb);
//
//         // eqn 10: kim2008
//         eangle *= (exp(-GammaAngle * (KAlpha*(theta-ThetaAlpha)*(theta-ThetaAlpha) + EpsilonAlpha)) +
//                     exp(-GammaAngle * (KBeta *(theta-ThetaBeta) *(theta-ThetaBeta))) );
//     }
//
//     return -GammaAngle*log(eangle);
// }

float Molecule::Eangle()
{
    double e_angle = 1.0;
    Vector3f ab;
    Vector3f cb;
    float theta;

    // TODO: iterate over flexible links only and check if there are enough residues at the beginning
    for (size_t i=1; i < linkCount; i++)
    {
        if (Links[i].update_e_angle)
        {
            ab = Residues[i-1].position - Residues[i].position;
            cb = Residues[i+1].position - Residues[i].position;
            theta = ab.angle(cb);
            // eqn 10: kim2008
            Links[i].e_angle = exp(-GammaAngle * (KAlpha * (theta - ThetaAlpha) * (theta - ThetaAlpha) + EpsilonAlpha)) +
                        exp(-GammaAngle * (KBeta *(theta - ThetaBeta) *(theta - ThetaBeta)));
            Links[i].update_e_angle = false;
        }
        e_angle *= Links[i].e_angle;
    }
    return -GammaAngleReciprocal * log(e_angle);
}

//torsional potential for the molecule
// float Molecule::Etorsion()
// {
//     float etorsion = 0;
//     uint i=0;
//     int n;
//     while (i<residueCount-4)
//     {
//         if (Links[i+3].terminal)
//         {
//             i+=4;
//         }
//         else
//         {
//             n = 1;
//             etorsion += (1+cos(n*Links[i].TorsionAngle - torsions.getSigma(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n)))*torsions.getV(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n);
//             n++;
//             etorsion += (1+cos(n*Links[i+1].TorsionAngle - torsions.getSigma(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n)))*torsions.getV(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n);
//             n++;
//             etorsion += (1+cos(n*Links[i+2].TorsionAngle - torsions.getSigma(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n)))*torsions.getV(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n);
//             n++;
//             etorsion += (1+cos(n*Links[i+3].TorsionAngle - torsions.getSigma(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n)))*torsions.getV(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n);
//         }
//         i++;
//     }
//     return etorsion;
// }

float Molecule::Etorsion()
{
    float e_torsion = 0.0;
    float torsion_angle;
    int r1;
    int r2;
    // TODO: iterate over flexible links only
    for (size_t i = 2; i < linkCount - 1; i++)
    {
        if (Links[i].update_e_torsion)
        {
            // TODO: calculate torsion angle here
            // eqn 11: kim2008
            Links[i].e_torsion = (1 + cos(torsion_angle - torsions.getSigma(r1, r2, 1))) * torsions.getV(r1, r2, 1) +
                                 (1 + cos(2 * torsion_angle - torsions.getSigma(r1, r2, 2))) * torsions.getV(r1, r2, 2) +
                                 (1 + cos(3 * torsion_angle - torsions.getSigma(r1, r2, 3))) * torsions.getV(r1, r2, 3) +
                                 (1 + cos(4 * torsion_angle - torsions.getSigma(r1, r2, 4))) * torsions.getV(r1, r2, 4);
            Links[i].update_e_torsion = false;
        }
        e_torsion += Links[i].e_torsion;
    }
    return e_torsion;
}

#endif

// float Molecule::Emembrane()
// {
//     float membranePotential = 0;
//     for (size_t i=0; i<residueCount; i++)
//     {
//         membranePotential += Residues[i].electrostaticCharge * EchargedMembrane(Residues[i].position.y,temperature) +
//                              0.05f*pow((Zreference/Residues[i].position.y),12.0f);
//     }
//     return membranePotential;
// }
//
// // interaction potential between negatively charged membrane and a residue
// // in: residue distance, temperature
// float Molecule::EchargedMembrane(float z, float T)
// {
//     float MembraneAlpha = (float)(exp((ephi0/(2.0f*K_b*T))-1.0f) / (exp(ephi0/(2.0f*K_b*T))+1.0f));
//     float Psi = 2.0f*K_b*T/E_charge*log((1.0f+MembraneAlpha*exp(-kappa*z))/(1.0f-MembraneAlpha*exp(-kappa*z)));
//     return Psi;
// }
