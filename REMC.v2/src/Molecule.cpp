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
    //vector<Link> vLinks;

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

    char line[512] = {' '};
    char tmpLine[512] = {' '};
#ifdef VERY_VERBOSE
    //cout << "\tParsing: " << pdbfilename << "\n";
#endif
    // for eachline in the file
    while (!input.eof())
    {
        input.getline(line,512);
        strncpy(tmpLine,line,strlen(line));
        char *token = strtok(tmpLine," ");

        bool quit = false;
        while (!quit && token != NULL)
        {
            if (strcmp(token,"ATOM")==0)
            {
                int serial(0);
                char name[5] = {' '};
                char altLoc = ' ';
                char resName[4] = {' '};
                char chainID = ' ';
                int resSeq(0);
                char iCode = ' ';
                float x(0);
                float y(0);
                float z(0);
                float occupancy(0);
                float tempFactor(0);
                char segID[5] = {' '};
                char element[3] = {' '};
                char charge[3] = {' '};

                char tmp[6] = {' '};

                // parse serial
                strncpy(tmp,line+6,4);
                sscanf (tmp,"%d",&serial);

                // parse name
                strncpy(tmp,line+12,4);
                sscanf (tmp,"%s",name);

                //parse altLoc, has to be this char if present
                altLoc = line[16];

                // parse resName
                strncpy(tmp,line+17,3);
                sscanf (tmp,"%s",resName);
                //cout << "resName" << resName << endl;

                strncpy(tmp,line+21,1);
                sscanf (tmp,"%c",&chainID);

                // parse resName
                strncpy(tmp,line+22,4);
                sscanf (tmp,"%d",&resSeq);

                strncpy(tmp,line+26,1);
                sscanf (tmp,"%c",&iCode);

                // parse x, y, z
                char tmpVals[36]  = {' '};

                //strncpy(tmpVals,line+30,35);
                //char *tmptokens = strtok(tmpVals," ");
                strncpy (tmpVals,line+30,8);
                //	sscanf (tmpVals,"%8.3f%8.3f%8.3f%6.2f%6.2f",&x,&y,&z,&occupancy,&tempFactor);
                sscanf (tmpVals,"%f",&x);
                //tmptokens = strtok(NULL," ");
                strncpy (tmpVals,line+38,8);
                sscanf (tmpVals,"%f",&y);
                strncpy (tmpVals,line+46,8);
                //tmptokens = strtok(NULL," ");
                sscanf (tmpVals,"%f",&z);
                strncpy (tmpVals,line+54,6);
                //tmptokens = strtok(NULL," ");
                sscanf (tmpVals,"%f",&occupancy);
                strncpy (tmpVals,line+60,6);
                //tmptokens = strtok(NULL," ");
                sscanf (tmpVals,"%f",&tempFactor);

                if (strlen(line)>=80)
                {
                    // parse segID
                    strncpy(tmp,line+72,4);
                    sscanf (tmp,"%s",segID);
                    // parse element
                    strncpy(tmp,line+76,2);
                    sscanf (tmp,"%s",element);
                    // parse charge
                    strncpy(tmp,line+78,2);
                    sscanf (tmp,"%s",charge);
                }

                // push the atom onto the vector of our structure if it is a CA
                if (strcmp(name,"CA") == 0) // it is a CA atom. center of our bead
                {
                    Residue R;
                    //Link L;
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
                    /*
                    		  			if (vResidues.size()>0)
                    		  			{
                    		  				//R.relativePosition= R.position-vResidues[vResidues.size()-1].position;
                    		  				L.Angle = (float) R.relativePosition.angle(vResidues[vResidues.size()-1].relativePosition);
                    		  				//cout << L.Angle << endl;
                    		  				L.TorsionAngle = 0; // twist relative to previous bead, not availiable from this atm
                    		  				L.BondLength = R.relativePosition.magnitude(); // distance between this and the previous bead
                    		  				//cout << "Link: angle=" << L.Angle << " terminal=" << L.terminal << endl;
                    		  			}
                    */

                    //R.sasa = 1.0f;
                    //R.vanderWaalRadius = float(AminoAcidsData.data[R.aminoAcidIndex].vanderWaalRadius);
                    R.moleculeID = index;
                    vResidues.push_back(R);
                    //vLinks.push_back(L);
                }

            }
            else if (strcmp(token,"TER")==0) // Optional,Mandatory if ATOM records exist.
            {
                chainCount++;
                //vLinks[vLinks.size()-1].terminal = true;
            }
            else if (strcmp(token,"END")==0)
            {
                quit = true;
            }

            token = NULL;
        }
    }
    input.close();


    // copy the resides and links read into the arrays.
    Residues = new Residue[vResidues.size()];
    residueCount = vResidues.size();
    //	Links = new Link[vLinks.size()];
    //	linkCount = vLinks.size();

    // set the residues and accumulate the positions to determine the center of the molecule
    for(size_t r=0; r<residueCount; r++)
    {
        Residue R = vResidues[r];
        memcpy (&Residues[r],&R,sizeof(R));
        center.x += R.position.x;
        center.y += R.position.y;
        center.z += R.position.z;

    }
    //calculate the center
    center = center / residueCount;

    // set the relative positions of all residues
    for(size_t r=0; r<residueCount; r++)
    {
        Residues[r].relativePosition.x = Residues[r].position.x - center.x;
        Residues[r].relativePosition.y = Residues[r].position.y - center.y;
        Residues[r].relativePosition.z = Residues[r].position.z - center.z;
    }

    calculateSASA();
    calculateVolume();

    /*
    for(size_t l=0;l<linkCount;l++)
    {
    	Link L = vLinks[l];
    	memcpy (&Links[l],&L,sizeof(L));
    }
    */

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



/*
// total bond energy for this molecule
float Molecule::Ebond()
{
	float ebond = 0;
	for (size_t i=1; i<residueCount; i++)
	{
		float rmag = float((Residues[i].position-Residues[i-1].position).magnitude());  // reduces the number of sqrts by 1
		// eqn 9: kim2008
		ebond += (rmag - R0*Angstrom)*(rmag - R0*Angstrom);
	}
	return 0.5 * K_spring * ebond;   // reduce ops by multiplying at the end by constants.
}

// total angular bond energy for this molecule
float Molecule::Eangle()
{
	double eangle = 0;
	Vector3f ab;
	Vector3f cb;
	float theta;

	for (size_t i=1; i<residueCount-1; i++)
	{
		ab = Residues[i-1].position - Residues[i].position;
		cb = Residues[i+1].position - Residues[i].position;
		theta = ab.angle(cb);

		// eqn 10: kim2008
		eangle *= (exp(-GammaAngle * (KAlpha*(theta-ThetaAlpha)*(theta-ThetaAlpha) + EpsilonAlpha)) +
					exp(-GammaAngle * (KBeta *(theta-ThetaBeta) *(theta-ThetaBeta))) );
	}

	return -GammaAngle*log(eangle);
}

//torsional potential for the molecule
float Molecule::Etorsion()
{
	float etorsion = 0;
	uint i=0;
	int n;
	while (i<residueCount-4)
	{
		if (Links[i+3].terminal)
		{
			i+=4;
		}
		else
		{
			n = 1;
			etorsion += (1+cos(n*Links[i].TorsionAngle - torsions.getSigma(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n)))*torsions.getV(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n);
			n++;
			etorsion += (1+cos(n*Links[i+1].TorsionAngle - torsions.getSigma(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n)))*torsions.getV(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n);
			n++;
			etorsion += (1+cos(n*Links[i+2].TorsionAngle - torsions.getSigma(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n)))*torsions.getV(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n);
			n++;
			etorsion += (1+cos(n*Links[i+3].TorsionAngle - torsions.getSigma(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n)))*torsions.getV(Residues[i+1].aminoAcidIndex,Residues[i+2].aminoAcidIndex,n);
		}
		i++;
	}
	return etorsion;
}

float Molecule::Emembrane()
{
	float membranePotential = 0;
	for (size_t i=0; i<residueCount; i++)
	{
		membranePotential += Residues[i].electrostaticCharge * EchargedMenbrane(Residues[i].position.y,temperature) +
							 0.05f*pow((Zreference/Residues[i].position.y),12.0f);
	}
	return membranePotential;
}

// interaction potential between negatively charged membrane and a residue
// in: residue distance, temperature
float Molecule::EchargedMenbrane(float z, float T)
{
	float MembraneAlpha = (float)(exp((ephi0/(2.0f*K_b*T))-1.0f) / (exp(ephi0/(2.0f*K_b*T))+1.0f));
	float Psi = 2.0f*K_b*T/E_charge*log((1.0f+MembraneAlpha*exp(-kappa*z))/(1.0f-MembraneAlpha*exp(-kappa*z)));
	return Psi;
}
*/
