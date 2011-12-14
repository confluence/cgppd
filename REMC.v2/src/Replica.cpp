#include "Replica.h"

using namespace std;

#if USING_CUDA || INCLUDE_TIMERS
#include <cutil.h>
#endif

Replica::Replica()
{
    temperature = 300.0f;
    label = -1;
    potential = 0.0f;
    moleculeCount = 0;
    moleculeArraySize = 0;
    residueCount = 0;
    maxMoleculeSize = 0;
    nonCrowderCount = 0;
    accept = 0;
    acceptA = 0;
    reject = 0;
    totalAcceptReject = 0;
    totalAccept = 0;
    boundSamples = 0;
    samplesSinceLastExchange = 0;
    totalBoundSamples = 0;
    totalSamples = 0;
    paircount = 0;
    nonCrowderResidues = 0;
    timersInit = false;
    translateStep = INITIAL_TRANSLATIONAL_STEP;
    rotateStep = INITIAL_ROTATIONAL_STEP;

}

#if INCLUDE_TIMERS
void Replica::initTimers()
{
    //timers for profiling the cuda functions
    replicaToGPUTimer = 0;
    replicaUpdateGPUTimer = 0;
    replicaECUDATimer = 0;
    replicaMoleculeUpdateTimer = 0;
    initGPUMemoryTimer = 0;

    CUT_SAFE_CALL( cutCreateTimer(&replicaToGPUTimer) );
    CUT_SAFE_CALL( cutCreateTimer(&replicaToHostTimer) );
    CUT_SAFE_CALL( cutCreateTimer(&replicaUpdateGPUTimer) );
    CUT_SAFE_CALL( cutCreateTimer(&replicaMoleculeUpdateTimer) );
    CUT_SAFE_CALL( cutCreateTimer(&replicaECUDATimer) );
    CUT_SAFE_CALL( cutCreateTimer(&initGPUMemoryTimer) );
    CUT_SAFE_CALL( cutCreateTimer(&replicaEHostTimer) );

    timersInit = true;
}
#endif

// Replica::Replica(const Replica& r)
// {
//     Replica();
//     label = r.label;
//     temperature = r.temperature;
//     moleculeCount = r.moleculeCount;
//     residueCount = r.residueCount;
//     aminoAcids = r.aminoAcids;
//     molecules = new Molecule[moleculeCount];
//     maxMoleculeSize = r.maxMoleculeSize;
//     nonCrowderCount = r.nonCrowderCount;
//     nonCrowderResidues = r.nonCrowderResidues;
//     for (size_t m=0; m<moleculeCount; m++)
//     {
//         molecules[m] = r.molecules[m];
//         cout << "molecule copy: " << &molecules[m] << " <- " << &r.molecules[m]<< endl;
//     }
// #if USING_CUDA
//     blockSize = r.blockSize;
//     sharedMemSize = r.sharedMemSize;
//     replicaIsOnDevice = false;
// #endif
//
//     paircount = r.paircount;
//     rotateStep = r.rotateStep;
//     translateStep = r.translateStep;
// }

void Replica::setAminoAcidData(AminoAcids a)
{
    aminoAcids = a;
}

Replica::~Replica()
{
    //aminoAcids = AminoAcids();
#if INCLUDE_TIMERS
    if (timersInit)
    {
        CUT_SAFE_CALL( cutDeleteTimer(replicaToGPUTimer) );
        CUT_SAFE_CALL( cutDeleteTimer(replicaToHostTimer) );
        CUT_SAFE_CALL( cutDeleteTimer(replicaUpdateGPUTimer) );
        CUT_SAFE_CALL( cutDeleteTimer(replicaMoleculeUpdateTimer) );
        CUT_SAFE_CALL( cutDeleteTimer(replicaECUDATimer) );
        CUT_SAFE_CALL( cutDeleteTimer(initGPUMemoryTimer) );
        CUT_SAFE_CALL( cutDeleteTimer(replicaEHostTimer) )
    }
#endif
    if (moleculeCount > 0)
        delete [] molecules;
}

int Replica::countpairs()
{
    paircount = 0;
    for (size_t mI=0; mI<moleculeCount; mI++)
        for (size_t mJ=mI+1; mJ<moleculeCount; mJ++)
            for (size_t mi=0; mi<molecules[mI].residueCount; mi++)
                for (size_t mj=0; mj<molecules[mJ].residueCount; mj++)
                    paircount++;
    return paircount;
}

void Replica::exchangeReplicas(Replica &r)
{
    swap(label,r.label);
    swap(temperature,r.temperature);

    // for statistics
    swap(totalAcceptReject,r.totalAcceptReject);
    swap(totalAccept,r.totalAccept);
    swap(boundSamples,r.boundSamples);
    swap(samplesSinceLastExchange,r.samplesSinceLastExchange);
    swap(totalBoundSamples,r.totalBoundSamples);
    swap(totalSamples,r.totalSamples);


}

void Replica::copy(const Replica &r)
{
    label = r.label;
    temperature = r.temperature;
    aminoAcids = r.aminoAcids;
    residueCount = r.residueCount;
    maxMoleculeSize = r.maxMoleculeSize;
    boundingValue = r.boundingValue;
    nonCrowderCount = r.nonCrowderCount;
    nonCrowderResidues = r.nonCrowderResidues;

//     if (molecules != NULL && moleculeCount != r.moleculeCount )
    if (moleculeCount != 0 && moleculeCount != r.moleculeCount )
    {
        try
        {
            delete [] molecules;
        }
        catch ( char * str )
        {
            cout << "Exception raised: delete [] molecules in Replica::copy()" << endl;
        }
    }

    molecules = new Molecule[r.moleculeCount];
    moleculeCount = r.moleculeCount;

    //the next bit is important because it makes sure that residues are contiguous in memory => better cpu performance
    contiguousResiduesSize = int(32.0f*ceil(float(residueCount)/32.0f));
    contiguousResidues = new Residue[contiguousResiduesSize];
    // contiguousResidues[residueCount..contiguousResiduesSize-1].aminoAcidIndex == PADDER_IDENTIFIER
    for (int i=residueCount; i<contiguousResiduesSize; i++)
    {
        contiguousResidues[i].aminoAcidIndex = int(PADDER_IDENTIFIER);
    }


    int rescount = 0;
    for (size_t m=0; m<moleculeCount; m++)
    {
        memcpy(contiguousResidues+rescount,r.molecules[m].Residues,r.molecules[m].residueCount*sizeof(Residue));

        molecules[m].Residues = contiguousResidues+rescount;
        molecules[m].residueCount = r.molecules[m].residueCount;
        molecules[m].translationalStep = r.molecules[m].translationalStep;
        molecules[m].rotationalStep = r.molecules[m].rotationalStep;
        molecules[m].AminoAcidsData = r.molecules[m].AminoAcidsData;
        molecules[m].center = r.molecules[m].center;
        molecules[m].position = r.molecules[m].position;
        molecules[m].rotation = r.molecules[m].rotation;
        molecules[m].moleculeRoleIdentifier = r.molecules[m].moleculeRoleIdentifier;
        rescount += r.molecules[m].residueCount;
    }
    potential = r.potential;
#if USING_CUDA
    blockSize = r.blockSize;
    sharedMemSize = r.sharedMemSize;
#endif
}

void Replica::reserveContiguousMoleculeArray(int size)
{
    if (moleculeCount != 0)
    {
        cout << " ** Molecules already exist. Cannot reserve size" << endl;
        return;
    }
    moleculeArraySize = size;
    molecules = new Molecule[size];
}

int Replica::loadMolecule(const char* pdbfilename)
{
    if (moleculeCount+1 > moleculeArraySize) //need a new array
    {
        moleculeArraySize += 3;
        Molecule *new_molecules = new Molecule[moleculeArraySize];
        if (moleculeCount > 0)
            memcpy(new_molecules,molecules,sizeof(Molecule)*moleculeCount);
        try	{
            delete [] molecules;    // free existing array before forgetting about it
        }
        catch ( char * str ) 	{
            cout << "Exception raised: delete [] molecules failed" << endl;
        }

        molecules = new Molecule[moleculeArraySize];
        memcpy(molecules,new_molecules,sizeof(Molecule)*moleculeCount);
        delete [] new_molecules;
    }

    molecules[moleculeCount].AminoAcidsData = aminoAcids;
    molecules[moleculeCount].index = moleculeCount;
    molecules[moleculeCount].initFromPDB(pdbfilename);
    residueCount += molecules[moleculeCount].residueCount;
    maxMoleculeSize = max(maxMoleculeSize,molecules[moleculeCount].residueCount);
    moleculeCount++;
    return moleculeCount-1;
}

int Replica::loadMolecule(const char* pdbfilename, Vector3f position, Vector3double rotationAxis, double rotationAmount)
{
    int i = loadMolecule(pdbfilename);
    rotationAxis.normalizeInPlace();
    molecules[i].setPosition(position);
    molecules[i].rotateQ(rotationAxis,rotationAmount);
    return i;
}

void Replica::initRNGs()
{
    unsigned long long seed = time (NULL);
    srand(time(NULL)+(label+1)*(label+1));

    rng_moleculeSelection = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (rng_moleculeSelection,random());

    rng_rotate = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (rng_rotate,random());

    rng_rotateAmount = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (rng_rotateAmount,random());

    rng_translate = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (rng_translate,random());

    rng_translateAmount = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (rng_translateAmount,random());

    MCRng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set (MCRng,random());

    MCKbRng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set (MCKbRng,random());
}

void Replica::freeRNGs()
{
    gsl_rng_free (rng_moleculeSelection);
    gsl_rng_free (rng_rotateAmount);
    gsl_rng_free (rng_translateAmount);
    gsl_rng_free (rng_rotate);
    gsl_rng_free (rng_translate);
    gsl_rng_free (MCRng);
    gsl_rng_free (MCKbRng);
}

Vector3f Replica::createNormalisedRandomVector(gsl_rng * r)
{
    //	float denominator = gsl_rng_max(r) - gsl_rng_min(r);
    //	Vector3f x = 2.0f*(Vector3f(gsl_rng_get(r)/denominator-0.5f,gsl_rng_get(r)/denominator-0.5f,gsl_rng_get(r)/denominator-0.5f));
    //Vector3f x = (Vector3f(gsl_rng_uniform(r)-0.5,0.5-gsl_rng_uniform(r),gsl_rng_uniform(r)-0.5)).normalize();
    //return x;
    return (Vector3f(gsl_rng_uniform(r)-0.5,0.5-gsl_rng_uniform(r),gsl_rng_uniform(r)-0.5)).normalize();
}

Vector3double Replica::createNormalisedRandomVectord(gsl_rng * r)
{
    //double denominator = gsl_rng_max(r) - gsl_rng_min(r);
    Vector3double x(gsl_rng_uniform(r)-0.5,0.5-gsl_rng_uniform(r),gsl_rng_uniform(r)-0.5);
    x.normalizeInPlace();
    return x;
}

#define _translate	0
#define _rotate 	1

// here for profiling
inline void Replica::rotate(const int m, const double rotateStep)
{
    molecules[m].rotateQ(createNormalisedRandomVectord(rng_rotate),rotateStep);
#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
    cout << "    Rotate: Replica "<< label << "/Molecule " << moleculeNo << endl;
#endif
}

// here for profiling
void Replica::translate(const int m, const float translateStep)
{
    // pick directions to rotate in, do by generating a random normalised vector of length equal to INITIAL_TRANSLATIONAL_STEP
    //Vector3f oldPosition = molecules[moleculeNo].center;
    Vector3f translateVector = translateStep * createNormalisedRandomVector(rng_translate);
#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
    cout << "    Translate: Replica "<< label << "/Molecule " << moleculeNo << endl;
#endif


    // bounding sphere conditions
#if BOUNDING_METHOD == BOUNDING_SPHERE
    // if the translate falls inside the bounding radius
    if ((molecules[m].center + translateVector).sumSquares() < boundingValue*boundingValue)
    {
        molecules[m].translate(translateVector);
    }
    // if the change causes the molecule to fall outside the bounding radius
    else
    {
#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
        cout << "  - Reject: Replica "<< label << "/Molecule " << moleculeNo <<  " : falls outside bounding sphere" << endl;
#endif
    }
#elif BOUNDING_METHOD == PERIODIC_BOUNDARY 	// periodic boundary conditions
    Vector3f newPosition = molecules[m].center + translateVector;
    newPosition.x = fmod(newPosition.x+boundingValue,boundingValue);
    newPosition.y = fmod(newPosition.y+boundingValue,boundingValue);
    newPosition.z = fmod(newPosition.z+boundingValue,boundingValue);
    molecules[m].setPosition(newPosition);
#endif
}



void Replica::MCSearch(int steps)
{
    //float translateStep = INITIAL_TRANSLATIONAL_STEP;
    //double rotateStep = INITIAL_ROTATIONAL_STEP;

    bool lastOperationWasRotate;

    float oldPotential = potential;

    Molecule savedMolecule;
    savedMolecule.reserveResidueSpace(maxMoleculeSize);

    for (int step=0; step<steps; step++)
    {
#if OUTOUTPUT_LEVEL >= PRINT_MC_STEP_COUNT
        cout << "Step: " << step << endl;
#endif
        uint moleculeNo = (int) gsl_rng_uniform_int(rng_moleculeSelection,moleculeCount);

        uint mutationType = gsl_ran_bernoulli (MCRng,translate_rotate_bernoulli_bias);

        // save the current state so we can roll back if it was not a good mutation.
        savedMolecule.saveBeforeStateChange(&molecules[moleculeNo]);

        switch (mutationType)
        {
        case _rotate:
        {
            rotate(moleculeNo, rotateStep);

#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
            cout << "    Rotate: Replica "<< label << "/Molecule " << moleculeNo << endl;
#endif
            break;
        }

        case _translate:
        {
            translate(moleculeNo, translateStep);

#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
            cout << "    Translate: Replica "<< label << "/Molecule " << moleculeNo << endl;
#endif

            break;
        }
        default:
            break;
        }

#if CUDA_E
        // copy host data to device. so we can do the calculations on it.
        MoleculeDataToDevice(moleculeNo);

        double newPotential(EonDevice());  // sequential cuda call
        //if (abs(temperature-300)<1) cout << newPotential << " " << EonDeviceNC() << endl;
#if PERFORM_GPU_AND_CPU_E
        float cpu_e(E());
        float err = abs(cpu_e-newPotential)/abs(cpu_e);
        printf("%24.20f %24.20f %24.20f\n",cpu_e,float(newPotential),err);
#endif
#else // only CPU calls
        float newPotential = E();
#endif

        float delta = newPotential - oldPotential;

        // accept change if its better.
        if (delta < 0.0)
        {
            potential = newPotential;
            oldPotential = potential;
            accept++;
#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
            cout << "  * Replace: Replica "<< label << "/Molecule " << moleculeNo <<  " : delta E = " << delta << " E = " << potential << endl;
#endif
        }
        // accept change if it meets the boltzmann criteria, must be (kJ/mol)/(RT), delta is in kcal/mol @ 294K
        else if (gsl_rng_uniform(MCKbRng)<exp(-(delta*4184.0f)/(Rgas*temperature)))
        {
            potential = newPotential;
            oldPotential = potential;
            acceptA++;
#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
            cout << "  **Replace: Replica "<< label << "/Molecule " << moleculeNo <<  " : delta E = " << delta << " U < " << exp(-delta * 4.184f/(Rgas*temperature)) << " E = " << potential << endl;
#endif
        }
        //reject
        else
        {
#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
            cout << "  - Reject: Replica "<< label << "/Molecule " << moleculeNo <<  " : delta E = " << delta << " E = " << oldPotential << endl;
#endif

            reject++;
            molecules[moleculeNo].undoStateChange(&savedMolecule);

#if CUDA_E
            MoleculeDataToDevice(moleculeNo); // you have to update the device again because the copy will be inconsistent
#endif
            potential = oldPotential;
        }
    }
    delete [] savedMolecule.Residues;

}

// bool Replica::savePDB(const char *filename) // saves multiple pdb files per replica
// {
//     char filenameExt[256];
//     char tmp[64];
//     for (size_t i=0; i<moleculeCount; i++)
//     {
//         strcpy (filenameExt,filename);
//         sprintf (tmp,"%02d",int(i));
//         strcat (filenameExt,tmp);
//         strcat (filenameExt,".pdb");
//         molecules[i].saveAsPDB(filenameExt);
//     }
//     return true;
// }
//
// void Replica::saveAsSinglePDB(const char *filename)
// {
//     FILE * output;
//     output = fopen (filename,"w");
//     fprintf(output,"REMARK %s \n",filename);
//     fprintf(output,"REMARK potential: %0.10f \n",float(potential));
//     fprintf(output,"REMARK temperature: %5.1f \n",temperature);
//
//     for (size_t i=0; i<moleculeCount; i++)
//     {
//         fprintf(output,"REMARK Molecule: %d\n",int(i));
//         fprintf(output,"REMARK Rotation relative to input Q(w,x,y,z): %f %f %f %f\n",molecules[i].rotation.w,molecules[i].rotation.x,molecules[i].rotation.y,molecules[i].rotation.z);
//         fprintf(output,"REMARK Centriod position P(x,y,z): %f %f %f\n",molecules[i].center.x,molecules[i].center.y,molecules[i].center.z);
//     }
//
//     char chainId = 'A';
//     int itemcount = 0;
//     int lastSeqNo = 0;
//     for (size_t m=0; m<moleculeCount; m++)
//     {
//         size_t i=0;
//         while (i<molecules[m].residueCount)
//         {
//             itemcount++;
//             fprintf(output,"ATOM  %5d %4s%C%3s %C%4d%C  %8.3f%8.3f%8.3f%6.2f%6.2f\n",itemcount,"CA",' ',aminoAcids.get(molecules[m].Residues[i].aminoAcidIndex).getSNAME(),chainId,molecules[m].Residues[i].resSeq,' ',molecules[m].Residues[i].position.x,molecules[m].Residues[i].position.y,molecules[m].Residues[i].position.z,1.0f,1.0f);
//             lastSeqNo = molecules[m].Residues[i].resSeq;
//             i++;
//         }
//         fprintf(output,"TER   %5d      %3s %C%4d \n",itemcount,aminoAcids.get(molecules[m].Residues[i-1].aminoAcidIndex).getSNAME(),chainId,lastSeqNo);
//         chainId++;
//         fflush(output);
//     }
//     fprintf(output,"END \n");
//     fflush(output);
//     fclose(output);
// }
// simulation evaluations

inline float crowderPairPotential(const float r)
{
    return powf(6.0f/r,12.0f);
}

inline float distance(Vector3f a,Vector3f b, float _boxdim)
{
    float Xab(a.x-b.x);
    float Yab(a.y-b.y);
    float Zab(a.z-b.z);

    Xab = Xab - _boxdim * round(Xab/_boxdim);
    Yab = Yab - _boxdim * round(Yab/_boxdim);
    Zab = Zab - _boxdim * round(Zab/_boxdim);

    return sqrtf(Xab*Xab+Yab*Yab+Zab*Zab);
}

inline void kahan_sum(double &potential, const double p_ij, double &c)
{
    double y(p_ij - c);
    double t(potential + y);
    c = (t - potential) - y;
    potential = t;
}

inline void sum(double &potential, const double p_ij, double &c)
{
#if COMPENSATE_KERNEL_SUM
    kahan_sum(potential, p_ij, c);
#else
    potential += p_ij;
#endif
}

double Replica::E()
{
#if INCLUDE_TIMERS
    CUT_SAFE_CALL(cutStartTimer(replicaEHostTimer));
#endif

    double epotential = 0.0f;
    double LJAccumulator = 0.0f;
    double DHAccumulator = 0.0f;
    double DH_constant_component =  DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f ;

    double c_lj(0.0f);
    double c_dh(0.0f);

#if FLEXIBLE_LINKS
    PotentialComponents mol_e;
    double bond_accumulator = 0.0f;
    double angle_accumulator = 1.0f;
    double torsion_accumulator = 0.0f;

    double c_b(0.0f);
    double c_t(0.0f);
#endif

#define iRes molecules[mI].Residues[mi]
#define jRes molecules[mJ].Residues[mj]

    for (size_t mI=0; mI < moleculeCount; mI++)
    {
#if REPULSIVE_CROWDING
        if (molecules[mI].moleculeRoleIdentifier == CROWDER_IDENTIFIER)
        {
            for (size_t mJ = mI + 1; mJ < moleculeCount; mJ++)
            {
                for (size_t mi = 0; mi < molecules[mI].residueCount; mi++)
                {
                    for (size_t mj = 0; mj < molecules[mJ].residueCount; mj++)
                    {
                        double r(distance(iRes.position, jRes.position, boundingValue) + EPS);
                        if (r < const_repulsive_cutoff)
                        {
                            sum(LJAccumulator, crowderPairPotential(r), c_lj);
                        }
                    }
                }
            }
        }
        else
        {
#endif
#if FLEXIBLE_LINKS
            mol_e = molecules[mI].E();
            angle_accumulator *= mol_e.angle;
            sum(LJAccumulator, mol_e.LJ, c_lj);
            sum(DHAccumulator, mol_e.DH, c_dh);
            sum(bond_accumulator, mol_e.bond, c_b);
            sum(torsion_accumulator, mol_e.torsion, c_t);
#endif
            for (size_t mJ = mI + 1; mJ < moleculeCount; mJ++)
            {
#if REPULSIVE_CROWDING
                if (molecules[mJ].moleculeRoleIdentifier == CROWDER_IDENTIFIER)
                {
                    for (size_t mi = 0; mi < molecules[mI].residueCount; mi++)
                    {
                        for (size_t mj = 0; mj < molecules[mJ].residueCount; mj++)
                        {
                            double r(distance(iRes.position, jRes.position, boundingValue) + EPS);
                            if (r < const_repulsive_cutoff)
                            {
                                sum(LJAccumulator, crowderPairPotential(r), c_lj);
                            }
                        }
                    }
                }
                else
                {
#endif
                    for (size_t mi = 0; mi < molecules[mI].residueCount; mi++)
                    {
                        for (size_t mj = 0; mj < molecules[mJ].residueCount; mj++)
                        {
                            double r(distance(iRes.position, jRes.position, boundingValue) + EPS);
                            double DH(iRes.DH_component(jRes, r));
                            double LJ(iRes.LJ_component(jRes, r, aminoAcids));
                            sum(DHAccumulator, DH, c_dh);
                            sum(LJAccumulator, LJ, c_lj);
                        }
                    }
#if REPULSIVE_CROWDING
                }
#endif
            }
#if REPULSIVE_CROWDING
        }
#endif
    }

#if FLEXIBLE_LINKS
    epotential = (DHAccumulator * DH_constant_component + LJAccumulator * LJ_CONVERSION_FACTOR + 0.5 * K_spring * bond_accumulator + -GammaAngleReciprocal * log(angle_accumulator) + torsion_accumulator) * KBTConversionFactor;
#else
    epotential = (LJAccumulator * LJ_CONVERSION_FACTOR + DHAccumulator * DH_constant_component) * KBTConversionFactor;
#endif

#if INCLUDE_TIMERS
    CUT_SAFE_CALL(cutStopTimer(replicaEHostTimer));
#endif
    return epotential;
}

double Replica::E(Molecule *a,Molecule *b)
{
    double epotential = 0.0f;
    double LJAccumulator = 0.0f;
    double DHAccumulator = 0.0f;
    double DH_constant_component =  1.602176487f * 1.602176487f * DH_CONVERSION_FACTOR;

    double c_lj(0.0f);
    double c_dh(0.0f);

#define aRes a->Residues[mi]
#define bRes b->Residues[mj]

    for (size_t mi = 0; mi < a->residueCount; mi++)
    {
        for (size_t mj = 0; mj < b->residueCount; mj++)
        {
            double r (distance(aRes.position, bRes.position, boundingValue) + EPS);
            double DH(aRes.DH_component(bRes, r));
            double LJ(aRes.LJ_component(bRes, r, aminoAcids));

            sum(DHAccumulator, DH, c_dh);
            sum(LJAccumulator, LJ, c_lj);
        }
    }

    epotential = (LJAccumulator * LJ_CONVERSION_FACTOR + DHAccumulator * DH_constant_component) * KBTConversionFactor;
    return epotential;
}

void Replica::printTimers()
{
#if INCLUDE_TIMERS

    cout << "Timer values (ms)" << endl;
    cout << "Total(ms)   \tAverage(ms)  \t Action" << endl;

    printf("%12.6f\t%12.6f\t Replica to GPU (init & transfer)\n",   cutGetTimerValue(replicaToGPUTimer),            cutGetAverageTimerValue(replicaToGPUTimer));
    printf("%12.6f\t%12.6f\t Replica Update on GPU (transter)\n",   cutGetTimerValue(replicaUpdateGPUTimer),        cutGetAverageTimerValue(replicaUpdateGPUTimer));
    printf("%12.6f\t%12.6f\t Kernel Timer (computation)\n",         cutGetTimerValue(replicaECUDATimer),            cutGetAverageTimerValue(replicaECUDATimer));
    printf("%12.6f\t%12.6f\t Host Timer (computation)\n",           cutGetTimerValue(replicaEHostTimer),            cutGetAverageTimerValue(replicaEHostTimer));
    printf("%12.6f\t%12.6f\t Update Molecule on GPU (transfer)\n",  cutGetTimerValue(replicaMoleculeUpdateTimer),   cutGetAverageTimerValue(replicaMoleculeUpdateTimer));
    printf("%12.6f\t%12.6f\t GPU Memory Initialisation (malloc)\n", cutGetTimerValue(initGPUMemoryTimer),           cutGetAverageTimerValue(initGPUMemoryTimer));
    printf("Kernel Speedup:  %0.1fx\n", cutGetAverageTimerValue(replicaEHostTimer)/cutGetAverageTimerValue(replicaECUDATimer));

#else

    cout << "Timers disabled: set INCLUDE_TIMERS 1 and USING_CUDA 1 in definitions.h" << endl;

#endif
}

void Replica::countNonCrowdingResidues()
{
    nonCrowderResidues = 0;
    for (int m=0; m<moleculeCount; m++)
    {
        if (molecules[m].moleculeRoleIdentifier >= 0)
        {
            nonCrowderResidues += molecules[m].residueCount;
        }
    }
}

#if USING_CUDA
// all functions following this line are dependent on CUDA

void Replica::setBlockSize(int blockSize)
{
    this->blockSize = blockSize;

    // also set the shared memory required per block as the 2 are linked
    sharedMemSize = sizeof(float);
#if POSITIONDATA_MEMORY == SHARED_MEM
    sharedMemSize += sizeof(float4);
#endif
#if METADATA_MEMORY == SHARED_MEM
    sharedMemSize += sizeof(float4);
#endif
    sharedMemSize *= blockSize;
}

#if CUDA_STREAMS
void Replica::ReserveSumSpace()
{
    if (!replicaIsOnDevice)
    {
        cout << "! Error: replica not on device implies paddedSize==0; cannot perform Replica::ReserveSumSpace()" << endl;
    }
    // result stored on the device
    // gridSize can be arbitrary
    gridSize = paddedSize/blockSize;
    resultSize = gridSize;
    //for a parallel sum each grid must have one cell in the array of results from all the threads
    cudaMallocHost((void**)&kernelResult, sizeof(float)*resultSize*resultSize);
    cudaMalloc((void **)&device_kernelResult,sizeof(float)*resultSize*resultSize);

    cudaError_t err;
    if ((err = cudaGetLastError()) != cudaSuccess)
        printf("CUDA error: %s\n", cudaGetErrorString(err));

    memset(kernelResult,0,sizeof(float)*resultSize*resultSize);
    //CUDA_memcpy_to_device_async(device_kernelResult,kernelResult,sizeof(float)*resultSize*resultSize,cudaStream);
    //CUDA_memcpy_to_device(device_kernelResult,kernelResult,sizeof(float)*resultSize*resultSize);
    cudaMemset(device_kernelResult, 0, sizeof(float)*resultSize*resultSize);
}

void Replica::FreeSumSpace()
{
    cudaFreeHost(kernelResult);
    cudaFree((void**)&device_kernelResult);
}

float Replica::SumGridResults()
{
    float potentialSum = 0.0f;
    cudaStreamSynchronize(cudaStream); // ensure that the results are written back in time.

    /*	float c(0.0f); //kahanSum
    for (int i=0;i<resultSize*resultSize;i++)
    {
    		float y(kernelResult[i] - c);
    		float t((potentialSum) + y);
    		c = (t-(potentialSum)) -y;
    		potentialSum = t;
    }
    */

    for (int i=0; i<resultSize*resultSize; i++)
        potentialSum += kernelResult[i];


    cudaMemset(device_kernelResult, 0, sizeof(float)*resultSize*resultSize);
    return potentialSum * KBTConversionFactor;
}
#endif

#if CUDA_STREAMS

// 1/3 of the above function, does the mutation on the gpu asynchronously
void Replica::MCSearchMutate()
{
    oldPotential = potential;

    uint moleculeNo = (int) gsl_rng_uniform_int(rng_moleculeSelection,moleculeCount);
    lastMutationIndex = moleculeNo;

    uint mutationType = gsl_ran_bernoulli (MCRng,translate_rotate_bernoulli_bias);

    // save the current state so we can roll back if it was not a good mutation.
    savedMolecule.saveBeforeStateChange(&molecules[moleculeNo]);

    switch (mutationType)
    {
    case _rotate:
    {
        rotate(moleculeNo, rotateStep);
        break;
    }

    case _translate:
    {
        translate(moleculeNo, translateStep);
        break;
    }

    default:
        break;
    }
    MoleculeDataToDevice(moleculeNo);
    newPotential = 0.0f;
}

void Replica::MCSearchEvaluate()
{
    EonDeviceAsync();
}

void Replica::MCSearchAcceptReject()
{
    //cudaStreamSynchronize(cudaStream);  // sync, newPotential needs to have been returned
    newPotential = SumGridResults();

    //cout << "new energy replica[" << temperature << "] = " << newPotential << endl;

#if PERFORM_GPU_AND_CPU_E
    float cpu_e(E());
    float err = abs(cpu_e-newPotential)/abs(cpu_e);
    printf("%24.20f %24.20f %24.20f\n",cpu_e,float(newPotential),err);
#endif

    float delta = (newPotential - oldPotential);  // needs to be in K_bT
    // if the change is bad then discard it.
    if (delta < 0.0f)
    {
        potential = newPotential;
        oldPotential = potential;
#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
        cout << "  * Replace: Replica "<< label << "/Molecule " << lastMutationIndex <<  " : delta E = " << delta << " E = " << potential << endl;
#endif
    }
    else if (gsl_rng_uniform(MCKbRng)<exp(-delta*4.184f/(Rgas*temperature)))
    {
        potential = newPotential;
        oldPotential = potential;
#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
        cout << "  **Replace: Replica "<< label << "/Molecule " << lastMutationIndex <<  " : delta E = " << delta << " U < " << exp(-delta*4.184f/(Rgas*temperature)) << " E = " << potential << endl;
#endif
    }
    else
    {
#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
        cout << "  - Reject: Replica "<< label << "/Molecule " << lastMutationIndex <<  " : delta E = " << delta << " E = " << oldPotential << endl;
#endif
        molecules[lastMutationIndex].undoStateChange(&savedMolecule);
        potential = oldPotential;
        MoleculeDataToDevice(lastMutationIndex); // you have to update the device again because the copy will be inconsistent
    }
}
#endif  // streams

// copy the replica to the device
void Replica::ReplicaDataToDevice()
{
#if INCLUDE_TIMERS
    CUT_SAFE_CALL(cutStartTimer(initGPUMemoryTimer));
#endif

    // allocate total size needed.

    paddedSize = int(ceil(float(residueCount)/float(blockSize)))*blockSize; // reserve a blocksize multiple because it allows for efficient summation
    dataSetSize = paddedSize;

#if CUDA_STREAMS
    cudaMallocHost((void**)&host_float4_residuePositions,sizeof(float4)*paddedSize);
#else
    host_float4_residuePositions = new float4[paddedSize];
#endif

    //host_float4_residuePositions = new float4[paddedSize];
    host_float4_residueMeta = new float4[paddedSize];
    host_moleculeStartPositions = new int[moleculeCount];

    cudaMalloc((void**)&device_float4_residuePositions,sizeof(float4)*paddedSize);
    cudaMalloc((void**)&device_float4_residueMeta,sizeof(float4)*paddedSize);
    cudaMalloc((void**)&device_moleculeStartPositions,sizeof(int)*moleculeCount);
    cudaMalloc((void**)&device_residueCount,sizeof(int));
    cudaMalloc((void**)&device_moleculeCount,sizeof(int));

#if INCLUDE_TIMERS
    CUT_SAFE_CALL(cutStopTimer(initGPUMemoryTimer));
    CUT_SAFE_CALL(cutStartTimer(replicaToGPUTimer));
#endif
    // keep a counter to calculate the offset from the beginning of the residues array
    int arrayIndex = 0;
    // pack the residues into an array and copy it to the correct place on the GPU
    for (int m=0; m<moleculeCount; m++)
    {
        host_moleculeStartPositions[m] = arrayIndex;   // start position of the first molecule in the host/gpu array

        float mf(m);

        for (int rc=0; rc<molecules[m].residueCount; rc++)
        {
            host_float4_residuePositions[arrayIndex].x = molecules[m].Residues[rc].position.x;
            host_float4_residuePositions[arrayIndex].y = molecules[m].Residues[rc].position.y;
            host_float4_residuePositions[arrayIndex].z = molecules[m].Residues[rc].position.z;
            host_float4_residuePositions[arrayIndex].w = mf;  // residue belongs to this molecule id
            host_float4_residueMeta[arrayIndex].x = molecules[m].Residues[rc].aminoAcidIndex;
            host_float4_residueMeta[arrayIndex].y = molecules[m].Residues[rc].electrostaticCharge;
            host_float4_residueMeta[arrayIndex].z = molecules[m].Residues[rc].vanderWaalRadius;
            host_float4_residueMeta[arrayIndex].w = molecules[m].moleculeRoleIdentifier;

            arrayIndex++;
        }
    }

    while ( arrayIndex < paddedSize )
    {
        host_float4_residuePositions[arrayIndex].x = 0;//float(rand()); // stops division by zero
        host_float4_residuePositions[arrayIndex].y = 0.0f;
        host_float4_residuePositions[arrayIndex].z = 0.0f;
        host_float4_residuePositions[arrayIndex].w = PADDER_IDENTIFIER;  // residue belongs to this molecule id

        //h_id[arrayIndex] = PADDER_IDENTIFIER;

        host_float4_residueMeta[arrayIndex].x = PADDER_AATYPE;  // amino acid index
        host_float4_residueMeta[arrayIndex].y = 0.0f;
        host_float4_residueMeta[arrayIndex].z = 0.0f;
        host_float4_residueMeta[arrayIndex].w = 0.0f;
        arrayIndex++;
    }


// cant stream as the amount of stuff is larger than pagable memory
#if 	CUDA_STREAMS
    cudaMemcpyAsync(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice,cudaStream);
    cudaStreamSynchronize(cudaStream);
#else
    // copy all the above to the device
    cudaMemcpy(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice);
#endif
    cudaMemcpy(device_float4_residueMeta, host_float4_residueMeta, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice);
    cudaMemcpy(device_moleculeStartPositions, host_moleculeStartPositions, sizeof(int)*moleculeCount, cudaMemcpyHostToDevice);
    cudaMemcpy(device_residueCount, &residueCount, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(device_moleculeCount, &moleculeCount, sizeof(int), cudaMemcpyHostToDevice);


#if METADATA_MEMORY == TEXTURE_MEM
    bindMetaDataToTexture(device_float4_residueMeta,sizeof(float4)*paddedSize);
#endif
#if POSITIONDATA_MEMORY == TEXTURE_MEM
    bindPositionDataToTexture(device_float4_residuePositions,sizeof(float4)*paddedSize);
#endif


#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(replicaToGPUTimer) );
    //cout << "Data to GPU Initialisation and Transfer time: " << (cutGetTimerValue(initGPUMemoryTimer)+cutGetTimerValue(replicaToGPUTimer)) << "ms"<< endl;
#endif

    replicaIsOnDevice = true;

    return;
}

void Replica::setLJpotentials(float *ljp)
{
    device_LJPotentials = ljp;
}

void Replica::ReplicaDataToHost()
{

#if INCLUDE_TIMERS
    CUT_SAFE_CALL(cutStartTimer(replicaToHostTimer));
#endif
    int residueCount = 0;

    for (int m=0; m<moleculeCount; m++)
    {
        residueCount += molecules[m].residueCount;
    }


//#ifdef 	CUDA_STREAMS ==== doesn't work because host code needs it to be synchronous for now
//	CUDA_memcpy_to_host_async(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*dataSetSize, cudaStream);
//	CUDA_memcpy_to_host_async(device_float4_residueMeta, host_float4_residueMeta, sizeof(float4)*dataSetSize, cudaStream);
//#else
    // copy all the above to the device
    cudaMemcpy(host_float4_residuePositions, device_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyDeviceToHost);
    //CUDA_memcpy_to_host(host_float4_residuePositions, device_float4_residuePositions, sizeof(float4)*paddedSize);

    // I dont think this next line is neccessary because it would not have changed on the device
    //CUDA_memcpy_to_host(device_float4_residueMeta, host_float4_residueMeta, sizeof(float4)*paddedSize);
//#endif


#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(replicaToHostTimer) );
#endif
    // position of the molecule in the gpu array
    int arrayIndex = 0;
    // unpack the residues from the host array
    for (int m=0; m<moleculeCount; m++)
    {
        arrayIndex = host_moleculeStartPositions[m];   // start position of the first molecule in the host/gpu array
        for (int offset=0; offset<molecules[m].residueCount; offset++)
        {
            molecules[m].Residues[offset].position.x = host_float4_residuePositions[arrayIndex+offset].x;
            molecules[m].Residues[offset].position.y = host_float4_residuePositions[arrayIndex+offset].y;
            molecules[m].Residues[offset].position.z = host_float4_residuePositions[arrayIndex+offset].z;
        }
    }

#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(replicaToHostTimer) );
#endif
    return;
}


void Replica::UpdateDeviceData()
{
    if (!replicaIsOnDevice)
    {
        cout << "ERROR: Replica::UpdateDeviceData() called without initialising device data"<< endl;
        ReplicaDataToDevice();
    }
#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStartTimer(replicaUpdateGPUTimer) );
#endif
    // works only if the replica is already on the device
    // keep a counter to calculate the offset from the beginning of the residues array
    int arrayIndex = 0;
    // pack the residues into an array and copy it to the correct place on the GPU
    for (int m=0; m<moleculeCount; m++)
    {
        host_moleculeStartPositions[m] = arrayIndex;   // start position of the first molecule in the host/gpu array
        for (int rc=0; rc<molecules[m].residueCount; rc++)
        {
            host_float4_residuePositions[arrayIndex].x = molecules[m].Residues[rc].position.x;
            host_float4_residuePositions[arrayIndex].y = molecules[m].Residues[rc].position.y;
            host_float4_residuePositions[arrayIndex].z = molecules[m].Residues[rc].position.z;
            host_float4_residuePositions[arrayIndex].w = float(m);  // residue belongs to this molecule id
            //	host_float4_residueMeta[arrayIndex].x = molecules[m].Residues[rc].aminoAcidIndex;
            //	host_float4_residueMeta[arrayIndex].y = molecules[m].Residues[rc].electrostaticCharge;
            //	host_float4_residueMeta[arrayIndex].z = molecules[m].Residues[rc].vanderWaalRadius;
            //	host_float4_residueMeta[arrayIndex].w = temperature;
            arrayIndex++;
        }
    }

    // copy all the above to the device
#if 	CUDA_STREAMS
    cudaMemcpyAsync(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice,cudaStream);
    //CUDA_memcpy_to_device_async(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize,cudaStream);

#else
    cudaMemcpy(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice);
    //CUDA_memcpy_to_device(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize);

#endif
#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(replicaUpdateGPUTimer) );
#endif
    return;
}

// update a molecule on the device
void Replica::MoleculeDataToDevice(int MoleculeID)
{
    if (!replicaIsOnDevice)
    {
        cout << "ERROR: Replica::MoleculeDataToDevice("<< MoleculeID << ") called without initialising device data."<< endl;
        ReplicaDataToDevice();
    }
#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStartTimer(replicaMoleculeUpdateTimer) );
#endif

    // TODO: change to use float4 for everything, will eliminate this copy

    /*
    1. get molecule on device's residues pointer
    2. make residue array
    3. copy to device at location found in 1
    */

    int residueIndex = host_moleculeStartPositions[MoleculeID];
    int memoryPosition = residueIndex;
    int moleculeSize = molecules[MoleculeID].residueCount;
    for (int rc=0; rc<molecules[MoleculeID].residueCount; rc++)
    {
        host_float4_residuePositions[residueIndex].x = molecules[MoleculeID].Residues[rc].position.x;
        host_float4_residuePositions[residueIndex].y = molecules[MoleculeID].Residues[rc].position.y;
        host_float4_residuePositions[residueIndex].z = molecules[MoleculeID].Residues[rc].position.z;
        host_float4_residuePositions[residueIndex].w = float(MoleculeID);  // residue belongs to this molecule id
        residueIndex++;
    }

    // copy to the device
#if CUDA_STREAMS
    cudaMemcpyAsync(&device_float4_residuePositions[memoryPosition], &host_float4_residuePositions[memoryPosition], sizeof(float4)*moleculeSize,cudaMemcpyHostToDevice,cudaStream);
    //CUDA_memcpy_to_device_async(&device_float4_residuePositions[memoryPosition], &host_float4_residuePositions[memoryPosition], sizeof(float4)*moleculeSize,cudaStream);
    //residueMeta will already be on the device
#else
    cudaMemcpy(&device_float4_residuePositions[memoryPosition], &host_float4_residuePositions[memoryPosition], sizeof(float4)*moleculeSize, cudaMemcpyHostToDevice);
    //CUDA_memcpy_to_device(&device_float4_residuePositions[memoryPosition], &host_float4_residuePositions[memoryPosition], sizeof(float4)*moleculeSize);
    //residueMeta will already be on the device
#endif
#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(replicaMoleculeUpdateTimer) );
#endif
    return;
}

void Replica::EonDeviceAsync()
{
#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStartTimer(replicaECUDATimer) );
#endif
#if CUDA_STREAMS
    // compute potential parts
    CUDA_EonDevice_async(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, device_kernelResult, resultSize, blockSize, dataSetSize, sharedMemSize, cudaStream);
    //CUDA_EonDevice(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, device_kernelResult,blockSize,dataSetSize);

    // write back to CPU
    //CUDA_memcpy_to_host_async(kernelResult,device_kernelResult,resultSize*resultSize*sizeof(float),cudaStream);
    cudaMemcpyAsync(kernelResult,device_kernelResult,resultSize*resultSize*sizeof(float),cudaMemcpyDeviceToHost,cudaStream);

#else
    cout << " ! Replica::EonDeviceAsync() can only be run using streams" << endl;
#endif
#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(replicaECUDATimer) );
#endif
}

double Replica::EonDevice()
{
#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStartTimer(replicaECUDATimer) );
#endif

    double result(0.0);

    CUDA_EonDevice(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, &result,blockSize,dataSetSize,sharedMemSize);

#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(replicaECUDATimer) );
#endif
    return result;
}

double Replica::EonDeviceNC()
{
    double result(0.0);
    int NCdataSetSize(ceil(float(nonCrowderResidues)/float(blockSize))*blockSize);
    //cout << "using " << NCdataSetSize <<  " as non NC size" << endl;

    CUDA_EonDeviceNC(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, &result,blockSize,NCdataSetSize,sharedMemSize);

    return result;
}

void Replica::FreeDevice()
{
    replicaIsOnDevice = false;  // no longer on device

#if CUDA_STREAMS
    cudaFreeHost(host_float4_residuePositions);
#else
    delete [] host_float4_residuePositions;
#endif
#if METADATA_MEMORY == TEXTURE_MEM
    freeMetaDataTexture();
#endif
#if POSITIONDATA_MEMORY == TEXTURE_MEM
    freePositionDataTexture();
#endif

    delete [] host_float4_residueMeta;
    delete [] host_moleculeStartPositions;

    cudaFree(device_float4_residuePositions);
    cudaFree(device_float4_residueMeta);
    cudaFree(device_moleculeStartPositions);
    cudaFree(device_residueCount);
    cudaFree(device_moleculeCount);

#if CUDA_MC
    cudaFree(device_moleculeCenters);
    cudaFree(device_moleculeLengths);
    cudaFree(device_translationVector);
    cudaFree(device_reverseTranslationVector);
    cudaFree(device_rotationVector);
    cudaFree(device_reverseRotationVector);
    cudaFreeHost(host_translationVector);
    cudaFreeHost(host_reverseTranslationVector);
    cudaFreeHost(host_rotationVector);
    cudaFreeHost(host_reverseRotationVector);
#endif
}

#if CUDA_MC
// rotate the molecule on device about a vector and amount
bool Replica::rotateOnDevice(int moleculeID, Vector3f rvector, float amount)
{
    //device_rotationVector;  // vector(x,y,z)|amount(w)

    host_rotationVector->x = rvector.x;
    host_rotationVector->y = rvector.y;
    host_rotationVector->z = rvector.z;
    host_rotationVector->w = amount;
    CUDA_memcpy_to_device_async(device_rotationVector,&host_rotationVector,sizeof(float4),cudaStream);

    CUDA_rotateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[moleculeID], &device_moleculeLengths[moleculeID], molecules[moleculeID].length, device_rotationVector, &device_moleculeCenters[moleculeID], cudaStream);

    // create the undo vector while its busy
    host_reverseRotationVector->x = rvector.x;
    host_reverseRotationVector->y = rvector.y;
    host_reverseRotationVector->z = rvector.z;
    host_reverseRotationVector->w = -amount;
    CUDA_memcpy_to_device_async(device_reverseTranslationVector,&host_reverseTranslationVector,sizeof(float4),cudaStream);


    return true;
}

bool Replica::translateOnDevice(int moleculeID, Vector3f translation)
{
    // send the vector to the device
    host_translationVector->x = translation.x;
    host_translationVector->y = translation.y;
    host_translationVector->z = translation.z;
    CUDA_memcpy_to_device_async(device_translationVector,&host_translationVector,sizeof(float4),cudaStream);

//  CUDA_translateMolecule (float4 *residuePositions, int *startPosition, int *moleculeLength, int moleculeLength, int *moleculeId, float4* translation, cudaStream_t stream)

    CUDA_translateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[moleculeID], &device_moleculeLengths[moleculeID], molecules[moleculeID].length, device_translationVector, &device_moleculeCenters[moleculeID], cudaStream);

    // create the undo vector while its busy
    host_reverseTranslationVector->x = -translation.x;
    host_reverseTranslationVector->y = -translation.y;
    host_reverseTranslationVector->z = -translation.z;
    CUDA_memcpy_to_device_async(device_reverseTranslationVector,&host_reverseTranslationVector,sizeof(float4),cudaStream);

    return true;
}

void Replica::cudaRollbackMutation()
{
    if (lastMutationWasATranslate)
    {
        CUDA_translateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[lastMutatedMolecule], &device_moleculeLengths[lastMutatedMolecule], molecules[lastMutatedMolecule].length, device_reverseTranslationVector,&device_moleculeCenters[lastMutatedMolecule], cudaStream);
    }
    else
    {
        CUDA_rotateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[lastMutatedMolecule], &device_moleculeLengths[lastMutatedMolecule], molecules[lastMutatedMolecule].length, device_reverseRotationVector,&device_moleculeCenters[lastMutatedMolecule], cudaStream);
    }
}
#endif // CUDA MC
#endif  // CUDA


// this code is a bit special...
bool Replica::sample(SimulationData *data, int current_step, float boundEnergyThreshHold, pthread_mutex_t *writeFileMutex)
{
    bool isBound = false;
    float nonCrowderPotential(0.0f);
    samplesSinceLastExchange++;

    // in the case of all molecules being of interest then just use the potential
    if (moleculeCount != nonCrowderCount)
    {
#if USING_CUDA
        nonCrowderPotential = EonDeviceNC(); // do on GPU, no stream support, can be implemented in 3mins if stream samlping is fixed.
#else
        // if there are crowders and molecules of interest the only use the energy of the interesting ones
        // TODO: does this work?! Are all the NC molecules at the front of the list?
        for (int i=0;i<nonCrowderCount;i++)
        {
            for (int j=i+1;j<nonCrowderCount;j++)
            {
                nonCrowderPotential += E(&molecules[i],&molecules[j]); // CPU failsafe
            }
        }
#endif
    }
    else
    {
        nonCrowderPotential = potential;
    }

    if (nonCrowderPotential < boundEnergyThreshHold)
    {
        isBound = true;
        boundSamples++;
        pthread_mutex_lock(writeFileMutex);
        fprintf(data->boundConformations, "%d; %0.5f (%0.5f); %0.1f\n", current_step,nonCrowderPotential,potential,temperature);

        for (int a=0; a<moleculeCount; a++)
        {
            if (a<nonCrowderCount)
                fprintf(data->boundConformations,"0 %2d: %f %f %f %f %f %f %f\n", a, molecules[a].rotation.w,molecules[a].rotation.x,molecules[a].rotation.y,molecules[a].rotation.z,molecules[a].center.x,molecules[a].center.y,molecules[a].center.z);
            else
                fprintf(data->boundConformations,"1 %2d: %f %f %f %f %f %f %f\n", a, molecules[a].rotation.w,molecules[a].rotation.x,molecules[a].rotation.y,molecules[a].rotation.z,molecules[a].center.x,molecules[a].center.y,molecules[a].center.z);
        }

        pthread_mutex_unlock(writeFileMutex);
    }

    // record the other samples, should be merged with the above block some time since we should record all samples anyway
    if (bool_recordAllSamples && !isBound)
    {
        pthread_mutex_lock(writeFileMutex);
        fprintf(data->boundConformations, "%d; %0.5f (%0.5f); %0.1f\n", current_step,nonCrowderPotential,potential,temperature);

        for (int a=0; a<moleculeCount; a++)
        {
            if (a<nonCrowderCount)
                fprintf(data->boundConformations,"0 %2d: %f %f %f %f %f %f %f\n", a, molecules[a].rotation.w,molecules[a].rotation.x,molecules[a].rotation.y,molecules[a].rotation.z,molecules[a].center.x,molecules[a].center.y,molecules[a].center.z);
            else
                fprintf(data->boundConformations,"1 %2d: %f %f %f %f %f %f %f\n", a, molecules[a].rotation.w,molecules[a].rotation.x,molecules[a].rotation.y,molecules[a].rotation.z,molecules[a].center.x,molecules[a].center.y,molecules[a].center.z);
        }

        pthread_mutex_unlock(writeFileMutex);
    }

    return isBound;
}

