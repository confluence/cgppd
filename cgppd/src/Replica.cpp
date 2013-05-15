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
    newPotential = 0.0f;
    moleculeCount = 0;
    moleculeArraySize = 0;
    residueCount = 0;
    max_residue_count = 0;
#if FLEXIBLE_LINKS
    max_segment_count = 0;
#endif
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
    translateStep = INITIAL_TRANSLATIONAL_STEP;
    rotateStep = INITIAL_ROTATIONAL_STEP;
#if USING_CUDA
    replicaIsOnDevice = false;
#endif
#if INCLUDE_TIMERS
    timersInit = false;
#endif
}


void Replica::init_first_replica(const vector<moldata> mdata, AminoAcids amino_acid_data, const float bounding_value, const int initial_molecule_array_size) // constructor for initial replica; not final parameter list
{
    aminoAcids = amino_acid_data;
    boundingValue = bounding_value;
    reserveContiguousMoleculeArray(initial_molecule_array_size);

    for (size_t s = 0; s < mdata.size(); s++)
    {
        int mi = loadMolecule(mdata[s].pdbfilename);

        if (mdata[s].translate)
        {
            molecules[mi].translate(Vector3f(mdata[s].px, mdata[s].py, mdata[s].pz));
        }
        else
        {
            molecules[mi].setPosition(Vector3f(mdata[s].px, mdata[s].py, mdata[s].pz));
        }

        if (mdata[s].ra >= 0.000)
        {
            Vector3double v = Vector3double(mdata[s].rx,mdata[s].ry,mdata[s].rz);
            v.normalizeInPlace();
            molecules[mi].rotate(v,mdata[s].ra);
        }

        if (mdata[s].crowder)
        {
            molecules[mi].setMoleculeRoleIdentifier(CROWDER_IDENTIFIER);
        }
        else
        {
            nonCrowderCount++;
            nonCrowderResidues += molecules[mi].residueCount;
        }
    }
}

void Replica::init_child_replica(const Replica& ir, const int label, const float temperature, const float rotate_step, const float translate_step, const int thread_count)// constructor for final replicas; not final parameter list
{
    // TODO: incorporate copy method into this one and try to clean it up
    copy(ir);
#if INCLUDE_TIMERS
    initTimers();
#endif
    this->label = label;
    this->temperature = temperature;
    rotateStep = rotate_step;
    translateStep = translate_step;
    // TODO: only one range for entire simulation (but is the library threadsafe?)
    initRNGs();
    threadCount = thread_count;
#if FLEXIBLE_LINKS
    savedMolecule.init_saved_molecule(max_residue_count, max_segment_count);
#else
    savedMolecule.init_saved_molecule(max_residue_count);
#endif
    // leaving out stuff inside thread function for now
}

#if INCLUDE_TIMERS
void Replica::initTimers()
{
    //timers for profiling the cuda functions
    replicaToGPUTimer = 0;
    replicaToHostTimer = 0;
    replicaUpdateGPUTimer = 0;
    replicaMoleculeUpdateTimer = 0;
    replicaECUDATimer = 0;
    initGPUMemoryTimer = 0;
    replicaEHostTimer = 0;

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

Replica::~Replica()
{
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
    {
        delete [] molecules;
    }

    if (contiguousResiduesSize > 0)
    {
        delete [] contiguousResidues;
    }
}

// TODO: ugly, and only used once ever.
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
    max_residue_count = r.max_residue_count;
#if FLEXIBLE_LINKS
    max_segment_count = r.max_segment_count;
#endif
    boundingValue = r.boundingValue;
    nonCrowderCount = r.nonCrowderCount;
    nonCrowderResidues = r.nonCrowderResidues;

    if (moleculeCount != 0)
    {
        delete [] molecules;
    }
    molecules = new Molecule[r.moleculeCount];
    moleculeCount = r.moleculeCount;

    //the next bit is important because it makes sure that residues are contiguous in memory => better cpu performance
    contiguousResiduesSize = int(32.0f*ceil(float(residueCount)/32.0f));
    contiguousResidues = new Residue[contiguousResiduesSize];
    for (int i=residueCount; i<contiguousResiduesSize; i++)
    {
        contiguousResidues[i].aminoAcidIndex = int(PADDER_IDENTIFIER);
    }

    int rescount = 0;
    for (size_t m = 0; m < moleculeCount; m++)
    {
        molecules[m].copy(r.molecules[m], contiguousResidues + rescount);
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
        LOG(ERROR, "Cannot reserve molecule size: molecules already exist.");
        return;
    }
    moleculeArraySize = size;
    molecules = new Molecule[size];
}

// TODO: also add translation / position / rotation to this?
int Replica::loadMolecule(const char* pdbfilename)
{
    if (moleculeCount+1 > moleculeArraySize) //need a new array
    {
        moleculeArraySize += 3;
        Molecule *new_molecules = new Molecule[moleculeArraySize];
        if (moleculeCount > 0)
        {
            memcpy(new_molecules, molecules, sizeof(Molecule)*moleculeCount);
            delete [] molecules;    // free existing array before forgetting about it
        }

        molecules = new Molecule[moleculeArraySize];
        memcpy(molecules, new_molecules, sizeof(Molecule)*moleculeCount);
        delete [] new_molecules;
    }

    molecules[moleculeCount].init(pdbfilename, aminoAcids, moleculeCount, boundingValue);

    residueCount += molecules[moleculeCount].residueCount;
    max_residue_count = max(max_residue_count, molecules[moleculeCount].residueCount);
#if FLEXIBLE_LINKS
    max_segment_count = max(max_segment_count, molecules[moleculeCount].segmentCount);
#endif
    moleculeCount++;
    return moleculeCount-1;
}

// TODO: all of this needs to be done once for the whole simulation
void Replica::initRNGs()
{
    srand(time(NULL)+(label+1)*(label+1));

    rng = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (rng, random());

#if FLEXIBLE_LINKS
    MC_move_weights = new double[4];
    MC_move_weights[0] =  WEIGHT_MC_TRANSLATE;
    MC_move_weights[1] = WEIGHT_MC_ROTATE;
    MC_move_weights[2] = WEIGHT_MC_ROTATE_DOMAIN;
    MC_move_weights[3] = WEIGHT_MC_LOCAL;

    MC_discrete_table = gsl_ran_discrete_preproc(4, MC_move_weights);
    for (size_t m = 0; m < moleculeCount; m++)
    {
        // TODO this is hacky, but it will be gone when this is per-simulation
        molecules[m].MC_discrete_table = MC_discrete_table;
    }
#endif
}

// TODO: why is this not called from the destructor?!
void Replica::freeRNGs()
{
    gsl_rng_free(rng);

#if FLEXIBLE_LINKS
    // TODO: will this be a problem for the molecules?
    gsl_ran_discrete_free(MC_discrete_table);
    delete [] MC_move_weights;
#endif
}

void Replica::MCSearch(int steps)
{
    for (int step=0; step<steps; step++)
    {
        LOG(DEBUG_MC, "Step %3d:\t", step);
        uint moleculeNo = (int) gsl_rng_uniform_int(rng, moleculeCount);
        // save the current state so we can roll back if it was not a good mutation.
//         savedMolecule.saveBeforeStateChange(&molecules[moleculeNo]);
        savedMolecule.MC_backup_restore(&molecules[moleculeNo]);
        LOG(DEBUG_MC, "Replica %d:\tMolecule %d:\t", label, moleculeNo);

        molecules[moleculeNo].make_MC_move(rng, rotateStep, translateStep);

#if CUDA_E
        // copy host data to device. so we can do the calculations on it.
        MoleculeDataToDevice(moleculeNo);
        newPotential = EonDevice();  // sequential cuda call
        //if (abs(temperature-300)<1) cout << newPotential << " " << EonDeviceNC() << endl;
#if PERFORM_GPU_AND_CPU_E
        float cpu_e(E());
        float err = abs(cpu_e - newPotential) / abs(cpu_e);
        printf("%24.20f %24.20f %24.20f\n", cpu_e, float(newPotential), err);
#endif
#else // only CPU calls
        newPotential = E();
#endif
        float delta = newPotential - potential;

        // accept change if its better.
        if (delta < 0.0)
        {
            potential = newPotential;
            accept++;
            LOG(DEBUG_MC, "* Replace:\tdelta E = %f;\tE = %f\n", delta, potential);
        }
        // accept change if it meets the boltzmann criteria, must be (kJ/mol)/(RT), delta is in kcal/mol @ 294K
        else if (gsl_rng_uniform(rng) < exp(-(delta*4184.0f)/(Rgas*temperature)))
        {
            potential = newPotential;
            acceptA++;
            LOG(DEBUG_MC, "**Replace:\tdelta E = %f;\tE = %f;\tU < %f\n", delta, potential, exp(-delta * 4.184f/(Rgas*temperature)));
        }
        //reject
        else
        {
            reject++;
//             molecules[moleculeNo].undoStateChange(&savedMolecule);
            molecules[moleculeNo].MC_backup_restore(&savedMolecule);
            LOG(DEBUG_MC, "- Reject:\tdelta E = %f;\tE = %f\n", delta, potential);
#if CUDA_E
            MoleculeDataToDevice(moleculeNo); // you have to update the device again because the copy will be inconsistent
#endif
        }
    }
}

inline float crowderPairPotential(const float r)
{
    return powf(6.0f/r,12.0f);
}

double Replica::E()
{
#if INCLUDE_TIMERS
    CUT_SAFE_CALL(cutStartTimer(replicaEHostTimer));
#endif

    Potential potential;

#define iRes molecules[mI].Residues[mi]
#define jRes molecules[mJ].Residues[mj]

    for (size_t mI = 0; mI < moleculeCount; mI++)
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
                        double r(iRes.distance(jRes, boundingValue) + EPS);
                        if (r < const_repulsive_cutoff)
                        {
                                potential.increment_LJ(crowderPairPotential(r));
                        }
                    }
                }
            }
        }
        else
        {
#endif
#if FLEXIBLE_LINKS
            potential.increment(molecules[mI].E());
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
                            double r(iRes.distance(jRes, boundingValue) + EPS);
                            if (r < const_repulsive_cutoff)
                            {
                                potential.increment_LJ(crowderPairPotential(r));
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
                            double r(iRes.distance(jRes, boundingValue) + EPS);
                            potential.increment_LJ(calculate_LJ(iRes, jRes, r, aminoAcids));
                            potential.increment_DH(calculate_DH(iRes, jRes, r));
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

    potential.print_log(INFO_POTENTIAL, "\nTotal potential");
    double total(potential.total());

#if INCLUDE_TIMERS
    CUT_SAFE_CALL(cutStopTimer(replicaEHostTimer));
#endif
    return total;
}

double Replica::E(Molecule *a,Molecule *b)
{
    Potential potential;

#define aRes a->Residues[mi]
#define bRes b->Residues[mj]

    for (size_t mi = 0; mi < a->residueCount; mi++)
    {
        for (size_t mj = 0; mj < b->residueCount; mj++)
        {
            double r(aRes.distance(bRes, boundingValue) + EPS);
            potential.increment_LJ(calculate_LJ(aRes, bRes, r, aminoAcids));
            potential.increment_DH(calculate_DH(aRes, bRes, r));
        }
    }
#if FLEXIBLE_LINKS
    potential.increment(a->E());
    potential.increment(b->E());
#endif
    return potential.total();
}


#if FLEXIBLE_LINKS
double Replica::internal_molecule_E() {
    Potential potential;

    for (size_t mI = 0; mI < moleculeCount; mI++)
    {
#if REPULSIVE_CROWDING
        if (molecules[mI].moleculeRoleIdentifier != CROWDER_IDENTIFIER)
        {
#endif
            potential.increment(molecules[mI].E());
#if REPULSIVE_CROWDING
        }
#endif
    }
    return potential.total();
}
#endif


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
    uint moleculeNo = (int) gsl_rng_uniform_int(rng, moleculeCount);
    lastMutationIndex = moleculeNo;
    // save the current state so we can roll back if it was not a good mutation.
//     savedMolecule.saveBeforeStateChange(&molecules[moleculeNo]);
    savedMolecule.MC_backup_restore(&molecules[moleculeNo]);
    LOG(DEBUG_MC, "Replica %d\tMolecule %d:\t", label, lastMutationIndex);
    molecules[moleculeNo].make_MC_move(rng, rotateStep, translateStep);
}

void Replica::MCSearchEvaluate()
{
    EonDeviceAsync();
}

void Replica::MCSearchAcceptReject()
{
    //cudaStreamSynchronize(cudaStream);  // sync, newPotential needs to have been returned
    newPotential = SumGridResults();
#if FLEXIBLE_LINKS
    newPotential += internal_molecule_E();
#endif

    //cout << "new energy replica[" << temperature << "] = " << newPotential << endl;

#if PERFORM_GPU_AND_CPU_E
    float cpu_e(E());
    float err = abs(cpu_e-newPotential)/abs(cpu_e);
    printf("%24.20f %24.20f %24.20f\n",cpu_e,float(newPotential),err);
#endif

    float delta = (newPotential - potential);  // needs to be in K_bT
    // if the change is bad then discard it.
    if (delta < 0.0f)
    {
        potential = newPotential;
        LOG(DEBUG_MC, "* Replace:\tdelta E = %f;\tE = %f\n", delta, potential);
    }
    else if (gsl_rng_uniform(rng) < exp(-delta*4.184f/(Rgas*temperature)))
    {
        potential = newPotential;
        LOG(DEBUG_MC, "* Replace: delta E = %f;\tE = %f;\tU < %f\n", delta, potential, exp(-delta*4.184f/(Rgas*temperature)));
    }
    else
    {
//         molecules[lastMutationIndex].undoStateChange(&savedMolecule);
        molecules[lastMutationIndex].MC_backup_restore(&savedMolecule);
        LOG(DEBUG_MC, "- Reject:\tdelta E = %f;\tE = %f\n", delta, potential);
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
    LOG(DEBUG_CUDA, "block size: %d, padded size: %d\n", blockSize, paddedSize);

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
    //TODO add new timer for this
#if FLEXIBLE_LINKS
    result += internal_molecule_E();
#endif
    return result;
}

double Replica::EonDeviceNC()
{
    double result(0.0);
    int NCdataSetSize(ceil(float(nonCrowderResidues)/float(blockSize))*blockSize);
    //cout << "using " << NCdataSetSize <<  " as non NC size" << endl;

    CUDA_EonDeviceNC(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, &result,blockSize,NCdataSetSize,sharedMemSize);

    //TODO add new timer for this
#if FLEXIBLE_LINKS
    result += internal_molecule_E();
#endif

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

// TODO: it's completely unnecessary to pass data in here; we only use the bound conformations. Eliminating this eliminates the circular dependency between SimulationData and Replica.
// this code is a bit special...
void Replica::sample(SimulationData *data, int current_step, float boundEnergyThreshHold, pthread_mutex_t *writeFileMutex)
{
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
        for (int i=0; i < nonCrowderCount; i++)
        {
            for (int j = i+1; j < nonCrowderCount; j++)
            {
                nonCrowderPotential += E(&molecules[i], &molecules[j]); // CPU failsafe
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
        boundSamples++;
    }

    if (bool_recordAllSamples || nonCrowderPotential < boundEnergyThreshHold)
    {
        pthread_mutex_lock(writeFileMutex);
        fprintf(data->boundConformations, "%d; %0.5f (%0.5f); %0.1f\n", current_step,nonCrowderPotential,potential,temperature);

        for (int a=0; a<moleculeCount; a++)
        {
            // TODO use ternary expression instead of this stupid if
            // TODO: we only use this file now to generate PDBs, so if we write the PDBs straight away we don't need it at all?
            // TODO: we need to omit crowders from the PDB
            if (a<nonCrowderCount)
                fprintf(data->boundConformations,"0 %2d: %f %f %f %f %f %f %f\n", a, molecules[a].rotation.w,molecules[a].rotation.x,molecules[a].rotation.y,molecules[a].rotation.z,molecules[a].center.x,molecules[a].center.y,molecules[a].center.z);
            else
                fprintf(data->boundConformations,"1 %2d: %f %f %f %f %f %f %f\n", a, molecules[a].rotation.w,molecules[a].rotation.x,molecules[a].rotation.y,molecules[a].rotation.z,molecules[a].center.x,molecules[a].center.y,molecules[a].center.z);
        }

        pthread_mutex_unlock(writeFileMutex);
    }
}

// TODO: option to omit crowders
bool Replica::savePDB(const char *filename)
{
    char filenameExt[256];
    char tmp[64];
    for (size_t i=0;i<moleculeCount;i++)
    {
        strcpy (filenameExt,filename);
        sprintf (tmp,"%02d",int(i));
        strcat (filenameExt,tmp);
        strcat (filenameExt,".pdb");
        molecules[i].saveAsPDB(filenameExt);
    }
    return true;
}

// TODO: option to omit crowders
void Replica::saveAsSinglePDB(const char *filename)
{
    FILE * output;
    output = fopen (filename,"w");
    fprintf(output,"REMARK %s \n",filename);
    fprintf(output,"REMARK potential: %0.10f \n",float(potential));
    fprintf(output,"REMARK temperature: %5.1f \n",temperature);

    for (size_t i=0;i<moleculeCount;i++)
    {
        fprintf(output,"REMARK Molecule: %d\n",int(i));
        fprintf(output,"REMARK Rotation relative to input Q(w,x,y,z): %f %f %f %f\n",molecules[i].rotation.w,molecules[i].rotation.x,molecules[i].rotation.y,molecules[i].rotation.z);
        fprintf(output,"REMARK Centroid position P(x,y,z): %f %f %f\n",molecules[i].center.x,molecules[i].center.y,molecules[i].center.z);
    }

    char chainId = 'A';
    int itemcount = 0;
    int lastSeqNo = 0;
    for (size_t m=0;m<moleculeCount;m++)
    {
        size_t i=0;
        while (i<molecules[m].residueCount)
        {
            itemcount++;
            fprintf(output,"ATOM  %5d %4s%C%3s %C%4d%C  %8.3f%8.3f%8.3f%6.2f%6.2f\n",itemcount,"CA",' ',aminoAcids.get(molecules[m].Residues[i].aminoAcidIndex).getSNAME(),chainId,molecules[m].Residues[i].resSeq,' ',molecules[m].Residues[i].position.x,molecules[m].Residues[i].position.y,molecules[m].Residues[i].position.z,1.0f,1.0f);
            lastSeqNo = molecules[m].Residues[i].resSeq;
            i++;
        }
        fprintf(output,"TER   %5d      %3s %C%4d \n",itemcount,aminoAcids.get(molecules[m].Residues[i-1].aminoAcidIndex).getSNAME(),chainId,lastSeqNo);
        chainId++;
        fflush(output);
    }
    fprintf(output,"END \n");
    fflush(output);
    fclose(output);
}