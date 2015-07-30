#include "Replica.h"

using namespace std;

Replica::Replica() : RNGs_initialised(false), nonCrowderPotential(0.0f), temperature(300.0f), label(-1), potential(0.0f), newPotential(0.0f), moleculeCount(0), moleculeArraySize(0), residueCount(0), max_residue_count(0), nonCrowderCount(0), accept(0), acceptA(0), reject(0), totalAcceptReject(0), totalAccept(0), boundSamples(0), samplesSinceLastExchange(0), totalBoundSamples(0), totalSamples(0), paircount(0), nonCrowderResidues(0), translateStep(INITIAL_TRANSLATIONAL_STEP), rotateStep(INITIAL_ROTATIONAL_STEP), contiguousResiduesSize(0),
#if FLEXIBLE_LINKS
    calculate_rigid_potential_only(false),
#endif
#if USING_CUDA
    replicaIsOnDevice(false),
#endif
#if INCLUDE_TIMERS
    timersInit(false),
#endif
    acceptanceRatio(0.0f), accumulativeAcceptanceRatio(0.0f), fractionBound(0.0f), accumulativeFractionBound(0.0f)
{
}

void Replica::init_first_replica(const argdata parameters, AminoAcids amino_acid_data, const int initial_molecule_array_size) // constructor for initial replica
{
    aminoAcids = amino_acid_data;
    boundingValue = parameters.bound;
    reserveContiguousMoleculeArray(initial_molecule_array_size);
    
    for (size_t s = 0; s < parameters.mdata.size(); s++)
    {
        loadMolecule(parameters.mdata[s]);
    }

    // ASSIGN UIDS (mostly needed for CUDA)
    
    int chain_offset = 1;
    int domain_offset = 1;
    int bond_offset = 1;

    for (int m = 0; m < moleculeCount; m++) {
        molecules[m].graph.assign_uids(molecules[m].Residues, chain_offset, domain_offset, bond_offset);
    }

    // COUNT PAIRS FOR NONBONDED POTENTIAL
    
    paircount = residueCount * (residueCount - 1) / 2; // start with the handshake algorithm
    for (size_t m = 0; m < moleculeCount; m++) {
        Molecule & mol =  molecules[m];
#if FLEXIBLE_LINKS
        paircount -= 3 * mol.residueCount; // exclude residues 3 or fewer apart on the backbone in each molecule
        paircount -= mol.graph.segment_bonds.size(); // exclude pairs which form a bond (e.g. between chains)
        paircount -= mol.graph.indirect_neighbours.size(); // exclude residues 3 or fewer apart along paths leading through bonds
        // exclude all pairs within the same rigid domain
        for (int d = 0; d < mol.graph.rigid_domains.size(); d++) {
            int d_size = mol.graph.rigid_domains[d].size();
            paircount -= d_size * (d_size - 1) / 2;
        }
#else
        paircount -= mol.residueCount * (mol.residueCount - 1) / 2; // exclude all pairs within the same molecule
#endif
    }
    
#if USING_CUDA
    if (parameters.auto_blockdim)
    {
        (residueCount < 1024) ? setBlockSize(32) : setBlockSize(64);
        VLOG(1) << "\tAutomatically calculated block dimension: " << blockSize;
    }
    else
    {
        setBlockSize(parameters.cuda_blockSize);
    }
#endif
    
    VLOG(1) << "\tLoaded: " << residueCount << " residues in " << moleculeCount << " molecules";
    for (int i = 0; i < moleculeCount; i++)
    {
        molecules[i].log_info(i);
    }
    VLOG(1) << "\tCounted : " << nonCrowderResidues << " complex residues and " << paircount << " residue interaction pairs.";
}

void Replica::init_child_replica(const Replica& ir, const int index, const double geometricTemperature, const double geometricRotation, const double geometricTranslate, const argdata parameters)
{
    // TODO: incorporate copy method into this one and try to clean it up
    // TODO: if it's not used anywhere else, remove setting of index and temperature -- it only confuses things
    copy(ir);

    label = index;
    temperature = parameters.temperatureMin * pow(geometricTemperature, index);
    rotateStep = MIN_ROTATION * pow(geometricRotation, index);
    translateStep = MIN_TRANSLATION * pow(geometricTranslate, index);
    threadCount = parameters.threads;

#if INCLUDE_TIMERS
    initTimers();
#endif
    // TODO: only one range for entire simulation (but is the library threadsafe?)
    initRNGs();
    savedMolecule.init_saved_molecule(max_residue_count);
    // leaving out stuff inside thread function for now
}

#if USING_CUDA

#if CUDA_STREAMS
void Replica::setup_CUDA(float * device_LJ_potentials, cudaStream_t * streams, int stream_index)
#else
void Replica::setup_CUDA(float * device_LJ_potentials)
#endif // CUDA_STREAMS
{
    device_LJPotentials = device_LJ_potentials;
#if CUDA_STREAMS
    cudaStream = streams[stream_index];  // use rotation to assign replicas to streams
#endif // CUDA_STREAMS
    ReplicaDataToDevice();
#if CUDA_STREAMS
    ReserveSumSpace(); // reserve a space for the potential summation to be stored
#endif // CUDA_STREAMS
}

void Replica::teardown_CUDA()
{
#if CUDA_STREAMS
    FreeSumSpace();
    getLastCudaError("Error freeing sum space");
#endif // CUDA_STREAMS
    FreeDevice();
    getLastCudaError("Error freeing device");
}

#endif // USING_CUDA

#if INCLUDE_TIMERS
void Replica::initTimers()
{
    //timers for profiling the cuda functions
    //replicaToGPUTimer = 0;
    //replicaToHostTimer = 0;
    //replicaUpdateGPUTimer = 0;
    //replicaMoleculeUpdateTimer = 0;
    //replicaECUDATimer = 0;
    //initGPUMemoryTimer = 0;
    //replicaEHostTimer = 0;
    replicaToGPUTimer = NULL;
    replicaToHostTimer = NULL;
    replicaUpdateGPUTimer = NULL;
    replicaMoleculeUpdateTimer = NULL;
    replicaECUDATimer = NULL;
    initGPUMemoryTimer = NULL;
    replicaEHostTimer = NULL;

    //CUT_SAFE_CALL( cutCreateTimer(&replicaToGPUTimer) );
    //CUT_SAFE_CALL( cutCreateTimer(&replicaToHostTimer) );
    //CUT_SAFE_CALL( cutCreateTimer(&replicaUpdateGPUTimer) );
    //CUT_SAFE_CALL( cutCreateTimer(&replicaMoleculeUpdateTimer) );
    //CUT_SAFE_CALL( cutCreateTimer(&replicaECUDATimer) );
    //CUT_SAFE_CALL( cutCreateTimer(&initGPUMemoryTimer) );
    //CUT_SAFE_CALL( cutCreateTimer(&replicaEHostTimer) );

    sdkCreateTimer(&replicaToGPUTimer);
    sdkCreateTimer(&replicaToHostTimer);
    sdkCreateTimer(&replicaUpdateGPUTimer);
    sdkCreateTimer(&replicaMoleculeUpdateTimer);
    sdkCreateTimer(&replicaECUDATimer);
    sdkCreateTimer(&initGPUMemoryTimer);
    sdkCreateTimer(&replicaEHostTimer);

    timersInit = true;
}
#endif

Replica::~Replica()
{
#if INCLUDE_TIMERS
    if (timersInit)
    {
        //CUT_SAFE_CALL( cutDeleteTimer(replicaToGPUTimer) );
        //CUT_SAFE_CALL( cutDeleteTimer(replicaToHostTimer) );
        //CUT_SAFE_CALL( cutDeleteTimer(replicaUpdateGPUTimer) );
        //CUT_SAFE_CALL( cutDeleteTimer(replicaMoleculeUpdateTimer) );
        //CUT_SAFE_CALL( cutDeleteTimer(replicaECUDATimer) );
        //CUT_SAFE_CALL( cutDeleteTimer(initGPUMemoryTimer) );
        //CUT_SAFE_CALL( cutDeleteTimer(replicaEHostTimer) )
        sdkDeleteTimer(&replicaToGPUTimer);
        sdkDeleteTimer(&replicaToHostTimer);
        sdkDeleteTimer(&replicaUpdateGPUTimer);
        sdkDeleteTimer(&replicaMoleculeUpdateTimer);
        sdkDeleteTimer(&replicaECUDATimer);
        sdkDeleteTimer(&initGPUMemoryTimer);
        sdkDeleteTimer(&replicaEHostTimer);
    }
#endif

    if (RNGs_initialised) {
        freeRNGs();
    }

    if (moleculeCount > 0)
    {
        delete [] molecules;
    }

    if (contiguousResiduesSize > 0)
    {
        delete [] contiguousResidues;
    }
}

void Replica::exchangeReplicas(Replica &r)
{
    swap(label, r.label);
    swap(temperature, r.temperature);

    // for statistics
    
    swap(accept, r.accept);
    swap(acceptA, r.acceptA);
    swap(reject, r.reject);
    
    swap(totalAcceptReject, r.totalAcceptReject);
    swap(totalAccept, r.totalAccept);
    swap(acceptanceRatio, r.acceptanceRatio);
    swap(accumulativeAcceptanceRatio, r.accumulativeAcceptanceRatio);

    swap(boundSamples, r.boundSamples);
    swap(samplesSinceLastExchange, r.samplesSinceLastExchange);
    swap(totalBoundSamples, r.totalBoundSamples);
    swap(totalSamples, r.totalSamples);
    swap(fractionBound, r.fractionBound);
    swap(accumulativeFractionBound, r.accumulativeFractionBound);
}

void Replica::copy(const Replica &r)
{
    label = r.label;
    temperature = r.temperature;
    aminoAcids = r.aminoAcids;
    residueCount = r.residueCount;
    max_residue_count = r.max_residue_count;
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
        LOG(ERROR) << "Cannot reserve molecule size: molecules already exist.";
        return;
    }
    moleculeArraySize = size;
    molecules = new Molecule[size];
}

void Replica::loadMolecule(const moldata mol)
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

    molecules[moleculeCount].init(mol, aminoAcids, moleculeCount, boundingValue);

    residueCount += molecules[moleculeCount].residueCount;
    max_residue_count = max(max_residue_count, molecules[moleculeCount].residueCount);

    if (!mol.crowder) {
        nonCrowderCount++;
        nonCrowderResidues += molecules[moleculeCount].residueCount;
    }

    moleculeCount++;
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
    MC_move_weights[2] = WEIGHT_MC_FLEX;
    MC_move_weights[3] = WEIGHT_MC_LOCAL;

    MC_discrete_table = gsl_ran_discrete_preproc(4, MC_move_weights);
    for (size_t m = 0; m < moleculeCount; m++)
    {
        // TODO this is hacky, but it will be gone when this is per-simulation
        molecules[m].MC_discrete_table = MC_discrete_table;
    }
#endif
    RNGs_initialised = true;
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

void Replica::MCSearch(int steps, int mcstep)
{
    for (int step = 0; step < steps; step++)
    {
        ostringstream debug_log;
        int moleculeNo = (int) gsl_rng_uniform_int(rng, moleculeCount);
        // save the current state so we can roll back if it was not a good mutation.
        savedMolecule.MC_backup_restore(&molecules[moleculeNo]);
        molecules[moleculeNo].make_MC_move(rng, rotateStep, translateStep);

        debug_log << "Step " << mcstep + step << ":\treplica " << label << "\tmolecule " << moleculeNo << ":\t" << molecules[moleculeNo].last_MC_move << "\t";
#if USING_CUDA
        // copy host data to device. so we can do the calculations on it.
        MoleculeDataToDevice(moleculeNo);
        newPotential = EonDevice();  // sequential cuda call
#if FLEXIBLE_LINKS
        if (!calculate_rigid_potential_only) {
            double bonded_potential = internal_molecule_E(false).total();
            debug_log << "new Eu: " << newPotential << ",\tnew Eb: " << bonded_potential << "\t";
            newPotential += bonded_potential;
        }
#endif // FLEXIBLE_LINKS
#else // only CPU calls
        Potential cpu_e = E();
#if FLEXIBLE_LINKS
        if (!calculate_rigid_potential_only) {
            cpu_e += internal_molecule_E(true);
            debug_log << "new Eu: " << cpu_e.total_LJ() + cpu_e.total_DH() << ",\tnew Eb: " << cpu_e.total_bond() + cpu_e.total_angle() + cpu_e.total_torsion() << "\t";
        }
#endif // FLEXIBLE_LINKS
        newPotential = cpu_e.total();
#endif // not USING_CUDA

        debug_log << "new E: " << newPotential << "\t";
    
        float delta = newPotential - potential;

        // accept change if its better.
        if (delta < 0.0f)
        {
            potential = newPotential;
            accept++;
            debug_log << "* Replace:\tdelta E = " << delta << ";\tE = " << potential;
        }
        // accept change if it meets the boltzmann criteria -- delta must be converted from kcal/mol to J/mol
        else if (gsl_rng_uniform(rng) < exp(-(delta*kcal)/(Rgas*temperature)))
        {
            potential = newPotential;
            acceptA++;
            debug_log << "**Replace:\tdelta E = " << delta << ";\tE = " << potential << ";\tU < " << exp(-delta * kcal/(Rgas*temperature));
        }
        //reject
        else
        {
            reject++;
            debug_log << "- Reject:\tdelta E = " << delta << ";\tE = " << potential;
            molecules[moleculeNo].MC_backup_restore(&savedMolecule);
#if USING_CUDA
            MoleculeDataToDevice(moleculeNo); // you have to update the device again because the copy will be inconsistent
#endif
        }

        // TODO: this is probably wrong
        DLOG(INFO) << debug_log.str();

    }
}

inline float crowderPairPotential(const float r)
{
    return powf(6.0f/r,12.0f);
}

Potential Replica::E()
{
#if INCLUDE_TIMERS
    //CUT_SAFE_CALL(cutStartTimer(replicaEHostTimer));
    sdkStartTimer(&replicaEHostTimer);
#endif

    Potential potential;

#define iRes molecules[mI].Residues[mi]
#define jRes molecules[mJ].Residues[mj]

    for (size_t mI = 0; mI < moleculeCount; mI++)
    {
#if REPULSIVE_CROWDING && !LJ_OFF
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
            for (size_t mJ = mI + 1; mJ < moleculeCount; mJ++)
            {
#if REPULSIVE_CROWDING && !LJ_OFF
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
#if !LJ_OFF
                            potential.increment_LJ(calculate_LJ(iRes, jRes, r, aminoAcids));
#endif
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

#if INCLUDE_TIMERS
    //CUT_SAFE_CALL(cutStopTimer(replicaEHostTimer));
    sdkStopTimer(&replicaEHostTimer);
#endif
    last_potential_obj = potential;
    return potential;
}

Potential Replica::E(Molecule *a, Molecule *b)
{
    Potential potential;

#define aRes a->Residues[mi]
#define bRes b->Residues[mj]

    for (size_t mi = 0; mi < a->residueCount; mi++)
    {
        for (size_t mj = 0; mj < b->residueCount; mj++)
        {
            double r(aRes.distance(bRes, boundingValue) + EPS);
#if !LJ_OFF
            potential.increment_LJ(calculate_LJ(aRes, bRes, r, aminoAcids));
#endif
            potential.increment_DH(calculate_DH(aRes, bRes, r));
        }
    }
#if FLEXIBLE_LINKS
    if (!calculate_rigid_potential_only) {
        potential += a->E(true);
        potential += b->E(true);
    }
#endif
    return potential;
}


#if FLEXIBLE_LINKS
Potential Replica::internal_molecule_E(bool include_LJ_and_DH) {
    Potential potential;

    for (size_t mI = 0; mI < moleculeCount; mI++)
    {
#if REPULSIVE_CROWDING
        if (molecules[mI].moleculeRoleIdentifier != CROWDER_IDENTIFIER)
        {
#endif
            potential += molecules[mI].E(include_LJ_and_DH);
#if REPULSIVE_CROWDING
        }
#endif
    }

    last_potential_obj = potential;
    return potential;
}
#endif


void Replica::printTimers()
{
#if INCLUDE_TIMERS
    LOG(INFO) << "Timer values (ms)";
    LOG(INFO) << "Total(ms)   \tAverage(ms)  \t Action";
    LOG(INFO) << sdkGetTimerValue(&replicaToGPUTimer) << sdkGetAverageTimerValue(&replicaToGPUTimer) << "Replica to GPU (init & transfer)";
    LOG(INFO) << sdkGetTimerValue(&replicaUpdateGPUTimer) << sdkGetAverageTimerValue(&replicaUpdateGPUTimer) << "Replica Update on GPU (transter)";
    LOG(INFO) << sdkGetTimerValue(&replicaECUDATimer) << sdkGetAverageTimerValue(&replicaECUDATimer) << "Kernel Timer (computation)";
    LOG(INFO) << sdkGetTimerValue(&replicaEHostTimer) << sdkGetAverageTimerValue(&replicaEHostTimer) << "Host Timer (computation)";
    LOG(INFO) << sdkGetTimerValue(&replicaMoleculeUpdateTimer) << sdkGetAverageTimerValue(&replicaMoleculeUpdateTimer) << "Update Molecule on GPU (transfer)";
    LOG(INFO) << sdkGetTimerValue(&initGPUMemoryTimer) << sdkGetAverageTimerValue(&initGPUMemoryTimer) << "GPU Memory Initialisation (malloc)";
    LOG(INFO) << "Kernel Speedup: " << sdkGetAverageTimerValue(&replicaEHostTimer)/sdkGetAverageTimerValue(&replicaECUDATimer);
#else
    LOG(INFO) << "Timers disabled: set INCLUDE_TIMERS 1 and USING_CUDA 1 in definitions.h";
#endif
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
        LOG(ERROR) << "! Error: replica not on device implies paddedSize==0; cannot perform Replica::ReserveSumSpace()";
    }
    // result stored on the device
    // gridSize can be arbitrary
    gridSize = paddedSize/blockSize;
    resultSize = gridSize;
    //for a parallel sum each grid must have one cell in the array of results from all the threads
    cudaMallocHost((void**)&kernelResult, sizeof(float)*resultSize*resultSize);
    getLastCudaError("Error allocating sum space on host");
    cudaMalloc((void **)&device_kernelResult,sizeof(float)*resultSize*resultSize);
    getLastCudaError("Error allocating sum space on device");

//     cudaError_t err;
//     if ((err = cudaGetLastError()) != cudaSuccess)
//         printf("CUDA error: %s\n", cudaGetErrorString(err));

    memset(kernelResult,0,sizeof(float)*resultSize*resultSize);
    getLastCudaError("Error zeroing sum space on host");
    //CUDA_memcpy_to_device_async(device_kernelResult,kernelResult,sizeof(float)*resultSize*resultSize,cudaStream);
    //CUDA_memcpy_to_device(device_kernelResult,kernelResult,sizeof(float)*resultSize*resultSize);
    cudaMemset(device_kernelResult, 0, sizeof(float)*resultSize*resultSize);
    getLastCudaError("Error zeroing sum space on device");
}

void Replica::FreeSumSpace()
{
    cudaFreeHost(kernelResult);
    getLastCudaError("Error freeing sum space on host");
    cudaFree(device_kernelResult);
    getLastCudaError("Error freeing sum space on device");
}

double Replica::SumGridResults()
{
    double potentialSum = 0.0f;
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
    getLastCudaError("Error zeroing sum space on device after summing");
    return potentialSum;
}
#endif

#if CUDA_STREAMS

// 1/3 of the above function, does the mutation on the gpu asynchronously
void Replica::MCSearchMutate(int mcstep)
{
    int moleculeNo = (int) gsl_rng_uniform_int(rng, moleculeCount);
    // save the current state so we can roll back if it was not a good mutation.
    savedMolecule.MC_backup_restore(&molecules[moleculeNo]);
    molecules[moleculeNo].make_MC_move(rng, rotateStep, translateStep);

    lastMutationIndex = moleculeNo;
}

void Replica::MCSearchEvaluate(int mcstep)
{
    MoleculeDataToDevice(lastMutationIndex);
    EonDeviceAsync();
}

void Replica::MCSearchAcceptReject(int mcstep)
{
    ostringstream debug_log;
    debug_log << "Step " << mcstep << ":\treplica " << label << "\tmolecule " << lastMutationIndex << ":\t" << molecules[lastMutationIndex].last_MC_move << "\t";

    newPotential = SumGridResults();

#if FLEXIBLE_LINKS
    if (!calculate_rigid_potential_only) {
        double bonded_potential = internal_molecule_E(false).total();
        debug_log << "new Eu: " << newPotential << ",\tnew Eb: " << bonded_potential << "\t";
        newPotential += bonded_potential;
    }
#endif

    debug_log << "new E: " << newPotential << "\t";

    float delta = (newPotential - potential);
    
    if (delta < 0.0f)
    {
        potential = newPotential;
        accept++;
        debug_log << "* Replace:\tdelta E = " << delta << ";\tE = " << potential;
    }
    // accept change if it meets the boltzmann criteria -- delta must be converted from kcal/mol to J
    else if (gsl_rng_uniform(rng) < exp(-delta*kcal/(Rgas*temperature)))
    {
        potential = newPotential;
        acceptA++;
        debug_log << "**Replace:\tdelta E = " << delta << ";\tE = " << potential << ";\tU < " << exp(-delta*kcal/(Rgas*temperature));
    }
    else
    // if the change is bad then discard it.
    {
        reject++;
        debug_log << "- Reject:\tdelta E = " << delta << ";\tE = " << potential;
        molecules[lastMutationIndex].MC_backup_restore(&savedMolecule);
        MoleculeDataToDevice(lastMutationIndex); // you have to update the device again because the copy will be inconsistent
    }
    
    DLOG(INFO) << debug_log.str();
}
#endif  // streams

// copy the replica to the device
void Replica::ReplicaDataToDevice()
{
#if INCLUDE_TIMERS
    //CUT_SAFE_CALL(cutStartTimer(initGPUMemoryTimer));
    sdkStartTimer(&initGPUMemoryTimer);
#endif

    // allocate total size needed.

    paddedSize = int(ceil(float(residueCount)/float(blockSize)))*blockSize; // reserve a blocksize multiple because it allows for efficient summation
    dataSetSize = paddedSize;
    DLOG(INFO) << "block size: " << blockSize << ", padded size: " << paddedSize;

#if CUDA_STREAMS
    cudaMallocHost((void**)&host_float4_residuePositions,sizeof(float4)*paddedSize); // pinned memory?
#else
    host_float4_residuePositions = new float4[paddedSize];
#endif

    //host_float4_residuePositions = new float4[paddedSize];
    host_float4_residueMeta = new float4[paddedSize]; // WTF? Why is this not different for streams?
    host_moleculeStartPositions = new int[moleculeCount];

    cudaMalloc((void**)&device_float4_residuePositions,sizeof(float4)*paddedSize);
    cudaMalloc((void**)&device_float4_residueMeta,sizeof(float4)*paddedSize);
    cudaMalloc((void**)&device_moleculeStartPositions,sizeof(int)*moleculeCount);
    cudaMalloc((void**)&device_residueCount,sizeof(int));
    cudaMalloc((void**)&device_moleculeCount,sizeof(int));

#if INCLUDE_TIMERS
    //CUT_SAFE_CALL(cutStopTimer(initGPUMemoryTimer));
    sdkStopTimer(&initGPUMemoryTimer);
    //CUT_SAFE_CALL(cutStartTimer(replicaToGPUTimer));
    sdkStartTimer(&replicaToGPUTimer);
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
            // POS
            host_float4_residuePositions[arrayIndex].x = molecules[m].Residues[rc].position.x;
            host_float4_residuePositions[arrayIndex].y = molecules[m].Residues[rc].position.y;
            host_float4_residuePositions[arrayIndex].z = molecules[m].Residues[rc].position.z;
#if FLEXIBLE_LINKS
            host_float4_residuePositions[arrayIndex].w = molecules[m].Residues[rc].pos_w; // DOMAIN_UID.BOND_UID
            if (calculate_rigid_potential_only) {
                // override the domain part, to put all residues in one molecule in the same rigid domain
                host_float4_residuePositions[arrayIndex].w -= int(host_float4_residuePositions[arrayIndex].w); // subtract the current integer part
                host_float4_residuePositions[arrayIndex].w += mf + 1.0f; // add the molecule ID + 1 (zero means "not in a rigid domain")
            }
#else
            host_float4_residuePositions[arrayIndex].w = mf;  // we only care if the residues are in the same molecule
#endif
            
            // META
            host_float4_residueMeta[arrayIndex].x = molecules[m].Residues[rc].aminoAcidIndex;
            host_float4_residueMeta[arrayIndex].y = molecules[m].Residues[rc].electrostaticCharge;
            host_float4_residueMeta[arrayIndex].z = molecules[m].Residues[rc].vanderWaalRadius;
            if (molecules[m].moleculeRoleIdentifier == CROWDER_IDENTIFIER)
            {
                // This molecule is a crowder
                host_float4_residueMeta[arrayIndex].w = CROWDER_IDENTIFIER;
            }
            else
            {
#if FLEXIBLE_LINKS
                host_float4_residueMeta[arrayIndex].w = molecules[m].Residues[rc].meta_w; // RESIDUE_ID.CHAIN_UID
#else
                host_float4_residueMeta[arrayIndex].w = 0.0f; // we don't use this 
#endif
            }

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
#if CUDA_STREAMS
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
    //CUT_SAFE_CALL( cutStopTimer(replicaToGPUTimer) );
    sdkStopTimer(&replicaToGPUTimer);
#endif

    replicaIsOnDevice = true;

    return;
}

// update a molecule on the device
void Replica::MoleculeDataToDevice(int moleculeId)
{
    if (!replicaIsOnDevice)
    {
        LOG(ERROR) << "ERROR: Replica::MoleculeDataToDevice("<< moleculeId << ") called without initialising device data.";
        ReplicaDataToDevice();
    }
#if INCLUDE_TIMERS
    //CUT_SAFE_CALL( cutStartTimer(replicaMoleculeUpdateTimer) );
    sdkStartTimer(&replicaMoleculeUpdateTimer);
#endif

    // TODO: change to use float4 for everything, will eliminate this copy

    /*
    1. get molecule on device's residues pointer
    2. make residue array
    3. copy to device at location found in 1
    */

    int residueIndex = host_moleculeStartPositions[moleculeId];
    int memoryPosition = residueIndex;
    int moleculeSize = molecules[moleculeId].residueCount;
    for (int rc=0; rc<molecules[moleculeId].residueCount; rc++)
    {
        host_float4_residuePositions[residueIndex].x = molecules[moleculeId].Residues[rc].position.x;
        host_float4_residuePositions[residueIndex].y = molecules[moleculeId].Residues[rc].position.y;
        host_float4_residuePositions[residueIndex].z = molecules[moleculeId].Residues[rc].position.z;
        
        // Only the positions change; nothing else should have to be updated.
        // If a unit test needs to change e.g. calculate_rigid_potential_only while reusing the same replica, it should call ReplicaDataToDevice explicitly afterwards.

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
    //CUT_SAFE_CALL( cutStopTimer(replicaMoleculeUpdateTimer) );
    sdkStopTimer(&replicaMoleculeUpdateTimer);
#endif
    return;
}

void Replica::EonDeviceAsync()
{
#if INCLUDE_TIMERS
    //CUT_SAFE_CALL( cutStartTimer(replicaECUDATimer) );
    sdkStartTimer(&replicaECUDATimer);
#endif
#if CUDA_STREAMS
    // compute potential parts
    CUDA_EonDevice_async(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, device_kernelResult, resultSize, blockSize, dataSetSize, sharedMemSize, cudaStream);
    //CUDA_EonDevice(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, device_kernelResult,blockSize,dataSetSize);

    // write back to CPU
    //CUDA_memcpy_to_host_async(kernelResult,device_kernelResult,resultSize*resultSize*sizeof(float),cudaStream);
    cudaMemcpyAsync(kernelResult,device_kernelResult,resultSize*resultSize*sizeof(float),cudaMemcpyDeviceToHost,cudaStream);

#else
    LOG(ERROR) << "Replica::EonDeviceAsync() can only be run using streams.";
#endif
#if INCLUDE_TIMERS
    //CUT_SAFE_CALL( cutStopTimer(replicaECUDATimer) );
    sdkStopTimer(&replicaECUDATimer);
#endif
}

double Replica::EonDevice()
{
#if INCLUDE_TIMERS
    //CUT_SAFE_CALL( cutStartTimer(replicaECUDATimer) );
    sdkStartTimer(&replicaECUDATimer);
#endif

    double result(0.0);

    CUDA_EonDevice(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, &result,blockSize,dataSetSize,sharedMemSize);

#if INCLUDE_TIMERS
    //CUT_SAFE_CALL( cutStopTimer(replicaECUDATimer) );
    sdkStopTimer(&replicaECUDATimer);
#endif

    return result;
}

double Replica::EonDeviceNC()
{
    double result(0.0);
    int NCdataSetSize(ceil(float(nonCrowderResidues)/float(blockSize))*blockSize);

    CUDA_EonDeviceNC(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, &result,blockSize,NCdataSetSize,sharedMemSize);

    //TODO add new timer for this
#if FLEXIBLE_LINKS
    if (!calculate_rigid_potential_only) {
        result += internal_molecule_E(false).total();
    }
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

//#if CUDA_MC
    //cudaFree(device_moleculeCenters);
    //cudaFree(device_moleculeLengths);
    //cudaFree(device_translationVector);
    //cudaFree(device_reverseTranslationVector);
    //cudaFree(device_rotationVector);
    //cudaFree(device_reverseRotationVector);
    //cudaFreeHost(host_translationVector);
    //cudaFreeHost(host_reverseTranslationVector);
    //cudaFreeHost(host_rotationVector);
    //cudaFreeHost(host_reverseRotationVector);
//#endif
}

//#if CUDA_MC
//// rotate the molecule on device about a vector and amount
//bool Replica::rotateOnDevice(int moleculeId, Vector3f rvector, float amount)
//{
    ////device_rotationVector;  // vector(x,y,z)|amount(w)

    //host_rotationVector->x = rvector.x;
    //host_rotationVector->y = rvector.y;
    //host_rotationVector->z = rvector.z;
    //host_rotationVector->w = amount;
    //CUDA_memcpy_to_device_async(device_rotationVector,&host_rotationVector,sizeof(float4),cudaStream);

    //CUDA_rotateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[moleculeId], &device_moleculeLengths[moleculeId], molecules[moleculeId].length, device_rotationVector, &device_moleculeCenters[moleculeId], cudaStream);

    //// create the undo vector while its busy
    //host_reverseRotationVector->x = rvector.x;
    //host_reverseRotationVector->y = rvector.y;
    //host_reverseRotationVector->z = rvector.z;
    //host_reverseRotationVector->w = -amount;
    //CUDA_memcpy_to_device_async(device_reverseTranslationVector,&host_reverseTranslationVector,sizeof(float4),cudaStream);


    //return true;
//}

//bool Replica::translateOnDevice(int moleculeId, Vector3f translation)
//{
    //// send the vector to the device
    //host_translationVector->x = translation.x;
    //host_translationVector->y = translation.y;
    //host_translationVector->z = translation.z;
    //CUDA_memcpy_to_device_async(device_translationVector,&host_translationVector,sizeof(float4),cudaStream);

////  CUDA_translateMolecule (float4 *residuePositions, int *startPosition, int *moleculeLength, int moleculeLength, int *moleculeId, float4* translation, cudaStream_t stream)

    //CUDA_translateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[moleculeId], &device_moleculeLengths[moleculeId], molecules[moleculeId].length, device_translationVector, &device_moleculeCenters[moleculeId], cudaStream);

    //// create the undo vector while its busy
    //host_reverseTranslationVector->x = -translation.x;
    //host_reverseTranslationVector->y = -translation.y;
    //host_reverseTranslationVector->z = -translation.z;
    //CUDA_memcpy_to_device_async(device_reverseTranslationVector,&host_reverseTranslationVector,sizeof(float4),cudaStream);

    //return true;
//}

//void Replica::cudaRollbackMutation()
//{
    //if (lastMutationWasATranslate)
    //{
        //CUDA_translateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[lastMutatedMolecule], &device_moleculeLengths[lastMutatedMolecule], molecules[lastMutatedMolecule].length, device_reverseTranslationVector,&device_moleculeCenters[lastMutatedMolecule], cudaStream);
    //}
    //else
    //{
        //CUDA_rotateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[lastMutatedMolecule], &device_moleculeLengths[lastMutatedMolecule], molecules[lastMutatedMolecule].length, device_reverseRotationVector,&device_moleculeCenters[lastMutatedMolecule], cudaStream);
    //}
//}
//#endif // CUDA MC
#endif  // CUDA

// this code is a bit special...
void Replica::sample(SimulationData * data, int current_step, float boundEnergyThreshHold)
{
    nonCrowderPotential = 0.0f;
    samplesSinceLastExchange++;

    // in the case of all molecules being of interest then just use the potential
    if (moleculeCount != nonCrowderCount)
    {
#if USING_CUDA
        nonCrowderPotential = EonDeviceNC(); // do on GPU, no stream support, can be implemented in 3mins if stream samlping is fixed.
        // TODO: arrrgh, what does this *mean*?
#else
        // if there are crowders and molecules of interest the only use the energy of the interesting ones
        // TODO: make absolutely sure all the NC molecules are at the front of the list
        for (int i = 0; i < nonCrowderCount; i++)
        {
            for (int j = i + 1; j < nonCrowderCount; j++)
            {
                nonCrowderPotential += E(&molecules[i], &molecules[j]).total(); // CPU failsafe
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
        char savename[256];
        memset(savename,0,256);
        sprintf(savename, "output/%s/pdb/sample_%d_%0.1fK_%5.2f.pdb", data->prefix, current_step, temperature, nonCrowderPotential);
        saveAsSinglePDB(savename, data->prefix, current_step);

        pthread_mutex_lock(data->writeFileMutex);
        fprintf(data->boundConformations, "%d; %0.5f (%0.5f); %0.1f\n%s\n", current_step, nonCrowderPotential, potential, temperature, savename);
        pthread_mutex_unlock(data->writeFileMutex);
    }
}

void Replica::fraction_bound(FILE * fractionBoundFile)
{
    totalBoundSamples += boundSamples;
    totalSamples += samplesSinceLastExchange;

    fractionBound = float(boundSamples) / max(1.0f, float(samplesSinceLastExchange));
    accumulativeFractionBound = float(totalBoundSamples)/max(1.0f,float(totalSamples));

    fprintf(fractionBoundFile,"| %6.4f %6.4f ", fractionBound, accumulativeFractionBound);

    boundSamples = 0;
    samplesSinceLastExchange = 0;
}

void Replica::acceptance_ratio(FILE * acceptanceRatioFile)
{
    totalAccept += acceptA + accept;
    totalAcceptReject += acceptA + accept + reject;

    acceptanceRatio = float(acceptA + accept) / float(acceptA + accept + reject);
    accumulativeAcceptanceRatio = float(totalAccept) / float(totalAcceptReject);

    fprintf(acceptanceRatioFile,"| %6.4f %6.4f ", acceptanceRatio, accumulativeAcceptanceRatio);

    acceptA = 0;
    accept = 0;
    reject = 0;
}

void Replica::saveAsSinglePDB(const char *filename, const char *prefix, int sample, bool skip_crowders)
{
    FILE * output;
    output = fopen (filename,"w");
    fprintf(output,"REMARK simulation: %s \n", prefix);
    fprintf(output,"REMARK temperature: %5.1f \n", temperature);
    fprintf(output,"REMARK sample: %d \n", sample);
    fprintf(output,"REMARK potential: %0.10f \n", float(potential));

    char out_chain_id = 0;
    int itemcount = 0;
    int lastSeqNo = 0;
    char last_chain = 0;
    size_t i;
    for (size_t m = 0; m < moleculeCount; m++)
    {
        if (skip_crowders && m >= nonCrowderCount)
        {
            break;
        }

        for (i = 0; i < molecules[m].residueCount; i++)
        {
            Residue r = molecules[m].Residues[i];

            out_chain_id = 64 + r.chain_UID;

            if (i > 0 && out_chain_id != last_chain) {
                fprintf(output,"TER   %5d      %3s %C%4d \n", itemcount, aminoAcids.get(molecules[m].Residues[i-1].aminoAcidIndex).getSNAME(), last_chain, lastSeqNo);
            }
            
            itemcount++;

            fprintf(output,"ATOM  %5d CA   %3s %C%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", itemcount, aminoAcids.get(r.aminoAcidIndex).getSNAME(), out_chain_id, r.resSeq, r.position.x, r.position.y, r.position.z, 0.0f, 0.0f);

            lastSeqNo = r.resSeq;
            last_chain = out_chain_id;
        }
        fprintf(output,"TER   %5d      %3s %C%4d \n", itemcount, aminoAcids.get(molecules[m].Residues[i-1].aminoAcidIndex).getSNAME(), out_chain_id, lastSeqNo);
        fflush(output);
    }
    fprintf(output,"END \n");
    fflush(output);
    fclose(output);
}
