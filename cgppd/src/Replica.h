#ifndef REPLICA_H_
#define REPLICA_H_
#include "vector3f.h"
#include "definitions.h"
#include "structures.h"
#include "Quaternion.h"
#include <fstream>
#include <vector>
#include "Molecule.h"
#include "Potential.h"
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "TorsionalLookupMatrix.h"
#include <pthread.h>
//#include "rngfunctions.h"


using namespace std;

#if USING_CUDA

#include <cutil.h>
#include <cuda.h>
#include <builtin_types.h>

#include "cudaExterns.h"

#endif  // USING_CUDA

class Replica
{
public:
    Replica();
    ~Replica();

    void init_first_replica(const vector<moldata> mdata, AminoAcids amino_acid_data, const float bounding_value, const int initial_molecule_array_size);
    void init_child_replica(const Replica& ir, const int label, const float temperature, const float rotate_step, const float translate_step, const int thread_count);

    void reserveContiguousMoleculeArray(int size);
    void initRNGs();
    void freeRNGs();

    // TODO: remove this once it's incorporated into final replica constructor
    void copy(const Replica& r); 						// assignment copy, duplicate it.

    void exchangeReplicas(Replica &r);

    double E();
    double E(Molecule *a,Molecule *b);   // use for fraction bound calcs

    // the search method that mutates the replicas molecules.
    void MCSearch(int steps);//, Replica * replaceReplica);

    int loadMolecule(const char *pdbfilename);

    AminoAcids aminoAcids;

#if FLEXIBLE_LINKS
    TorsionalLookupMatrix torsions;
    double internal_molecule_E();
#endif

    Residue * contiguousResidues;
    int contiguousResiduesSize;
    Molecule * molecules;
    int moleculeCount;
    int moleculeArraySize;
    int residueCount;

    float boundingValue;

    float temperature;
    short label;
    double potential;
// #if CUDA_E && !CUDA_STREAMS // TODO: really?
    double newPotential;
// #else
//     float newPotential;
// #endif

    float translateStep;
    double rotateStep;

    // rngs for this object
    gsl_rng * rng; // we only need one Mersenne Twister for ALL THE THINGS, really.
#if FLEXIBLE_LINKS
    double * MC_move_weights; // weightings for discrete MC move probabilities
    gsl_ran_discrete_t * MC_discrete_table; // preprocessed lookup table for MC move generator
#endif

    int threadCount;
    int nonCrowderCount;
//     int maxMoleculeSize;
    int max_residue_count;
    int max_segment_count;
    Molecule savedMolecule;

    int accept;
    int acceptA;
    int reject;
    int totalAcceptReject;
    int totalAccept;

    int boundSamples;
    int samplesSinceLastExchange;

    int totalBoundSamples;
    int totalSamples;

    bool sample(SimulationData *data, int current_step, float boundEnergyThreshHold, pthread_mutex_t *writeFileMutex);

    int countpairs();
    int paircount;

    int nonCrowderResidues; // the number of non crowder residues for fast bound configurations test
    void countNonCrowdingResidues();

#if USING_CUDA

    float4 *device_float4_residuePositions; // x,y,z,w == position.x,position.y,position.z,id
    float4 *device_float4_residueMeta;		// x,y,z,w == index,charge,vdwr,temperature

    float4 *host_float4_residuePositions;
    float4 *host_float4_residueMeta;

    float *device_LJPotentials;

    int *device_residueCount;
    int *device_moleculeCount;
    int *device_moleculeStartPositions;		// position offset where each molecule starts in the above arrays
    int *host_moleculeStartPositions; 		// tracker for the above on the hose side so that individual molecules can be updated without additional calculation

    int blockSize;
    int dataSetSize;	// size of the data on gpu including padding residues
    int paddedSize;		// size taken in gpu memory because of alignment padding
    int gridSize;
    int resultSize;
    int sharedMemSize;

    bool replicaIsOnDevice;


#if CUDA_MC
    int *device_moleculeLengths;
    float4 *device_moleculeCenters;				// rotational centers of each molecule

    float4 *device_translationVector;
    float4 *device_reverseTranslationVector;
    float4 *device_rotationVector;  // vector(x,y,z)|amount(w)
    float4 *device_reverseRotationVector;
    float4 *host_translationVector;
    float4 *host_reverseTranslationVector;
    float4 *host_rotationVector;  // vector(x,y,z)|amount(w)
    float4 *host_reverseRotationVector;

    bool rotateOnDevice(int moleculeID, Vector3f vector, float amount);
    bool translateOnDevice(int moleculeID, Vector3f translation);
    bool lastMutationWasATranslate;
    int lastMutatedMolecule;
    void cudaRollbackMutation();
#endif

    void ReplicaDataToDevice();	// copy the replica to the device
    void ReplicaDataToHost(); // TODO: Never used?
    void UpdateDeviceData(); // TODO: Never used?
    void MoleculeDataToDevice(int MoleculeID); // update a molecule on the device
    double EonDevice();
    double EonDeviceNC();
    void EonDeviceAsync();
    void setBlockSize(int blockSize); // this setter actually does stuff and is necessary

    void FreeDevice();

#if CUDA_STREAMS
    cudaStream_t cudaStream;
    float *device_kernelResult;  	// holds the result of a kernel on the device
    float *kernelResult;			// as above but on host, sync cudaStream to make them the same
    void ReserveSumSpace();	// space on host and device to allow for the sum of kernel returns
    void FreeSumSpace();		// frees the above
    float SumGridResults();	//sums the grid returned by the potential kernel
    uint lastMutationIndex;
    //split the functions to allow for latency hiding and overlapping calls
    void MCSearchMutate();
    void MCSearchEvaluate();
    void MCSearchAcceptReject();
#endif

#endif

#if INCLUDE_TIMERS
    //timers for profiling the cuda functions
    uint replicaToGPUTimer;
    uint replicaToHostTimer;
    uint replicaUpdateGPUTimer;
    uint replicaECUDATimer;
    uint replicaMoleculeUpdateTimer;
    uint initGPUMemoryTimer;
    uint replicaEHostTimer;
    bool timersInit;
#endif

    void initTimers();
    void printTimers();

    Replica * GLreplica;
};

#endif /*REPLICA_H_*/