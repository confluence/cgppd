#ifndef REPLICA_H_
#define REPLICA_H_
#include "vector3f.h"
#include "definitions.h"
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

    void setAminoAcidData(AminoAcids a);
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

    // TODO: are both of these used?
    int loadMolecule(const char *pdbfilename);
    int loadMolecule(const char *pdbfilename, Vector3f position, Vector3double rotationAxis, double rotationAmount);

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
    void setBoundingValue(float value) { // TODO: eliminate this setter
        boundingValue = value;
    };

    float temperature;
    short label;
    double potential;

    float translateStep;
    double rotateStep;

    void setTranslateStep(float t) { // TODO: eliminate this setter
        translateStep = t;
    }
    void setRotateStep(double r) { // TODO: eliminate this setter
        rotateStep = r;
    }
    void setTemperature(float t) { // TODO: eliminate this setter
        temperature = t;
    }
    void setLabel(short l) { // TODO: eliminate this setter
        label = l;
    }

    // rngs for this object
    gsl_rng * rng; // we only need one Mersenne Twister for ALL THE THINGS, really.
#if FLEXIBLE_LINKS
    double * MC_move_weights; // weightings for discrete MC move probabilities
    gsl_ran_discrete_t * MC_discrete_table; // preprocessed lookup table for MC move generator
#endif

    int threadCount;
    int nonCrowderCount;
    int maxMoleculeSize;
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

    uint get_MC_mutation_type(const Molecule* m); // randomly select type of MC move

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

    void setDeviceLJPotentials(float * ljp) { // TODO: eliminate this setter
        device_LJPotentials = ljp;
    }
    void ReplicaDataToHost(); // TODO: Never used? Remove?
    void UpdateDeviceData(); // TODO: Never used? Remove?
    void MoleculeDataToDevice(int MoleculeID); // update a molecule on the device
    double EonDevice();
    double EonDeviceNC();
    void EonDeviceAsync();
    void setLJpotentials(float *ljp); // TODO: eliminate this setter
    void setBlockSize(int blockSize);

    void FreeDevice();

#if CUDA_STREAMS
    cudaStream_t cudaStream;
    float *device_kernelResult;  	// holds the result of a kernel on the device
    float *kernelResult;			// as above but on host, sync cudaStream to make them the same
    void ReserveSumSpace();	// space on host and device to allow for the sum of kernel returns
    void FreeSumSpace();		// frees the above
    float SumGridResults();	//sums the grid returned by the potential kernel
    uint lastMutationIndex;
    float oldPotential;
    float newPotential;
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

    // TODO: remove these once the ones in molecule are used everywhere
    Vector3f createNormalisedRandomVector(gsl_rng * r);
    Vector3double createNormalisedRandomVectord(gsl_rng * r);

    Replica * GLreplica;

// TODO: actually necessary for these to be private?
private:
//     void rotate(const int m, const double step);
//     void translate(const int m, const float step);
// #if FLEXIBLE_LINKS
//     void rotate_domain(const int m, const float rotateStep, const int ri, const bool before);
//     void local(const int m, const float step, const int num_moves);
// #endif

};

#endif /*REPLICA_H_*/
