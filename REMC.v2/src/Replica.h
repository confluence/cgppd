#ifndef REPLICA_H_
#define REPLICA_H_
#include "vector3f.h"
#include "definitions.h"
#include "Quaternion.h"
#include <fstream>
#include <vector>
#include "Molecule.h"
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
//     Replica(const Replica& r);
    void setAminoAcidData(AminoAcids a);
    void reserveContiguousMoleculeArray(int size);
    void initRNGs();
    void freeRNGs();
    void copy(const Replica& r); 						// assignment copy, duplicate it.

    void exchangeReplicas(Replica &r);

    void kahan_sum(double potential, double p_ij, double &c);
    double E();

    float E(Molecule *a,Molecule *b);   // use for fraction bound calcs

    // the search method that mutates the replicas molecules.
    void MCSearch(int steps);//, Replica * replaceReplica);

//     bool savePDB(const char *filename); // save multiple pdb files per replica
//     void saveAsSinglePDB(const char *filename); //save one pdb, each molecule is given a chain id from A upward

    int loadMolecule(const char *pdbfilename);
    int loadMolecule(const char *pdbfilename, Vector3f position, Vector3double rotationAxis, double rotationAmount);

    AminoAcids aminoAcids;

#ifdef FLEXIBLE_LINKS
    TorsionalLookupMatrix torsions;
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
    gsl_rng * rng_translate;	// the translation vector rngs
    gsl_rng * rng_translateAmount;	// the translation vector rngs
    gsl_rng * rng_rotate;		// the rotation vector rng
    gsl_rng * rng_rotateAmount;	// the translation vector rngs
    gsl_rng * rng_moleculeSelection; // molecule selection rng
    gsl_rng * MCRng;  // used to determine what change to make
    gsl_rng * MCKbRng;	// boltzmann acceptance rng in MC

    int threadCount;
    int nonCrowderCount;
    int maxMoleculeSize;

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
    Molecule savedMolecule;
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

    Vector3f createNormalisedRandomVector(gsl_rng * r);
    Vector3double createNormalisedRandomVectord(gsl_rng * r);

    Replica * GLreplica;


private:
    void rotate(const int m, const double step);
    void translate(const int m, const float step);
};

#endif /*REPLICA_H_*/
