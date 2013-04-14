#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <pthread.h>
#include <sys/mman.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "definitions.h"
#include "AminoAcid.h"
#include "TorsionalLookupMatrix.h"
#include "Replica.h"
#include "vector3f.h"
#include "Quaternion.h"
#include <map>
#include <unistd.h>  // include for profiling, TAU cannot find getpid() as it is an externed call

#if USING_CUDA
#include <cuda.h>
#include <cutil.h>  // CUDA c util package
#include <cuda_runtime_api.h>
#include <cutil_inline.h>
#include <cutil_inline_runtime.h>
#include "cudaExterns.h"
#endif

#if GLVIS
#include "openglvis.h"
#endif

class Simulation
{
public:
    #if GLVIS
    Replica *GLreplica;
    #endif

    AminoAcids aminoAcidData;

    float e = 2.71828182845904523536f;
    size_t lowestEnergy;
    Replica initialReplica;
    Replica replica[REPLICA_COUNT];   // container for all replicas in the simulation
    bool replicasInitialised = false;

    // TODO: create one rng; pass into replicas! Thesis assumes one rng for all replicas.
    gsl_rng  * REMCRng;         // random numbers for the REMC method to swap labels

    // required objects for synchronisation
    pthread_mutex_t waitingThreadMutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t waitingCounterMutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t reMutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t writeFileMutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t waitingThreadCond;
    pthread_cond_t waitingReplicaExchangeCond;
    int waitingThreadCount;

    Simulation(argdata parameters);
    ~Simulation();

    void run_check();
    void init();

    void *MCthreadableFunction(void *arg);
    void run();

    void initSamplingFiles (FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile);
    void closeSamplingFiles (FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile);
};

#endif /*SIMULATION_H_*/