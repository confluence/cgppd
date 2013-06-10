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
#include <map>
#include <unistd.h>  // include for profiling, TAU cannot find getpid() as it is an externed call
#include <getopt.h>

#include "definitions.h"

#if USING_CUDA
#include <cuda.h>
#include <cutil.h>  // CUDA c util package
#include <cuda_runtime_api.h>
#include <cutil_inline.h>
#include <cutil_inline_runtime.h>
#endif

#include "structures.h"
#include "AminoAcid.h"
#include "TorsionalLookupMatrix.h"
#include "Replica.h"
#include "vector3f.h"
#include "Quaternion.h"

#if USING_CUDA
#include "cudaExterns.h"
#endif

#if GLVIS
#include "openglvis.h"
#endif

using namespace std;

class Simulation
{
public:
    argdata parameters;

    AminoAcids aminoAcidData;

    FILE *boundConformationsFile;
    FILE *fractionBoundFile;
    FILE *acceptanceRatioFile;
    FILE *exchangeFrequencyFile;

    Replica initialReplica;
    Replica replica[REPLICA_COUNT];   // container for all replicas in the simulation
    Replica * _300kReplica;

    int waitingThreads;
    int exchanges;  // number of exchanges performed
    int tests;      // number of replica exchange tests
    int totalExchanges; // total number of exchanges
    int totalTests; // total number of exchanges
    int offset;  // RE offset
    int steps;   // # re steps performed

    double geometricTemperature;
    double geometricTranslate;
    double geometricRotation;

    pthread_t *thread;
    SimulationData *data;
    bool thread_created;
    pthread_attr_t attr;
    pthread_mutex_t waitingCounterMutex;
    pthread_mutex_t writeFileMutex;
    pthread_cond_t waitingThreadCond;
    pthread_cond_t waitingReplicaExchangeCond;

    gsl_rng  * REMCRng;         // random numbers for the REMC method to swap labels

    map<float, int> replicaTemperatureMap;
    map<int, float> temperatureMap;

#if INCLUDE_TIMERS
    uint RELoopTimer;
    uint MCLoopTimer;
#endif

    Simulation();
    ~Simulation();


    void init(int argc, char **argv, int pid);

    void printHelp();
    void getArgs(int argc, char **argv);
    void loadArgsFromFile();
    void check_and_modify_parameters();
    void writeFileIndex();
    void printArgs();

    void run();

    void initSamplingFiles();
    void closeSamplingFiles();
};

void *MCthreadableFunction(void *arg);

#endif /*SIMULATION_H_*/