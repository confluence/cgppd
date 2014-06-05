#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <string.h>
#include <pthread.h>
#include <sys/mman.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <map>
#include <unistd.h>  // include for profiling, TAU cannot find getpid() as it is an externed call
#include <getopt.h>

#include <glog/logging.h>

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
    float near_300k; // temperature closest to 300k in this simulation
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
#ifdef VERBOSE_THREAD_LOGGING
    pthread_mutex_t logMutex;
#endif
    pthread_cond_t waitingThreadCond;
    pthread_cond_t waitingReplicaExchangeCond;

    gsl_rng  * REMCRng;         // random numbers for the REMC method to swap labels

    map<float, int> position_of_temperature;
    map<int, float> temperature;

#if INCLUDE_TIMERS
    uint RELoopTimer;
    uint MCLoopTimer;
#endif

    Simulation();
    ~Simulation();


    void init(int argc, char **argv, int pid);

    void printHelp();
    void getArgs(int argc, char **argv, bool first_pass);
    void loadArgsFromFile();
    void check_and_modify_parameters();
    void writeFileIndex();
    void printArgs();

    void run();

    void initSamplingFile(const char * name, FILE ** file_addr);
    void closeSamplingFile(const char * name, FILE ** file_addr);
    void initSamplingFiles();
    void writeSamplingFileHeaders();
    void closeSamplingFiles();

    float exchangeFraction;
    float accumulativeExchangeFraction;
    void exchange_frequency();
};

#if USING_CUDA
#if CUDA_STREAMS
void setup_CUDA(int device_id, float box_dimension, float * device_LJ_potentials, AminoAcids * amino_acid_data, cudaStream_t * streams, int streams_per_thread);
void teardown_CUDA(float * device_LJ_potentials, cudaStream_t * streams, int streams_per_thread);
#else
void setup_CUDA(int device_id, float box_dimension, float * device_LJ_potentials, AminoAcids * amino_acid_data);
void teardown_CUDA(float * device_LJ_potentials);
#endif // CUDA_STREAMS
#endif // USING_CUDA

void *MCthreadableFunction(void *arg);

#endif /*SIMULATION_H_*/
