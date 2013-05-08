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

//     float e = 2.71828182845904523536f;
    size_t lowestEnergy;
    Replica initialReplica;
    Replica replica[REPLICA_COUNT];   // container for all replicas in the simulation
    bool replicasInitialised;

    // TODO: create one rng; pass into replicas! Thesis assumes one rng for all replicas.
    gsl_rng  * REMCRng;         // random numbers for the REMC method to swap labels

    Simulation();
    ~Simulation();

    void printHelp();
    void getArgs(int argc, char **argv);
    void loadArgsFromFile();
    void checkParameterSanity();
    void writeFileIndex();
    void printArgs();
    void init(int argc, char **argv, int pid);

#if USING_CUDA
    void run_check();
#endif

    void run();

    void initSamplingFiles (FILE ** fractionBoundFile, FILE ** boundConformationsFile, FILE ** acceptanceRatioFile,  FILE ** exchangeFrequencyFile);
    void closeSamplingFiles (FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile, FILE * exchangeFrequencyFile);
};

void *MCthreadableFunction(void *arg);

#endif /*SIMULATION_H_*/