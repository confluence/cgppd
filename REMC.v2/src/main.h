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

//all cuda things
#if USING_CUDA
#include <cuda.h>
#include <cutil.h>  // CUDA c util package
#include <cuda_runtime_api.h>
#include <cutil_inline.h>
#include <cutil_inline_runtime.h>
#include "cudaExterns.h"

int cuda_blockSize = TILE_DIM;
bool auto_blockdim = true;

#endif

using namespace std;

float e = 2.71828182845904523536f;
size_t lowestEnergy;
Replica replica[REPLICA_COUNT];   // container for all replicas in the simulation
vector<uint> nonCrowderMolecules; // contains a list of each molecule/protein we want, lets us quickly determine which ones to save

// TODO: create one rng; pass into replicas! Thesis assumes one rng for all replicas.
gsl_rng  * REMCRng;         // random numbers for the REMC method to swap labels
bool exitCondition = false;
bool viewConditions = false;
bool skipsimulation = false;
int mcstepcount = 0;
int threads = THREAD_COUNT;   // 2 threads per core i think
int streams = STREAM_COUNT;

void REMCSimulation(Replica *initialReplica, argdata *parameters);
bool getArgs(argdata * d, int argc, char **argv);
void printHelp(bool badArg); // TODO: wrapper for bad args + help; rename to usage
void loadArgsFromFile(argdata * parameters);
void *MCthreadableFunction(void *arg);
void initSamplingFiles (argdata *args, FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile);
void closeSamplingFiles (argdata *args, FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile);

// All the opengl stuff
#if GLVIS
Replica *GLreplica;
#include "openglvis.h"
#endif
