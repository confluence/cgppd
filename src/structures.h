#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <vector>
#include <map>
#include <iostream>
#include <string.h>

#include "definitions.h"

using namespace std;

class Replica;

struct SimulationData
{
    long index;
    Replica *replica;
    int replicaCount;
    int GPUID;
    int threads;
    int streams;
    int replicas_in_this_thread;
    int max_replicas_per_thread;
    int replicas_per_stream;
    int streams_per_thread;
    int MCsteps;
    int REsteps; // TODO: will this ever be used in the threadable function if we precalculate mcstepsPerRE?
    int MC_steps_per_RE;
    int sampleFrequency;
    int sampleStartsAfter;
    int checkpointFrequency;
    float bound;
    int * waitingThreadCount;
    char prefix[256];
    char title[256];
    FILE * fractionBound;
    FILE * boundConformations;
    pthread_mutex_t * waitingCounterMutex;
    pthread_mutex_t * writeFileMutex;
#ifdef VERBOSE_THREAD_LOGGING
    pthread_mutex_t * logMutex;
#endif
    pthread_cond_t * waitingThreadCond;
    pthread_cond_t * waitingReplicaExchangeCond;
};

struct segdata
{
    vector<int> residue_indices;
};

struct moldata
{
    char pdbfilename[256];
    char name[256];
    float px;
    float py;
    float pz;
    bool translate;
    float rx;
    float ry;
    float rz;
    float ra;
    bool crowder;
    bool all_flexible;

    vector<segdata> segments;

    moldata(): px(0), py(0), pz(0), translate(false), rx(0), ry(0), rz(0), ra(0), crowder(false), all_flexible(false)
    {
        memset(pdbfilename, 0, 256);
        memset(name, 0, 256);
    }
};

struct argdata
{
    bool viewConditions;
    bool skipsimulation;
#if USING_CUDA
    int cuda_blockSize;
    bool auto_blockdim;
#endif
    int threads;
    int streams;
    int max_replicas_per_thread;
    int gpus;
    int MCsteps;
    int REsteps;
    int replicas;
    float bound;
    int sampleFrequency;
    int sampleStartsAfter;
    int checkpointFrequency;
    char prefix[256];
    char title[256];
    char file[256];
    bool inputFile;
//         TODO: for some reason I wanted to remove this
    int nonCrowders;
    float temperatureMin;
    float temperatureMax;
    int pid;
    // TODO: this is gone until checkpointing goes back in
//         char checkpointfilename[256];
    bool resume;
    vector<moldata> mdata;
    map<string, int> mdata_map;

    argdata() :
        nonCrowders(0), viewConditions(false), skipsimulation(false),
#if USING_CUDA
        cuda_blockSize(TILE_DIM), auto_blockdim(true),
#endif
        gpus(1), threads(THREAD_COUNT), streams(STREAM_COUNT), MCsteps(MC_STEPS), REsteps(REMC_STEPS),
        sampleFrequency(SAMPLE_FREQ), sampleStartsAfter(STEPS_BEFORE_SAMPLE), inputFile(false),
        replicas(REPLICA_COUNT), bound(BOUNDING_VALUE), temperatureMin(LOWTEMP), temperatureMax(HIGHTEMP),
        pid(0), resume(false), checkpointFrequency(CHECKPOINTFREQUENCY), max_replicas_per_thread(0)
    {
        memset(prefix, 0, 256);
        memset(file, 0, 256);
        memset(title, 0, 256);
    }
};

#endif
