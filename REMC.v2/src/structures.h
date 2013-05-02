#ifndef STRUCTURES_H_
    #define STRUCTURES_H_

    #include <vector>
    #include <iostream>
    #include <string.h>

    #include "definitions.h"

    class Replica;

    struct SimulationData
    {
        long index;
        Replica *replica;
        int replicaCount;
        int GPUID;
        int threads;
        int streams;
        int MCsteps;
        int REsteps;
        int sampleFrequency;
        int sampleStartsAfter;
        float bound;
        int * waitingThreadCount;
        int * conformationsBound;
        FILE * fractionBound;
        FILE * boundConformations;
    };

    struct moldata
    {
        char pdbfilename[256];
        float px;
        float py;
        float pz;
        bool translate;
        float rx;
        float ry;
        float rz;
        float ra;
        bool crowder;
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
        int gpus;
        int MCsteps;
        int REsteps;
        int replicas;
        float bound;
        int sampleFrequency;
        int sampleStartsAfter;
        char prependageString[256];
        char file[256];
        char logfile[256];
        bool inputFile;
//         TODO: for some reason I wanted to remove this
        int nonCrowders;
        float temperatureMin;
        float temperatureMax;
        int pid;
        // TODO: this is gone until checkpointing goes back in
//         char checkpointfilename[256];
        bool resume;
        int currentStep;
        std::vector<moldata> mdata;

        argdata() :
            resume(false), nonCrowders(0), viewConditions(false), skipsimulation(false),
#if USING_CUDA
            cuda_blockSize(TILE_DIM), auto_blockdim(true),
#endif
            gpus(1), threads(THREAD_COUNT), streams(STREAM_COUNT), MCsteps(MC_STEPS), REsteps(REMC_STEPS),
            sampleFrequency(SAMPLE_FREQ), sampleStartsAfter(STEPS_BEFORE_SAMPLE), inputFile(false),
            replicas(REPLICA_COUNT), bound(BOUNDING_VALUE), temperatureMin(LOWTEMP), temperatureMax(HIGHTEMP), pid(0)
        {
            memset(prependageString, 0, 256);
            memset(file, 0, 256);
            memset(logfile, 0, 256);
        }
    };

#endif