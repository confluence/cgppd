#ifndef STRUCTURES_H_
    #define STRUCTURES_H_

    #include <vector>
    #include <iostream>

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
        char checkpointfilename[256];
        bool resume;
        int currentStep;
        std::vector<moldata> mdata;
    };

#endif