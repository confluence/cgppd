#include <UnitTest++.h>
#include <Replica.h>
#include <cutil.h>
#include <cutil_inline.h>

struct ReplicaFixture
{
    ReplicaFixture()
    {
        cuInit(0);

        aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
        aminoAcidData.loadLJPotentialData(LJPDSOURCE);

        // get values for the example conformations
        egnames = new char[60];

        cudaMalloc((void**)&ljp_t,LJArraySize);
        cutilCheckMsg("Failed to cudaMalloc");
        copyLJPotentialDataToDevice(ljp_t,&aminoAcidData);

        // set box dimensions
        testboxdim = 118.4f;
        CUDA_setBoxDimension(testboxdim);

#if LJ_LOOKUP_METHOD == TEXTURE_MEM
        bindLJTexture(ljp_t);
#endif
    }

    ~ReplicaFixture()
    {
        cudaFree(ljp_t);
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
        unbindLJTexture();
#endif
        cout.flush();
    }

    AminoAcids aminoAcidData;
    Replica replicas[10];
    char *egnames;
    float * ljp_t;
    float testboxdim;
};

TEST_FIXTURE(ReplicaFixture, TestReplica)
{
    struct ExpectedResult {
        float cpu;
        float cpu_nc;
        float gpu;
        float gpu_nc;
    };

    ExpectedResult results[10] = {
        { -0.293705,  -0.293705,  -0.293705,  -0.293705},
        { -1.056291,  -1.056291,  -1.056291,  -1.056291},
        {-10.277430, -10.277433, -10.277432, -10.277432},
        { -7.580391,  -7.580381,  -7.580392,  -7.580392},
        {  0.000104,   0.000104,   0.000104,   0.000104},
        { -5.562237,  -5.562236,  -5.562238,  -5.562238},
        { -5.480217,  -5.480216,  -5.480219,  -5.480219},
        {-10.711965, -10.711964, -10.711965, -10.711965},
        { -9.900360,  -9.900360,  -9.900360,  -9.900360},
        { -8.527747,  -8.527744,  -8.527748,  -8.527748}
    };

    for (int i = 0; i<=9; i++)
    {
        replicas[i].aminoAcids = aminoAcidData;
        replicas[i].label = i+1;
        replicas[i].setBoundingValue(testboxdim);

        sprintf(egnames,"data/conf%d/1a.pdb",i+1);
        replicas[i].loadMolecule(egnames);
        sprintf(egnames,"data/conf%d/1b.pdb",i+1);
        replicas[i].loadMolecule(egnames);

        replicas[i].initTimers();
        replicas[i].countNonCrowdingResidues();
        replicas[i].E();

        replicas[i].setDeviceLJPotentials(ljp_t);
        replicas[i].setBlockSize(TILE_DIM);
        replicas[i].ReplicaDataToDevice();

        double gpu = replicas[i].EonDevice();
        double gpu_nc = replicas[i].EonDeviceNC();

        double cpu = replicas[i].E();
        double cpu_nc = replicas[i].E(&replicas[i].molecules[0],&replicas[i].molecules[1]);

        float e = 0.000001;
        CHECK_CLOSE(cpu, results[i].cpu, e);
        CHECK_CLOSE(cpu_nc, results[i].cpu_nc, e);
        CHECK_CLOSE(gpu, results[i].gpu, e);
        CHECK_CLOSE(gpu_nc, results[i].gpu_nc, e);

        replicas[i].FreeDevice();
    }
}

// TODO: check timers -- decouple from printing
// TODO: stop printing random stuff in replica?
//         replicas[i].printTimers();
