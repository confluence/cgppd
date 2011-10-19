#include <UnitTest++.h>
#include <Replica.h>
#include <cutil.h>
#include <cutil_inline.h>
#include "definitions.h"

struct TenReplicasFixture
{
    AminoAcids aminoAcidData;
    Replica replicas[10];
    char *test_molecule_file;
    float * ljp_t;
    float testboxdim;

    TenReplicasFixture()
    {
        cuInit(0);

        aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
        aminoAcidData.loadLJPotentialData(LJPDSOURCE);

        cudaMalloc((void**)&ljp_t,LJArraySize);
        cutilCheckMsg("Failed to cudaMalloc");
        copyLJPotentialDataToDevice(ljp_t,&aminoAcidData);

        // set box dimensions
        testboxdim = 118.4f;
        CUDA_setBoxDimension(testboxdim);

        test_molecule_file = new char[60];

#if LJ_LOOKUP_METHOD == TEXTURE_MEM
        bindLJTexture(ljp_t);
#endif
    }

    ~TenReplicasFixture()
    {
        cudaFree(ljp_t);
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
        unbindLJTexture();
#endif
        cout.flush();
    }
};

TEST_FIXTURE(TenReplicasFixture, TestTenReplicas)
{
//     Replica replicas[10];
//     char *test_molecule_file = new char[60];

    struct ExpectedResult
    {
        float cpu;
        float cpu_nc;
        float gpu;
        float gpu_nc;
    };

    ExpectedResult expected_results[10] = {
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

    float expected_averages [7] = {0.029f, 0.0f, 0.271f, 0.981f, 0.0f, 0.0f, 0.09f};
    int exceeded_averages [7] = {0, 0, 0, 0, 0, 0, 0};

    for (int i = 0; i < 10; i++)
    {
        replicas[i].aminoAcids = aminoAcidData;
        replicas[i].label = i+1;
        replicas[i].setBoundingValue(testboxdim);

        replicas[i].reserveContiguousMoleculeArray(2);

        sprintf(test_molecule_file,"data/conf%d/1a.pdb",i+1);
        replicas[i].loadMolecule(test_molecule_file);
        sprintf(test_molecule_file,"data/conf%d/1b.pdb",i+1);
        replicas[i].loadMolecule(test_molecule_file);

        replicas[i].initTimers();
        replicas[i].countNonCrowdingResidues();

        replicas[i].setDeviceLJPotentials(ljp_t);
        replicas[i].setBlockSize(TILE_DIM);
        replicas[i].ReplicaDataToDevice();

        double gpu = replicas[i].EonDevice();
        double gpu_nc = replicas[i].EonDeviceNC();

        double cpu = replicas[i].E();
        double cpu_nc = replicas[i].E(&replicas[i].molecules[0],&replicas[i].molecules[1]);

        float e = 0.000001;
        CHECK_CLOSE(expected_results[i].cpu, cpu, e);
        CHECK_CLOSE(expected_results[i].cpu_nc, cpu_nc, e);
        CHECK_CLOSE(expected_results[i].gpu, gpu, e);
        CHECK_CLOSE(expected_results[i].gpu_nc, gpu_nc, e);

        float averages[7];

        averages[0] = cutGetAverageTimerValue(replicas[i].replicaToGPUTimer);
        averages[1] = cutGetAverageTimerValue(replicas[i].replicaUpdateGPUTimer);
        averages[2] = cutGetAverageTimerValue(replicas[i].replicaECUDATimer);
        averages[3] = cutGetAverageTimerValue(replicas[i].replicaEHostTimer);
        averages[4] = cutGetAverageTimerValue(replicas[i].replicaMoleculeUpdateTimer);
        averages[5] = cutGetAverageTimerValue(replicas[i].replicaDeviceMCTimer);
        averages[6] = cutGetAverageTimerValue(replicas[i].initGPUMemoryTimer);

        for (int j = 0; j < 7; j++)
        {
            if (averages[j] > expected_averages[j] + 0.01)
            {
                exceeded_averages[j]++;
            }
        }

        replicas[i].FreeDevice();
    }

    for (int j = 0; j < 7; j++)
    {
        // ignore one or two outliers
        CHECK(exceeded_averages[j] <= 2);
    }
}
