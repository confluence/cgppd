#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>
#include <Replica.h>
#include <cutil.h>
#include <cutil_inline.h>
#include "definitions.h"

class TestTenReplicas : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestTenReplicas);
    CPPUNIT_TEST(testSanity);
    CPPUNIT_TEST_SUITE_END();

private:
    AminoAcids aminoAcidData;
    Replica replicas[10];
    char *test_molecule_file;
    float * ljp_t;
    float testboxdim;

public:
    void setUp();
    void testSanity();
    void tearDown();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestTenReplicas);

void TestTenReplicas::setUp()
{
    cuInit(0);

    aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
    aminoAcidData.loadLJPotentialData(LJPDSOURCE);
    testboxdim = 118.4f;
    test_molecule_file = new char[60];

#if USING_CUDA
    cudaMalloc((void**)&ljp_t,LJArraySize);
    cutilCheckMsg("Failed to cudaMalloc");
    copyLJPotentialDataToDevice(ljp_t,&aminoAcidData);

    // set box dimensions
    CUDA_setBoxDimension(testboxdim);
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    bindLJTexture(ljp_t);
#endif
#endif
}

void TestTenReplicas::tearDown()
{
#if USING_CUDA
    cudaFree(ljp_t);
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    unbindLJTexture();
#endif
#endif
    cout.flush();
}

void TestTenReplicas::testSanity()
{
    struct ExpectedResult
    {
        float cpu;
        float cpu_nc;
        float gpu;
        float gpu_nc;
    };

    static const ExpectedResult expected_results[10] = {
        { -0.293705,  -0.293705,  -0.293705,  -0.293705},
        { -1.056291,  -1.056291,  -1.056291,  -1.056291},
        {-10.277430, -10.277433, -10.277432, -10.277432},
        { -7.577059,  -7.577051,  -7.577060,  -7.577060},
        {  0.000106,   0.000106,   0.000106,   0.000106},
        { -5.559506,  -5.559506,  -5.559508,  -5.559508},
        { -5.441210,  -5.441209,  -5.441211,  -5.441211},
        {-10.657519, -10.657518, -10.657518, -10.657518},
        { -9.891660,  -9.891660,  -9.891659,  -9.891659},
        { -8.511853,  -8.511852,  -8.511855,  -8.511855}
    };

    static const float expected_averages [6] = {0.029f, 0.0f, 0.271f, 0.981f, 0.0f, 0.09f};
    static int exceeded_averages [6] = {0, 0, 0, 0, 0, 0};

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

#if USING_CUDA
        replicas[i].setDeviceLJPotentials(ljp_t);
        replicas[i].setBlockSize(TILE_DIM);
        replicas[i].ReplicaDataToDevice();
#endif

        float e = 0.000001;
        double cpu = replicas[i].E();
        double cpu_nc = replicas[i].E(&replicas[i].molecules[0],&replicas[i].molecules[1]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].cpu, cpu, e);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].cpu_nc, cpu_nc, e);

#if USING_CUDA
        double gpu = replicas[i].EonDevice();
        double gpu_nc = replicas[i].EonDeviceNC();
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].gpu, gpu, e);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].gpu_nc, gpu_nc, e);
#endif

        float averages[6];

        averages[0] = cutGetAverageTimerValue(replicas[i].replicaToGPUTimer);
        averages[1] = cutGetAverageTimerValue(replicas[i].replicaUpdateGPUTimer);
        averages[2] = cutGetAverageTimerValue(replicas[i].replicaECUDATimer);
        averages[3] = cutGetAverageTimerValue(replicas[i].replicaEHostTimer);
        averages[4] = cutGetAverageTimerValue(replicas[i].replicaMoleculeUpdateTimer);
        averages[5] = cutGetAverageTimerValue(replicas[i].initGPUMemoryTimer);

        for (int j = 0; j < 6; j++)
        {
            if (averages[j] > expected_averages[j] + 0.01)
            {
                exceeded_averages[j]++;
            }
        }

#if USING_CUDA
        replicas[i].FreeDevice();
#endif
    }

    for (int j = 0; j < 6; j++)
    {
        // ignore one or two outliers
        CPPUNIT_ASSERT(exceeded_averages[j] <= 2);
    }
}
