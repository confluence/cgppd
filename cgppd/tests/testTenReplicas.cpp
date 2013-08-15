#include <cppunit/extensions/TestFactoryRegistry.h>

#include <cppunit/extensions/HelperMacros.h>
#include "definitions.h"
#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include <Replica.h>

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
#if USING_CUDA
    cuInit(0);
#endif

    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);
    testboxdim = 118.4f;

#if USING_CUDA
    cudaMalloc((void**)&ljp_t, LJArraySize);
    cutilCheckMsg("Failed to cudaMalloc");
    copyLJPotentialDataToDevice(ljp_t, &aminoAcidData);

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
        float gpu;
    };

    static const ExpectedResult expected_results[10] = {
        { -0.293705,  -0.293705},
        { -1.056291,  -1.056291},
        {-10.277430, -10.277432},
        { -7.577059,  -7.577060},
        {  0.000106,   0.000106},
        { -5.559506,  -5.559508},
        { -5.441210,  -5.441211},
        {-10.657519, -10.657518},
        { -9.891660,  -9.891659},
        { -8.511853,  -8.511855}
    };

    for (int i = 0; i < 10; i++)
    {
        argdata parameters;
        parameters.bound = testboxdim;
#if USING_CUDA
        parameters.auto_blockdim = false;
#endif

        moldata m1;
        m1.translate = true;
        sprintf(m1.pdbfilename, "data/conf%d/1a.pdb", i+1);

        moldata m2;
        m2.translate = true;
        sprintf(m2.pdbfilename, "data/conf%d/1b.pdb", i+1);

        parameters.mdata.push_back(m1);
        parameters.mdata.push_back(m2);

        replicas[i].init_first_replica(parameters, aminoAcidData, 2);
        replicas[i].label = i+1;

#if INCLUDE_TIMERS
        replicas[i].initTimers();
#endif

#if USING_CUDA
        replicas[i].setup_CUDA(ljp_t);
#endif

        float e = 0.000001;
        double cpu = replicas[i].E().total();
        double cpu_nc = replicas[i].E(&replicas[i].molecules[0], &replicas[i].molecules[1]).total();
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].cpu, cpu, e);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].cpu, cpu_nc, e);

#if USING_CUDA
        double gpu = replicas[i].EonDevice();
        double gpu_nc = replicas[i].EonDeviceNC();
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].gpu, gpu, e);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].gpu, gpu_nc, e);

        replicas[i].FreeDevice();
#endif
    }

}
