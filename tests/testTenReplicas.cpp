#include "definitions.h"
#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include <Simulation.h>
#include <Replica.h>
#include <testCommon.h>

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

#if !LJ_REPULSIVE && !LJ_OFF  // ignore this whole test unless we're using the normal LJ constants
CPPUNIT_TEST_SUITE_REGISTRATION(TestTenReplicas);
#endif // !LJ_REPULSIVE && !LJ_OFF

void TestTenReplicas::setUp()
{
    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);
    testboxdim = 1000.0f;
#if USING_CUDA
    setup_CUDA(0, testboxdim, ljp_t, &aminoAcidData);
#endif // USING_CUDA
}

void TestTenReplicas::tearDown()
{
#if USING_CUDA
    teardown_CUDA(ljp_t);
#endif // USING_CUDA
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
        { -0.000079,  -0.000079},
        { -5.559506,  -5.559508},
        { -5.441210,  -5.441211},
        {-10.657519, -10.657518},
        { -9.891660,  -9.891659},
        { -8.511853,  -8.511855}
    };

#if PRINT_REFERENCE_CONFORMATIONS_FROM_TEST
    struct Subtotals
    {
        float total;
        float LJ;
        float DH;
    };

    static const Subtotals CHARMM_results[10] = {
        { -0.294,  -0.081,  -0.213},
        { -1.056,  -1.322,  0.266},
        { -10.277,  -9.095,  -1.182},
        { -7.580,  -5.903,  -1.678},
        {-0.000079, -0.000021, -0.000058},
        { -5.562,  -4.810,  -0.752},
        { -5.480,  -4.213,  -1.267},
        { -10.712,  -9.266,  -1.446},
        { -9.900,  -7.951,  -1.949},
        { -8.528,  -7.459,  -1.069}
    };

    float mean_relative_error_total(0);
    float mean_relative_error_LJ(0);
    float mean_relative_error_DH(0);

    float max_relative_error_total(0);
    float max_relative_error_LJ(0);
    float max_relative_error_DH(0);
    printf("\tTOTAL\tLJ\tDH\n");
#endif

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
        Potential cpu = replicas[i].E();
        Potential cpu_nc = replicas[i].E(&replicas[i].molecules[0], &replicas[i].molecules[1]);

#if PRINT_REFERENCE_CONFORMATIONS_FROM_TEST
        if (i == 4)
        {
            printf("%2d CPU: %10.3e %10.3e %10.3e\n", i + 1, cpu.total(), cpu.total_LJ(), cpu.total_DH());
        }
        else
        {
            printf("%2d CPU: %10.3f %10.3f %10.3f\n", i + 1, cpu.total(), cpu.total_LJ(), cpu.total_DH());
        }
#endif

//         CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].cpu, cpu.total(), e);
//         CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].cpu, cpu_nc.total(), e);

#if PRINT_REFERENCE_CONFORMATIONS_FROM_TEST
        float relative_error_total = abs(CHARMM_results[i].total - cpu.total())/abs(CHARMM_results[i].total);
        float relative_error_LJ = abs(CHARMM_results[i].LJ - cpu.total_LJ())/abs(CHARMM_results[i].LJ);
        float relative_error_DH = abs(CHARMM_results[i].DH - cpu.total_DH())/abs(CHARMM_results[i].DH);

        printf("ERRORS: %f %f %f\n", relative_error_total, relative_error_LJ, relative_error_DH);

        mean_relative_error_total += relative_error_total/ 10;
        mean_relative_error_LJ += relative_error_LJ / 10;
        mean_relative_error_DH += relative_error_DH / 10;

        max_relative_error_total = max(max_relative_error_total, relative_error_total);
        max_relative_error_LJ = max(max_relative_error_LJ, relative_error_LJ);
        max_relative_error_DH = max(max_relative_error_DH, relative_error_DH);
#endif

#if USING_CUDA
        double gpu = replicas[i].EonDevice();
        double gpu_nc = replicas[i].EonDeviceNC();

        #if PRINT_REFERENCE_CONFORMATIONS_FROM_TEST
        if (i == 4)
        {
            printf("%2d GPU: %10.3e\n\n", i + 1, gpu);
        }
        else
        {
            printf("%2d GPU: %10.3f\n\n", i + 1, gpu);
        }

#endif
//         CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].gpu, gpu, e);
//         CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_results[i].gpu, gpu_nc, e);

        replicas[i].FreeDevice();
#endif
    }

#if PRINT_REFERENCE_CONFORMATIONS_FROM_TEST
    printf("\tMean relative error\tMax relative error\n");
    printf("Tot\t%f\t%f\n", mean_relative_error_total, max_relative_error_total);
    printf("LJ\t%f\t%f\n", mean_relative_error_LJ, max_relative_error_LJ);
    printf("DH\t%f\t%f\n", mean_relative_error_DH, max_relative_error_DH);
#endif
}
