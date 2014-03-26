#include "definitions.h"
#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include <Simulation.h>
#include <Replica.h>
#include <testCommon.h>

class TestReferenceConformations : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestReferenceConformations);
    CPPUNIT_TEST(testExactResults);
    CPPUNIT_TEST(testVsKimResults);
    CPPUNIT_TEST(testVsCharmmResults);
    CPPUNIT_TEST_SUITE_END();

private:
    AminoAcids aminoAcidData;
    Replica replicas[10];
    float * ljp_t;
    float testboxdim;
    Potential cpu[10];
    Potential cpu_nc[10];
#if USING_CUDA
    double gpu[10];
    double gpu_nc[10];
#endif

public:
    void setUp();
    void testExactResults();
    void testVsKimResults();
    void testVsCharmmResults();
    void tearDown();
};

#if !LJ_REPULSIVE && !LJ_OFF  // ignore this whole test unless we're using the normal LJ constants
CPPUNIT_TEST_SUITE_REGISTRATION(TestReferenceConformations);
#endif // !LJ_REPULSIVE && !LJ_OFF

void TestReferenceConformations::setUp()
{
    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);
    testboxdim = 1000.0f;
#if USING_CUDA
    setup_CUDA(0, testboxdim, ljp_t, &aminoAcidData);
#endif // USING_CUDA
    
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

        cpu[i] = replicas[i].E();
        cpu_nc[i] = replicas[i].E(&replicas[i].molecules[0], &replicas[i].molecules[1]);
#if USING_CUDA
        gpu[i] = replicas[i].EonDevice();
        gpu_nc[i] = replicas[i].EonDeviceNC();
#endif
    }
}

void TestReferenceConformations::tearDown()
{
#if USING_CUDA
    teardown_CUDA(ljp_t);

    for (int i = 0; i < 10; i++)
    {
        replicas[i].FreeDevice();
    }
#endif // USING_CUDA
}

void TestReferenceConformations::testExactResults()
{
    // This test checks that nothing changes unexpectedly in our implementation of the potential calculation.
    // This should be tested more thoroughly (including the bonded components) in the Molecule and Replica unit tests.

    vector<Potential> expected_exact_results = {
        Potential(-0.0796731, -0.212202),
        Potential(-1.30422, 0.264889),
        Potential(-8.97201, -1.17831),
        Potential(-5.81976, -1.67169),
        Potential(-2.08631e-05, -5.76825e-05),
        Potential(-4.74252, -0.749355),
        Potential(-4.11734, -1.26302),
        Potential(-9.08724, -1.44068),
        Potential(-7.83475, -1.94273),
        Potential(-7.34204, -1.06544)
    };
    
    for (int i = 0; i < 10; i++) {
        ASSERT_POTENTIALS_EQUAL(expected_exact_results[i], cpu[i]);
        ASSERT_POTENTIALS_EQUAL(expected_exact_results[i], cpu_nc[i]);
#if USING_CUDA
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_exact_results[i].total(), gpu[i], 0.00001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_exact_results[i].total(), gpu_nc[i], 0.00001);
#endif
    }
}

void TestReferenceConformations::testVsKimResults()
{
    // This test compares our implementation to Young Kim's original results.
    // We need these fudge factors to compensate for different conversion factors, etc..
    
    double LJ_fudge_factor(0.6/RT_to_kcalmol); // Kim used 0.6 for the RT to kcal/mol conversion
    double DH_fudge_factor(1.00465); // we have no idea; maybe truncated constants?
    
    vector<Potential> kim_results = {
        Potential(-0.08092694, -0.21315828),
        Potential(-1.32250603,  0.26608889),
        Potential(-9.09456258, -1.18374210),
        Potential(-5.90453581, -1.67963588),
        Potential(-0.00002116, -0.00005795),
        Potential(-4.81153585, -0.75302819),
        Potential(-4.18365757, -1.26891119),
        Potential(-9.22300158, -1.44730149),
        Potential(-7.95237776, -1.95173336),
        Potential(-7.44768961, -1.07043531)
    };
    
    double mean_e_LJ(0);
    double mean_e_DH(0);
    double mean_e_total(0);
    
    double max_e_LJ(0);
    double max_e_DH(0);
    double max_e_total(0);
    
    for (int i = 0; i < 10; i++) {
        const double & kim_LJ = kim_results[i].total_LJ();
        const double & kim_DH = kim_results[i].total_DH();
        const double & kim_total = kim_results[i].total();
        
        double LJ = cpu[i].total_LJ() * LJ_fudge_factor;
        double DH = cpu[i].total_DH() * DH_fudge_factor;
        double total = LJ + DH;
        
        // relative error, as described in Tunbridge 2011, p.113
        
        double e_LJ = fabs(kim_LJ - LJ)/fabs(kim_LJ);
        double e_DH = fabs(kim_DH - DH)/fabs(kim_DH);
        double e_total = fabs(kim_total - total)/fabs(kim_total);
        
        cout << "Error " << i + 1 << ": LJ: " << e_LJ << " DH: " << e_DH << " total: " << e_total << endl;
        
        mean_e_LJ += e_LJ/10;
        mean_e_DH += e_DH/10;
        mean_e_total += e_total/10;
        
        max_e_LJ = max(max_e_LJ, e_LJ);
        max_e_DH = max(max_e_DH, e_DH);
        max_e_total = max(max_e_total, e_total);
    }
    
    cout << "Mean error: LJ: " << mean_e_LJ << " DH: " << mean_e_DH << " total: " << mean_e_total << endl;
    cout << "Max error: LJ: " << max_e_LJ << " DH: " << max_e_DH << " total: " << max_e_total << endl;
    
    // Apart from the LJ in conformation 6, we seem to be doing quite well. We should probably investigate that more closely.
    
    // TODO TODO TODO add the actual test, and clean up the printing so that it only happens if the flag is set. Is there a better way to do this?
    
}

void TestReferenceConformations::testVsCharmmResults()
{
    // TODO TODO TODO do this comparison too?
}
