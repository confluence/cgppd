#include "definitions.h"
#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include <Simulation.h>
#include <testCommon.h>

class TestReplica : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestReplica);
    CPPUNIT_TEST(testCPUandGPUPotential);
    CPPUNIT_TEST_SUITE_END();

private:
    AminoAcids aminoAcidData;
    Replica initial_replica;
    Replica replica;
    float * ljp_t;
    float testboxdim;

public:
    void setUp();
    void testCPUandGPUPotential();
    void tearDown();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestReplica);

void TestReplica::setUp()
{
    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);
    testboxdim = 118.4f;

#if USING_CUDA
#if CUDA_STREAMS
    setup_CUDA(0, testboxdim, ljp_t, &aminoAcidData, &initial_replica.cudaStream, 1);
#else
    setup_CUDA(0, testboxdim, ljp_t, &aminoAcidData);
#endif // CUDA_STREAMS
#endif // USING_CUDA

    argdata parameters;
    parameters.bound = testboxdim;
#if USING_CUDA
    parameters.auto_blockdim = false;
#endif

    moldata m1;
    m1.translate = true;
    sprintf(m1.pdbfilename, "data/conf1/1a.pdb");

    moldata m2;
    m2.translate = true;
    sprintf(m2.pdbfilename, "data/conf1/1b.pdb");
    
    segdata seg1;
    seg1.residue_indices = {72, 73, 74, 75};
    m2.segments.push_back(seg1);
    
    parameters.mdata.push_back(m1);
    parameters.mdata.push_back(m2);

    initial_replica.init_first_replica(parameters, aminoAcidData, 2);
    initial_replica.label = 1;
}

void TestReplica::tearDown()
{
#if USING_CUDA
#if CUDA_STREAMS
    teardown_CUDA(ljp_t, &initial_replica.cudaStream, 1);
#else
    teardown_CUDA(ljp_t);
#endif // CUDA_STREAMS
#endif // USING_CUDA
    cout.flush();
}

void TestReplica::testCPUandGPUPotential()
{
    Replica replica;
    argdata parameters;
    replica.init_child_replica(initial_replica, 1, 300.0f, 0.2f, 0.5f, parameters);
#if USING_CUDA
#if CUDA_STREAMS
    replica.setup_CUDA(ljp_t, &initial_replica.cudaStream, 0);
#else
    replica.setup_CUDA(ljp_t);
#endif // CUDA_STREAMS
#endif // USING_CUDA

    // TODO: why is the internal LJ so huge?
#if !LJ_REPULSIVE && !LJ_OFF
    double flex_lj(296.795363);
#elif LJ_OFF
    double flex_lj(0);
#elif LJ_REPULSIVE
    double flex_lj(609.503791);
# endif // !LJ_REPULSIVE && !LJ_OFF
    
    Potential expected_flexible_potential(flex_lj, 0.335944, 0.180426, 1.555213, 5.420859);
    ASSERT_POTENTIALS_EQUAL(expected_flexible_potential, replica.E());

#if USING_CUDA
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_flexible_potential.total(), replica.EonDevice(), 0.001);
#endif // USING_CUDA

    replica.calculate_rigid_potential_only = true;
#if USING_CUDA
    replica.ReplicaDataToDevice();
#endif // USING_CUDA

#if !LJ_REPULSIVE && !LJ_OFF
    double rigid_lj(-0.079673);
#elif LJ_OFF
    double rigid_lj(0);
#elif LJ_REPULSIVE
    double rigid_lj(-0.275268);
# endif // !LJ_REPULSIVE && !LJ_OFF

    Potential expected_rigid_potential(rigid_lj, -0.212202, 0.000000, -0.000000, 0.000000);
    ASSERT_POTENTIALS_EQUAL(expected_rigid_potential, replica.E());

#if USING_CUDA
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_rigid_potential.total(), replica.EonDevice(), 0.001);
#endif // USING_CUDA

#if USING_CUDA
    replica.teardown_CUDA();
#endif
}

