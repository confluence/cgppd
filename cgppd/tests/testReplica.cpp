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
    setup_CUDA(0, testboxdim, ljp_t, &aminoAcidData);
#endif

    replica.aminoAcids = aminoAcidData;
    replica.label = 1;
    replica.boundingValue = testboxdim;

    replica.reserveContiguousMoleculeArray(2);

    replica.loadMolecule("data/conf1_flexible/1a.pdb");
    replica.loadMolecule("data/conf1_flexible/1b.pdb");

#if USING_CUDA
    int cuda_blockSize;
    (replica.residueCount < 1024) ? cuda_blockSize = 32 : cuda_blockSize = 64;
    replica.setBlockSize(cuda_blockSize);
#endif

}

void TestReplica::tearDown()
{
#if USING_CUDA
    cudaFree(ljp_t);
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    unbindLJTexture();
#endif
#endif
    cout.flush();
}

void TestReplica::testCPUandGPUPotential()
{
    Replica child_replica;
    argdata parameters;
    child_replica.init_child_replica(replica, 1, 300.0f, 0.2f, 0.5f, parameters);
#if USING_CUDA
    child_replica.device_LJPotentials = ljp_t;
    child_replica.ReplicaDataToDevice();
#endif
    double expected_flexible_potential[6] = {300.864756, 0.337112, 0.180426, 1.555213, 3.495782, 306.433290};
    ASSERT_POTENTIAL_EQUALS(expected_flexible_potential, child_replica.E());
#if USING_CUDA
    CPPUNIT_ASSERT_DOUBLES_EQUAL(306.433290, child_replica.EonDevice(), 0.00001);
#endif

    child_replica.calculate_rigid_potential_only = true;
    double expected_rigid_potential[6] = {-0.080765, -0.212940, 0.000000, -0.000000, 0.000000, -0.293705};
    ASSERT_POTENTIAL_EQUALS(expected_rigid_potential, child_replica.E());
#if USING_CUDA
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.293705, child_replica.EonDevice(), 0.00001);
#endif

#if USING_CUDA
    child_replica.FreeDevice();
#endif
}

