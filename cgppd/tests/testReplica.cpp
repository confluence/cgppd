#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>

#include "definitions.h"
#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include <Simulation.h>

class TestReplica : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestReplica);
    CPPUNIT_TEST(testMC);
    CPPUNIT_TEST(testCPUandGPUPotential);
    CPPUNIT_TEST_SUITE_END();

private:
    AminoAcids aminoAcidData;
    Replica replica;
    float * ljp_t;
    float testboxdim;

public:
    void setUp();
    void testMC();
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

    replica.loadMolecule("data/conf1/1a.pdb");
    replica.loadMolecule("data/conf1/1b.pdb");

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

void TestReplica::testMC()
{
    Replica child_replica;
    argdata parameters;
    child_replica.init_child_replica(replica, 1, 300.0f, 0.2f, 0.5f, parameters);
#if USING_CUDA
    child_replica.device_LJPotentials = ljp_t;
    child_replica.ReplicaDataToDevice();
#endif
    child_replica.potential = child_replica.E().total();
    child_replica.MCSearch(20, 0);
#if USING_CUDA
    child_replica.FreeDevice();
#endif
}

void TestReplica::testCPUandGPUPotential()
{
    //TODO Test CUDA+rigid molecules vs CPU+rigid molecules vs reference conformation 1
    //TODO Test CUDA+flexible molecules vs CPU+flexible molecules vs hand-calculated reference
}