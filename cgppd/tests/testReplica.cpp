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
    cout << "START" << endl;
    Potential potential_flexible = child_replica.E();
    cout << "END" << endl;
    potential_flexible.print_log(ALWAYS, "Flexible");

    child_replica.calculate_rigid_potential_only = true;
    Potential potential_rigid = child_replica.E();
    potential_rigid.print_log(ALWAYS, "Rigid");

#if USING_CUDA
    child_replica.FreeDevice();
#endif
}

