#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>
#include <Replica.h>
#include <cutil.h>
#include <cutil_inline.h>
#include "definitions.h"

class TestReplica : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestReplica);
    CPPUNIT_TEST(testMC);
    CPPUNIT_TEST_SUITE_END();

private:
    AminoAcids aminoAcidData;
    Replica replica;
    float * ljp_t;
    float testboxdim;

public:
    void setUp();
    void testMC();
    void tearDown();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestReplica);

void TestReplica::setUp()
{
    aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
    aminoAcidData.loadLJPotentialData(LJPDSOURCE);
    testboxdim = 118.4f;

#if USING_CUDA
    cuInit(0);
    cudaMalloc((void**)&ljp_t,LJArraySize);
    cutilCheckMsg("Failed to cudaMalloc");
    copyLJPotentialDataToDevice(ljp_t,&aminoAcidData);

    // set box dimensions
    CUDA_setBoxDimension(testboxdim);

#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    bindLJTexture(ljp_t);
#endif
#endif

    replica.aminoAcids = aminoAcidData;
    replica.label = 1;
    replica.setBoundingValue(testboxdim);

    replica.reserveContiguousMoleculeArray(2);

    replica.loadMolecule("data/conf1/1a.pdb");
//     replica.loadMolecule("data/conf1/1b.pdb");
    replica.loadMolecule("tests/1UBQ.pdb");
}

void TestReplica::tearDown()
{
#if USING_CUDA
    if (replica.replicaIsOnDevice)
    {
        replica.FreeDevice();
    }
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
    child_replica.init_child_replica(replica, 1, 300.0f, 0.2f, 0.5f, 1);
    child_replica.potential = child_replica.E();
    child_replica.MCSearch(20);
}
