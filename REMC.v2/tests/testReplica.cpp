#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>
#include <Replica.h>
#include <cutil.h>
#include <cutil_inline.h>
#include "definitions.h"

class TestReplica : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestReplica);
    CPPUNIT_TEST(testCopy);
    CPPUNIT_TEST_SUITE_END();

private:
    AminoAcids aminoAcidData;
    Replica replica;
    float * ljp_t;
    float testboxdim;

public:
    void setUp();
    void testCopy();
    void tearDown();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestReplica);

void TestReplica::setUp()
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

#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    bindLJTexture(ljp_t);
#endif
    replica.aminoAcids = aminoAcidData;
    replica.label = 1;
    replica.setBoundingValue(testboxdim);

    replica.reserveContiguousMoleculeArray(2);

    replica.loadMolecule("data/conf1/1a.pdb");
    replica.loadMolecule("data/conf1/1b.pdb");
}

void TestReplica::tearDown()
{
    if (replica.replicaIsOnDevice)
    {
        replica.FreeDevice();
    }
    cudaFree(ljp_t);
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    unbindLJTexture();
#endif
    cout.flush();
}

void TestReplica::testCopy()
{
//     replica.setDeviceLJPotentials(ljp_t);
//     replica.setBlockSize(TILE_DIM);
//     replica.ReplicaDataToDevice();
//
//     Replica replica_copy;
//     replica_copy.copy(*replica);
//
//     double original_E = replica.EonDevice();
//     double copy_E = replica_copy.EonDevice();
//
//     CPPUNIT_ASSERT_DOUBLES_EQUAL(original_E, copy_E, 0.000001);
}