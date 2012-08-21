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
    replica.loadMolecule("data/conf1/1b.pdb");
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

// TODO: this is pointless; remove it
void TestReplica::testCopy()
{
    Replica replica_copy;
    replica_copy.copy(replica);

    CPPUNIT_ASSERT_EQUAL(replica.label, replica_copy.label);
    CPPUNIT_ASSERT_EQUAL(replica.temperature, replica_copy.temperature);
    CPPUNIT_ASSERT_EQUAL(replica.maxMoleculeSize, replica_copy.maxMoleculeSize);
    CPPUNIT_ASSERT_EQUAL(replica.boundingValue, replica_copy.boundingValue);
    CPPUNIT_ASSERT_EQUAL(replica.residueCount, replica_copy.residueCount);
    CPPUNIT_ASSERT_EQUAL(replica.moleculeCount, replica_copy.moleculeCount);
    CPPUNIT_ASSERT_EQUAL(replica.nonCrowderCount, replica_copy.nonCrowderCount);
    CPPUNIT_ASSERT_EQUAL(replica.nonCrowderResidues, replica_copy.nonCrowderResidues);
    CPPUNIT_ASSERT_EQUAL(replica.potential, replica_copy.potential);
#if USING_CUDA
    CPPUNIT_ASSERT_EQUAL(replica.blockSize, replica_copy.blockSize);
    CPPUNIT_ASSERT_EQUAL(replica.sharedMemSize, replica_copy.sharedMemSize);
#endif
}

void TestReplica::testMC()
{
    Replica child_replica;
    child_replica.init_child_replica(replica, 1, 300.0f, 0.2f, 0.5f, 1);
#define LOGLEVEL DEBUG
    child_replica.MCSearch(20);
}