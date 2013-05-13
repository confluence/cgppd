#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>

#include "definitions.h"
#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include <Replica.h>

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
    cutilCheckMsg("Failed to initialise CUDA runtime.");
    cudaMalloc((void**)&ljp_t,LJArraySize);
    cutilCheckMsg("Failed to allocate contact potential memory on the GPU");
    CUDA_setBoxDimension(testboxdim);
    cutilCheckMsg("Failed to copy box dimensions to GPU");
    copyLJPotentialDataToDevice(ljp_t,&aminoAcidData);
    cutilCheckMsg("Failed to code contact potentials to device.");
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    bindLJTexture(ljp_t);
#endif
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
    child_replica.init_child_replica(replica, 1, 300.0f, 0.2f, 0.5f, 1);
#if USING_CUDA
    child_replica.device_LJPotentials = ljp_t;
    child_replica.ReplicaDataToDevice();
#endif
    child_replica.potential = child_replica.E();
    child_replica.MCSearch(20);
#if USING_CUDA
    child_replica.FreeDevice();
#endif
}