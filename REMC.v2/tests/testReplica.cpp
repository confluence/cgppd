#include <UnitTest++.h>
#include <Replica.h>
#include <cutil.h>
#include <cutil_inline.h>
#include "definitions.h"

struct ReplicaFixture
{
    AminoAcids aminoAcidData;
    Replica replica;
    float * ljp_t;
    float testboxdim;

    ReplicaFixture()
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

    ~ReplicaFixture()
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
};

// TEST_FIXTURE(ReplicaFixture, TestReplicaCopy)
// {
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
//     CHECK_CLOSE(original_E, copy_E, 0.000001);
// }
//
//
