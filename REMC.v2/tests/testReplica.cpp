#include <UnitTest++.h>
#include <Replica.h>
#include <cutil.h>
#include <cutil_inline.h>

struct ReplicaFixture
{
    ReplicaFixture()
    {
        cuInit(0);
        cout.precision(8);

        aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
        aminoAcidData.loadLJPotentialData(LJPDSOURCE);

        // get values for the example conformations
        egnames = new char[60];

        cudaMalloc((void**)&ljp_t,LJArraySize);
        cutilCheckMsg("Failed to cudaMalloc");
        copyLJPotentialDataToDevice(ljp_t,&aminoAcidData);

        // set box dimensions
        testboxdim = 118.4f;
        CUDA_setBoxDimension(testboxdim);

#if LJ_LOOKUP_METHOD == TEXTURE_MEM
        bindLJTexture(ljp_t);
#endif
    }

    ~ReplicaFixture()
    {
        cudaFree(ljp_t);
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
        unbindLJTexture();
#endif
        cout.flush();
    }

    AminoAcids aminoAcidData;
    Replica exampleReplicas[10];
    char *egnames;
    float * ljp_t;
    float testboxdim;
};

TEST_FIXTURE(ReplicaFixture, TestReplica)
{
    for (int i = 1; i<=10; i++)
    {
        exampleReplicas[i-1].aminoAcids = aminoAcidData;
        exampleReplicas[i-1].label = i;
        exampleReplicas[i-1].setBoundingValue(testboxdim);
        sprintf(egnames,"data/conf%d/1a.pdb",i);
        exampleReplicas[i-1].loadMolecule(egnames);
        sprintf(egnames,"data/conf%d/1b.pdb",i);
        exampleReplicas[i-1].loadMolecule(egnames);
        exampleReplicas[i-1].initTimers();
        uint cpu_E_timer;
        exampleReplicas[i-1].countNonCrowdingResidues();
        exampleReplicas[i-1].E();

        exampleReplicas[i-1].setDeviceLJPotentials(ljp_t);
        exampleReplicas[i-1].setBlockSize(TILE_DIM);
        exampleReplicas[i-1].ReplicaDataToDevice();
        ///cout << "CPU execution time: " << cutGetTimerValue(cpu_E_timer)/100.0f << "ms"<<  endl;
        uint gpu_E_timer;
        double r = exampleReplicas[i-1].EonDevice();
        double r2 = exampleReplicas[i-1].EonDeviceNC();

        double cpu_ee = exampleReplicas[i-1].E();
        double cpu_ee2 = exampleReplicas[i-1].E(&exampleReplicas[i-1].molecules[0],&exampleReplicas[i-1].molecules[1]);

        printf ("CPU (CPUnc), GPU (NC), diff(diff): %13f (%13f), %13f (%13f)    \n",cpu_ee,cpu_ee2,r,r2,abs(cpu_ee-r),abs(cpu_ee2-r2));
        //exampleReplicas[i-1].printTimers();
        exampleReplicas[i-1].FreeDevice();
    }
}
