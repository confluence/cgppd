// Replica Kernel functions to call cuda kernels for energy calculations in a replica
#ifndef CUDA_FUNCTIONS_CU
#define CUDA_FUNCTIONS_CU

using namespace std;

#include "CudaFunctions.h"

// Block index
#define bx blockIdx.x
#define by blockIdx.y

// Thread index
#define tx threadIdx.x
#define ty threadIdx.y

#if USING_CUDA

// special sauce function for throwing an exception when there's a CUDA error instead of exiting
void __throwLastCudaError(const char *errorMessage, const char *file, const int line)
{
    cudaError_t err = cudaGetLastError();

    if (cudaSuccess != err)
    {
        fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
                file, line, errorMessage, (int)err, cudaGetErrorString(err));
        DEVICE_RESET
        throw errorMessage;
        //exit(EXIT_FAILURE);
    }
}

#if LJ_LOOKUP_METHOD == TEXTURE_MEM

void bindLJTexture(float * ljp)
{
    size_t offset;
    cudaBindTexture(&offset, LJTexture, ljp, LJArraySize);
    if (int(offset)>0)
        printf("!!! Texture offset %d\n",int(offset));
}

void unbindLJTexture()
{
    cudaUnbindTexture(LJTexture);
}

#endif

#if METADATA_MEMORY == TEXTURE_MEM
int bindMetaDataToTexture(void* deviceMemory, size_t size)
{
    return cudaBindTexture(0, residueMetaTex, deviceMemory, size);
}

int freeMetaDataTexture()
{
    return cudaUnbindTexture(residueMetaTex);
}
#endif

#if POSITIONDATA_MEMORY == TEXTURE_MEM
int bindPositionDataToTexture(void* deviceMemory, size_t size)
{
    return cudaBindTexture(0, residuePositionTex, deviceMemory, size);
}

int freePositionDataTexture()
{
    return cudaUnbindTexture(residuePositionTex);
}
#endif

float * LJPotentialDataToDevice (AminoAcids *a)
{
    long tableSize(sizeof(float[AA_COUNT*AA_COUNT]));
    long containerSize(tableSize);  // pad the container such that shared memory tables align nicely when preloading

#if LJ_LOOKUP_METHOD == SHARED_MEM
    containerSize = TILE_DIM * ceil(float(tableSize)/float(TILE_DIM));
#endif

    float * safeData = new float[AA_COUNT*AA_COUNT];

    for (int i=0; i<AA_COUNT; i++)
    {
        memcpy (&safeData[i*AA_COUNT],&a->LJpotentials[i][0],sizeof(float)*AA_COUNT);
    }

#if LJ_LOOKUP_METHOD == CONST_MEM
    cudaMemcpyToSymbol(const_LJPotentialData, &safeData, sizeof(tableSize));
    delete [] safeData;
    return LJPotentialData;
#else
    float * _LJPotentialData;  // global device memory
    checkCudaErrors(cudaMalloc((void**)&_LJPotentialData,containerSize));
    checkCudaErrors(cudaMemset(_LJPotentialData,0,containerSize));
    // copy host memory to device
    checkCudaErrors(cudaMemcpy(_LJPotentialData, safeData, tableSize, cudaMemcpyHostToDevice) );
    delete [] safeData;
    return _LJPotentialData;
#endif
}

void copyLJPotentialDataToDevice (float * dev_LJPotentialData, AminoAcids *a)
{
#define mem_size int(sizeof(float)*AA_COUNT*AA_COUNT)
#define table_size (AA_COUNT*AA_COUNT)
    float * safeData = new float[table_size];

    for (int i=0; i<table_size; i++)
        safeData[i] = 0.0f;

    for (int i=0; i<AA_COUNT; i++)
    {
        memcpy (&safeData[i*AA_COUNT],a->LJpotentials[i],sizeof(float)*AA_COUNT);
    }

    //CUDA_SAFE_CALL(cudaMalloc((void**)&dev_LJPotentialData,mem_size));
    // copy host memory to device
    cudaMemcpy(dev_LJPotentialData, safeData, mem_size, cudaMemcpyHostToDevice);
    getLastCudaError("Failed to copy contact potentials");
    delete [] safeData;
}

//void cudaInfo()
//{
    //int deviceCount;
    //checkCudaErrors(cudaGetDeviceCount(&deviceCount)); // this used to be CUDA_SAFE_CALL_NO_SYNC -- what's the difference?
    //if (deviceCount == 0)
    //{
        //printf("ERROR: no devices supporting CUDA.\n");
        //exit(EXIT_FAILURE);
    //}
//#if OUTPUT_LEVEL > 0
    //struct cudaDeviceProp devInfo;
    //printf("\n------------------------------------\n");
    //for (int i=0; i<deviceCount; i++)
    //{
        //cudaGetDeviceProperties( &devInfo, 0 );
        //printf ("CUDA Device %d Info:\n",i);
        //printf("Name: %s\n", devInfo.name );
        //printf("totalGlobalMem: %i\n",(int)devInfo.totalGlobalMem );
        //printf("sharedMemPerBlock: %i\n",(int)devInfo.sharedMemPerBlock );
        //printf("regsPerBlock: %i\n",devInfo.regsPerBlock );
        //printf("warpSize: %i\n",devInfo.warpSize );
        //printf("memPitch %i\n",(int)devInfo.memPitch );
        //printf("maxThreadsPerBlock: %i\n",devInfo.maxThreadsPerBlock);
        //if (i < deviceCount-1)
            //printf("\n");
    //}
    //printf("------------------------------------\n");
//#endif
//}

void CUDA_setBoxDimension(float value)
{
    cudaMemcpyToSymbol(const_boxdim, &value, sizeof(value));
    getLastCudaError ("Set boxdim");

};
void CUDA_freeBoxDimension()
{
    cudaFree(&const_boxdim);
};



//void MCSearchOnDevice()
//{
    //printf("MCSearchOnDevice() does nothing\n");
//};

// stream kernel for overlapping calls
void CUDA_EonDevice_async(float4 *residuePositions, float4 *residueMeta, int * residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, float* kernelResult, int resultSize, int blockSize, int datasetSize, int sm_size, cudaStream_t stream)
{
    dim3 dimBlock(blockSize,1,1);   // threads
    dim3 dimGrid(datasetSize/blockSize,datasetSize/blockSize,1);    	// blocks

    // Launch the device computation
    // result must become an array of dimentions the same as grid so that the multiple thread blocks can sum in it properly
    //cudaStreamSynchronize(stream);

    E_TiledKernel<<< dimGrid,dimBlock,sm_size,stream >>>(residuePositions, residueMeta, residueCount, moleculePositions, moleculeCount, LJPotentials, kernelResult);

    getLastCudaError("Kernel execution failed");
    //cudaMemcpyAsync(hostResult,kernelResult, sizeof(float)*resultSize*resultSize, cudaMemcpyDeviceToHost,stream);
    return;
};

//void CUDA_Esum_async(float* result, float *d_resultMatrix, int resultSize, int datasetSize, cudaStream_t stream)
//{
    //// parallel sum on gpu if required
//}

// alternative summation algorithm will work with up to 33554432 residues
void CUDA_EonDevice(float4 *residuePositions, float4 *residueMeta, int * residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, double* result, int blockSize, int datasetSize, int sm_size)
{

    // result stored on the device
    // gridSize can be arbitrary
    int gridSize = datasetSize/blockSize;
    // the parallel sum of the results must have an array with size that is a power of 2
    int resultSize = gridSize; //int(pow(2,ceil(log(float(gridSize))/log(2.0))));
    float *tmpSums = new float[resultSize*resultSize];

    //for a parallel sum each grid must have one cell in the array of results from all the threads
    float *d_result;
    cudaMalloc((void **)&d_result,sizeof(float)*resultSize*resultSize);
    cudaMemset(d_result,0,sizeof(float)*resultSize*resultSize);

    //float4 testmem;
    //cudaMemcpy(&testmem,residuePositions, sizeof(float4), cudaMemcpyDeviceToHost);


    // Launch the device computation
    // result must become an array of dimensions the same as grid so that the multiple thread blocks can sum in it properly

    dim3 dimBlock(blockSize,1,1);   		// threads
    dim3 dimGrid(gridSize,gridSize,1);    	// blocks


    E_TiledKernel<<< dimGrid,dimBlock,sm_size,0 >>>(residuePositions, residueMeta, residueCount, moleculePositions, moleculeCount, LJPotentials, d_result);
    // check if kernel invocation generated an error

#ifdef _EMU

    printf("After kernel:\n");
    for (int i=0; i<resultSize; i++)
    {
        for (int j=0; j<resultSize; j++)
            printf (" %10.6f",d_result[i*resultSize+j]);
        printf("\n");
    }
#endif

#if PARALLEL_SUM
    parallelSum_kernel<<<1,gridSize>>>(d_result);
    cudaMemcpy(result,d_result, sizeof(float), cudaMemcpyDeviceToHost);
#else
    cudaMemcpy(tmpSums,d_result, sizeof(float)*resultSize*resultSize, cudaMemcpyDeviceToHost);
    result[0] = 0.0f;
#if COMPENSATE_KERNEL_SUM
    float c(0.0f);
    /*
    		function kahanSum(input)
    		 var potential = 0
    		 var c = 0
    		 for i = 0 to blockdim-1
    		  y = lj+dh - c
    		  t = potential + y
    		  c = (t - potential) - y
    		  potential = t
    		 next i
    		return sum
    		*/
#endif
    for (int i=0; i<resultSize*resultSize; i++)
    {
#if COMPENSATE_KERNEL_SUM
        float y(tmpSums[i] -c);
        float t((*result) + y);
        c = (t-(*result)) -y;
        *result = t;
#else
    result[0] += tmpSums[i];
#endif
    }

    delete [] tmpSums;
#endif
    // free the memory assigned for this iteration.
    cudaFree(d_result);
    getLastCudaError("Kernel execution failed");
    return;
};

//void CUDA_EonDeviceTest(float *d_x, float *d_y,float *d_z, int *d_id, float4 *residueMeta, float* LJPotentials, float* result, int blockSize, int datasetSize)
//{

    //// result stored on the device
    //// gridSize can be arbitrary
    //int gridSize(datasetSize/blockSize);
    //float *tmpSums = new float[gridSize*gridSize];

    ////for a parallel sum each grid must have one cell in the array of results from all the threads
    //float *d_result;
    //cudaMalloc((void **)&d_result,sizeof(float)*gridSize*gridSize);
    //cudaMemset(d_result,0,sizeof(float)*gridSize*gridSize);

    ////float4 testmem;
    ////cudaMemcpy(&testmem,residuePositions, sizeof(float4), cudaMemcpyDeviceToHost);


    //// Launch the device computation
    //// result must become an array of dimensions the same as grid so that the multiple thread blocks can sum in it properly

    //dim3 dimBlock(blockSize,1,1);   		// threads
    //dim3 dimGrid(gridSize,gridSize,1);    	// blocks

    //printf ("+++++++++++++++++++ NO KERNEL EXECUTED +++++++++++++++++++++++++++++");
    ////E_TestTiledKernel<<< dimGrid,dimBlock >>>(d_x,d_y,d_z,d_id, residueMeta, LJPotentials, d_result);

    //cudaMemcpy(tmpSums,d_result, sizeof(float)*gridSize*gridSize, cudaMemcpyDeviceToHost);
    //result[0] = 0.0f;
    //for (int i=0; i<gridSize*gridSize; i++)
    //{
        //result[0] += tmpSums[i];
    //}
    //delete [] tmpSums;

    //// free the memory assigned for this iteration.
    //cudaFree(d_result);
    //// check if kernel invocation generated an error
    //getLastCudaError("Kernel execution failed");

    //return;
//};

__device__ float length(const float4 a, const float4 b)
{
    float Xab(a.x-b.x);
    float Yab(a.y-b.y);
    float Zab(a.z-b.z);

    Xab = Xab - const_boxdim * rintf(Xab/const_boxdim);
    Yab = Yab - const_boxdim * rintf(Yab/const_boxdim);
    Zab = Zab - const_boxdim * rintf(Zab/const_boxdim);

    return sqrtf(Xab*Xab+Yab*Yab+Zab*Zab);

};

__device__ float crowderPotential(float r)
{
    return __powf(6.0f/r,12.0f);
};

__device__ float dhPotential(float q1, float q2, float r)
{
    return  q1 * q2 *__expf(-r/Xi) / r ;
};

///* full kernel with optional texture containers
__global__ void E_TiledKernel(float4 * residuePositions, float4 * residueMeta, int * residueCount, int * moleculePositions, int * moleculeCount, float* LJPotentialData, float* result)
{
    extern __shared__ float sharedmem_results[];	// array of temp results to sum, one per thread run.
#if POSITIONDATA_MEMORY == SHARED_MEM
    const int positionOfPosition = blockDim.x/4;
    extern __shared__ float4 X_tile_residuePositions[];
#endif
#if METADATA_MEMORY == SHARED_MEM
#if POSITIONDATA_MEMORY == SHARED_MEM
    const int positionOfMeta = __mul24(5,blockDim.x)/4;
#else
    const int positionOfMeta = blockDim.x/4;
#endif
    extern __shared__ float4 X_tile_residueMeta[];
#endif

#if LJ_LOOKUP_METHOD == SHARED_MEM_LJ
    __shared__ float sharedLJ[AA_COUNT*AA_COUNT];

    if(tx==0)
    {
        for (int i=0; i<AA_COUNT*AA_COUNT; i+=AA_COUNT)
            memcpy(&(sharedLJ[i]),&(LJPotentialData[i]),sizeof(float)*AA_COUNT);
    }
#endif

    sharedmem_results[tx] = 0.0f;//EPS;
#if  CULL_LOWER_TRIANGULAR_SUM  // only works if we throw away the entire block, otherwise results cannot be guaranteed correct due to the parallel sum
    if (bx < by)
        return;
#endif

#if POSITIONDATA_MEMORY == TEXTURE_MEM
    float4 yresiduep = tex1Dfetch(residuePositionTex, by*blockDim.x+tx);
#else
    float4 yresiduep = residuePositions[by*blockDim.x+tx];
#endif
#if METADATA_MEMORY == TEXTURE_MEM
    float4 yresiduem = tex1Dfetch(residueMetaTex, by*blockDim.x+tx);
#else
    float4 yresiduem = residueMeta[by*blockDim.x+tx];
#endif
    // think of the position as a column of residues so 0 == by:0,tx:0, N = by:gridDimY*bDX+tx

    // copy residues from GRAM to shared memory tiles

#if POSITIONDATA_MEMORY == SHARED_MEM
    X_tile_residuePositions[positionOfPosition + tx] = residuePositions[bx*blockDim.x+tx];
#endif
#if METADATA_MEMORY == SHARED_MEM
    X_tile_residueMeta[positionOfMeta + tx] = residueMeta[bx*blockDim.x+tx];
#endif

    __syncthreads();  // very important so that the entire shared mem array is assigned before any threads use it

    if (yresiduep.w > PADDER_IDENTIFIER) // the residue in question is not a padding element
    {
        float lj_subtotal(0.0f);
        float dh_subtotal(0.0f);
#if COMPENSATE_KERNEL_SUM
        float c_lj(0.0f);   // Kahan summation variable
        float c_dh(0.0f);
        float y;
        float t;

        //		function kahanSum(input)
        //		 var potential = 0
        //		 var c = 0
        //		 for i = 0 to blockdim-1
        //		  y = lj+dh - c
        //		  t = potential + y
        //		  c = (t - potential) - y
        //		  potential = t
        //		 next i
        //		return sum

#endif
        //#pragma unroll 4 - cant unroll dynamic variable
        for (int i = 0; i < blockDim.x; i++)
        {
#if POSITIONDATA_MEMORY == TEXTURE_MEM
            float4 pos = tex1Dfetch(residuePositionTex,(bx*blockDim.x+i));
#elif POSITIONDATA_MEMORY == SHARED_MEM
            float4 pos = X_tile_residuePositions[positionOfPosition + i];
#endif
            // pos.x = residue X coordinate
            // pos.y = residue Y coordinate
            // pos.z = residue Z coordinate
            // FLEXIBLE: pos.w = PADDER_IDENTIFIER or DOMAIN_UID.BOND_UID
            // RIGID: pos.w = PADDER_IDENTIFIER or molecule id

#if METADATA_MEMORY == TEXTURE_MEM
            float4 meta = tex1Dfetch(residueMetaTex,(bx*blockDim.x+i));
#elif METADATA_MEMORY == SHARED_MEM
            float4 meta = X_tile_residueMeta[positionOfMeta + i];
#endif
            // meta.x = aminoacid index
            // meta.y = charge
            // meta.z = vdw radius
            // FLEXIBLE: meta.w = CROWDER_IDENTIFIER or MOLECULE_ID.CHAIN_UID
            // RIGID: meta.w = CROWDER_IDENTIFIER or unused
#if FLEXIBLE_LINKS && !ASSUME_POLYMER_FOLDING_TEST
            // padder or same rigid domain or same bond or residues close on the backbone
            // currently ignoring case where residues are close because of a bond; subtracting this component on the CPU instead
            
            // we unpack the floats

            float xdomain(0.0f);
            float ydomain(0.0f);

            float xbond = modff(pos.w, &xdomain);
            float ybond = modff(yresiduep.w, &ydomain);
            
            float xmol(0.0f); // NOT USED IN THIS KERNEL; USED IN NC KERNEL
            float ymol(0.0f); // NOT USED IN THIS KERNEL; USED IN NC KERNEL
            
            float xchain = modff(meta.w, &xmol);
            float ychain = modff(yresiduem.w, &ymol);

            float xresid = bx * blockDim.x + i; // TODO: IS THIS RIGHT?
            float yresid = by * blockDim.x + tx; // TODO: IS THIS RIGHT?
            
            if (pos.w == PADDER_IDENTIFIER || (xdomain && xdomain == ydomain) || (xbond && xbond == ybond) || (xchain == ychain && fabs(xresid - yresid) < 4))
#elif ASSUME_POLYMER_FOLDING_TEST
            // More efficient handling of this special case.  Assume we are folding a single polymer.
            // There's only one chain; we calculate all pairs except close neighbours on the backbone.
            if (fabs(meta.w - yresiduem.w) < 4)
#else
             // same molecule or padder
            if (yresiduep.w == pos.w || pos.w == PADDER_IDENTIFIER)
#endif
            {
                // DO NOTHING
            }
            else
            {

                float r(length(yresiduep,pos) + EPS);  // add eps so that r is never 0, can happen in a collision

#if USE_POTENTIAL_CUTOFF
                if (r>interaction_cutoff) continue;
#endif

                float LJ(0.0f);
                float DH(0.0f);

#if REPULSIVE_CROWDING && !LJ_OFF
                if (yresiduem.w == CROWDER_IDENTIFIER || meta.w == CROWDER_IDENTIFIER) // repulsive crowder interaction
                {
                    if (r<const_repulsive_cutoff)
                        LJ = crowderPotential(r);
                }
                else  // normal LJ interaction
                {
#endif
#if !LJ_OFF
                    int ijX(rint( AA_COUNT*yresiduem.x + meta.x));
                    //do the texture fetch first
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
                    float Eij(LJ_lambda*(tex1Dfetch(LJTexture,ijX) - e0));
#elif LJ_LOOKUP_METHOD == SHARED_MEM
                    float Eij(LJ_lambda*(sharedLJ[ijX] - e0));
#elif LJ_LOOKUP_METHOD == CONST_MEM
                    float Eij(LJ_lambda*(const_LJPotentialData[ijX] - e0));
#else  // __global__ or __constant__
                    float Eij(LJ_lambda*(LJPotentialData[ijX] - e0));
#endif // LJ_LOOKUP_METHOD
#endif // !LJ_OFF
                    DH = dhPotential(yresiduem.y,meta.y,r);
                    dh_subtotal += DH;
#if !LJ_OFF
                    // sigmaij is the average atomic radius determined by the van der waals radius in kim2008
                    float sigmaij((yresiduem.z + meta.z) * 0.5f);

                    float LJtmp(__powf(sigmaij/r,6.0f));

                    LJ = -4.0f*Eij*LJtmp*(LJtmp-1.0f);
                    //LJ = __fmul_rn(-4.0f,__fmul_rn(Eij,__fmul_rn(LJtmp,__fadd_rn(LJtmp,-1.0f))));
                    if (Eij>0.0f && r < (sigmaij*r0_constant))
                    {
                        LJ = -LJ + 2.0f*Eij;
                    }
#endif // !LJ_OFF
#if REPULSIVE_CROWDING && !LJ_OFF
                }  // end conditional branch for LJ or repulsive short-range energy
#endif
#if !LJ_OFF
#if COMPENSATE_KERNEL_SUM
                y = LJ - c_lj;
                t = lj_subtotal + y;
                c_lj = (t-lj_subtotal)-y;
                lj_subtotal = t;
#else
                lj_subtotal += LJ;
#endif
#endif // !LJ_OFF
            } // if !(X_tile_residuePositions.w == Y_tile_residuePositions.w || X_tile_residuePositions.w < CROWDER_IDENTIFIER )
        } // for i = 0..Bdx
        sharedmem_results[tx] = (lj_subtotal * RT_to_kcalmol) + (dh_subtotal * DH_CONVERSION_FACTOR);

#if CULL_LOWER_TRIANGULAR_SUM
        if (bx==by)
            sharedmem_results[tx] = sharedmem_results[tx]*0.5f;
#else
        sharedmem_results[tx] *= 0.5f
#endif //  CULL_LOWER_TRIANGULAR_SUM
    } // if !padding

    __syncthreads();  // all threads must have computed values so that they can be read concurrently from shared memory

    // do reduction, see SDK reduction example
    // FULLY unrolled reduction

    // do reduction in shared mem
    if (blockDim.x >= 512) {
        if (tx < 256) {
            sharedmem_results[tx] += sharedmem_results[tx + 256];
        }
        __syncthreads();
    }
    if (blockDim.x >= 256) {
        if (tx < 128) {
            sharedmem_results[tx] += sharedmem_results[tx + 128];
        }
        __syncthreads();
    }
    if (blockDim.x >= 128) {
        if (tx <  64) {
            sharedmem_results[tx] += sharedmem_results[tx + 64];
        }
        __syncthreads();
    }

    if (tx < 32)
    {
        if (blockDim.x >=  64) {
            sharedmem_results[tx] += sharedmem_results[tx + 32];
        }
        if (blockDim.x >=  32) {
            sharedmem_results[tx] += sharedmem_results[tx + 16];
        }
        if (blockDim.x >=  16) {
            sharedmem_results[tx] += sharedmem_results[tx +  8];
        }
        if (blockDim.x >=   8) {
            sharedmem_results[tx] += sharedmem_results[tx +  4];
        }
        if (blockDim.x >=   4) {
            sharedmem_results[tx] += sharedmem_results[tx +  2];
        }
        if (blockDim.x >=   2) {
            sharedmem_results[tx] += sharedmem_results[tx +  1];
        }
    }

    // write results of the calculation back to global memory, only thread 0 has the correct answer
    // only one warp processes < 32, therefore no syncthreads required here
    if (tx == 0)
        result[bx*gridDim.x+by] = sharedmem_results[tx];

    return;
};


//void CUDA_setMoleculeBlength(int length)
//{
    //cudaMemcpyToSymbol(const_lenght_of_b, &length, sizeof(length));
    //getLastCudaError ("CUDA_setMoleculeBlength: Set b size");
//}

//// kernel to evaluate the potential between two molecules.
//// assigns a tread to each residue in molecule A, and iterates through molecule B
//// TODO: check it works
//__global__ void E_MoleculeKernel(float4 * residuePositionsA, float4 * residueMetaA, float4 * residuePositionsB, float4 * residueMetaB, float* LJPotentialData, float* result )
//{
    //extern __shared__ float sharedmem_results[];	// array of temp results to sum, one per thread run.

    //const int positionOfPosition = blockDim.x/4;
    //const int positionOfMeta = __mul24(5,blockDim.x)/4;
    //extern __shared__ float4 B_Positions[];
    //extern __shared__ float4 B_Meta[];

    //sharedmem_results[tx] = 0.0f;

    //float4 A_pos = residuePositionsA[bx*blockDim.x+tx];
    //float4 A_meta = residueMetaA[bx*blockDim.x+tx];

    //// copy residues from GRAM to shared memory tiles

    //for( int tile(0); tile < const_lenght_of_b/blockDim.x; tile++)
    //{
        //B_Positions[positionOfPosition + tx] = residuePositionsB[tile*blockDim.x+tx];
        //B_Meta[positionOfMeta + tx] 		 = residueMetaB[tile*blockDim.x+tx];

        //__syncthreads();  // very important so that the entire shared mem array is assigned before any threads use it

        //if (A_pos.w > PADDER_IDENTIFIER && B_Positions[positionOfPosition + tx].w > PADDER_IDENTIFIER) // the residue in question is not a padding element
        //{
            //float lj_subtotal(0.0f);
            //float dh_subtotal(0.0f);

            //for (int i = 0; i < blockDim.x; ++i)
            //{
                //int ijX(rint( AA_COUNT*A_meta.x + B_Meta[positionOfMeta + i].x));
                //float r(length(A_pos,B_Positions[positionOfPosition + i]) + EPS);  // add eps so that r is never 0, can happen in a collision
                //// get contacts for short range, do long while it loads
                //float Eij(LJ_lambda*(tex1Dfetch(LJTexture,ijX) - e0));
                //// do long range
                //float DH(dhPotential(A_meta.y,B_Meta[positionOfMeta + i].y,r));
                //dh_subtotal += DH;
                ////continue doing short range
                //float sigmaij((A_meta.z + B_Meta[positionOfMeta + i].z) * 0.5f);
                //float LJtmp(__powf(sigmaij/r,6.0f));
                //float LJ (-4.0f*Eij*LJtmp*(LJtmp-1.0f));
                //if (Eij>0.0f && r < (sigmaij*r0_constant))
                    //LJ = -LJ + 2.0f*Eij;
                //lj_subtotal += LJ;

            //} // for i = 0..Bdx
            //sharedmem_results[tx] = (lj_subtotal * RT_to_kcalmol) + (dh_subtotal * DH_CONVERSION_FACTOR);
        //} // if !padding

        //__syncthreads();  // make sure all threads are done before moving to the enxt tile
    //}
    //// do reduction, see SDK reduction example
    //// FULLY unrolled reduction

    //// do reduction in shared mem
    //if (blockDim.x >= 512) {
        //if (tx < 256) {
            //sharedmem_results[tx] += sharedmem_results[tx + 256];
        //}
        //__syncthreads();
    //}
    //if (blockDim.x >= 256) {
        //if (tx < 128) {
            //sharedmem_results[tx] += sharedmem_results[tx + 128];
        //}
        //__syncthreads();
    //}
    //if (blockDim.x >= 128) {
        //if (tx <  64) {
            //sharedmem_results[tx] += sharedmem_results[tx + 64];
        //}
        //__syncthreads();
    //}

    //if (tx < 32)
    //{
        //if (blockDim.x >=  64) {
            //sharedmem_results[tx] += sharedmem_results[tx + 32];
        //}
        //if (blockDim.x >=  32) {
            //sharedmem_results[tx] += sharedmem_results[tx + 16];
        //}
        //if (blockDim.x >=  16) {
            //sharedmem_results[tx] += sharedmem_results[tx +  8];
        //}
        //if (blockDim.x >=   8) {
            //sharedmem_results[tx] += sharedmem_results[tx +  4];
        //}
        //if (blockDim.x >=   4) {
            //sharedmem_results[tx] += sharedmem_results[tx +  2];
        //}
        //if (blockDim.x >=   2) {
            //sharedmem_results[tx] += sharedmem_results[tx +  1];
        //}
    //}

    //// write results of the calculation back to global memory, only thread 0 has the correct answer
    //// only one warp processes < 32, therefore no syncthreads required here
    //if (tx == 0)
        //result[bx*gridDim.x+by] = sharedmem_results[tx];

    //return;
//};


void CUDA_EonDeviceNC(float4 *residuePositions, float4 *residueMeta, int * residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, double* result, int blockSize, int datasetSize, int sm_size)
{
    int gridSize = datasetSize/blockSize;
    // the parallel sum of the results must have an array with size that is a power of 2
    int resultSize = gridSize;
    float *tmpSums = new float[resultSize*resultSize];

    //for a parallel sum each grid must have one cell in the array of results from all the threads
    float *d_result;
    cudaMalloc((void **)&d_result,sizeof(float)*resultSize*resultSize);
    cudaMemset(d_result,0,sizeof(float)*resultSize*resultSize);

    dim3 dimBlock(blockSize,1,1);   		// threads
    dim3 dimGrid(gridSize,gridSize,1);    	// blocks

    E_TiledKernelNC<<< dimGrid,dimBlock,sm_size,0 >>>(residuePositions, residueMeta, residueCount, moleculePositions, moleculeCount, LJPotentials, d_result);

    cudaMemcpy(tmpSums,d_result, sizeof(float)*resultSize*resultSize, cudaMemcpyDeviceToHost);
    result[0] = 0.0f;
    for (int i=0; i<resultSize*resultSize; i++)
    {
        result[0] += tmpSums[i];
    }

    delete [] tmpSums;

    // free the memory assigned for this iteration.
    cudaFree(d_result);
    getLastCudaError("Kernel execution failed EonDeviceNC");
    return;
};

///* full kernel with optional texture containers
__global__ void E_TiledKernelNC(float4 * residuePositions, float4 * residueMeta, int * residueCount, int * moleculePositions, int * moleculeCount, float* LJPotentialData, float* result)
{
    extern __shared__ float sharedmem_results[];	// array of temp results to sum, one per thread run.
#if POSITIONDATA_MEMORY == SHARED_MEM
    const int positionOfPosition = blockDim.x/4;
    extern __shared__ float4 X_tile_residuePositions[];
#endif
#if METADATA_MEMORY == SHARED_MEM
#if POSITIONDATA_MEMORY == SHARED_MEM
    const int positionOfMeta = __mul24(5,blockDim.x)/4;
#else
    const int positionOfMeta = blockDim.x/4;
#endif
    extern __shared__ float4 X_tile_residueMeta[];
#endif

#if LJ_LOOKUP_METHOD == SHARED_MEM_LJ
    __shared__ float sharedLJ[AA_COUNT*AA_COUNT];

    if(tx==0)
    {
        for (int i=0; i<AA_COUNT*AA_COUNT; i+=AA_COUNT)
            memcpy(&(sharedLJ[i]),&(LJPotentialData[i]),sizeof(float)*AA_COUNT);
    }
#endif

    sharedmem_results[tx] = 0.0f;//EPS;
#if  CULL_LOWER_TRIANGULAR_SUM  // only works if we throw away the entire block, otherwise results cannot be guaranteed correct due to the parallel sum
    if (bx < by)
        return;
#endif

#if POSITIONDATA_MEMORY == TEXTURE_MEM
    float4 yresiduep = tex1Dfetch(residuePositionTex, by*blockDim.x+tx);
#else
    float4 yresiduep = residuePositions[by*blockDim.x+tx];
#endif
#if METADATA_MEMORY == TEXTURE_MEM
    float4 yresiduem = tex1Dfetch(residueMetaTex, by*blockDim.x+tx);
#else
    float4 yresiduem = residueMeta[by*blockDim.x+tx];
#endif

#if POSITIONDATA_MEMORY == SHARED_MEM
    X_tile_residuePositions[positionOfPosition + tx] = residuePositions[bx*blockDim.x+tx];
#endif
#if METADATA_MEMORY == SHARED_MEM
    X_tile_residueMeta[positionOfMeta + tx] = residueMeta[bx*blockDim.x+tx];
#endif

    __syncthreads();  // very important so that the entire shared mem array is assigned before any threads use it

    if (yresiduep.w > PADDER_IDENTIFIER && yresiduem.w != CROWDER_IDENTIFIER ) // the residue in question is not a padding element or a crowder
    {
        float lj_subtotal(0.0f);
        float dh_subtotal(0.0f);

        //#pragma unroll 4 - cant unroll dynamic variable
        for (int i = 0; i < blockDim.x; i++)
        {
#if POSITIONDATA_MEMORY == TEXTURE_MEM
            float4 pos = tex1Dfetch(residuePositionTex,(bx*blockDim.x+i));
#elif POSITIONDATA_MEMORY == SHARED_MEM
            float4 pos = X_tile_residuePositions[positionOfPosition + i];
#endif

#if METADATA_MEMORY == TEXTURE_MEM
            float4 meta = tex1Dfetch(residueMetaTex,(bx*blockDim.x+i));
#elif METADATA_MEMORY == SHARED_MEM
            float4 meta = X_tile_residueMeta[positionOfMeta + i];
#endif
            
#if FLEXIBLE_LINKS && !ASSUME_POLYMER_FOLDING_TEST
            // IN THIS KERNEL WE ONLY CALCULATE THE NONBONDED POTENTIAL
            // but we have to use the same data as the other kernel
            
            // we unpack the floats
            
            //float xmol(0.0f);
            //float ymol(0.0f);
            
            //float xchain = modff(meta.w, &xmol); // NOT USED HERE; ONLY USED IN NORMAL KERNEL
            //float ychain = modff(yresiduem.w, &ymol); // NOT USED HERE; ONLY USED IN NORMAL KERNEL
            
            // Is there a more efficient function we could use to get just the integer part?

            float xmol = truncf(meta.w);
            float ymol = truncf(yresiduem.w);
                        
            if (pos.w == PADDER_IDENTIFIER || meta.w == CROWDER_IDENTIFIER || xmol == ymol)
#elif ASSUME_POLYMER_FOLDING_TEST
            // More efficient handling of this special case.  Assume we are folding a single polymer.
            // There's only one chain; we calculate all pairs except close neighbours on the backbone.
            // This is actually totally irrelevant and we shouldn't be calling this kernel!
            if (fabs(meta.w - yresiduem.w) < 4)
#else
             // same molecule or padder or crowder
            if (yresiduep.w == pos.w || pos.w == PADDER_IDENTIFIER || meta.w == CROWDER_IDENTIFIER)
#endif  
            {
                // DO NOTHING
            }
            else
            {
                //m.x = aminoacid index
                //m.y = charge
                //m.z = vdw radius
                //m.w = crowder if == CROWDER_IDENTIFIER
                float r(length(yresiduep,pos) + EPS);  // add eps so that r is never 0, can happen in a collision

                float LJ(0.0f);
                float DH(0.0f);
#if !LJ_OFF
                int ijX(rint( AA_COUNT*yresiduem.x + meta.x));

                //do the texture fetch first
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
                float Eij(LJ_lambda*(tex1Dfetch(LJTexture,ijX) - e0));
#elif LJ_LOOKUP_METHOD == SHARED_MEM
                float Eij(LJ_lambda*(sharedLJ[ijX] - e0));
#elif LJ_LOOKUP_METHOD == CONST_MEM
                float Eij(LJ_lambda*(const_LJPotentialData[ijX] - e0));
#else  // __global__ or __constant__
                float Eij(LJ_lambda*(LJPotentialData[ijX] - e0));
#endif // LJ_LOOKUP_METHOD
#endif // !LJ_OFF
                DH = dhPotential(yresiduem.y,meta.y,r);
                dh_subtotal += DH;
#if !LJ_OFF
                // sigmaij is the average atomic radius determined by the van der waals radius in kim2008
                float sigmaij((yresiduem.z + meta.z) * 0.5f);

                float LJtmp(__powf(sigmaij/r,6.0f));

                LJ = -4.0f*Eij*LJtmp*(LJtmp-1.0f);
                if (Eij>0.0f && r < (sigmaij*r0_constant))
                {
                    LJ = -LJ + 2.0f*Eij;
                }

                lj_subtotal += LJ;
#endif // !LJ_OFF
            } // if pair not excluded
        } // for i = 0..Bdx
        sharedmem_results[tx] = (lj_subtotal * RT_to_kcalmol) + (dh_subtotal * DH_CONVERSION_FACTOR);

#if CULL_LOWER_TRIANGULAR_SUM
        if (bx==by)
            sharedmem_results[tx] = sharedmem_results[tx]*0.5f;
#else
        sharedmem_results[tx] *= 0.5f
#endif //  CULL_LOWER_TRIANGULAR_SUM
    } // if !padding

    __syncthreads();  // all threads must have computed values so that they can be read concurrently from shared memory

    // do reduction in shared mem
    if (blockDim.x >= 512) {
        if (tx < 256) {
            sharedmem_results[tx] += sharedmem_results[tx + 256];
        }
        __syncthreads();
    }
    if (blockDim.x >= 256) {
        if (tx < 128) {
            sharedmem_results[tx] += sharedmem_results[tx + 128];
        }
        __syncthreads();
    }
    if (blockDim.x >= 128) {
        if (tx <  64) {
            sharedmem_results[tx] += sharedmem_results[tx + 64];
        }
        __syncthreads();
    }

    if (tx < 32)
    {
        if (blockDim.x >=  64) {
            sharedmem_results[tx] += sharedmem_results[tx + 32];
        }
        if (blockDim.x >=  32) {
            sharedmem_results[tx] += sharedmem_results[tx + 16];
        }
        if (blockDim.x >=  16) {
            sharedmem_results[tx] += sharedmem_results[tx +  8];
        }
        if (blockDim.x >=   8) {
            sharedmem_results[tx] += sharedmem_results[tx +  4];
        }
        if (blockDim.x >=   4) {
            sharedmem_results[tx] += sharedmem_results[tx +  2];
        }
        if (blockDim.x >=   2) {
            sharedmem_results[tx] += sharedmem_results[tx +  1];
        }
    }

    // write results of the calculation back to global memory, only thread 0 has the correct answer
    // only one warp processes < 32, therefore no syncthreads required here
    if (tx == 0)
        result[bx*gridDim.x+by] = sharedmem_results[tx];

    return;
};












// performs a parallel sum on the data in values on the gpu with dimentions blockDim.x*blockDim.y
__global__ void parallelSum_kernel(float * values)
{
#if REDUCTION_BLOCK_DIM > 256
    if (REDUCTION_BLOCK_DIM == 512) {
        if (tx < 256) {
            values[tx] += values[tx + 256];
        }
        __syncthreads();
    }
#endif
#if REDUCTION_BLOCK_DIM > 128
    if (REDUCTION_BLOCK_DIM >= 256) {
        if (tx < 128) {
            values[tx] += values[tx + 128];
        }
        __syncthreads();
    }
#endif
#if REDUCTION_BLOCK_DIM > 64
    if (REDUCTION_BLOCK_DIM >= 128) {
        if (tx <  64) {
            values[tx] += values[tx +  64];
        }
        __syncthreads();
    }
#endif
    if (tx < 32)
    {
        if (REDUCTION_BLOCK_DIM >=  64) {
            values[tx] += values[tx + 32];
        }
        if (REDUCTION_BLOCK_DIM >=  32) {
            values[tx] += values[tx + 16];
        }
        if (REDUCTION_BLOCK_DIM >=  16) {
            values[tx] += values[tx +  8];
        }
        if (REDUCTION_BLOCK_DIM >=   8) {
            values[tx] += values[tx +  4];
        }
        if (REDUCTION_BLOCK_DIM >=   4) {
            values[tx] += values[tx +  2];
        }
        if (REDUCTION_BLOCK_DIM >=   2) {
            values[tx] += values[tx +  1];
        }
    }
};

//__global__ void rotateMolecule_kernel (float4 *residuePositions, int *startPosition, int *length, float4* rotationVector, float4 *center)
//{
    //// cache frequently used variables
    //__shared__ float4 centerV;
    //__shared__ float4 q;

    //centerV.x = center->x;
    //centerV.y = center->y;
    //centerV.z = center->z;

    //// create rotation quaternion
    //q.w = cos(rotationVector->w*0.5f);
    //q.x = rotationVector->x*sin(rotationVector->w*0.5f);
    //q.y = rotationVector->y*sin(rotationVector->w*0.5f);
    //q.z = rotationVector->z*sin(rotationVector->w*0.5f);

    //__syncthreads();

    //// made it a unit quaternion
    //float rq_len = rsqrt(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);

    //__syncthreads();

    //q.w = q.w*rq_len;
    //q.x = q.x*rq_len;
    //q.y = q.y*rq_len;
    //q.z = q.z*rq_len;

    //__syncthreads();

    ////rotate the molecule using the quaternion made earlier

//#define vx (residuePositions[bx*blockDim.x+tx+startPosition[0]].x - centerV.x)
//#define vy (residuePositions[bx*blockDim.x+tx+startPosition[0]].y - centerV.y)
//#define vz (residuePositions[bx*blockDim.x+tx+startPosition[0]].z - centerV.z)

    //if (bx*blockDim.x+tx < length[0])  // not addressing out of bounds
    //{
        //residuePositions[bx*blockDim.x+startPosition[0]+tx].x = ( vx*(q.w*q.w+q.x*q.x-q.y*q.y-q.z*q.z) + vy*(q.x*q.y-q.w*q.z)*2.0f  		+ vz*(q.w*q.y+q.x*q.z)*2.0f )        	 + centerV.x;
        //residuePositions[bx*blockDim.x+startPosition[0]+tx].y = ( vx*(q.w*q.z+q.x*q.y)*2.0f   		 + vy*(q.w*q.w-q.x*q.x+q.y*q.y-q.z*q.z) + vz*(q.y*q.z-q.w*q.x)*2.0f ) 		     + centerV.y;
        //residuePositions[bx*blockDim.x+startPosition[0]+tx].z = ( vx*(q.x*q.z-q.w*q.y)*2.0f   		 + vy*(q.y*q.z+q.w*q.x)*2.0f   		    + vz*(q.w*q.w-q.x*q.x-q.y*q.y+q.z*q.z) ) + centerV.z;
    //}

//};

//__global__ void translateMolecule_kernel (float4 *residuePositions, int *startPosition, int *length, float4* translationVector, float4 *center)
//{
    //if (bx*blockDim.x+tx < length[0])  // not addressing out of bounds
    //{
        //residuePositions[bx*blockDim.x+startPosition[0]+tx].x = (residuePositions[bx*blockDim.x+startPosition[0]+tx].x + translationVector->x);
        //residuePositions[bx*blockDim.x+startPosition[0]+tx].y = (residuePositions[bx*blockDim.x+startPosition[0]+tx].y + translationVector->y);
        //residuePositions[bx*blockDim.x+startPosition[0]+tx].z = (residuePositions[bx*blockDim.x+startPosition[0]+tx].z + translationVector->z);
    //}

    //float3 tmpCenter;
    //tmpCenter.x = center->x + translationVector->x;
    //tmpCenter.y = center->y + translationVector->y;
    //tmpCenter.z = center->z + translationVector->z;

    //__syncthreads();

    //center->x = tmpCenter.x;
    //center->y = tmpCenter.y;
    //center->z = tmpCenter.z;
//};





#endif //USING CUDA
#endif
