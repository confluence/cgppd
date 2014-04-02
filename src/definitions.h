#ifndef DEFINITIONS_H_
    #define DEFINITIONS_H_

    // Constant Data Input Files

    #define AMINOACIDDATASOURCE 	"data/AminoAcids"
    #define TORSIONALPAIRDATA 		"data/torsional_pair_potentials"

    // Simulation Values
    #define EPS 1e-38f
    #define BOUNDING_RADIUS		100.0f   // angstrom. size of sphere for spherical boundary conditions
    #define BOUNDING_VALUE		150.0f   // size of box edge for periodic conditions
    #define BOUNDING_SPHERE  	0
    #define PERIODIC_BOUNDARY 	1
    #define BOUNDING_METHOD	  	PERIODIC_BOUNDARY
    #define MEMBRANE_PRESENT 0 // use membrane calculations, molecule 0 is the embedded molecule in this case.

    #define REPULSIVE_CROWDING 	0

    #define REPLICA_COUNT 						20
    #define MC_STEPS 							5000
    #define HIGHTEMP 							500.0f
    #define LOWTEMP 							250.0f
    #define REMC_STEPS							10
    #define STEPS_BEFORE_SAMPLE 				MC_STEPS
    #define SAMPLE_FREQ							MC_STEPS
    #define INITIAL_TRANSLATIONAL_STEP			0.5f // angstrom
    #define MIN_TRANSLATION						0.5f
    #define MAX_TRANSLATION						2.0f
    #define INITIAL_ROTATIONAL_STEP				0.2f // rads
    #define MIN_ROTATION						0.2f
    #define MAX_ROTATION						0.4f
    #define BOUND_ENERGY_VALUE					-1.1844f	// 2k_bT in kcal/mol at 294K TODO: this is wrong! Should use 298K
    #define translate_rotate_bernoulli_bias 	0.5f // translations / mutations (rotations+translations) ratio
    #define GEOMETRIC_CUTOFF					8.0f
    #define bool_recordAllSamples				true //false

    #define const_repulsive_cutoff	15.0f //repulsive crowder cutoff, angstrom
    #define USE_POTENTIAL_CUTOFF        0 //use the following cutoff, 0:=no, 1:=yes (GPU only)
    #define interaction_cutoff      45.0f //interaction cutoff, angstrom

    #define PADDER_IDENTIFIER	-2.0f	//flag for identifying a padding residue
    #define CROWDER_IDENTIFIER	-1.0f 	//flag for identifying a crowder residue
    #define RIGID_IDENTIFIER    -3.0f   //flag for identifying a residue in a rigid molecule
    #define PADDER_AATYPE 		20.0f

    // code things
//     #define VERBOSE_THREAD_LOGGING 1

    #define ALWAYS   1

    #define ERROR   1
    #define WARN   1
    #define INFO    1
    #define DEBUG   1

    #define LOG(level, ...) \
            do { if (level) printf(__VA_ARGS__); } while (0)

//     #define DEBUGLOG(fmt, ...) \
//             do { printf("%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)

    #ifdef EnableOPENGL			//set here to override the makefile
        #define GLVIS 1
    #else
        #define GLVIS 0
    #endif

    #define GL_AXES 			1
    #define THREADING 			1
    #define THREAD_COUNT		1
    #define STREAM_COUNT		1
    //#define CHECKPOINTING		1
    #define CHECKPOINTFREQUENCY 10000

    #ifdef EnableCUDA			//set here to override the makefile
        #define USING_CUDA 1
    #else
        #define USING_CUDA 0
    #endif

    #define COMPENSATE_KERNEL_SUM 0

    #ifdef EnableFlexibleLinkers
        #define FLEXIBLE_LINKS  1
    #else
        #define FLEXIBLE_LINKS  0
    #endif

    #define MC_TRANSLATE 0
    #define MC_ROTATE 1

    #if FLEXIBLE_LINKS
        #define MC_FLEX 2
        #define MC_LOCAL 3

        #define MC_LOCAL_TRANSLATE 0
        #define MC_LOCAL_CRANKSHAFT 1

        #define WEIGHT_MC_TRANSLATE 1.0f
        #define WEIGHT_MC_ROTATE 1.0f
        #define WEIGHT_MC_FLEX 1.0f
        #define WEIGHT_MC_LOCAL 1.0f

        #define NUM_LOCAL_MOVES 1
        #define LOCAL_TRANSLATE_BIAS 0.5

        #define LOCAL_TRANSLATE_STEP_SCALING_FACTOR 0.1f
    #endif

    #if USING_CUDA
        #define GPU_COUNT          		1
        #define USING_TILE_KERNEL  		1
        #define PARALLEL_SUM			0

        #define SHARED_MEM				0
        #define CONST_MEM				1
        #define GLOBAL_MEM				2
        #define TEXTURE_MEM				3

        #define LJ_LOOKUP_METHOD        TEXTURE_MEM // type of memory to use for LJ lookups
        #define METADATA_MEMORY			SHARED_MEM  // type of memory to use for residue metadata
        #define POSITIONDATA_MEMORY		SHARED_MEM  // type of memory to use for residue position data

        #define USE_LOOP_UNROLL	          0
        #define CULL_LOWER_TRIANGULAR_SUM 1		// only computer the upper triangular sum for the summations, saves 0.5*n(n-1) work
        #ifdef EnableStreams
            #define CUDA_STREAMS 1
        #else
            #define CUDA_STREAMS 0
        #endif


        #define CUDA_MC		               0	// use cuda to do the MC operations

        #define TILE_DIM 			64 // block size for looping kernel, must be a power of 2 <= 512
        #define REDUCTION_BLOCK_DIM		128   // block size for reductions

        #define BLOCK_SIZE 				22  // N*N <= 512
        #if BLOCK_SIZE <= 4
            #define SHARED_MEM_RESULT_SIZE	16 		// [nearest power of 2 > BLOCK_SIZE]^2
        #elif BLOCK_SIZE <= 5
            #define SHARED_MEM_RESULT_SIZE	32
        #elif BLOCK_SIZE <= 8
            #define SHARED_MEM_RESULT_SIZE	64
        #elif BLOCK_SIZE <= 11
            #define SHARED_MEM_RESULT_SIZE	128
        #elif BLOCK_SIZE <= 16
            #define SHARED_MEM_RESULT_SIZE	256
        #elif BLOCK_SIZE <= 22						//sqrt 512 == 22.something
            #define SHARED_MEM_RESULT_SIZE	512
        #endif
        #define INCLUDE_TIMERS      1   //include timers
    #else
        #define INCLUDE_TIMERS 0 // the timers use cutil -- no CUDA; no timers.
    #endif // USING_CUDA 1

    #if GLVIS
        void GlutDisplay();
    #endif

    #define PRINT_REFERENCE_CONFORMATIONS_FROM_TEST 0
    #define PRINT_STARTING_PDB 1

#endif /*DEFINITIONS_H_*/
