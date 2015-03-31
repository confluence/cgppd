#include "Simulation.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv)
{
    //FLAGS_log_dir = "logs";
    //FLAGS_alsologtostderr = 1;

    //google::InitGoogleLogging(argv[0]);

    VLOG(0) << "CGPPD version: " << HGVERSION;
    VLOG(0) << "Compiled with:";

#ifdef GLVIS
    VLOG(0) << "\tOpenGL support";
#endif

#if USING_CUDA
    VLOG(0) << "\tCUDA support";
#if CUDA_STREAMS
    VLOG(0) << "\t\tAsynchronous GPU calls (CUDA capability 1.1+ required)";
#endif // CUDA_STREAMS
    VLOG(0) << "\t\tTile size: " << TILE_DIM;
    string mem_type;
#if LJ_LOOKUP_METHOD == SHARED_MEM
    mem_type = "Shared";
#elif LJ_LOOKUP_METHOD == CONST_MEM
    mem_type = "Constant";
#elif LJ_LOOKUP_METHOD == GLOBAL_MEM
    mem_type = "Global";
#elif LJ_LOOKUP_METHOD == TEXTURE_MEM
    mem_type = "Texture";
#endif // LJ_LOOKUP_METHOD
    VLOG(0) << "\t\tLJ lookup memory type: " << mem_type;
#endif // USING_CUDA

    VLOG_IF(COMPENSATE_KERNEL_SUM, 0) << "\tKahan summation in kernels";
    VLOG_IF(FLEXIBLE_LINKS, 0) << "\tFlexible linkers";
    VLOG_IF(LJ_REPULSIVE, 0) << "\tLennard-Jones potentials always repulsive";
    VLOG_IF(LJ_OFF, 0) << "\tLennard-Jones potentials off";

    int pid = int( getpid() );

    Simulation simulation;
    simulation.init(argc, argv, pid);

#if GLVIS
    if (simulation.parameters.viewConditions)
    {
        init_glvis(&simulation.replica, &simulation.initialReplica, argc, argv, simulation.parameters.replicas, int(simulation.parameters.bound));
    }
#endif

#if PRINT_STARTING_PDB
    char start_pdb_filename[256];
    memset(start_pdb_filename, 0, 256);
    sprintf(start_pdb_filename, "%s_start.pdb", simulation.parameters.prefix);
    simulation.initialReplica.saveAsSinglePDB(start_pdb_filename, simulation.parameters.title);
#endif

    if (!simulation.parameters.skipsimulation)
    {
        simulation.run();
    }

#if GLVIS
    if (simulation.parameters.viewConditions)
    {
        enter_viewing_mode(&simulation.initialReplica);
    }
#endif

    VLOG(0) << "Finished.";
    return 0;
}
