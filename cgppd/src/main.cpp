#include "Simulation.h"

int main(int argc, char **argv)
{
    LOG(ALWAYS, "CGPPD version %s\n", HGVERSION);
    LOG(ALWAYS, "Compiled with:\n");

#ifdef GLVIS
    LOG(ALWAYS, "\tOpenGL support\n");
#endif

#if USING_CUDA
    LOG(ALWAYS, "\tCUDA support\n");
#if CUDA_STREAMS
    LOG(ALWAYS, "\t\tAsynchronous GPU calls (CUDA capability 1.1+ required)\n");
#endif // CUDA_STREAMS
    LOG(ALWAYS, "\t\tTile Size: %d\n", TILE_DIM);
    LOG(ALWAYS, "\t\tLJ lookup memory type: ");
#if LJ_LOOKUP_METHOD == SHARED_MEM
    LOG(ALWAYS, "Shared\n");
#elif LJ_LOOKUP_METHOD == CONST_MEM
    LOG(ALWAYS, "Constant\n");
#elif LJ_LOOKUP_METHOD == GLOBAL_MEM
    LOG(ALWAYS, "Global\n");
#elif LJ_LOOKUP_METHOD == TEXTURE_MEM
    LOG(ALWAYS, "Texture\n");
#endif // LJ_LOOKUP_METHOD
#endif // USING_CUDA

    LOG(COMPENSATE_KERNEL_SUM, "\tKahan summation in kernels\n");
    LOG(FLEXIBLE_LINKS, "\tFlexible linkers\n");
    LOG(LJ_REPULSIVE, "\tLennard-Jones potentials always repulsive\n");
    LOG(LJ_OFF, "\tLennard-Jones potentials off\n");

    int pid = int( getpid() );

    Simulation simulation;
    simulation.init(argc, argv, pid);

#if GLVIS
    if (simulation.parameters.viewConditions)
    {
        init_glvis(&simulation.replica, &simulation.initialReplica, argc, argv, simulation.parameters.replicas, int(simulation.parameters.bound));
    }
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

    LOG(ALWAYS, "Finished.\n");
    return 0;
}
