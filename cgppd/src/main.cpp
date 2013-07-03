#include "Simulation.h"

int main(int argc, char **argv)
{
    cout.precision(8);

    cout << "Version: " << HGVERSION << endl;
    cout << "Compiled with:" << endl;
#ifdef GLVIS
    cout << "  OpenGL support" << endl;
#endif
#if USING_CUDA
    cout << "  CUDA support" << endl;
#endif
#if CUDA_STREAMS
    cout << "  Asynchronous GPU calls (CUDA capability 1.1+ required)" << endl;
#endif
#if COMPENSATE_KERNEL_SUM
    cout << "  Kahan summation in kernels" << endl;
#endif

#if USING_CUDA
#if LOGLEVEL >= 0
    cout << "CUDA parameters and options for this run:\n" ;
    cout << "-----------------------------------------\n" ;
    cout << "Tile Size " << TILE_DIM << endl;
    cout << "LJ lookup memory type: ";
#if LJ_LOOKUP_METHOD == SHARED_MEM
    cout << "Shared" << endl;
#elif LJ_LOOKUP_METHOD == CONST_MEM
    cout << "Constant" << endl;
#elif LJ_LOOKUP_METHOD == GLOBAL_MEM
    cout << "Global" << endl;
#elif LJ_LOOKUP_METHOD == TEXTURE_MEM
    cout << "Texture" << endl;
#endif // LJ_LOOKUP_METHOD
    cout << "-----------------------------------------\n" ;
#endif // LOGLEVEL
#endif // USING_CUDA

    int pid = int( getpid() );

    Simulation simulation;
    simulation.init(argc, argv, pid);

// TODO: this is pretty gross; hide it in a function.
#if GLVIS
    if (simulation.parameters.viewConditions)
    {
        replica = &simulation.replica;
        glutInit(&argc, argv);
        camera.setPosition(Vector3f(-15,15,15),Vector3f(1,-1,-1),Vector3f(0,1,0));
        char windowName[64] = {"Coarse-Grained Protein-Protein Docker (cgppd)"};
        GlutInit(WIDTH,HEIGHT,windowName);
        gl_replicaCount = simulation.parameters.replicas;
        gl_boundingValue = int(simulation.parameters.bound);
        GLreplica = &simulation.initialReplica;
    }
#endif

    if (!simulation.parameters.skipsimulation)
    {
        simulation.run();
    }

#if GLVIS
    if (simulation.parameters.viewConditions)
    {
        // TODO: this is overwritten inside run. Why do we set it back here?
        GLreplica = &simulation.initialReplica;
        cout << "Entering free gl viewing mode." << endl;
        glutMainLoop();
    }
#endif

    cout << "Finished." << endl;
    return 0;
}
