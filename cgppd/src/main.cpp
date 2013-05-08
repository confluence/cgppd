#include "Simulation.h"

int main(int argc, char **argv)
{
    cout.precision(8);
    int sysreturn;
    // TODO: move this inside Simulation init?
    sysreturn = system("mkdir -p checkpoints");
    sysreturn = system("mkdir -p output");
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

    // get global options for the simulation

    int pid = int( getpid() );


    Simulation simulation;

    // TODO: make all the parameter functions into methods; fold everything into init:
    simulation.init(argc, argv, pid);

    // TODO: only if verbose output
    simulation.printArgs();

#if USING_CUDA
    if (simulation.parameters.runcheck)
    {
        simulation.run_check();
    }
    else
    {
#endif

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

#if USING_CUDA
    }
#endif
    cout << "Finished." << endl;
    return 0;
}
