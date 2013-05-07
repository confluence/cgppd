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
    argdata parameters;

    int pid = int( getpid() );
    getArgs(&parameters, argc, argv, pid);

    Simulation simulation;

    loadArgsFromFile(&parameters);
    checkParameterSanity(&parameters);
    writeFileIndex(&parameters);
    getArgs(&parameters, argc, argv, pid); // second pass to override any variables if doing performance tests

    // TODO: only if verbose output
    printArgs(&parameters);

#if USING_CUDA
    if (parameters.runcheck)
    {
        simulation.run_check();
    }
    else
    {
#endif
        // TODO: make all the parameter functions into methods; fold everything into init:
        simulation.init(parameters);

#if GLVIS
        if (parameters.viewConditions)
        {
            replica = &simulation.replica;
            glutInit(&argc, argv);
            camera.setPosition(Vector3f(-15,15,15),Vector3f(1,-1,-1),Vector3f(0,1,0));
            char windowName[64] = {"Coarse-Grained Protein-Protein Docker (cgppd)"};
            GlutInit(WIDTH,HEIGHT,windowName);
            gl_replicaCount = parameters.replicas;
            gl_boundingValue = int(parameters.bound);
            GLreplica = &simulation.initialReplica;
        }
#endif

        if (!parameters.skipsimulation)
        {
            simulation.run();
        }

#if GLVIS
        if (parameters.viewConditions)
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
