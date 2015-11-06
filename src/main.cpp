#include "Simulation.h"

int main(int argc, char **argv)
{
    LOGOG_INITIALIZE();
    
    {
        logog::Cout log_output;
        
        try {
            Simulation simulation;
            simulation.init(argc, argv, int(getpid()));

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

            LOGOG_INFO("Finished.");
        } catch (const char* error_msg) {
            if (strcmp(error_msg, "help") != 0) {
                LOGOG_ERROR("Aborting: %s", error_msg);
            }
        }
    } // We kill the log target(s) here

    LOGOG_SHUTDOWN();
    return 0;
}
