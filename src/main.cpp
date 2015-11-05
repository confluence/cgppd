#include "Simulation.h"

//INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv)
{
    LOGOG_INITIALIZE();
    logog::Cout log_output;
    
    INFO("Hello, world!");

    //START_EASYLOGGINGPP(argc, argv);
    //el::Configurations conf("logs/log.conf");
    //el::Loggers::reconfigureLogger("default", conf);

    //VLOG(0) << "------------------------------------------------------------------------";
    //VLOG(0) << "LOGGING TO " << conf.get(el::Level::Info, el::ConfigurationType::Filename)->value();
    //VLOG(0) << "------------------------------------------------------------------------";

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

    //VLOG(0) << "Finished.";
    LOGOG_SHUTDOWN();
    return 0;
}
