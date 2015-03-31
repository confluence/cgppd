#define CATCH_CONFIG_RUNNER
#include <testCommon.h>

INITIALIZE_EASYLOGGINGPP

int main( int argc, char **argv)
{
    START_EASYLOGGINGPP(argc, argv);
    el::Configurations conf("logs/testlog.conf");
    el::Loggers::reconfigureLogger("default", conf);

    VLOG(0) << "------------------------------------------------------------------------";
    VLOG(0) << "LOGGING TO " << conf.get(el::Level::Info, el::ConfigurationType::Filename)->value();
    VLOG(0) << "------------------------------------------------------------------------";

    int result = Catch::Session().run( argc, argv );
    return result;
}
