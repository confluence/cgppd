#define CATCH_CONFIG_RUNNER
#include <testCommon.h>

int main( int argc, char **argv)
{
    LOGOG_INITIALIZE();

    int result;
    
    {
        logog::Cout log_output;
        LOGOG_SET_LEVEL(0);

        result = Catch::Session().run( argc, argv );
    } // We kill the log target(s) here

    LOGOG_SHUTDOWN();
    
    return result;
}
