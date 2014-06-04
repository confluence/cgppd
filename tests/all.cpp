#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <string.h>
#include <testCommon.h>

using namespace std;

int main( int argc, char **argv)
{
    FLAGS_log_dir = "logs";
    FLAGS_alsologtostderr = 1;
    google::InitGoogleLogging(argv[0]);

    CppUnit::TextUi::TestRunner runner;
    CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
    runner.addTest( registry.makeTest() );
    
    string testname = "";
    
    if (argc > 1) {
        testname = string(argv[1]);
    }
    
    bool wasSuccessful = runner.run(testname, false);
    return wasSuccessful;
}
