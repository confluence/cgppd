#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <string.h>
#include <testCommon.h>

INITIALIZE_EASYLOGGINGPP

using namespace std;

int main( int argc, char **argv)
{
    START_EASYLOGGINGPP(argc, argv);
    
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
