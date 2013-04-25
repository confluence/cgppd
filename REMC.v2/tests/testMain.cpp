#include <iostream>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>

#include "definitions.h"

#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include "main.h"

// extern bool getArgs(argdata * d, int argc, char **argv);

class TestMain : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestMain);
    CPPUNIT_TEST(testGetArgs);
//     CPPUNIT_TEST(testLoadArgsFromFile);
    CPPUNIT_TEST_SUITE_END();

private:


public:
    void setUp();
    void testGetArgs();
//     void testLoadArgsFromFile();
    void tearDown();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestMain);

void TestMain::setUp()
{
}

void TestMain::testGetArgs()
{
    argdata a;
    const char * c_argv[29] = {"my_programname", "-f", "my_filename", "-c", "-h", "-v", "-0", "-t", "23", "-s", "23", "-g", "23", "-m", "23", "-e", "23", "-r", "23", "-o", "my_suffix", "-b", "23", "-t0", "23", "-t1", "23", "-bx", "23"};
    char * argv[29];
    for(int i = 0; i < 30; i++) {
        argv[i] = strdup(c_argv[i]);
    }
    int argc(30);

    getArgs(&a, argc, argv);

    cout << "viewConditions: " << a.viewConditions << endl;
    cout << "skipsimulation: " << a.skipsimulation << endl;
#if USING_CUDA
    cout << "cuda_blockSize: " << a.cuda_blockSize << endl;
    cout << "auto_blockdim: " << a.auto_blockdim << endl;
#endif
    cout << "threads: " << a.threads << endl;
    cout << "streams: " << a.streams << endl;
    cout << "gpus: " << a.gpus << endl;
    cout << "MCsteps: " << a.MCsteps << endl;
    cout << "REsteps: " << a.REsteps << endl;
    cout << "replicas: " << a.replicas << endl;
    cout << "bound: " << a.bound << endl;
    cout << "sampleFrequency: " << a.sampleFrequency << endl;
    cout << "sampleStartsAfter: " << a.sampleStartsAfter << endl;
    cout << "prependageString: " << a.prependageString << endl;
    cout << "file: " << a.file << endl;
    cout << "logfile: " << a.logfile << endl;
    cout << "inputFile: " << a.inputFile << endl;
    cout << "nonCrowders: " << a.nonCrowders << endl;
    cout << "temperatureMin: " << a.temperatureMin << endl;
    cout << "temperatureMax: " << a.temperatureMax << endl;
    cout << "pid: " << a.pid << endl;
    cout << "checkpointfilename: " << a.checkpointfilename << endl;
    cout << "resume: " << a.resume << endl;
    cout << "currentStep: " << a.currentStep << endl;
}

void TestMain::tearDown()
{
}