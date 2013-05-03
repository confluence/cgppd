#include <iostream>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>

#include "definitions.h"

#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include "Simulation.h"

class TestSimulation : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestSimulation);
    CPPUNIT_TEST(testGetArgs);
    CPPUNIT_TEST(testLoadArgsFromFile);
//     CPPUNIT_TEST(testPrintArgs);
    CPPUNIT_TEST_SUITE_END();

private:


public:
    void setUp();
    void testGetArgs();
    void testLoadArgsFromFile();
//     void testPrintArgs();
    void tearDown();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestSimulation);

void TestSimulation::setUp()
{
}

void TestSimulation::testGetArgs()
{
#define LEN_TEST_ARGV 28
    argdata a;
    const char * c_argv[LEN_TEST_ARGV] = {"my_programname", "-f", "my_filename", "-c", "-v", "-q", "-t", "23", "-s", "23", "-g", "23", "-m", "23", "-e", "23", "-r", "23", "-o", "my_suffix", "-b", "23", "-n", "23", "-x", "23", "-d", "23"};
    char * argv[LEN_TEST_ARGV];
    for(int i = 0; i < LEN_TEST_ARGV; i++) {
        argv[i] = strdup(c_argv[i]);
    }
    int argc(LEN_TEST_ARGV);

    getArgs(&a, argc, argv, 5);

    CPPUNIT_ASSERT_EQUAL(true, a.viewConditions);
    CPPUNIT_ASSERT_EQUAL(true, a.skipsimulation);
#if USING_CUDA
    CPPUNIT_ASSERT_EQUAL(23, a.cuda_blockSize);
    CPPUNIT_ASSERT_EQUAL(false, a.auto_blockdim);
#endif
    CPPUNIT_ASSERT_EQUAL(23, a.threads);
    CPPUNIT_ASSERT_EQUAL(23, a.streams);
    CPPUNIT_ASSERT_EQUAL(23, a.gpus);
    CPPUNIT_ASSERT_EQUAL(23, a.MCsteps);
    CPPUNIT_ASSERT_EQUAL(23, a.REsteps);
    CPPUNIT_ASSERT_EQUAL(23, a.replicas);
    CPPUNIT_ASSERT_EQUAL(5000, a.sampleFrequency);
    CPPUNIT_ASSERT_EQUAL(5000, a.sampleStartsAfter);
    CPPUNIT_ASSERT_EQUAL(true, a.inputFile);
    CPPUNIT_ASSERT_EQUAL(0, a.nonCrowders);
    CPPUNIT_ASSERT_EQUAL(5, a.pid);
    CPPUNIT_ASSERT_EQUAL(false, a.resume);
//     CPPUNIT_ASSERT_EQUAL(0, a.currentStep); // TODO: fix this

    CPPUNIT_ASSERT_DOUBLES_EQUAL(a.bound, 23.0, 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(a.temperatureMin, 23.0, 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(a.temperatureMax, 23.0, 0);

//     CPPUNIT_ASSERT(strcmp(a.prependageString, "") == 0); // TODO: fix this
    CPPUNIT_ASSERT(strcmp(a.file, "my_filename") == 0);
    CPPUNIT_ASSERT(strcmp(a.logfile, "my_suffix") == 0);
}

void TestSimulation::testLoadArgsFromFile()
{
    argdata a;
    strcpy(a.file, "tests/argtestfile");
    a.inputFile = true;

    loadArgsFromFile(&a);

    cout << a.MCsteps;
}

void TestSimulation::tearDown()
{
}