#include <iostream>
#include "definitions.h"

#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include "Simulation.h"
#include <testCommon.h>

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
#define LEN_TEST_ARGV 27
    Simulation s;
    const char * c_argv[LEN_TEST_ARGV] = {"my_programname", "-f", "my_filename", "-p", "-q", "-t", "23", "-s", "23", "-g", "23", "-m", "23", "-e", "23", "-r", "23", "-o", "my_suffix", "-b", "23", "-n", "23", "-x", "23", "-d", "23"};
    char * argv[LEN_TEST_ARGV];
    for(int i = 0; i < LEN_TEST_ARGV; i++) {
        argv[i] = strdup(c_argv[i]);
    }
    int argc(LEN_TEST_ARGV);

    s.getFileArg(argc, argv);
    s.getArgs(argc, argv);

    CPPUNIT_ASSERT_EQUAL(true, s.parameters.viewConditions);
    CPPUNIT_ASSERT_EQUAL(true, s.parameters.skipsimulation);
#if USING_CUDA
    CPPUNIT_ASSERT_EQUAL(23, s.parameters.cuda_blockSize);
    CPPUNIT_ASSERT_EQUAL(false, s.parameters.auto_blockdim);
#endif
    CPPUNIT_ASSERT_EQUAL(23, s.parameters.threads);
    CPPUNIT_ASSERT_EQUAL(23, s.parameters.streams);
    CPPUNIT_ASSERT_EQUAL(23, s.parameters.gpus);
    CPPUNIT_ASSERT_EQUAL(23, s.parameters.MCsteps);
    CPPUNIT_ASSERT_EQUAL(23, s.parameters.REsteps);
    CPPUNIT_ASSERT_EQUAL(23, s.parameters.replicas);
    CPPUNIT_ASSERT_EQUAL(5000, s.parameters.sampleFrequency);
    CPPUNIT_ASSERT_EQUAL(5000, s.parameters.sampleStartsAfter);
    CPPUNIT_ASSERT_EQUAL(true, s.parameters.inputFile);
    CPPUNIT_ASSERT_EQUAL(0, s.parameters.nonCrowders);
    CPPUNIT_ASSERT_EQUAL(0, s.parameters.pid);
    CPPUNIT_ASSERT_EQUAL(false, s.parameters.resume);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(23.0, s.parameters.bound, 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(23.0, s.parameters.temperatureMin, 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(23.0, s.parameters.temperatureMax, 0);

    CPPUNIT_ASSERT(strcmp(s.parameters.prefix, "") == 0);
    CPPUNIT_ASSERT(strcmp(s.parameters.file, "my_filename") == 0);
    CPPUNIT_ASSERT(strcmp(s.parameters.logfile, "my_suffix") == 0);

    CPPUNIT_ASSERT_EQUAL(0, s.parameters.verbosity);
}

void TestSimulation::testLoadArgsFromFile()
{
    Simulation s;
    strcpy(s.parameters.file, "tests/argtestfile");
    s.parameters.inputFile = true;

    s.loadArgsFromFile();

    CPPUNIT_ASSERT_EQUAL(false, s.parameters.viewConditions);
    CPPUNIT_ASSERT_EQUAL(false, s.parameters.skipsimulation);
#if USING_CUDA
    CPPUNIT_ASSERT_EQUAL(TILE_DIM, s.parameters.cuda_blockSize);
    CPPUNIT_ASSERT_EQUAL(true, s.parameters.auto_blockdim);
#endif
    CPPUNIT_ASSERT_EQUAL(1, s.parameters.threads);
    CPPUNIT_ASSERT_EQUAL(1, s.parameters.streams);
    CPPUNIT_ASSERT_EQUAL(1, s.parameters.gpus);
    CPPUNIT_ASSERT_EQUAL(10000, s.parameters.MCsteps);
    CPPUNIT_ASSERT_EQUAL(10000, s.parameters.REsteps);
    CPPUNIT_ASSERT_EQUAL(1, s.parameters.replicas);
    CPPUNIT_ASSERT_EQUAL(1000, s.parameters.sampleFrequency);
    CPPUNIT_ASSERT_EQUAL(1000, s.parameters.sampleStartsAfter);
    CPPUNIT_ASSERT_EQUAL(true, s.parameters.inputFile);
    CPPUNIT_ASSERT_EQUAL(2, s.parameters.nonCrowders);
    CPPUNIT_ASSERT_EQUAL(0, s.parameters.pid);
    CPPUNIT_ASSERT_EQUAL(false, s.parameters.resume);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(100.0, s.parameters.bound, 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(288.0, s.parameters.temperatureMin, 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(300.0, s.parameters.temperatureMax, 0);

    CPPUNIT_ASSERT(strcmp(s.parameters.prefix, "") == 0);
    CPPUNIT_ASSERT(strcmp(s.parameters.file, "tests/argtestfile") == 0);
    CPPUNIT_ASSERT(strcmp(s.parameters.logfile, "") == 0);

    CPPUNIT_ASSERT_EQUAL(7, (int)s.parameters.mdata.size());
    CPPUNIT_ASSERT_EQUAL(2, (int)s.parameters.mdata_map.size());
    
    moldata uim_data = s.parameters.mdata[s.parameters.mdata_map["UIM1"]];
    moldata ubq_data = s.parameters.mdata[s.parameters.mdata_map["ubiquitin"]];
    
    CPPUNIT_ASSERT_EQUAL(true, uim_data.all_flexible);
    CPPUNIT_ASSERT_EQUAL(0, (int)uim_data.segments.size());
    
    CPPUNIT_ASSERT_EQUAL(false, ubq_data.all_flexible);
    CPPUNIT_ASSERT_EQUAL(1, (int)ubq_data.segments.size());
    
    CPPUNIT_ASSERT_EQUAL(72, ubq_data.segments[0].residue_indices[0]);
    CPPUNIT_ASSERT_EQUAL(73, ubq_data.segments[0].residue_indices[1]);
    CPPUNIT_ASSERT_EQUAL(74, ubq_data.segments[0].residue_indices[2]);
    CPPUNIT_ASSERT_EQUAL(75, ubq_data.segments[0].residue_indices[3]);
}

void TestSimulation::tearDown()
{
}
