#include "Simulation.h"
#include <testCommon.h>

TEST_CASE("Simulation", "[simulation]") {

    Simulation s;

    SECTION ("Get arguments from command line") {

        char * argv[] = {"my_programname", "-f", "tests/argtestfile", "-p", "-q", "-t", "23", "-s", "23", "-g", "23", "-m", "23", "-e", "23", "-r", "23", "-b", "23", "-n", "23", "-x", "23", "-d", "23", NULL};
        int argc(25);

        s.getArgs(argc, argv, true);
        s.getArgs(argc, argv, false);

        REQUIRE(s.parameters.viewConditions);
        REQUIRE(s.parameters.skipsimulation);
    #if USING_CUDA
        REQUIRE(s.parameters.cuda_blockSize == 23);
        REQUIRE(s.parameters.auto_blockdim == false);
    #endif
        REQUIRE(s.parameters.threads == 23);
        REQUIRE(s.parameters.streams == 23);
        REQUIRE(s.parameters.gpus == 23);
        REQUIRE(s.parameters.MCsteps == 23);
        REQUIRE(s.parameters.REsteps == 23);
        REQUIRE(s.parameters.replicas == 23);
        REQUIRE(s.parameters.sampleFrequency == 5000);
        REQUIRE(s.parameters.sampleStartsAfter == 5000);
        REQUIRE(s.parameters.inputFile);
        REQUIRE(s.parameters.nonCrowders == 0);
        REQUIRE(s.parameters.pid == 0);
        REQUIRE(s.parameters.resume == false);

        REQUIRE(s.parameters.bound == Approx(23.0));
        REQUIRE(s.parameters.temperatureMin == Approx(23.0));
        REQUIRE(s.parameters.temperatureMax == Approx(23.0));

        REQUIRE(strcmp(s.parameters.prefix, "") == 0);
        REQUIRE(strcmp(s.parameters.file, "tests/argtestfile") == 0);
    }

    SECTION("Load configuration from file") {
        strcpy(s.parameters.file, "tests/argtestfile");
        s.parameters.inputFile = true;

        s.loadArgsFromFile();

        REQUIRE(s.parameters.viewConditions == false);
        REQUIRE(s.parameters.skipsimulation == false);
        #if USING_CUDA
        REQUIRE(s.parameters.cuda_blockSize == TILE_DIM);
        REQUIRE(s.parameters.auto_blockdim);
        #endif
        REQUIRE(s.parameters.threads == 1);
        REQUIRE(s.parameters.streams == 1);
        REQUIRE(s.parameters.gpus == 1);
        REQUIRE(s.parameters.MCsteps == 10000);
        REQUIRE(s.parameters.REsteps == 10000);
        REQUIRE(s.parameters.replicas == 1);
        REQUIRE(s.parameters.sampleFrequency == 1000);
        REQUIRE(s.parameters.sampleStartsAfter == 1000);
        REQUIRE(s.parameters.inputFile);
        REQUIRE(s.parameters.nonCrowders == 2);
        REQUIRE(s.parameters.pid == 0);
        REQUIRE(s.parameters.resume == false);

        REQUIRE(s.parameters.bound == Approx(100.0));
        REQUIRE(s.parameters.temperatureMin == Approx(288.0));
        REQUIRE(s.parameters.temperatureMax == Approx(300.0));

        REQUIRE(strcmp(s.parameters.prefix, "") == 0);
        REQUIRE(strcmp(s.parameters.file, "tests/argtestfile") == 0);

        REQUIRE(s.parameters.mdata.size() == 7);
        REQUIRE(s.parameters.mdata_map.size() == 2);

        moldata uim_data = s.parameters.mdata[s.parameters.mdata_map["UIM1"]];
        moldata ubq_data = s.parameters.mdata[s.parameters.mdata_map["ubiquitin"]];

        REQUIRE(uim_data.all_flexible);
        REQUIRE(uim_data.segments.size() == 0);

        REQUIRE(ubq_data.all_flexible == false);
        REQUIRE(ubq_data.segments.size() == 1);

        REQUIRE(ubq_data.segments[0].residue_indices[0] == 72);
        REQUIRE(ubq_data.segments[0].residue_indices[1] == 73);
        REQUIRE(ubq_data.segments[0].residue_indices[2] == 74);
        REQUIRE(ubq_data.segments[0].residue_indices[3] == 75);
    }

    SECTION("REMC sanity check") {

        char * argv[] = {"program_name_goes_here", "-f", "tests/remctestfile", NULL};

        // TODO: use a constant instead of a pid
        s.init(3, argv, int(getpid()));
        
        // seed the simulation range with a constant
        // seed each replica range with a constant
        
        gsl_rng_set(s.REMCRng, 0);
        
        for (int i = 0; i < s.parameters.replicas; i++) {
            gsl_rng_set(s.replica[i].rng, i);
        }

        // run the simulation
        
        s.run();

        // TODO: compare output directory to reference directory with diff -- exclude log file
        // TODO: delete test output directory

    }
}
