#include <Simulation.h>
#include <testCommon.h>

#if !LJ_REPULSIVE && !LJ_OFF  // ignore this whole test unless we're using the normal LJ constants

TEST_CASE("Reference conformations" "[reference]") {
    AminoAcids aminoAcidData;
    Replica replicas[10];
    float bounding_value = 1000.0f;
    Potential cpu[10];
    Potential cpu_nc[10];
#if USING_CUDA
#if CUDA_STREAMS
    cudaStream_t streams[10];
#endif
    float * ljp_t;
    double gpu[10];
    double gpu_nc[10];
#endif

    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);

#if USING_CUDA
#if CUDA_STREAMS
    setup_CUDA(0, bounding_value, ljp_t, &aminoAcidData, streams, 10);
#else
    setup_CUDA(0, bounding_value, ljp_t, &aminoAcidData);
#endif // CUDA_STREAMS
#endif // USING_CUDA
    
    for (int i = 0; i < 10; i++)
    {
        argdata parameters;
        parameters.bound = bounding_value;
#if USING_CUDA
        parameters.auto_blockdim = false;
#endif

        moldata m1;
        m1.translate = true;
        sprintf(m1.pdbfilename, "data/conf%d/1a.pdb", i+1);

        moldata m2;
        m2.translate = true;
        sprintf(m2.pdbfilename, "data/conf%d/1b.pdb", i+1);

        parameters.mdata.push_back(m1);
        parameters.mdata.push_back(m2);

        replicas[i].init_first_replica(parameters, aminoAcidData, 2);
        replicas[i].label = i+1;

#if INCLUDE_TIMERS
        replicas[i].initTimers();
#endif

#if USING_CUDA
#if CUDA_STREAMS
        replicas[i].setup_CUDA(ljp_t, streams, i);
#else
        replicas[i].setup_CUDA(ljp_t);
#endif
#endif

        cpu[i] = replicas[i].E();
        cpu_nc[i] = replicas[i].E(&replicas[i].molecules[0], &replicas[i].molecules[1]);
#if USING_CUDA
        gpu[i] = replicas[i].EonDevice();
        gpu_nc[i] = replicas[i].EonDeviceNC();
#endif
    }

    SECTION("Test that exact results are unchanged") {
        // This test checks that nothing changes unexpectedly in our implementation of the potential calculation.
        // It also checks that the GPU, CPU and NC functions produce the same totals.
        // This should be tested more thoroughly (including the bonded components) in the Molecule and Replica unit tests.

        vector<Potential> expected_exact_results = {
            Potential(-0.0796731, -0.212202),
            Potential(-1.30422, 0.264889),
            Potential(-8.97201, -1.17831),
            Potential(-5.81976, -1.67169),
            Potential(-2.08631e-05, -5.76825e-05),
            Potential(-4.74252, -0.749355),
            Potential(-4.11734, -1.26302),
            Potential(-9.08724, -1.44068),
            Potential(-7.83475, -1.94273),
            Potential(-7.34204, -1.06544)
        };
        
        for (int i = 0; i < 10; i++) {
            REQUIRE(cpu[i].almost_equal(expected_exact_results[i]));
            REQUIRE(cpu_nc[i].almost_equal(expected_exact_results[i]));
    #if USING_CUDA
            REQUIRE(gpu[i] == Approx(expected_exact_results[i].total()));
            REQUIRE(gpu_nc[i] == Approx(expected_exact_results[i].total()));
    #endif
        }
    }

    SECTION("Test vs Kim's results") {
        // This test compares our implementation to Young Kim's original results.
        // We need these fudge factors to compensate for different conversion factors, etc..
        
        double LJ_fudge_factor(0.6/RT_to_kcalmol); // Kim used 0.6 for the RT to kcal/mol conversion
        double DH_fudge_factor(1.00465); // we have no idea; maybe truncated constants?
        
        vector<Potential> kim_results = {
            Potential(-0.08092694, -0.21315828),
            Potential(-1.32250603,  0.26608889),
            Potential(-9.09456258, -1.18374210),
            Potential(-5.90453581, -1.67963588),
            Potential(-0.00002116, -0.00005795),
            Potential(-4.81153585, -0.75302819),
            Potential(-4.18365757, -1.26891119),
            Potential(-9.22300158, -1.44730149),
            Potential(-7.95237776, -1.95173336),
            Potential(-7.44768961, -1.07043531)
        };
        
        double mean_e_LJ(0);
        double mean_e_DH(0);
        double mean_e_total(0);
        
        double max_e_LJ(0);
        double max_e_DH(0);
        double max_e_total(0);
        
        LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << "Comparison to Kim's results";
        LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << "Relative errors:";
        
        for (int i = 0; i < 10; i++) {
            const double & kim_LJ = kim_results[i].total_LJ();
            const double & kim_DH = kim_results[i].total_DH();
            const double & kim_total = kim_results[i].total();
            
            double LJ = cpu[i].total_LJ() * LJ_fudge_factor;
            double DH = cpu[i].total_DH() * DH_fudge_factor;
            double total = LJ + DH;
            
            // relative error, as described in Tunbridge 2011, p.113
            
            double e_LJ = fabs(kim_LJ - LJ)/fabs(kim_LJ);
            double e_DH = fabs(kim_DH - DH)/fabs(kim_DH);
            double e_total = fabs(kim_total - total)/fabs(kim_total);
            
            LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << i + 1 << ": LJ: " << e_LJ << " DH: " << e_DH << " total: " << e_total;

            mean_e_LJ += e_LJ/10;
            mean_e_DH += e_DH/10;
            mean_e_total += e_total/10;
            
            max_e_LJ = max(max_e_LJ, e_LJ);
            max_e_DH = max(max_e_DH, e_DH);
            max_e_total = max(max_e_total, e_total);
        }
        
        LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << "Mean error: LJ: " << mean_e_LJ << " DH: " << mean_e_DH << " total " << mean_e_total;
        LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << "Max error: LJ: " << max_e_LJ << " DH: " << max_e_DH << " total " << max_e_total;
        
        // Apart from the LJ in conformations 1 and 6, we seem to be doing quite well. We should probably investigate that more closely.
        
        REQUIRE(mean_e_LJ == Approx(0.001505));
        REQUIRE(mean_e_DH == Approx(0.000078));
        REQUIRE(mean_e_total == Approx(0.001073));
        
        REQUIRE(max_e_LJ == Approx(0.002866));
        REQUIRE(max_e_DH == Approx(0.000250));
        REQUIRE(max_e_total == Approx(0.002202));
    }

    SECTION("Test vs CHARMM implementation") {
        // This test compares our implementation to R. Best's CHARMM patch.
        
        double LJ_fudge_factor(1.0153); // we have no idea
        double DH_fudge_factor(1.0037); // CHARMM uses a less accurate CCELEC constant. That doesn't account for most of this, though.

        vector<Potential> charmm_results = {
            Potential(-0.081, -0.213),
            Potential(-1.322, 0.266),
            Potential(-9.095, -1.182),
            Potential(-5.903, -1.678),
            Potential(-0.000021, -0.000058),
            Potential(-4.810, -0.752),
            Potential(-4.213, -1.267),
            Potential(-9.266, -1.446),
            Potential(-7.951, -1.949),
            Potential(-7.459, -1.069)
        };
        
        double mean_e_LJ(0);
        double mean_e_DH(0);
        double mean_e_total(0);
        
        double max_e_LJ(0);
        double max_e_DH(0);
        double max_e_total(0);
        
        LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << "Comparison to CHARMM results";
        LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << "Relative errors:";
        
        for (int i = 0; i < 10; i++) {
            const double & charmm_LJ = charmm_results[i].total_LJ();
            const double & charmm_DH = charmm_results[i].total_DH();
            const double & charmm_total = charmm_results[i].total();
            
            double LJ = cpu[i].total_LJ() * LJ_fudge_factor;
            double DH = cpu[i].total_DH() * DH_fudge_factor;
            double total = LJ + DH;
            
            // relative error, as described in Tunbridge 2011, p.113
            
            double e_LJ = fabs(charmm_LJ - LJ)/fabs(charmm_LJ);
            double e_DH = fabs(charmm_DH - DH)/fabs(charmm_DH);
            double e_total = fabs(charmm_total - total)/fabs(charmm_total);
            
            LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << i + 1 << ": LJ: " << e_LJ << " DH: " << e_DH << " total: " << e_total;

            mean_e_LJ += e_LJ/10;
            mean_e_DH += e_DH/10;
            mean_e_total += e_total/10;
            
            max_e_LJ = max(max_e_LJ, e_LJ);
            max_e_DH = max(max_e_DH, e_DH);
            max_e_total = max(max_e_total, e_total);
        }
        
        LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << "Mean error: LJ: " << mean_e_LJ << " DH: " << mean_e_DH << " total " << mean_e_total;
        LOG_IF(PRINT_REFERENCE_CONFORMATIONS_FROM_TEST, INFO) << "Max error: LJ: " << max_e_LJ << " DH: " << max_e_DH << " total " << max_e_total;
        
        // This is less good, but also less important than the distance from the original results.
            
        REQUIRE(mean_e_LJ == Approx(0.002838));
        REQUIRE(mean_e_DH == Approx(0.000455));
        REQUIRE(mean_e_total == Approx(0.001722));
        
        REQUIRE(max_e_LJ == Approx(0.008679));
        REQUIRE(max_e_DH == Approx(0.001795));
        REQUIRE(max_e_total == Approx(0.005833));
    }

#if USING_CUDA
    for (int i = 0; i < 10; i++)
    {
        replicas[i].teardown_CUDA();
    }

#if CUDA_STREAMS
    teardown_CUDA(ljp_t, streams, 1);
#else
    teardown_CUDA(ljp_t);
#endif // CUDA_STREAMS
#endif // USING_CUDA
}

#endif // !LJ_REPULSIVE && !LJ_OFF
