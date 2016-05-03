#include <Simulation.h>
#include <testCommon.h>

TEST_CASE("Replica", "[replica]") {

    //GLOBAL

    AminoAcids aminoAcidData;
    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);

    float bounding_value = 118.4f;
#if USING_CUDA
#if CUDA_STREAMS
    cudaStream_t streams[1];
#endif // CUDA_STREAMS
    float * ljp_t;
#endif // USING_CUDA
    
#if USING_CUDA
#if CUDA_STREAMS
    setup_CUDA(0, bounding_value, ljp_t, &aminoAcidData, streams, 1);
#else
    setup_CUDA(0, bounding_value, ljp_t, &aminoAcidData);
#endif // CUDA_STREAMS
#endif // USING_CUDA

    // PARAMETERS

    argdata parameters;
    parameters.bound = bounding_value;
#if USING_CUDA
    parameters.auto_blockdim = false;
#endif

    moldata m1;
    m1.translate = true;
    sprintf(m1.pdbfilename, "data/conf1/1a.pdb");

    moldata m2;
    m2.translate = true;
    sprintf(m2.pdbfilename, "data/conf1/1b.pdb");
    
    segdata seg1;
    seg1.residue_indices = {72, 73, 74, 75};
    m2.segments.push_back(seg1);
    
    parameters.mdata.push_back(m1);
    parameters.mdata.push_back(m2);

    // INITIAL REPLICA
    Replica initial_replica;
    initial_replica.init_first_replica(parameters, aminoAcidData, 2);
    initial_replica.label = 1;

    // ONE CHILD REPLICA
    Replica replica;
    replica.init_child_replica(initial_replica, 1, 300.0f, 0.2f, 0.5f, parameters);
    
#if USING_CUDA
#if CUDA_STREAMS
    replica.setup_CUDA(ljp_t, streams, 0);
#else
    replica.setup_CUDA(ljp_t);
#endif // CUDA_STREAMS
#endif // USING_CUDA

#if !ASSUME_POLYMER_FOLDING_TEST // ignore this whole test if we're using the shortcut conditions for a single polymer
    SECTION("Potential") {
        
#if !LJ_REPULSIVE && !LJ_OFF
        double flex_lj(-0.068205);
        double rigid_lj(-0.079673);
#elif LJ_OFF
        double flex_lj(0);
        double rigid_lj(0);
#elif LJ_REPULSIVE
        double flex_lj(-0.527692);
        double rigid_lj(-0.275268);
#endif // !LJ_REPULSIVE && !LJ_OFF

        Potential expected_flexible_potential(flex_lj, -0.240383, 0.180426, 1.555213, 5.420859);
        Potential expected_rigid_potential(rigid_lj, -0.212202, 0.000000, -0.000000, 0.000000);

        REQUIRE((replica.E() + replica.internal_molecule_E(true)).almost_equal(expected_flexible_potential));
#if USING_CUDA
        REQUIRE((replica.EonDevice() + replica.internal_molecule_E(false).total()) == Approx(expected_flexible_potential.total()));
        REQUIRE(replica.EonDeviceNC() == Approx(expected_rigid_potential.total())); // THIS SHOULD BE RIGID!
#if CUDA_STREAMS
        replica.EonDeviceAsync();
        REQUIRE((replica.SumGridResults() + replica.internal_molecule_E(false).total()) == Approx(expected_flexible_potential.total()));
#endif // CUDA_STREAMS
#else // check that the CPU-only NC potential works
        REQUIRE(replica.E(&replica.molecules[0], &replica.molecules[1]).total() == Approx(expected_rigid_potential.total()));
#endif // USING_CUDA

        replica.calculate_rigid_potential_only = true;
#if USING_CUDA
        replica.ReplicaDataToDevice();
#endif // USING_CUDA

        REQUIRE(replica.E().almost_equal(expected_rigid_potential));
#if USING_CUDA
        REQUIRE(replica.EonDevice() == Approx(expected_rigid_potential.total()));
        REQUIRE(replica.EonDeviceNC() == Approx(expected_rigid_potential.total()));
#if CUDA_STREAMS
        replica.EonDeviceAsync();
        REQUIRE(replica.SumGridResults() == Approx(expected_rigid_potential.total()));
#endif // CUDA_STREAMS
#else // check that the CPU-only NC potential works
        REQUIRE(replica.E(&replica.molecules[0], &replica.molecules[1]).total() == Approx(expected_rigid_potential.total()));
#endif // USING_CUDA
    }
#endif // ASSUME_POLYMER_FOLDING_TEST

    SECTION("Geometry") {
        Molecule & m1 = replica.molecules[0];
        Molecule & m2 = replica.molecules[1];
        
        // Test counting of chains, rigid domains and segment bonds
        
        REQUIRE((int)m1.graph.rigid_domains.size() == 1);
        REQUIRE((int)m1.graph.segment_bonds.size() == 0);
        REQUIRE(m1.graph.num_chains == 1);
        
        REQUIRE((int)m2.graph.rigid_domains.size() == 1);
        REQUIRE((int)m2.graph.segment_bonds.size() == 0);
        REQUIRE(m2.graph.num_chains == 1);
        
        // Test assignment of UIDs to chains, rigid domains and segment bonds
        
        REQUIRE((int)m1.graph.domain_uid.size() == 1);
        REQUIRE((int)m1.graph.bond_uid.size() == 0);
        REQUIRE((int)m1.graph.chain_uid.size() == 1);
        
        REQUIRE((int)m2.graph.domain_uid.size() == 1);
        REQUIRE((int)m2.graph.bond_uid.size() == 0);
        REQUIRE((int)m2.graph.chain_uid.size() == 1);
        
        REQUIRE(m1.graph.domain_uid[0] == 1);
        REQUIRE(m2.graph.domain_uid[0] == 2);
        
        REQUIRE(m1.graph.chain_uid[0] == 1);
        REQUIRE(m2.graph.chain_uid[0] == 2);
        
        // Test application of UIDs to residues and packing of residue floats for GPU
        
        for (int i = 0; i < 24; i++) {
            REQUIRE(m1.Residues[i].rigid_domain_UID == 1);
            REQUIRE(m1.Residues[i].chain_UID == 1);
            REQUIRE(m1.Residues[i].segment_bond_UID == 0);
#if USING_CUDA
            REQUIRE(m1.Residues[i].pos_w == Approx(1.0f));
            REQUIRE(m1.Residues[i].meta_w == Approx(0.25f));
#endif
        }
        
        for (int i = 0; i < 73; i++) {
            REQUIRE(m2.Residues[i].rigid_domain_UID == 2);
            REQUIRE(m2.Residues[i].chain_UID == 2);
            REQUIRE(m2.Residues[i].segment_bond_UID == 0);
#if USING_CUDA
            REQUIRE(m2.Residues[i].pos_w == Approx(2.0f));
            REQUIRE(m2.Residues[i].meta_w == Approx(1.5f));
#endif
        }
        
        for (int i = 73; i < 76; i++) {
            REQUIRE(m2.Residues[i].rigid_domain_UID == 0);
            REQUIRE(m2.Residues[i].chain_UID == 2);
            REQUIRE(m2.Residues[i].segment_bond_UID == 0);
#if USING_CUDA
            REQUIRE(m2.Residues[i].pos_w == Approx(0.0f));
            REQUIRE(m2.Residues[i].meta_w == Approx(1.5f));
#endif
        }
    }

#if USING_CUDA
    replica.teardown_CUDA();
#if CUDA_STREAMS
    teardown_CUDA(ljp_t, streams, 1);
#else
    teardown_CUDA(ljp_t);
#endif // CUDA_STREAMS
#endif // USING_CUDA
    cout.flush();
}

