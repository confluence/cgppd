#include "definitions.h"
#if USING_CUDA
#include <cutil.h>
#include <cutil_inline.h>
#endif

#include <Simulation.h>
#include <testCommon.h>

class TestReplica : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestReplica);
    CPPUNIT_TEST(testCPUandGPUPotential);
    CPPUNIT_TEST(testPotentialGeometry);
    CPPUNIT_TEST_SUITE_END();

private:
    AminoAcids aminoAcidData;
    Replica initial_replica;
    Replica replica;
    float bounding_value;
#if USING_CUDA
#if CUDA_STREAMS
    cudaStream_t streams[1];
#endif // CUDA_STREAMS
    float * ljp_t;
#endif // USING_CUDA

public:
    void setUp();
    void testCPUandGPUPotential();
    void testPotentialGeometry();
    void tearDown();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestReplica);

void TestReplica::setUp()
{
    // GLOBAL STUFF
    
    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);
    
    bounding_value = 118.4f;
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

    initial_replica.init_first_replica(parameters, aminoAcidData, 2);
    initial_replica.label = 1;

    // ONE CHILD REPLICA
    
    replica.init_child_replica(initial_replica, 1, 300.0f, 0.2f, 0.5f, parameters);
    
#if USING_CUDA
#if CUDA_STREAMS
    replica.setup_CUDA(ljp_t, streams, 0);
#else
    replica.setup_CUDA(ljp_t);
#endif // CUDA_STREAMS
#endif // USING_CUDA
}

void TestReplica::tearDown()
{
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

void TestReplica::testCPUandGPUPotential()
{
#if !LJ_REPULSIVE && !LJ_OFF
    double flex_lj(-0.068205);
#elif LJ_OFF
    double flex_lj(0);
#elif LJ_REPULSIVE
    double flex_lj(-0.527692);
# endif // !LJ_REPULSIVE && !LJ_OFF
    Potential expected_flexible_potential(flex_lj, -0.240383, 0.180426, 1.555213, 5.420859);
    ASSERT_POTENTIALS_EQUAL(expected_flexible_potential, replica.E() + replica.internal_molecule_E(true));

#if USING_CUDA
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_flexible_potential.total(), replica.EonDevice() + replica.internal_molecule_E(false).total(), 0.00001);
#endif // USING_CUDA

    replica.calculate_rigid_potential_only = true;
#if USING_CUDA
    replica.ReplicaDataToDevice();
#endif // USING_CUDA

#if !LJ_REPULSIVE && !LJ_OFF
    double rigid_lj(-0.079673);
#elif LJ_OFF
    double rigid_lj(0);
#elif LJ_REPULSIVE
    double rigid_lj(-0.275268);
# endif // !LJ_REPULSIVE && !LJ_OFF

    Potential expected_rigid_potential(rigid_lj, -0.212202, 0.000000, -0.000000, 0.000000);
    ASSERT_POTENTIALS_EQUAL(expected_rigid_potential, replica.E());

#if USING_CUDA
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_rigid_potential.total(), replica.EonDevice(), 0.00001);
#endif // USING_CUDA
}

void TestReplica::testPotentialGeometry()
{
    Molecule & m1 = replica.molecules[0];
    Molecule & m2 = replica.molecules[1];
    
    // Test counting of chains, rigid domains and segment bonds
    
    CPPUNIT_ASSERT_EQUAL(1, (int)m1.graph.rigid_domains.size());
    CPPUNIT_ASSERT_EQUAL(0, (int)m1.graph.segment_bonds.size());
    CPPUNIT_ASSERT_EQUAL(1, m1.graph.num_chains);
    
    CPPUNIT_ASSERT_EQUAL(1, (int)m2.graph.rigid_domains.size());
    CPPUNIT_ASSERT_EQUAL(0, (int)m2.graph.segment_bonds.size());
    CPPUNIT_ASSERT_EQUAL(1, m2.graph.num_chains);
    
    // Test assignment of UIDs to chains, rigid domains and segment bonds
    
    CPPUNIT_ASSERT_EQUAL(1, (int)m1.graph.domain_uid.size());
    CPPUNIT_ASSERT_EQUAL(0, (int)m1.graph.bond_uid.size());
    CPPUNIT_ASSERT_EQUAL(1, (int)m1.graph.chain_uid.size());
    
    CPPUNIT_ASSERT_EQUAL(1, (int)m2.graph.domain_uid.size());
    CPPUNIT_ASSERT_EQUAL(0, (int)m2.graph.bond_uid.size());
    CPPUNIT_ASSERT_EQUAL(1, (int)m2.graph.chain_uid.size());
    
    CPPUNIT_ASSERT_EQUAL(1, m1.graph.domain_uid[0]);
    CPPUNIT_ASSERT_EQUAL(2, m2.graph.domain_uid[0]);
    
    CPPUNIT_ASSERT_EQUAL(1, m1.graph.chain_uid[0]);
    CPPUNIT_ASSERT_EQUAL(2, m2.graph.chain_uid[0]);
    
    // Test application of UIDs to residues and packing of residue floats for GPU
    
    for (int i = 0; i < 24; i++) {
        CPPUNIT_ASSERT_EQUAL(1, m1.Residues[i].rigid_domain_UID);
        CPPUNIT_ASSERT_EQUAL(1, m1.Residues[i].chain_UID);
        CPPUNIT_ASSERT_EQUAL(0, m1.Residues[i].segment_bond_UID);
#if USING_CUDA
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0f, m1.Residues[i].pos_w, 0.00001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(float(i + 1) + 0.5f, m1.Residues[i].meta_w, 0.00001);
#endif
    }
    
    for (int i = 0; i < 73; i++) {
        CPPUNIT_ASSERT_EQUAL(2, m2.Residues[i].rigid_domain_UID);
        CPPUNIT_ASSERT_EQUAL(2, m2.Residues[i].chain_UID);
        CPPUNIT_ASSERT_EQUAL(0, m2.Residues[i].segment_bond_UID);
#if USING_CUDA
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0f, m2.Residues[i].pos_w, 0.00001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(float(i + 1) + 0.25f, m2.Residues[i].meta_w, 0.00001);
#endif
    }
    
    for (int i = 73; i < 76; i++) {
        CPPUNIT_ASSERT_EQUAL(0, m2.Residues[i].rigid_domain_UID);
        CPPUNIT_ASSERT_EQUAL(2, m2.Residues[i].chain_UID);
        CPPUNIT_ASSERT_EQUAL(0, m2.Residues[i].segment_bond_UID);
#if USING_CUDA
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0f, m2.Residues[i].pos_w, 0.00001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(float(i + 1) + 0.25f, m2.Residues[i].meta_w, 0.00001);
#endif
    }
}
