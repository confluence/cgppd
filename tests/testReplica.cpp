// #include "definitions.h"
// #if USING_CUDA
// #include <cutil.h>
// #include <cutil_inline.h>
// #endif
//
// #include <Simulation.h>
// #include <testCommon.h>
//
// class TestReplica : public CppUnit::TestFixture
// {
//     CPPUNIT_TEST_SUITE(TestReplica);
//     CPPUNIT_TEST(testCPUandGPUPotential);
//     CPPUNIT_TEST_SUITE_END();
//
// private:
//     AminoAcids aminoAcidData;
//     Replica replica;
//     float * ljp_t;
//     float testboxdim;
//
// public:
//     void setUp();
//     void testCPUandGPUPotential();
//     void tearDown();
// };
//
// CPPUNIT_TEST_SUITE_REGISTRATION(TestReplica);
//
// void TestReplica::setUp()
// {
//     aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);
//     testboxdim = 118.4f;
//
// #if USING_CUDA
// #if CUDA_STREAMS
//     setup_CUDA(0, testboxdim, ljp_t, &aminoAcidData, &replica.cudaStream, 1);
// #else
//     setup_CUDA(0, testboxdim, ljp_t, &aminoAcidData);
// #endif // CUDA_STREAMS
// #endif // USING_CUDA
//
//     replica.aminoAcids = aminoAcidData;
//     replica.label = 1;
//     replica.boundingValue = testboxdim;
//
//     replica.reserveContiguousMoleculeArray(2);
//
//     replica.loadMolecule("data/conf1_flexible/1a.pdb");
//     replica.loadMolecule("data/conf1_flexible/1b.pdb");
//
// #if USING_CUDA
//     int cuda_blockSize;
//     (replica.residueCount < 1024) ? cuda_blockSize = 32 : cuda_blockSize = 64;
//     replica.setBlockSize(cuda_blockSize);
// #endif
//
// }
//
// void TestReplica::tearDown()
// {
// #if USING_CUDA
// #if CUDA_STREAMS
//     teardown_CUDA(ljp_t, &replica.cudaStream, 1);
// #else
//     teardown_CUDA(ljp_t);
// #endif // CUDA_STREAMS
// #endif // USING_CUDA
//     cout.flush();
// }
//
// void TestReplica::testCPUandGPUPotential()
// {
//     Replica child_replica;
//     argdata parameters;
//     child_replica.init_child_replica(replica, 1, 300.0f, 0.2f, 0.5f, parameters);
// #if USING_CUDA
// #if CUDA_STREAMS
//     child_replica.setup_CUDA(ljp_t, &replica.cudaStream, 0);
// #else
//     child_replica.setup_CUDA(ljp_t);
// #endif // CUDA_STREAMS
// #endif // USING_CUDA
//
// #if !LJ_REPULSIVE && !LJ_OFF
//     Potential expected_flexible_potential(300.864756, 0.337112, 0.180426, 1.555213, 3.495782);
// #elif LJ_OFF
//     Potential expected_flexible_potential(0.0, 0.337112, 0.180426, 1.555213, 3.495782);
// #elif LJ_REPULSIVE
//     Potential expected_flexible_potential(617.860549, 0.337112, 0.180426, 1.555213, 3.495782); // check value
// # endif // !LJ_REPULSIVE && !LJ_OFF
//
//     ASSERT_POTENTIALS_EQUAL(expected_flexible_potential, child_replica.E());
// #if USING_CUDA
//     CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_flexible_potential[5], child_replica.EonDevice(), 0.001);
// #endif // USING_CUDA
//
//     child_replica.calculate_rigid_potential_only = true;
// #if USING_CUDA
//     child_replica.ReplicaDataToDevice();
// #endif // USING_CUDA
//
// #if !LJ_REPULSIVE && !LJ_OFF
//     Potential expected_rigid_potential(-0.080765, -0.212940, 0.000000, -0.000000, 0.000000);
// #elif LJ_OFF
//     Potential expected_rigid_potential(0.0, -0.212940, 0.000000, -0.000000, 0.000000);
// #elif LJ_REPULSIVE
//     Potential expected_rigid_potential(-0.279042, -0.212940, 0.000000, -0.000000, 0.000000); // check value
// # endif // !LJ_REPULSIVE && !LJ_OFF
//
//     ASSERT_POTENTIALS_EQUAL(expected_rigid_potential, child_replica.E());
// #if USING_CUDA
//     CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_rigid_potential[5], child_replica.EonDevice(), 0.001);
// #endif // USING_CUDA
//
// #if USING_CUDA
//     child_replica.teardown_CUDA();
// #endif
// }
//
