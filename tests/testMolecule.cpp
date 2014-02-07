// #include <string>
// #include <ctype.h>
// #include <Molecule.h>
// #include <testCommon.h>
//
// class TestMolecule : public CppUnit::TestFixture
// {
//     CPPUNIT_TEST_SUITE(TestMolecule);
//     CPPUNIT_TEST(testPDBSequence);
//     CPPUNIT_TEST(testPDBCentre);
//     CPPUNIT_TEST(testPositions);
//     CPPUNIT_TEST(testRotate);
//     CPPUNIT_TEST(testTranslate);
//     CPPUNIT_TEST(testPotential);
// #if FLEXIBLE_LINKS
//     CPPUNIT_TEST(testSegments);
//     CPPUNIT_TEST(testGeometry);
//     CPPUNIT_TEST(testFlex);
//     CPPUNIT_TEST(testLocalTranslate);
//     CPPUNIT_TEST(testCrankshaft);
// #endif // FLEXIBLE_LINKS
//     CPPUNIT_TEST_SUITE_END();
//
// private:
//     Molecule ubiquitin;
//     Molecule angiotensin;
//     Molecule polyalanine8;
//     AminoAcids aminoAcidData;
//
// public:
//     void setUp();
//     void testPDBSequence();
//     void testPDBCentre();
//     void testPositions();
//     void testRotate();
//     void testTranslate();
//     void testPotential();
// #if FLEXIBLE_LINKS
//     void testSegments();
//     void testGeometry();
//     void testFlex();
//     void testLocalTranslate();
//     void testCrankshaft();
// #endif // FLEXIBLE_LINKS
// };
//
// CPPUNIT_TEST_SUITE_REGISTRATION(TestMolecule);
//
// void TestMolecule::setUp()
// {
//     aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);
//
//     ubiquitin.init("tests/1UBQ.pdb", aminoAcidData, 0, 1000.0f);
//     angiotensin.init("tests/angiotensin.pdb", aminoAcidData, 0, 1000.0f);
//
//     // small bounding value, so we can test wrapping over the boundary
//     polyalanine8.init("tests/ala8.pdb", aminoAcidData, 0, 15.0f);
// }
//
// void TestMolecule::testPDBSequence()
// {
//     const int ubiquitin_size = 76;
//     char ubiquitin_names[ubiquitin_size][4] = {"MET", "GLN", "ILE", "PHE", "VAL", "LYS", "THR", "LEU", "THR", "GLY", "LYS", "THR", "ILE", "THR", "LEU", "GLU", "VAL", "GLU", "PRO", "SER", "ASP", "THR", "ILE", "GLU", "ASN", "VAL", "LYS", "ALA", "LYS", "ILE", "GLN", "ASP", "LYS", "GLU", "GLY", "ILE", "PRO", "PRO", "ASP", "GLN", "GLN", "ARG", "LEU", "ILE", "PHE", "ALA", "GLY", "LYS", "GLN", "LEU", "GLU", "ASP", "GLY", "ARG", "THR", "LEU", "SER", "ASP", "TYR", "ASN", "ILE", "GLN", "LYS", "GLU", "SER", "THR", "LEU", "HIS", "LEU", "VAL", "LEU", "ARG", "LEU", "ARG", "GLY", "GLY"};
//
//     const int angiotensin_size = 10;
//     char angiotensin_names[angiotensin_size][4] = {"ASP","ARG","VAL","TYR","ILE","HIS","PRO","PHE","HIS","LEU"};
//
//     CPPUNIT_ASSERT_EQUAL(ubiquitin_size, ubiquitin.residueCount);
//     CPPUNIT_ASSERT_EQUAL(angiotensin_size, angiotensin.residueCount);
// #if FLEXIBLE_LINKS
//     CPPUNIT_ASSERT_EQUAL(angiotensin_size - 1, angiotensin.linkCount);
// #endif // FLEXIBLE_LINKS
//
//     int expected;
//     int got;
//
//     // TODO: rewrite these so they compare strings instead
//     for (int i = 0; i < ubiquitin_size; i++)
//     {
//         expected = aminoAcidData.getAminoAcidIndex(ubiquitin_names[i]);
//         got = ubiquitin.Residues[i].aminoAcidIndex;
//         CPPUNIT_ASSERT_EQUAL(expected, got);
//     }
//
//     for (int i = 0; i < angiotensin_size; i++)
//     {
//         expected = aminoAcidData.getAminoAcidIndex(angiotensin_names[i]);
//         got = angiotensin.Residues[i].aminoAcidIndex;
//         CPPUNIT_ASSERT_EQUAL(expected, got);
//     }
// }
//
// void TestMolecule::testPDBCentre()
// {
//     ASSERT_VECTOR3F_EQUALS(Vector3f(30.442894736842103, 29.00898684210527, 15.5163947368421), ubiquitin.center);
// }
//
// void TestMolecule::testPositions()
// {
//     ASSERT_VECTOR3F_EQUALS(Vector3f(1.458, 0, 0), polyalanine8.Residues[0].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(3.368, 3.081, 1.154), polyalanine8.Residues[1].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(6.993, 3.658, 0.156), polyalanine8.Residues[2].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(8.955, 6.649, 1.45), polyalanine8.Residues[3].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(12.524, 7.321, 0.318), polyalanine8.Residues[4].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(14.546, 10.212, 1.741), polyalanine8.Residues[5].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(18.052, 10.989, 0.485), polyalanine8.Residues[6].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(20.141, 13.77, 2.026), polyalanine8.Residues[7].position);
//
//     ASSERT_VECTOR3F_EQUALS(Vector3f(10.754625, 6.96, 0.91625), polyalanine8.center);
//
//     ASSERT_VECTOR3F_EQUALS(Vector3f(-9.29663, -6.96, -0.91625), polyalanine8.Residues[0].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(-7.38663, -3.879, 0.23775), polyalanine8.Residues[1].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(-3.76163, -3.302, -0.76025), polyalanine8.Residues[2].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(-1.79963, -0.311, 0.53375), polyalanine8.Residues[3].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(1.76937, 0.361, -0.59825), polyalanine8.Residues[4].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(3.79137, 3.252, 0.82475), polyalanine8.Residues[5].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(7.29737, 4.029, -0.43125), polyalanine8.Residues[6].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(9.38637, 6.81, 1.10975), polyalanine8.Residues[7].relativePosition);
// }
//
// void TestMolecule::testRotate()
// {
//     // TODO
//     // Rotate
//     // Check that absolute positions have changed
//     // Check that centre has not changed
//     // Check that relative positions have not changed
// }
//
// void TestMolecule::testTranslate()
// {
//     // Translate
//
//     Vector3f translation_vector(1, 1, 1);
//     polyalanine8.translate(translation_vector);
//
//     // Check that absolute positions have changed
//
//     ASSERT_VECTOR3F_EQUALS(Vector3f(2.458, 1, 1), polyalanine8.Residues[0].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(4.368, 4.081, 2.154), polyalanine8.Residues[1].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(7.993, 4.658, 1.156), polyalanine8.Residues[2].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(9.955, 7.649, 2.45), polyalanine8.Residues[3].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(13.524, 8.321, 1.318), polyalanine8.Residues[4].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(15.546, 11.212, 2.741), polyalanine8.Residues[5].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(19.052, 11.989, 1.485), polyalanine8.Residues[6].position);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(21.141, 14.77, 3.026), polyalanine8.Residues[7].position);
//
//     // Check that centre has changed
//
//     ASSERT_VECTOR3F_EQUALS(Vector3f(11.754625, 7.96, 1.91625), polyalanine8.center);
//
//     // Check that relative positions have not changed
//
//     ASSERT_VECTOR3F_EQUALS(Vector3f(-9.29663, -6.96, -0.91625), polyalanine8.Residues[0].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(-7.38663, -3.879, 0.23775), polyalanine8.Residues[1].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(-3.76163, -3.302, -0.76025), polyalanine8.Residues[2].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(-1.79963, -0.311, 0.53375), polyalanine8.Residues[3].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(1.76937, 0.361, -0.59825), polyalanine8.Residues[4].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(3.79137, 3.252, 0.82475), polyalanine8.Residues[5].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(7.29737, 4.029, -0.43125), polyalanine8.Residues[6].relativePosition);
//     ASSERT_VECTOR3F_EQUALS(Vector3f(9.38637, 6.81, 1.10975), polyalanine8.Residues[7].relativePosition);
// }
//
// void TestMolecule::testPotential()
// {
//     Potential potential_ala = polyalanine8.E(true); // with internal LJ and DH
//     Potential potential_ang = angiotensin.E(true); // with internal LJ and DH
// #if !LJ_REPULSIVE && !LJ_OFF
//     double expected_ala_potential[6] = {0.945580, 0.000000, 0.045626, 1.272977, 10.872785, 13.136967};
//     ASSERT_POTENTIAL_EQUALS(expected_ala_potential, potential_ala);
//
//     double expected_ang_potential[6] = {-0.253939, 0.052889, 0.642773, 4.887786, 15.871156, 21.200665};
//     ASSERT_POTENTIAL_EQUALS(expected_ang_potential, potential_ang);
// #elif LJ_OFF
//     CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, potential_ala.total_LJ(), 0.00001);
//     CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, potential_ang.total_LJ(), 0.00001);
// #endif // !LJ_REPULSIVE && !LJ_OFF
// }
//
// #if FLEXIBLE_LINKS
// void TestMolecule::testSegments()
// {
//     CPPUNIT_ASSERT_EQUAL(1, angiotensin.segmentCount);
// //     CPPUNIT_ASSERT_EQUAL(1, angiotensin.linkerCount);
//
//     CPPUNIT_ASSERT_EQUAL(0, angiotensin.Segments[0].start);
//     CPPUNIT_ASSERT_EQUAL(9, angiotensin.Segments[0].end);
//     CPPUNIT_ASSERT_EQUAL(true, angiotensin.Segments[0].flexible);
//
// //     CPPUNIT_ASSERT_EQUAL(0, angiotensin.Linkers[0]->start);
// //     CPPUNIT_ASSERT_EQUAL(9, angiotensin.Linkers[0]->end);
// }
//
// void TestMolecule::testGeometry()
// {
//     float expected_bond_lengths[9] = {
//         3.821749687194824,
//         3.8273043632507324,
//         3.8032095432281494,
//         3.8131093978881836,
//         3.76322865486145,
//         3.8027637004852295,
//         3.7888474464416504,
//         3.8237955570220947,
//         3.8045248985290527
//     };
//
//     float expected_angles[9] = {
//         0.0, // padding for simplicity
//         123.56746673583984,
//         86.7016830444336,
//         126.89990234375,
//         126.8516845703125,
//         129.04217529296875,
//         84.29325103759766,
//         130.4828338623047,
//         97.1257095336914,
//     };
//
//     float expected_torsion_angles[9] = {
//         0.0, // padding for simplicity
//         -165.99693298339844,
//         26.120737075805664,
//         -171.85838317871094,
//         -170.0393524169922,
//         -156.48416137695313,
//         4.634008884429932,
//         150.31800842285156,
//         0.0 // padding for simplicity
//     };
//
//     Potential e = angiotensin.E();
//
//     for (size_t i = 0; i < angiotensin.linkCount; i++)
//     {
//         CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_bond_lengths[i], angiotensin.Links[i].pseudo_bond, 0.00001);
//
//         // convert degrees (from VMD) to radians
//         CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_angles[i]/57.29577951308232, angiotensin.Residues[i].pseudo_angle, 0.00001);
//
//         // convert degrees (from VMD) to radians
//         CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_torsion_angles[i]/57.29577951308232, angiotensin.Links[i].pseudo_torsion, 0.00001);
//     }
// }
//
// void TestMolecule::testFlex()
// {
//     // TODO
// }
//
// void TestMolecule::testLocalTranslate()
// {
//     // TODO
// }
//
// void TestMolecule::testCrankshaft()
// {
//     // TODO
// }
//
// #endif // FLEXIBLE_LINKS
//
