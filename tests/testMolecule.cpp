#include <string>
#include <ctype.h>
#include <structures.h>
#include <Molecule.h>
#include <testCommon.h>

class TestMolecule : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestMolecule);
    CPPUNIT_TEST(testPDBSequence);
    CPPUNIT_TEST(testStartingPositions);
    CPPUNIT_TEST(testRotate);
    CPPUNIT_TEST(testTranslate);
    CPPUNIT_TEST(testCopy);
//     CPPUNIT_TEST(testPotential);
#if FLEXIBLE_LINKS
    CPPUNIT_TEST(testSegments);
    CPPUNIT_TEST(testGeometry);
    CPPUNIT_TEST(testFlex);
    CPPUNIT_TEST(testLocalTranslate);
    CPPUNIT_TEST(testCrankshaft);
#endif // FLEXIBLE_LINKS
    CPPUNIT_TEST_SUITE_END();

private:
    moldata ubq_data;
    moldata ang_data;
    moldata ala8_data;

    Molecule ubiquitin;
    Molecule angiotensin;
    Molecule polyalanine8;

    AminoAcids aminoAcidData;
    
    vector<Vector3f> ala8_start;

public:
    void setUp();
    void testPDBSequence();
    void testStartingPositions();
    void testRotate();
    void testTranslate();
    void testCopy();
//     void testPotential();
#if FLEXIBLE_LINKS
    void testSegments();
    void testGeometry();
    void testFlex();
    void testLocalTranslate();
    void testCrankshaft();
#endif // FLEXIBLE_LINKS
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestMolecule);

void TestMolecule::setUp()
{
    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);

    strcpy(ubq_data.pdbfilename, "tests/1UBQ.pdb");
    segdata seg;
    seg.residue_indices.push_back(72);
    seg.residue_indices.push_back(73);
    seg.residue_indices.push_back(74);
    seg.residue_indices.push_back(75);
    ubq_data.segments.push_back(seg);
    ubiquitin.init(ubq_data, aminoAcidData, 0, 1000.0f);

    strcpy(ang_data.pdbfilename, "tests/angiotensin.pdb");
    ang_data.all_flexible = true;
    angiotensin.init(ang_data, aminoAcidData, 0, 1000.0f);

    strcpy(ala8_data.pdbfilename, "tests/ala8.pdb");
    ala8_data.all_flexible = true;
    // small bounding value, so we can test wrapping over the boundary
    polyalanine8.init(ala8_data, aminoAcidData, 0, 15.0f);
    
    ala8_start = {
        Vector3f(-9.29663, -6.96, -0.91625),
        Vector3f(-7.38663, -3.879, 0.23775),
        Vector3f(-3.76163, -3.302, -0.76025),
        Vector3f(-1.79963, -0.311, 0.53375),
        Vector3f(1.76937, 0.361, -0.59825),
        Vector3f(3.79137, 3.252, 0.82475),
        Vector3f(7.29737, 4.029, -0.43125),
        Vector3f(9.38637, 6.81, 1.10975)
    };
}

void TestMolecule::testPDBSequence()
{
    vector<string> ubiquitin_names = {"MET", "GLN", "ILE", "PHE", "VAL", "LYS", "THR", "LEU", "THR", "GLY", "LYS", "THR", "ILE", "THR", "LEU", "GLU", "VAL", "GLU", "PRO", "SER", "ASP", "THR", "ILE", "GLU", "ASN", "VAL", "LYS", "ALA", "LYS", "ILE", "GLN", "ASP", "LYS", "GLU", "GLY", "ILE", "PRO", "PRO", "ASP", "GLN", "GLN", "ARG", "LEU", "ILE", "PHE", "ALA", "GLY", "LYS", "GLN", "LEU", "GLU", "ASP", "GLY", "ARG", "THR", "LEU", "SER", "ASP", "TYR", "ASN", "ILE", "GLN", "LYS", "GLU", "SER", "THR", "LEU", "HIS", "LEU", "VAL", "LEU", "ARG", "LEU", "ARG", "GLY", "GLY"};

    vector<string> angiotensin_names = {"ASP","ARG","VAL","TYR","ILE","HIS","PRO","PHE","HIS","LEU"};

    CPPUNIT_ASSERT_EQUAL((int)ubiquitin_names.size(), ubiquitin.residueCount);
    CPPUNIT_ASSERT_EQUAL((int)angiotensin_names.size(), angiotensin.residueCount);

    for (int i = 0; i < ubiquitin_names.size(); i++)
    {
        char name[4];
        strcpy(name, ubiquitin_names[i].c_str());
        CPPUNIT_ASSERT_EQUAL(aminoAcidData.getAminoAcidIndex(name), ubiquitin.Residues[i].aminoAcidIndex);
    }
    
    for (int i = 0; i < angiotensin_names.size(); i++)
    {
        char name[4];
        strcpy(name, angiotensin_names[i].c_str());
        CPPUNIT_ASSERT_EQUAL(aminoAcidData.getAminoAcidIndex(name), angiotensin.Residues[i].aminoAcidIndex);
    }
}

void TestMolecule::testStartingPositions()
{
    // We set the position to zero by default
    ASSERT_VECTOR3F_EQUALS(Vector3f(0, 0, 0), polyalanine8.center);
    
    // TODO: make a test for using a translation by 0 rather than setting the position
    
    for (int i = 0; i < 8; i++) {
        ASSERT_VECTOR3F_EQUALS(ala8_start[i], polyalanine8.Residues[i].relativePosition);
    }

    // These are the same because it's centered on zero
    for (int i = 0; i < 8; i++) {
        ASSERT_VECTOR3F_EQUALS(ala8_start[i], polyalanine8.Residues[i].position);
    }
}

void TestMolecule::testRotate()
{
    // Rotate 180 degrees about the x axis
    Vector3f rotation_vector(1, 0, 0);
    polyalanine8.rotate(rotation_vector, M_PIl);

    // Check that absolute positions have changed -- x coordinate should be the same; y and z should be flipped
    for (int i = 0; i < polyalanine8.residueCount; i++) {
        Vector3f expected_position = Vector3f(ala8_start[i].x, -ala8_start[i].y, -ala8_start[i].z);
        ASSERT_VECTOR3F_EQUALS(expected_position, polyalanine8.Residues[i].position);
    }
    
    // Check that relative positions have changed -- should be the same
    for (int i = 0; i < polyalanine8.residueCount; i++) {
        Vector3f expected_position = Vector3f(ala8_start[i].x, -ala8_start[i].y, -ala8_start[i].z);
        ASSERT_VECTOR3F_EQUALS(expected_position, polyalanine8.Residues[i].relativePosition);
    }

    // Check that centre has not changed
    ASSERT_VECTOR3F_EQUALS(Vector3f(0, 0, 0), polyalanine8.center);
}

void TestMolecule::testTranslate()
{
    // Translate

    Vector3f translation_vector(1, 1, 1);
    polyalanine8.translate(translation_vector);

    // Check that absolute positions have changed
    
    for (int i = 0; i < 8; i++) {
        ASSERT_VECTOR3F_EQUALS(ala8_start[i] + translation_vector, polyalanine8.Residues[i].position);
    }

    // Check that centre has changed

    ASSERT_VECTOR3F_EQUALS(Vector3f(1, 1, 1), polyalanine8.center);

    // Check that relative positions have not changed
    for (int i = 0; i < 8; i++) {
        ASSERT_VECTOR3F_EQUALS(ala8_start[i], polyalanine8.Residues[i].relativePosition);
    }
}

void TestMolecule::testCopy()
{
    // TODO
}

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

#if FLEXIBLE_LINKS
void TestMolecule::testSegments()
{
   CPPUNIT_ASSERT_EQUAL(3, (int)ubiquitin.graph.bonds.size());
   CPPUNIT_ASSERT_EQUAL(3, (int)ubiquitin.graph.angles.size());
   CPPUNIT_ASSERT_EQUAL(3, (int)ubiquitin.graph.torsions.size());

   CPPUNIT_ASSERT_EQUAL(7, (int)polyalanine8.graph.bonds.size());
   CPPUNIT_ASSERT_EQUAL(6, (int)polyalanine8.graph.angles.size());
   CPPUNIT_ASSERT_EQUAL(5, (int)polyalanine8.graph.torsions.size());

   CPPUNIT_ASSERT_EQUAL(9, (int)angiotensin.graph.bonds.size());
   CPPUNIT_ASSERT_EQUAL(8, (int)angiotensin.graph.angles.size());
   CPPUNIT_ASSERT_EQUAL(7, (int)angiotensin.graph.torsions.size());
    
}

void TestMolecule::testGeometry()
{
    vector<float> expected_bond_lengths = {
        3.821749687194824,
        3.8273043632507324,
        3.8032095432281494,
        3.8131093978881836,
        3.76322865486145,
        3.8027637004852295,
        3.7888474464416504,
        3.8237955570220947,
        3.8045248985290527
    };

    vector<float> expected_angles = {
        123.56746673583984,
        86.7016830444336,
        126.89990234375,
        126.8516845703125,
        129.04217529296875,
        84.29325103759766,
        130.4828338623047,
        97.1257095336914,
    };

    vector<float> expected_torsion_angles = {
        -165.99693298339844,
        26.120737075805664,
        -171.85838317871094,
        -170.0393524169922,
        -156.48416137695313,
        4.634008884429932,
        150.31800842285156,
    };

    Potential e = angiotensin.E();
    
    // TODO TODO TODO update this
    
    for (int i = 0; i < expected_bond_lengths.size(); i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_bond_lengths[i], angiotensin.graph.bonds[i].length, 0.00001);
    }
    
    for (int i = 0; i < expected_angles.size(); i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_angles[i]/57.29577951308232, angiotensin.graph.angles[i].theta, 0.00001);
    }
    
    for (int i = 0; i < expected_torsion_angles.size(); i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_torsion_angles[i]/57.29577951308232, angiotensin.graph.torsions[i].phi, 0.00001);
    }
}

void TestMolecule::testFlex()
{
    // TODO
}

void TestMolecule::testLocalTranslate()
{
    // TODO
}

void TestMolecule::testCrankshaft()
{
    // TODO
}

#endif // FLEXIBLE_LINKS

