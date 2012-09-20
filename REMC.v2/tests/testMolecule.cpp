#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>
#include <ctype.h>
#include <Molecule.h>

class TestMolecule : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestMolecule);
    CPPUNIT_TEST(testPDBSequence);
    CPPUNIT_TEST(testPDBCentre);
    CPPUNIT_TEST_SUITE_END();

private:
    Molecule molecule;
    AminoAcids aminoAcidData;

public:
    void setUp();
    void testPDBSequence();
    void testPDBCentre();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestMolecule);

void TestMolecule::setUp()
{
    aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
    aminoAcidData.loadLJPotentialData(LJPDSOURCE);

    molecule.init_amino_acid_data(aminoAcidData);
    molecule.index = 0;
    float testboxdim(118.4f);
    molecule.bounding_value = testboxdim;
    molecule.initFromPDB("tests/1UBQ.pdb");
}

void TestMolecule::testPDBSequence()
{
    const int expected_size = 76;
    CPPUNIT_ASSERT_EQUAL(expected_size, molecule.residueCount);

    char expected_sequence_names[expected_size][4] = {"MET", "GLN", "ILE", "PHE", "VAL", "LYS", "THR", "LEU", "THR", "GLY", "LYS", "THR", "ILE", "THR", "LEU", "GLU", "VAL", "GLU", "PRO", "SER", "ASP", "THR", "ILE", "GLU", "ASN", "VAL", "LYS", "ALA", "LYS", "ILE", "GLN", "ASP", "LYS", "GLU", "GLY", "ILE", "PRO", "PRO", "ASP", "GLN", "GLN", "ARG", "LEU", "ILE", "PHE", "ALA", "GLY", "LYS", "GLN", "LEU", "GLU", "ASP", "GLY", "ARG", "THR", "LEU", "SER", "ASP", "TYR", "ASN", "ILE", "GLN", "LYS", "GLU", "SER", "THR", "LEU", "HIS", "LEU", "VAL", "LEU", "ARG", "LEU", "ARG", "GLY", "GLY"};
    int expected;
    int got;

    // TODO: rewrite this so it compares strings instead
    for (int i = 0; i < expected_size; i++)
    {
        expected = aminoAcidData.getAminoAcidIndex(expected_sequence_names[i]);
        got = molecule.Residues[i].aminoAcidIndex;
        CPPUNIT_ASSERT_EQUAL(expected, got);
    }
}

void TestMolecule::testPDBCentre()
{
    static const float expected_centre[3] = {30.442894736842103, 29.00898684210527, 15.5163947368421};

    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_centre[0], molecule.center.x, 0.00001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_centre[1], molecule.center.y, 0.00001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_centre[2], molecule.center.z, 0.00001);
}

#if FLEXIBLE_LINKS
class TestLinkerPotentials : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestLinkerPotentials);
    CPPUNIT_TEST(testPDBSequence);
    CPPUNIT_TEST(testSegments);
    CPPUNIT_TEST(testGeometry);
    CPPUNIT_TEST_SUITE_END();

private:
    Molecule molecule;
    AminoAcids aminoAcidData;

public:
    void setUp();
    void testPDBSequence();
    void testSegments();
    void testGeometry();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestLinkerPotentials);

void TestLinkerPotentials::setUp()
{
    aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
    aminoAcidData.loadLJPotentialData(LJPDSOURCE);

    molecule.init_amino_acid_data(aminoAcidData);
    molecule.index = 0;
    float testboxdim(118.4f);
    molecule.bounding_value = testboxdim;
    molecule.initFromPDB("tests/angiotensin.pdb");
}

void TestLinkerPotentials::testPDBSequence()
{
    const int expected_size = 10;
    CPPUNIT_ASSERT_EQUAL(expected_size, molecule.residueCount);
    CPPUNIT_ASSERT_EQUAL(expected_size - 1, molecule.linkCount);

    char expected_sequence_names[expected_size][4] = {"ASP","ARG","VAL","TYR","ILE","HIS","PRO","PHE","HIS","LEU"};
    int expected;
    int got;

    // TODO: rewrite this so it compares strings instead
    for (int i = 0; i < expected_size; i++)
    {
        expected = aminoAcidData.getAminoAcidIndex(expected_sequence_names[i]);
        got = molecule.Residues[i].aminoAcidIndex;
        CPPUNIT_ASSERT_EQUAL(expected, got);
    }
}

void TestLinkerPotentials::testSegments()
{
    CPPUNIT_ASSERT_EQUAL(1, molecule.segmentCount);
//     CPPUNIT_ASSERT_EQUAL(1, molecule.linkerCount);

    CPPUNIT_ASSERT_EQUAL(0, molecule.Segments[0].start);
    CPPUNIT_ASSERT_EQUAL(9, molecule.Segments[0].end);
    CPPUNIT_ASSERT_EQUAL(true, molecule.Segments[0].flexible);

//     CPPUNIT_ASSERT_EQUAL(0, molecule.Linkers[0]->start);
//     CPPUNIT_ASSERT_EQUAL(9, molecule.Linkers[0]->end);
}

void TestLinkerPotentials::testGeometry()
{
    float expected_bond_lengths[9] = {
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

    float expected_angles[9] = {
        0.0, // padding for simplicity
        123.56746673583984,
        86.7016830444336,
        126.89990234375,
        126.8516845703125,
        129.04217529296875,
        84.29325103759766,
        130.4828338623047,
        97.1257095336914,
    };

    float expected_torsion_angles[9] = {
        0.0, // padding for simplicity
        -165.99693298339844,
        26.120737075805664,
        -171.85838317871094,
        -170.0393524169922,
        -156.48416137695313,
        4.634008884429932,
        150.31800842285156,
        0.0 // padding for simplicity
    };

    Potential e = molecule.E();
//     cout << e.total_LJ() << endl;
//     cout << e.total_DH() << endl;
//     cout << e.total_bond() << endl;
//     cout << e.total_angle() << endl;
//     cout << e.total_torsion() << endl;

    for (size_t i = 0; i < molecule.linkCount; i++)
    {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_bond_lengths[i], molecule.Links[i].pseudo_bond, 0.00001);

        // convert degrees (from VMD) to radians
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_angles[i]/57.29577951308232, molecule.Residues[i].pseudo_angle, 0.00001);

        // convert degrees (from VMD) to radians
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_torsion_angles[i]/57.29577951308232, molecule.Links[i].pseudo_torsion, 0.00001);
    }
}
#endif // FLEXIBLE_LINKS