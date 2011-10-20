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

    molecule.AminoAcidsData = aminoAcidData;
    molecule.index = 0;
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