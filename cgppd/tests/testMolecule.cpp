#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>
#include <ctype.h>
#include <Molecule.h>

class TestMolecule : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestMolecule);
    CPPUNIT_TEST(testPDBSequence);
    CPPUNIT_TEST(testPDBCentre);
#if FLEXIBLE_LINKS
    CPPUNIT_TEST(testSegments);
    CPPUNIT_TEST(testGeometry);
#endif // FLEXIBLE_LINKS
    CPPUNIT_TEST(testStartingPotential);
    CPPUNIT_TEST_SUITE_END();

private:
    Molecule ubiquitin;
    Molecule angiotensin;
    AminoAcids aminoAcidData;
    float testboxdim;

public:
    void setUp();
    void testPDBSequence();
    void testPDBCentre();
#if FLEXIBLE_LINKS
    void testSegments();
    void testGeometry();
#endif // FLEXIBLE_LINKS
    void testStartingPotential();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestMolecule);

void TestMolecule::setUp()
{
    testboxdim = 118.4f;

    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);

    ubiquitin.init("tests/1UBQ.pdb", aminoAcidData, 0, testboxdim);
    angiotensin.init("tests/angiotensin.pdb", aminoAcidData, 0, testboxdim);
}

void TestMolecule::testPDBSequence()
{
    const int ubiquitin_size = 76;
    char ubiquitin_names[ubiquitin_size][4] = {"MET", "GLN", "ILE", "PHE", "VAL", "LYS", "THR", "LEU", "THR", "GLY", "LYS", "THR", "ILE", "THR", "LEU", "GLU", "VAL", "GLU", "PRO", "SER", "ASP", "THR", "ILE", "GLU", "ASN", "VAL", "LYS", "ALA", "LYS", "ILE", "GLN", "ASP", "LYS", "GLU", "GLY", "ILE", "PRO", "PRO", "ASP", "GLN", "GLN", "ARG", "LEU", "ILE", "PHE", "ALA", "GLY", "LYS", "GLN", "LEU", "GLU", "ASP", "GLY", "ARG", "THR", "LEU", "SER", "ASP", "TYR", "ASN", "ILE", "GLN", "LYS", "GLU", "SER", "THR", "LEU", "HIS", "LEU", "VAL", "LEU", "ARG", "LEU", "ARG", "GLY", "GLY"};

    const int angiotensin_size = 10;
    char angiotensin_names[angiotensin_size][4] = {"ASP","ARG","VAL","TYR","ILE","HIS","PRO","PHE","HIS","LEU"};

    CPPUNIT_ASSERT_EQUAL(ubiquitin_size, ubiquitin.residueCount);
    CPPUNIT_ASSERT_EQUAL(angiotensin_size, angiotensin.residueCount);
#if FLEXIBLE_LINKS
    CPPUNIT_ASSERT_EQUAL(angiotensin_size - 1, angiotensin.linkCount);
#endif // FLEXIBLE_LINKS

    int expected;
    int got;

    // TODO: rewrite these so they compare strings instead
    for (int i = 0; i < ubiquitin_size; i++)
    {
        expected = aminoAcidData.getAminoAcidIndex(ubiquitin_names[i]);
        got = ubiquitin.Residues[i].aminoAcidIndex;
        CPPUNIT_ASSERT_EQUAL(expected, got);
    }

    for (int i = 0; i < angiotensin_size; i++)
    {
        expected = aminoAcidData.getAminoAcidIndex(angiotensin_names[i]);
        got = angiotensin.Residues[i].aminoAcidIndex;
        CPPUNIT_ASSERT_EQUAL(expected, got);
    }
}

void TestMolecule::testPDBCentre()
{
    static const float expected_centre[3] = {30.442894736842103, 29.00898684210527, 15.5163947368421};

    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_centre[0], ubiquitin.center.x, 0.00001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_centre[1], ubiquitin.center.y, 0.00001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_centre[2], ubiquitin.center.z, 0.00001);
}

#if FLEXIBLE_LINKS
void TestMolecule::testSegments()
{
    CPPUNIT_ASSERT_EQUAL(1, angiotensin.segmentCount);
//     CPPUNIT_ASSERT_EQUAL(1, angiotensin.linkerCount);

    CPPUNIT_ASSERT_EQUAL(0, angiotensin.Segments[0].start);
    CPPUNIT_ASSERT_EQUAL(9, angiotensin.Segments[0].end);
    CPPUNIT_ASSERT_EQUAL(true, angiotensin.Segments[0].flexible);

//     CPPUNIT_ASSERT_EQUAL(0, angiotensin.Linkers[0]->start);
//     CPPUNIT_ASSERT_EQUAL(9, angiotensin.Linkers[0]->end);
}

void TestMolecule::testGeometry()
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

    Potential e = angiotensin.E();

    for (size_t i = 0; i < angiotensin.linkCount; i++)
    {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_bond_lengths[i], angiotensin.Links[i].pseudo_bond, 0.00001);

        // convert degrees (from VMD) to radians
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_angles[i]/57.29577951308232, angiotensin.Residues[i].pseudo_angle, 0.00001);

        // convert degrees (from VMD) to radians
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_torsion_angles[i]/57.29577951308232, angiotensin.Links[i].pseudo_torsion, 0.00001);
    }
}
#endif // FLEXIBLE_LINKS

void TestMolecule::testStartingPotential()
{
}
