#include <string>
#include <ctype.h>
#include <Geometry.h>
#include <testCommon.h>

class TestGeometry : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestGeometry);
    CPPUNIT_TEST(testInit);
    CPPUNIT_TEST_SUITE_END();

private:
    Graph graph;

public:
    void testInit();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestGeometry);

void TestGeometry::testInit()
{
    vector<Residue> residues; // we're going to fill this with dummy residues for one of the diubiquitins, with just a label and a chain ID
    
    char chains[2] = {'A', 'B'};
    
    for (int c = 0; c < 2; c++) { // chain
        for (int r = 0; r < 76; r++) {
            Residue res;
            res.resSeq = 1 + c + r;
            res.chainId = chains[c];
            residues.push_back(res);
       }
    }
    
    segdata seg1;
    seg1.residue_indices = {72, 73, 74, 75, 123};

    segdata seg2;
    seg2.residue_indices = {148, 149, 150, 151};
    
    vector<segdata> segments = {seg1, seg2};

    graph.init(residues, false, segments);
    
    //--------------------------------------------------------------------------
    
    vector<Bond> expected_bonds = { Bond(72, 73), Bond(73, 74), Bond(74, 75), Bond(75, 123), Bond(148, 149), Bond(149, 150), Bond(150, 151) };

    vector<Angle> expected_angles = { Angle(71, 72, 73), Angle(72, 73, 74), Angle(73, 74, 75), Angle(74, 75, 123), Angle(75, 123, 122), Angle(75, 123, 124), Angle(147, 148, 149), Angle(148, 149, 150), Angle(149, 150, 151) };

    // TODO TODO TODO this is wrong; 75-123-122-121 is missing!
    vector<Torsion> expected_torsions = { Torsion(70, 71, 72, 73), Torsion(71, 72, 73, 74), Torsion(72, 73, 74, 75), Torsion(73, 74, 75, 123), Torsion(74, 75, 123, 122), Torsion(74, 75, 123, 124), Torsion(75, 123, 124, 125), Torsion(121, 122, 123, 75), Torsion(146, 147, 148, 149), Torsion(147, 148, 149, 150), Torsion(148, 149, 150, 151) };

    ASSERT_ITERABLE_EQUALS(expected_bonds, graph.bonds);
    ASSERT_ITERABLE_EQUALS(expected_angles, graph.angles);
    ASSERT_ITERABLE_EQUALS(expected_torsions, graph.torsions);
    
    //--------------------------------------------------------------------------
    
    vector<int> expected_crankshaft_residues = {73, 74, 75, 149, 150};
    vector<int> expected_local_residues = {73, 74, 75, 149, 150, 151};
    vector<int> expected_flex_residues = {72, 73, 74, 75, 123, 148, 149, 150, 151};

    ASSERT_ITERABLE_EQUALS(expected_crankshaft_residues, graph.MC_crankshaft_residues);
    ASSERT_ITERABLE_EQUALS(expected_local_residues, graph.MC_local_residues);
    ASSERT_ITERABLE_EQUALS(expected_flex_residues, graph.MC_flex_residues);
    
    //--------------------------------------------------------------------------

    set<int> expected_residues_branch_123_122 = to_set(make_range(76, 123));
    set<int> expected_residues_branch_75_123 = to_set(make_range(75, 151));
    set<int> expected_residues_branch_148_149 = to_set(make_range(148, 151));
    
    ASSERT_ITERABLE_EQUALS(expected_residues_branch_123_122, graph.branch(123, 122));
    ASSERT_ITERABLE_EQUALS(expected_residues_branch_75_123, graph.branch(75, 123));
    ASSERT_ITERABLE_EQUALS(expected_residues_branch_148_149, graph.branch(148, 149));

    //--------------------------------------------------------------------------
    
    set<int> expected_bonds_72({0});
    set<int> expected_bonds_73 = {0, 1};
    
    set<int> expected_angles_71 = {0};
    set<int> expected_angles_73 = {0, 1, 2};
    set<int> expected_angles_123 = {3, 4, 5};
    
    set<int> expected_torsions_70 = {0};
    set<int> expected_torsions_72 = {0, 1, 10};
    set<int> expected_torsions_73 = {0, 1, 2, 3};
//     set<int> expected_torsions_123 = { Torsion(73, 74, 75, 123), Torsion(74, 75, 123, 122), Torsion(74, 75, 123, 124), Torsion(75, 123, 122, 121), Torsion(75, 123, 124, 125) };
    
    ASSERT_ITERABLE_EQUALS(expected_bonds_72, graph.bonds_for_residue[72]);
}
