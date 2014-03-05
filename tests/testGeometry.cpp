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
    
    vector<segdata> segments;
    
    segdata seg1;
    seg1.residue_indices.push_back(72);
    seg1.residue_indices.push_back(73);
    seg1.residue_indices.push_back(74);
    seg1.residue_indices.push_back(75);
    seg1.residue_indices.push_back(123);

    segdata seg2;
    seg2.residue_indices.push_back(148);
    seg2.residue_indices.push_back(149);
    seg2.residue_indices.push_back(150);
    seg2.residue_indices.push_back(151);
    
    segments.push_back(seg1);
    segments.push_back(seg2);

    graph.init(residues, false, segments);
    
    set<Bond> expected_bonds;
    expected_bonds.insert(Bond(72, 73));
    expected_bonds.insert(Bond(73, 74));
    expected_bonds.insert(Bond(74, 75));
    expected_bonds.insert(Bond(75, 123));
    expected_bonds.insert(Bond(148, 149));
    expected_bonds.insert(Bond(149, 150));
    expected_bonds.insert(Bond(150, 151));
    
    set<Angle> expected_angles;
    expected_angles.insert(Angle(71, 72, 73));
    expected_angles.insert(Angle(72, 73, 74));
    expected_angles.insert(Angle(73, 74, 75));
    expected_angles.insert(Angle(74, 75, 123));
    expected_angles.insert(Angle(75, 123, 122));
    expected_angles.insert(Angle(75, 123, 124));
    expected_angles.insert(Angle(147, 148, 149));
    expected_angles.insert(Angle(148, 149, 150));
    expected_angles.insert(Angle(149, 150, 151));
    
    set<Torsion> expected_torsions;
    expected_torsions.insert(Torsion(70, 71, 72, 73));
    expected_torsions.insert(Torsion(71, 72, 73, 74));
    expected_torsions.insert(Torsion(72, 73, 74, 75));
    expected_torsions.insert(Torsion(73, 74, 75, 123));
    expected_torsions.insert(Torsion(74, 75, 123, 122));
    expected_torsions.insert(Torsion(74, 75, 123, 124));
    expected_torsions.insert(Torsion(75, 123, 124, 125));
    expected_torsions.insert(Torsion(121, 122, 123, 75));
    expected_torsions.insert(Torsion(146, 147, 148, 149));
    expected_torsions.insert(Torsion(147, 148, 149, 150));
    expected_torsions.insert(Torsion(148, 149, 150, 151));
    
    ASSERT_SET_EQUALS(expected_bonds, graph.bonds);
    ASSERT_SET_EQUALS(expected_angles, graph.angles);
    ASSERT_SET_EQUALS(expected_torsions, graph.torsions);
    
    set<int> expected_residues;
    expected_residues.insert(73);
    expected_residues.insert(74);
    expected_residues.insert(75);
    expected_residues.insert(149);
    expected_residues.insert(150);
    
    ASSERT_SET_EQUALS(expected_residues, graph.MC_crankshaft_residues);
    
    expected_residues.insert(151);
    
    ASSERT_SET_EQUALS(expected_residues, graph.MC_local_residues);
    
    expected_residues.insert(72);
    expected_residues.insert(123);
    expected_residues.insert(148);
    
    ASSERT_SET_EQUALS(expected_residues, graph.MC_flex_residues);
    
    set<int> expected_residues_branch_123_122;
    for (int i = 76; i < 124; i++) {
        expected_residues_branch_123_122.insert(i);
    }
    
    ASSERT_SET_EQUALS(expected_residues_branch_123_122, graph.branch(123, 122));

    set<int> expected_residues_branch_75_123;
    for (int i = 75; i < 152; i++) {
        expected_residues_branch_75_123.insert(i);
    }
    
    ASSERT_SET_EQUALS(expected_residues_branch_75_123, graph.branch(75, 123));
    
    set<int> expected_residues_branch_148_149;
    for (int i = 148; i < 152; i++) {
        expected_residues_branch_148_149.insert(i);
    }
    
    ASSERT_SET_EQUALS(expected_residues_branch_148_149, graph.branch(148, 149));

    // TODO TODO TODO add a test for cached bonds / angles / torsions for each residue
}
