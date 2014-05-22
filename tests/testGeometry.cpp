#include <string>
#include <ctype.h>
#include <Geometry.h>
#include <testCommon.h>

class TestGeometry : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestGeometry);
    CPPUNIT_TEST(testGeometricFeatures);
    CPPUNIT_TEST(testMCResidues);
    CPPUNIT_TEST(testFlexBranch);
    CPPUNIT_TEST(testFeaturesForPotential);
    CPPUNIT_TEST_SUITE_END();

private:
    Graph graph;

public:
    void setUp();
    void testGeometricFeatures();
    void testMCResidues();
    void testFlexBranch();
    void testFeaturesForPotential();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestGeometry);

void TestGeometry::setUp()
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

    graph.init(residues, false, segments, 2);
}

void TestGeometry::testGeometricFeatures()
{
    set<Bond> expected_bonds = { Bond(72, 73), Bond(73, 74), Bond(74, 75), Bond(75, 123), Bond(148, 149), Bond(149, 150), Bond(150, 151) };

    set<Angle> expected_angles = { Angle(71, 72, 73), Angle(72, 73, 74), Angle(73, 74, 75), Angle(74, 75, 123), Angle(75, 123, 122), Angle(75, 123, 124), Angle(147, 148, 149), Angle(148, 149, 150), Angle(149, 150, 151) };

    set<Torsion> expected_torsions = { Torsion(70, 71, 72, 73), Torsion(71, 72, 73, 74), Torsion(72, 73, 74, 75), Torsion(73, 74, 75, 123), Torsion(74, 75, 123, 122), Torsion(74, 75, 123, 124), Torsion(75, 123, 124, 125), Torsion(121, 122, 123, 75), Torsion(146, 147, 148, 149), Torsion(147, 148, 149, 150), Torsion(148, 149, 150, 151) };

    ASSERT_ITERABLE_EQUALS(expected_bonds, graph.bonds);
    ASSERT_ITERABLE_EQUALS(expected_angles, graph.angles);
    ASSERT_ITERABLE_EQUALS(expected_torsions, graph.torsions);
}

void TestGeometry::testMCResidues()
{
    vector<int> expected_crankshaft_residues = {73, 74, 75, 149, 150};
    vector<int> expected_local_residues = {73, 74, 75, 149, 150, 151};
    vector<int> expected_flex_residues = {72, 73, 74, 75, 123, 148, 149, 150, 151};

    ASSERT_ITERABLE_EQUALS(expected_crankshaft_residues, graph.MC_crankshaft_residues);
    ASSERT_ITERABLE_EQUALS(expected_local_residues, graph.MC_local_residues);
    ASSERT_ITERABLE_EQUALS(expected_flex_residues, graph.MC_flex_residues);
}

void TestGeometry::testFlexBranch()
{
    set<int> expected_residues_branch_123_122 = to_set(make_range(76, 123));
    set<int> expected_residues_branch_75_123 = to_set(make_range(75, 151));
    set<int> expected_residues_branch_148_149 = to_set(make_range(148, 151));
    
    ASSERT_ITERABLE_EQUALS(expected_residues_branch_123_122, graph.branch(123, 122));
    ASSERT_ITERABLE_EQUALS(expected_residues_branch_75_123, graph.branch(75, 123));
    ASSERT_ITERABLE_EQUALS(expected_residues_branch_148_149, graph.branch(148, 149));
}

void TestGeometry::testFeaturesForPotential()
{
    // test segment bonds
    
    vector<Pair> expected_segment_bonds = {{75, 123}};
    ASSERT_ITERABLE_EQUALS(expected_segment_bonds, graph.segment_bonds);
    
    // test indirect neighbours
    
    set<Pair> expected_indirect_neighbours = {{73, 123}, {74, 123}, {74, 122}, {74, 124}, {75, 122}, {75, 121}, {75, 124}, {75, 125}};
    ASSERT_ITERABLE_EQUALS(expected_indirect_neighbours, graph.indirect_neighbours);
    
    // test rigid_domain_around
    
    set<int> first_domain = to_set(make_range(0, 72));
    set<int> second_domain = to_set(make_range(76, 148));
    
    ASSERT_ITERABLE_EQUALS(first_domain, graph.rigid_domain_around(0));
    ASSERT_ITERABLE_EQUALS(first_domain, graph.rigid_domain_around(50));
    ASSERT_ITERABLE_EQUALS(first_domain, graph.rigid_domain_around(72));
    
    ASSERT_ITERABLE_EQUALS(second_domain, graph.rigid_domain_around(76));
    ASSERT_ITERABLE_EQUALS(second_domain, graph.rigid_domain_around(100));
    ASSERT_ITERABLE_EQUALS(second_domain, graph.rigid_domain_around(123));
    ASSERT_ITERABLE_EQUALS(second_domain, graph.rigid_domain_around(148));
    
    // test rigid domains
    
    CPPUNIT_ASSERT_EQUAL(2, (int)graph.rigid_domains.size());
    ASSERT_ITERABLE_EQUALS(first_domain, graph.rigid_domains[0]);
    ASSERT_ITERABLE_EQUALS(second_domain, graph.rigid_domains[1]);
}
