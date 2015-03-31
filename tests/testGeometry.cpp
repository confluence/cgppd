//#include <string>
//#include <ctype.h>
#include <Geometry.h>
#include <testCommon.h>

#if FLEXIBLE_LINKS

TEST_CASE("Geometry", "[geometry]") {
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
    
    Graph graph;
    graph.init(residues, false, segments, 2);

    SECTION ("Bonds") {
        set<Bond> expected_bonds = { Bond(72, 73), Bond(73, 74), Bond(74, 75), Bond(75, 123), Bond(148, 149), Bond(149, 150), Bond(150, 151) };
        REQUIRE(graph.bonds == expected_bonds);
    }
    
    SECTION ("Angles") {
        set<Angle> expected_angles = { Angle(71, 72, 73), Angle(72, 73, 74), Angle(73, 74, 75), Angle(74, 75, 123), Angle(75, 123, 122), Angle(75, 123, 124), Angle(147, 148, 149), Angle(148, 149, 150), Angle(149, 150, 151) };
        REQUIRE(graph.angles == expected_angles);
    }
    
    SECTION ("Torsions") {
        set<Torsion> expected_torsions = { Torsion(70, 71, 72, 73), Torsion(71, 72, 73, 74), Torsion(72, 73, 74, 75), Torsion(73, 74, 75, 123), Torsion(74, 75, 123, 122), Torsion(74, 75, 123, 124), Torsion(75, 123, 124, 125), Torsion(121, 122, 123, 75), Torsion(146, 147, 148, 149), Torsion(147, 148, 149, 150), Torsion(148, 149, 150, 151) };
        REQUIRE(graph.torsions == expected_torsions);
    }
    
    SECTION ("Crankshaft residues") {
        vector<int> expected_crankshaft_residues = {73, 74, 75, 149, 150};
        REQUIRE(graph.MC_crankshaft_residues == expected_crankshaft_residues);
    }
    
    SECTION ("Local translation residues") {
        vector<int> expected_local_residues = {73, 74, 75, 149, 150, 151};
        REQUIRE(graph.MC_local_residues == expected_local_residues);
    }
    
    SECTION ("Flex residues") {
        vector<int> expected_flex_residues = {72, 73, 74, 75, 123, 148, 149, 150, 151};
        REQUIRE(graph.MC_flex_residues == expected_flex_residues);
    }

    SECTION("Flex branch") {
        REQUIRE(graph.branch(123, 122) == to_set(make_range(76, 123)));
        REQUIRE(graph.branch(75, 123) == to_set(make_range(75, 151)));
        REQUIRE(graph.branch(148, 149) == to_set(make_range(148, 151)));
    }

    SECTION("Segment bonds") {
        vector<Pair> expected_segment_bonds = {{75, 123}};
        REQUIRE(graph.segment_bonds == expected_segment_bonds);
    }

    SECTION("Indirect neighbours") {

        set<Pair> expected_indirect_neighbours = {{73, 123}, {74, 123}, {74, 122}, {74, 124}, {75, 122}, {75, 121}, {75, 124}, {75, 125}};
        REQUIRE(graph.indirect_neighbours == expected_indirect_neighbours);
    }

    SECTION("Rigid domains") {
        set<int> first_domain = to_set(make_range(0, 72));
        set<int> second_domain = to_set(make_range(76, 148));
        
        REQUIRE(graph.rigid_domains.size() == 2);
        REQUIRE(graph.rigid_domains[0] == first_domain);
        REQUIRE(graph.rigid_domains[1] == second_domain);

        SECTION("Rigid domain around") {
            REQUIRE(graph.rigid_domain_around(0) == first_domain);
            REQUIRE(graph.rigid_domain_around(50) == first_domain);
            REQUIRE(graph.rigid_domain_around(72) == first_domain);
            
            REQUIRE(graph.rigid_domain_around(76) == second_domain);
            REQUIRE(graph.rigid_domain_around(100) == second_domain);
            REQUIRE(graph.rigid_domain_around(123) == second_domain);
            REQUIRE(graph.rigid_domain_around(148) == second_domain);
        }
    }
}

#endif // FLEXIBLE_LINKS
