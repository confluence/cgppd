//#include <string>
//#include <ctype.h>
//#include <structures.h>
#include <Molecule.h>
#include <testCommon.h>

TEST_CASE("Molecule", "[molecule]") {
    AminoAcids aminoAcidData;
    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);

    moldata ubq_data;
    Molecule ubiquitin;
    strcpy(ubq_data.pdbfilename, "tests/1UBQ.pdb");
    segdata seg;
    seg.residue_indices.push_back(72);
    seg.residue_indices.push_back(73);
    seg.residue_indices.push_back(74);
    seg.residue_indices.push_back(75);
    ubq_data.segments.push_back(seg);
    ubiquitin.init(ubq_data, aminoAcidData, 0, 1000.0f);

    moldata ang_data;
    Molecule angiotensin;
    strcpy(ang_data.pdbfilename, "tests/angiotensin.pdb");
    ang_data.all_flexible = true;
    angiotensin.init(ang_data, aminoAcidData, 0, 1000.0f);
    
    moldata ala8_data;
    Molecule polyalanine8;
    strcpy(ala8_data.pdbfilename, "tests/ala8.pdb");
    ala8_data.all_flexible = true;
    // small bounding value, so we can test wrapping over the boundary
    polyalanine8.init(ala8_data, aminoAcidData, 0, 15.0f);

    SECTION("Residue sequence") {
        vector<string> expected_ubiquitin_sequence = {"MET", "GLN", "ILE", "PHE", "VAL", "LYS", "THR", "LEU", "THR", "GLY", "LYS", "THR", "ILE", "THR", "LEU", "GLU", "VAL", "GLU", "PRO", "SER", "ASP", "THR", "ILE", "GLU", "ASN", "VAL", "LYS", "ALA", "LYS", "ILE", "GLN", "ASP", "LYS", "GLU", "GLY", "ILE", "PRO", "PRO", "ASP", "GLN", "GLN", "ARG", "LEU", "ILE", "PHE", "ALA", "GLY", "LYS", "GLN", "LEU", "GLU", "ASP", "GLY", "ARG", "THR", "LEU", "SER", "ASP", "TYR", "ASN", "ILE", "GLN", "LYS", "GLU", "SER", "THR", "LEU", "HIS", "LEU", "VAL", "LEU", "ARG", "LEU", "ARG", "GLY", "GLY"};
        REQUIRE(ubiquitin.residue_sequence() == expected_ubiquitin_sequence);

        vector<string> expected_angiotensin_sequence = {"ASP","ARG","VAL","TYR","ILE","HIS","PRO","PHE","HIS","LEU"};
        REQUIRE(angiotensin.residue_sequence() == expected_angiotensin_sequence);

    }

    SECTION("Position") {

        vector<Vector3f> ala8_start = {
            Vector3f(-9.29663, -6.96, -0.91625),
            Vector3f(-7.38663, -3.879, 0.23775),
            Vector3f(-3.76163, -3.302, -0.76025),
            Vector3f(-1.79963, -0.311, 0.53375),
            Vector3f(1.76937, 0.361, -0.59825),
            Vector3f(3.79137, 3.252, 0.82475),
            Vector3f(7.29737, 4.029, -0.43125),
            Vector3f(9.38637, 6.81, 1.10975)
        }; 

        SECTION("Starting positions") {
            // We set the position to zero by default
            REQUIRE(polyalanine8.center.almost_equal(Vector3f(0, 0, 0)));
            
            // TODO: make a test for using a translation by 0 rather than setting the position
            
            for (int i = 0; i < 8; i++) {
                REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(ala8_start[i]));
            }

            // These are the same because it's centered on zero
            for (int i = 0; i < 8; i++) {
                REQUIRE(polyalanine8.Residues[i].position.almost_equal(ala8_start[i]));
            }
        }

        SECTION("Rotate") {
            // Rotate 180 degrees about the x axis
            Vector3f rotation_vector(1, 0, 0);
            polyalanine8.rotate(rotation_vector, M_PIl);

            // Check that absolute positions have changed -- x coordinate should be the same; y and z should be flipped
            for (int i = 0; i < polyalanine8.residueCount; i++) {
                Vector3f expected_position = Vector3f(ala8_start[i].x, -ala8_start[i].y, -ala8_start[i].z);
                REQUIRE(polyalanine8.Residues[i].position.almost_equal(expected_position));
            }
            
            // Check that relative positions have changed -- should be the same
            for (int i = 0; i < polyalanine8.residueCount; i++) {
                Vector3f expected_position = Vector3f(ala8_start[i].x, -ala8_start[i].y, -ala8_start[i].z);
                REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(expected_position));
            }

            // Check that centre has not changed
            REQUIRE(polyalanine8.center.almost_equal(Vector3f(0, 0, 0)));
        }

        SECTION("Translate") {
            Vector3f translation_vector(1, 1, 1);
            polyalanine8.translate(translation_vector);

            // Check that absolute positions have changed

            for (int i = 0; i < 8; i++) {
                REQUIRE(polyalanine8.Residues[i].position.almost_equal(ala8_start[i] + translation_vector));
            }

            // Check that centre has changed

            REQUIRE(polyalanine8.center.almost_equal(Vector3f(1, 1, 1)));

            // Check that relative positions have not changed
            for (int i = 0; i < 8; i++) {
                REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(ala8_start[i]));
            }
        }

        SECTION ("Wrap") {
            polyalanine8.translate(Vector3f(16, 0, 0));
    
            // check that we wrapped right around
            REQUIRE(polyalanine8.center.almost_equal(Vector3f(1, 0, 0)));
            
            // check that absolute positions are wrapped
            for (int i = 0; i < 8; i++) {
                REQUIRE(polyalanine8.Residues[i].position.almost_equal(ala8_start[i] + Vector3f(1, 0, 0)));
            }

            // Check that relative positions have not changed
            for (int i = 0; i < 8; i++) {
                REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(ala8_start[i]));
            }
        }

        
#if FLEXIBLE_LINKS

        SECTION("Flex") {
            // Rotate branch starting with edge 4->5 180 degrees about axis passing through residue 4 and parallel to the x axis
            Vector3f rotation_vector(1, 0, 0);
            polyalanine8.flex(rotation_vector, M_PIl, 3, 4);
            
            Vector3f rotation_centre = polyalanine8.Residues[3].position;
            
            for (int i = 0; i < 4; i++) {
                REQUIRE(polyalanine8.Residues[i].position.almost_equal(ala8_start[i]));
            }

            Vector3f v;
            
            // rotated residues should have y and z flipped relative to residue 4
            for (int i = 4; i < 8; i++) {
                REQUIRE(polyalanine8.Residues[i].position.almost_equal(Vector3f(ala8_start[i].x, - ala8_start[i].y + 2 * (rotation_centre.y), - ala8_start[i].z + 2 * (rotation_centre.z))));
                v += polyalanine8.Residues[i].position - ala8_start[i];
            }
            
            // check that centre is updated
            
            REQUIRE(polyalanine8.center.almost_equal(v/8));
            
            // check that all the relative positions are correct
            
            for (int i = 0; i < 8; i++) {
                REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(polyalanine8.Residues[i].position - polyalanine8.center));
            }
        }

        SECTION("Local translate") {
            // Translate third residue

            Vector3f translation_vector(1, 1, 1);
            polyalanine8.translate(translation_vector, 2);
            
            // Check that absolute positions have changed

            for (int i = 0; i < 8; i++) {
                if (i == 2) {
                    REQUIRE(polyalanine8.Residues[i].position.almost_equal(ala8_start[i] + translation_vector));
                } else {
                    REQUIRE(polyalanine8.Residues[i].position.almost_equal(ala8_start[i]));
                }
            }

            // Check that centre has changed

            REQUIRE(polyalanine8.center.almost_equal(translation_vector/8.0));

            // Check that relative positions have changed
            for (int i = 0; i < 8; i++) {
                if (i == 2) {
                    REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(ala8_start[i] + translation_vector - translation_vector/8.0));
                } else {
                    REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(ala8_start[i] - translation_vector/8.0));
                }
            }
        }

        SECTION("Crankshaft") {
            Vector3f old_i = polyalanine8.Residues[1].position;
            Vector3f old_j = polyalanine8.Residues[2].position;
            Vector3f old_k = polyalanine8.Residues[3].position;
            
            Vector3f old_ij = old_j - old_i;
            Vector3f old_jk = old_k - old_j;

            polyalanine8.crankshaft(M_PIl, false, 2);

            Vector3f new_i = polyalanine8.Residues[1].position;
            Vector3f new_j = polyalanine8.Residues[2].position;
            Vector3f new_k = polyalanine8.Residues[3].position;
            
            Vector3f new_ij = new_j - new_i;
            Vector3f new_jk = new_k - new_j;
            
            // check that all absolute positions are the same except for j
            
            Vector3f translation_vector = new_j - old_j;
            
            REQUIRE(translation_vector.magnitude());
            
            for (int i = 0; i < 8; i++) {
                if (i == 2) {
                    REQUIRE(polyalanine8.Residues[i].position.almost_equal(ala8_start[i] + translation_vector));
                } else {
                    REQUIRE(polyalanine8.Residues[i].position.almost_equal(ala8_start[i]));
                }
            }
            
            // Check that centre has changed
            REQUIRE(polyalanine8.center.almost_equal(translation_vector/8.0));

            // Check that relative positions have changed
            for (int i = 0; i < 8; i++) {
                if (i == 2) {
                    REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(ala8_start[i] + translation_vector - translation_vector/8.0));
                } else {
                    REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(ala8_start[i] - translation_vector/8.0));
                }
            }
            
            // check that ij and jk are still the same length

            REQUIRE(old_ij.magnitude() == Approx(new_ij.magnitude()).epsilon(0.00001));
            REQUIRE(old_jk.magnitude() == Approx(new_jk.magnitude()).epsilon(0.00001));
            
            // check that j has actually been flipped by 180 degrees:
            // since old j != new j, this must be true if i, j, k and new j are all in the same plane
            // which is true if the volume of the tetrahedron they form is zero
            
            Vector3f & a = new_i;
            Vector3f & b = old_j;
            Vector3f & c = new_k;
            Vector3f & d = new_j;
            
            // omit absolute value and division by 6
            REQUIRE((a - d).dot((b - d).cross(c - d)) == Approx(0.0).epsilon(0.00001));
        }

        SECTION("Local translate wrap") {
            // test wrapping from a local move, because that uses a helper function which is slightly different
    
            // translate whole molecule so that centre is right next to the boundary
            
            polyalanine8.translate(Vector3f(14.99, 0, 0));

            // translate one residue so far over the boundary that the centre has to wrap
            polyalanine8.translate(Vector3f(8, 0, 0), 7);
            
            // check that the centre has wrapped
            REQUIRE(polyalanine8.center.almost_equal(Vector3f(14.99 + 1 - 15, 0, 0)));

            // check that absolute positions have wrapped
            
            for (int i = 0; i < 7; i++) {
                REQUIRE(polyalanine8.Residues[i].position.almost_equal(Vector3f(ala8_start[i].x + 14.99 - 15, ala8_start[i].y, ala8_start[i].z)));
            }
            
            REQUIRE(polyalanine8.Residues[7].position.almost_equal(Vector3f(ala8_start[7].x + 14.99 + 8 - 15, ala8_start[7].y, ala8_start[7].z)));
            
            // check that relative positions have been updated, and reflect only a move of the last residue by (8,0,0)
            
            for (int i = 0; i < 7; i++) {
                REQUIRE(polyalanine8.Residues[i].relativePosition.almost_equal(Vector3f(ala8_start[i].x - 1, ala8_start[i].y, ala8_start[i].z)));
            }
            
            REQUIRE(polyalanine8.Residues[7].relativePosition.almost_equal(Vector3f(ala8_start[7].x - 1 + 8, ala8_start[7].y, ala8_start[7].z)));
        }

#endif // FLEXIBLE_LINKS

    }

    SECTION("Potential") {
        Potential actual_ala = polyalanine8.E(true); // with internal LJ and DH
        Potential actual_ang = angiotensin.E(true); // with internal LJ and DH
#if !LJ_REPULSIVE && !LJ_OFF
        double ala_LJ = 0.932791;
        double ang_LJ = -0.250504;
#elif LJ_REPULSIVE
        double ala_LJ = 5.638410;
        double ang_LJ = -0.628958;
#elif LJ_OFF
        double ala_LJ = 0;
        double ang_LJ = 0;
#endif // !LJ_REPULSIVE && !LJ_OFF
        Potential expected_ala = Potential(ala_LJ, 0.000000, 0.045626, 1.272977, 10.872785);
        Potential expected_ang = Potential(ang_LJ, 0.052706, 0.642775, 4.887787, 15.871157);

        REQUIRE(actual_ala.almost_equal(expected_ala));
        REQUIRE(actual_ang.almost_equal(expected_ang));
    }

#if FLEXIBLE_LINKS

    SECTION("Segments") {
        REQUIRE(ubiquitin.graph.bonds.size() == 3);
        REQUIRE(ubiquitin.graph.angles.size() == 3);
        REQUIRE(ubiquitin.graph.torsions.size() == 3);

        REQUIRE(polyalanine8.graph.bonds.size() == 7);
        REQUIRE(polyalanine8.graph.angles.size() == 6);
        REQUIRE(polyalanine8.graph.torsions.size() == 5);

        REQUIRE(angiotensin.graph.bonds.size() == 9);
        REQUIRE(angiotensin.graph.angles.size() == 8);
        REQUIRE(angiotensin.graph.torsions.size() == 7);
    }

#endif // FLEXIBLE_LINKS
}
