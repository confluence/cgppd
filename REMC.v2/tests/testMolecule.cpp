#include <UnitTest++.h>
#include <ctype.h>
#include <Molecule.h>

struct MoleculeFixture
{
    Molecule molecule;
    AminoAcids aminoAcidData;

    MoleculeFixture()
    {
        aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
        aminoAcidData.loadLJPotentialData(LJPDSOURCE);

        molecule.AminoAcidsData = aminoAcidData;
        molecule.index = 0;
        molecule.initFromPDB("tests/1UBQ.pdb");
    }

    ~MoleculeFixture()
    {
    }
};

TEST_FIXTURE(MoleculeFixture, TestInitFromPDB) {
    const int expected_size = 76;
    CHECK_EQUAL(expected_size, molecule.residueCount);

    char expected_sequence_names[expected_size][4] = {"MET", "GLN", "ILE", "PHE", "VAL", "LYS", "THR", "LEU", "THR", "GLY", "LYS", "THR", "ILE", "THR", "LEU", "GLU", "VAL", "GLU", "PRO", "SER", "ASP", "THR", "ILE", "GLU", "ASN", "VAL", "LYS", "ALA", "LYS", "ILE", "GLN", "ASP", "LYS", "GLU", "GLY", "ILE", "PRO", "PRO", "ASP", "GLN", "GLN", "ARG", "LEU", "ILE", "PHE", "ALA", "GLY", "LYS", "GLN", "LEU", "GLU", "ASP", "GLY", "ARG", "THR", "LEU", "SER", "ASP", "TYR", "ASN", "ILE", "GLN", "LYS", "GLU", "SER", "THR", "LEU", "HIS", "LEU", "VAL", "LEU", "ARG", "LEU", "ARG", "GLY", "GLY"};
    int expected_sequence[expected_size];
    int got_sequence[expected_size];

    for (int i = 0; i < expected_size; i++)
    {
        expected_sequence[i] = aminoAcidData.getAminoAcidIndex(expected_sequence_names[i]);
        got_sequence[i] = molecule.Residues[i].aminoAcidIndex;
    }

    CHECK_ARRAY_EQUAL(expected_sequence, got_sequence, expected_size);

    static const float expected_centre[3] = {30.442894736842103, 29.00898684210527, 15.5163947368421};
    float got_centre[3] = {molecule.center.x, molecule.center.y, molecule.center.z};

    CHECK_ARRAY_CLOSE(expected_centre, got_centre, 3, 0.00001);
}
