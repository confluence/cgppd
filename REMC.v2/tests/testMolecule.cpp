#include <UnitTest++.h>
#include <Molecule.h>

struct MoleculeFixture
{
    Molecule molecule;
//     char* pdbfilename;
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

TEST_FIXTURE(MoleculeFixture, TestMolecule) {
    CHECK_EQUAL(76, molecule.residueCount);
    cout << aminoAcidData.data[molecule.Residues[0].aminoAcidIndex].sname << endl;
    cout << aminoAcidData.data[molecule.Residues[1].aminoAcidIndex].sname << endl;
    cout << aminoAcidData.data[molecule.Residues[2].aminoAcidIndex].sname << endl;

    char expected_sequence[76][4] = {"MET", "GLN", "ILE", "PHE", "VAL", "LYS", "THR", "LEU", "THR", "GLY", "LYS", "THR", "ILE", "THR", "LEU", "GLU", "VAL", "GLU", "PRO", "SER", "ASP", "THR", "ILE", "GLU", "ASN", "VAL", "LYS", "ALA", "LYS", "ILE", "GLN", "ASP", "LYS", "GLU", "GLY", "ILE", "PRO", "PRO", "ASP", "GLN", "GLN", "ARG", "LEU", "ILE", "PHE", "ALA", "GLY", "LYS", "GLN", "LEU", "GLU", "ASP", "GLY", "ARG", "THR", "LEU", "SER", "ASP", "TYR", "ASN", "ILE", "GLN", "LYS", "GLU", "SER", "THR", "LEU", "HIS", "LEU", "VAL", "LEU", "ARG", "LEU", "ARG", "GLY", "GLY"}
}
