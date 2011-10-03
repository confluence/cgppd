#include <UnitTest++.h>
#include <Molecule.h>

struct MoleculeFixture
{
    Molecule molecule;
//     char* pdbfilename;
    AminoAcids aminoAcidData;

    MoleculeFixture()
    {
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
}
