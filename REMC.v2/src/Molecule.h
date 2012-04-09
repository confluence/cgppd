#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <fstream>
#include <cstring>
#include <vector>
#include "definitions.h"
#include "AminoAcids.h"
#include "vector3f.h"
#include "Quaternion.h"
#include "TorsionalLookupMatrix.h"
#include "Residue.h"
#include "Link.h"
#include "Segment.h"
#include "Potential.h"

using namespace std;

class Molecule
{
public:
    // TODO: clean this up
    Molecule();
    ~Molecule();
    Molecule(const Molecule& m);
    Molecule operator = (const Molecule& m) {
        return Molecule(m);
    };
    void init(AminoAcid residues[], char *description, char *data);
    bool initFromPDB(const char* pdbfilename);
    void copy(const Molecule& m); // TODO: eliminate?
    void saveBeforeStateChange(const Molecule* m);
    void undoStateChange(const Molecule* m);
    void reserveResidueSpace(int size);
    char *print(); // TODO: eliminate?
    void saveAsPDB(const char *filename);
    void recalculate_relative_positions();
    void recalculate_center(); // TODO: actually use this
    void setPosition(Vector3f v);
    bool translate(Vector3f v);
    void setRotation(Quaternion q);
    bool rotateQ(const Vector3double Raxis, const double angle);

    void rotate();
    void translate();

#if FLEXIBLE_LINKS
    void recalculate_center(Vector3f difference);
    void mark_cached_potentials_for_update(const int ri);

    // TODO: eliminate these; merge into wrappers
    bool translate(Vector3f v, const int ri);
    bool crankshaft(double angle, const bool flip_angle, const int ri);
    bool rotate_domain(const Vector3double raxis, const double angle, const int ri, const bool before);

    void make_local_moves(gsl_rng * rng_linker, gsl_rng * rng_residue, rng_flip);
    void translate_residue();
    void crankshaft();
    void rotate_domain(gsl_rng * rng_linker, gsl_rng * rng_residue, rng_flip);
#endif

    void setMoleculeRoleIdentifier(float moleculeRoleIdentifier);
    bool amIACrowder; // TODO: REMOVE

    float getVolume();// TODO: REMOVE
    float calculateVolume();
    float calculateSASA();

    AminoAcids AminoAcidsData;
    void init_amino_acid_data(AminoAcids &a); // TODO: this is a hack; use proper constructors

    Residue *Residues;
    int residueCount;

    uint random_residue_index(gsl_rng * r); // return random residue index

#if FLEXIBLE_LINKS
    int linkCount;
    int segmentCount;
    int linkerCount;
    Link *Links;
    Segment *Segments; // Boundaries of all domains and linkers
    Segment **Linkers; // Array of pointers to all flexible segments -- TODO: will this do what I think this will do?
    TorsionalLookupMatrix torsions;
    Potential E(const float bounding_value);

    double LJ; // cached inter-segment LJ component of this molecule
    double DH; // cached inter-segment DH component of this molecule
    bool update_LJ_and_DH; // LJ and DH values need to be updated

    uint random_linker_index(gsl_rng * r); // return random linker index
    uint random_residue_index(gsl_rng * r, int li); // return random residue index from given linker
#endif
    size_t chainCount;
    Quaternion rotation;
    Vector3f center;
    int index;
    float volume;
    float moleculeRoleIdentifier;

    float translationalStep;// TODO: REMOVE
    float rotationalStep;// TODO: REMOVE
    char *filename;
    bool hasFilename;
    Vector3f xAxis;
    Vector3f yAxis;
    Vector3f zAxis;
    Vector3f translation;
};

#endif /*MOLECULE_H_*/



