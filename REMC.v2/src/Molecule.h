#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <fstream>
#include <cstring>
#include <vector>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
    Molecule(const char* pdbfilename); // constructor for molecule loaded from file; not final parameter list
    Molecule(const int size); // constructor for blank molecule of a particular size (for saving before MC move); not final parameter list

    ~Molecule();
    Molecule(const Molecule& m);// TODO: is this used?
    // TODO: is this used?
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
    void translate(Vector3f v);
    void setRotation(Quaternion q);
    void rotate(const Vector3double Raxis, const double angle);

    Vector3f normalised_random_vector_f(gsl_rng * rng);
    Vector3double normalised_random_vector_d(gsl_rng * rng);
    void rotate(gsl_rng * rng);
    // TODO: add boundary conditions to everything?
    void translate(gsl_rng * rng, const float bounding_value);

#if FLEXIBLE_LINKS
    uint random_linker_index(gsl_rng * rng); // return random linker index
    uint random_residue_index(gsl_rng * rng, int li); // return random residue index from given linker

    void recalculate_center(Vector3f difference);
    void mark_cached_potentials_for_update(const int ri);

    void translate(Vector3f v, const int ri);
    void crankshaft(double angle, const bool flip_angle, const int ri);
    void rotate_domain(const Vector3double raxis, const double angle, const int ri, const bool before);

    void rotate_domain(gsl_rng * rng);
    void make_local_moves(gsl_rng * rng);
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

    uint random_residue_index(gsl_rng * rng); // return random residue index

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



