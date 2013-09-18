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
    ~Molecule();

    // TODO: write a proper init function
    void init(const char* pdbfilename, AminoAcids &a, int index, const float bounding_value);
    bool initFromPDB(const char* pdbfilename);
    void copy(const Molecule& m, Residue * contiguous_residue_offset);
    void log_info(int index, int verbosity);

    void MC_backup_restore(const Molecule* m);

#if FLEXIBLE_LINKS
    void init_saved_molecule(int max_residue_count, int max_segment_count);
#else
    void init_saved_molecule(int max_residue_count);
#endif
    char * amino_acid_name(const int ri);
    void saveAsPDB(const char *filename);
    void recalculate_relative_positions();
    Vector3f recalculate_center(); // TODO: actually use this
    void setPosition(Vector3f v);
    void translate(Vector3f v);
    void setRotation(Quaternion q);
    void rotate(const Vector3double Raxis, const double angle);

    float bounding_value;

    Vector3f normalised_random_vector_f(gsl_rng * rng);
    Vector3double normalised_random_vector_d(gsl_rng * rng);
    Vector3f apply_boundary_conditions(Vector3f old_center, Vector3f new_center);

    void rotate(gsl_rng * rng, const double rotate_step);
    void translate(gsl_rng * rng, const double translate_step);

#if FLEXIBLE_LINKS
    uint random_linker_index(gsl_rng * rng); // return random linker index
    uint random_residue_index(gsl_rng * rng, int li); // return random residue index from given linker
    uint random_residue_index_middle(gsl_rng * rng, int li); // return random residue index from given linker which is not on a linker edge

    Vector3f recalculate_center(Vector3f difference);
    void mark_cached_potentials_for_update(const int ri, const bool crankshaft);

    void translate(Vector3f v, const int ri);
    void crankshaft(double angle, const bool flip_angle, const int ri);
    void rotate_domain(const Vector3double raxis, const double angle, const int ri, const bool before);

    void rotate_domain(gsl_rng * rng, const double rotate_step);
    void make_local_moves(gsl_rng * rng, const double rotate_step, const double translate_step);
    bool local_move_successful;

    gsl_ran_discrete_t * MC_discrete_table;
#endif

    uint get_MC_mutation_type(gsl_rng * rng); // randomly select type of MC move
    void make_MC_move(gsl_rng * rng, const double rotate_step, const double translate_step);

    void setMoleculeRoleIdentifier(float moleculeRoleIdentifier); // this setter is actually useful

    void calculateVolume();
    void calculate_length();
    void calculate_radius();

    AminoAcids AminoAcidsData;

    Residue *Residues;
    int residueCount;
    bool contiguous; // whether the array is pointing to a contiguous array created from the replica

    uint random_residue_index(gsl_rng * rng); // return random residue index

#if FLEXIBLE_LINKS
    int linkCount;
    int segmentCount;
    int linkerCount;
    Link *Links;
    Segment *Segments; // Boundaries of all domains and linkers
    int *Linkers; // Array of indices of all flexible segments
    TorsionalLookupMatrix torsions;
    Potential E();

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

    float length;
    float radius;

    char *filename;
    bool hasFilename;
    Vector3f xAxis;
    Vector3f yAxis;
    Vector3f zAxis;
    Vector3f translation;
};

#endif /*MOLECULE_H_*/



