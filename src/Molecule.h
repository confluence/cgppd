#ifndef MOLECULE_H_
#define MOLECULE_H_

#define ELPP_NO_DEFAULT_LOG_FILE 1
#include "easylogging++.h"

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
#include "Potential.h"
#include "Geometry.h"

using namespace std;

class Molecule
{
public:
    // TODO: clean this up
    Molecule();
    ~Molecule();

    void init(const moldata mol, AminoAcids &a, int index, const float bounding_value);
    vector<Residue> initFromPDB(const char* pdbfilename);
    void copy(const Molecule& m, Residue * contiguous_residue_offset);
    void log_info(int index);

    void MC_backup_restore(const Molecule* m);
    void init_saved_molecule(int max_residue_count);

    void recalculate_relative_positions();
    Vector3f recalculate_center(); // TODO: actually use this
    void setPosition(Vector3f v);
    void translate(Vector3f v);
    void setRotation(Quaternion q);
    void rotate(const Vector3double Raxis, const double angle);

    float bounding_value;

    Vector3f normalised_random_vector_f(gsl_rng * rng);
    Vector3double normalised_random_vector_d(gsl_rng * rng);

    bool centre_outside_boundary(Vector3f c);
    Vector3f wrapped_centre(Vector3f c);
    void wrap_if_necessary();

    void rotate(gsl_rng * rng, const double rotate_step);
    void translate(gsl_rng * rng, const double translate_step);

#if FLEXIBLE_LINKS
    Vector3f recalculate_center(Vector3f difference);

    void translate(Vector3f v, const int ri);
    void crankshaft(double angle, const bool flip_angle, const int ri);
    void flex(const Vector3double raxis, const double angle, const int ri, const int neighbour);

    void flex(gsl_rng * rng, const double rotate_step);
    void make_local_moves(gsl_rng * rng, const double rotate_step, const double translate_step);

    gsl_ran_discrete_t * MC_discrete_table;
#endif

    int get_MC_mutation_type(gsl_rng * rng); // randomly select type of MC move
    void make_MC_move(gsl_rng * rng, const double rotate_step, const double translate_step);

    void setMoleculeRoleIdentifier(float moleculeRoleIdentifier); // this setter is actually useful

    void calculateVolume();
    void calculate_length();

    AminoAcids AminoAcidsData;

    Residue *Residues;
    int residueCount;
    bool contiguous; // whether the array is pointing to a contiguous array created from the replica
    vector<string> residue_sequence();

    uint random_residue_index(gsl_rng * rng); // return random residue index

#if FLEXIBLE_LINKS
    bool is_flexible;

    Graph graph;
    TorsionalLookupMatrix torsion_data;

    Potential E(bool include_LJ_and_DH=true); // TODO TODO TODO refactor

#endif
    int chainCount;
    Quaternion rotation; // TODO remove

    Vector3f center;

    int index;
    float volume;
    float moleculeRoleIdentifier;

    float length;

    char filename[256];
    char name[256];
    char last_MC_move[256];

    Vector3f translation; // TODO remove TODO is this just a copy of the centre, or does it ignore wrapping?
};

#endif /*MOLECULE_H_*/



