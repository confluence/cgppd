#ifndef RESIDUE_H_
#define RESIDUE_H_
#include "vector3f.h"
#include "AminoAcids.h"
#include "definitions.h"
#include "Quaternion.h"

class Residue
{

public:
    Residue();
    Residue(const Residue& r);
    ~Residue();
    int aminoAcidIndex;  // stores the index of the amino acid for lookup
    float electrostaticCharge; // TODO: we could look this up in AminoAcids
    float vanderWaalRadius; // TODO: we could look this up in AminoAcids

    Vector3f position; // used for collisions
    Vector3f relativePosition; // vs center

    char chainId;
    int chain_int_id;
    int resSeq;
    int moleculeId; // TODO: this is probably unused
    
    // global identifiers for chain and optionally rigid domain and segment bond
    int chain_UID;
    int rigid_domain_UID;
    int segment_bond_UID;
    
#if USING_CUDA
    // cached values to be copied to the GPU
    // these store a pregenerated combination of chain uid, residue id, rigid domain id and bond id
    float pos_w; // domain.bond
    float meta_w; // residue.chain
    
    void pack_GPU_floats(int num_segment_bonds, int num_chains);
#endif

    double distance(const Vector3f p, const float bounding_value=0.0f);
    double distance(const Residue& r, const float bounding_value=0.0f);

    void set_rotation(Quaternion q, Vector3f origin);
    void set_rotation_about_center(Quaternion q, Vector3f center);
    
    // These are properties of a pair of residues within the same molecule which cause them to be excluded from the unbonded potential total
    bool chain_neighbour_of(Residue &other);
    bool same_rigid_domain_as(Residue &other);
    bool bonded_to(Residue &other);
};

#endif
