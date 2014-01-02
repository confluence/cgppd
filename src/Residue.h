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
    float electrostaticCharge;
    float vanderWaalRadius;

    Vector3f position; // used for collisions
    Vector3f relativePosition; // vs center
    Vector3f new_position; // placeholder for provisional position

    bool isCrowder;

    // upper and lower refer to the numerical index of the residues in the chain, ie: upper = current+1, lower = current-1

    //float sasa; // surface accessable to solvent area
    char chainId;
    int resSeq;
    int moleculeID;
    //AminoAcid aa;

    double distance(const Vector3f p, const float bounding_value=0.0f);
    double distance(const Residue& r, const float bounding_value=0.0f);

    void set_rotation(Quaternion q, Vector3f origin);
    void set_rotation_about_center(Quaternion q, Vector3f center);

#if FLEXIBLE_LINKS
    /* This value is only cached so that it can be checked by unit tests. */
    float pseudo_angle;
    /* Angle potential for this residue, before natural log and division
     * by gamma angle, which are performed on the final product over all links.*/
    float e_angle;
    // This is set to 1 if local MC moves necessitate an update of the cached value
    bool update_e_angle;
#endif
};

#endif