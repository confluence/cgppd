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
};

#endif
