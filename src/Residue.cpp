#include "Residue.h"

Residue::Residue() : aminoAcidIndex(0.0), electrostaticCharge(0.0), vanderWaalRadius(0.0), chainId(' '), resSeq(-1), moleculeId(-1),
#if USING_CUDA
pos_w(0.0f), meta_w(0.0f),
#endif
chain_UID(0), rigid_domain_UID(0), segment_bond_UID(0)
{
    aminoAcidIndex = int(PADDER_IDENTIFIER);
}

Residue::~Residue()
{
}

Residue::Residue(const Residue & r)
{
    aminoAcidIndex = r.aminoAcidIndex;
    electrostaticCharge = r.electrostaticCharge;
    vanderWaalRadius = r.vanderWaalRadius;
    position = r.position;
    relativePosition = r.relativePosition;
    chainId = r.chainId;
    resSeq = r.resSeq;
    moleculeId = r.moleculeId;
}

#if USING_CUDA
void Residue::pack_GPU_floats(int num_segment_bonds, int num_chains)
{
    // pos.w = DOMAIN_UID.BOND_UID
    // meta.w = RESIDUE_ID.CHAIN_UID
    
    // Because we can only compare residue IDs within the same chain, we can use meta.w for the comparison without stripping off the fractional part (which will be the same).
    // We do have to extract both parts of pos.w

   float bond_fraction(0);
   if (segment_bond_UID) {
       //bond_fraction = pow(2.0, -segment_bond_UID);
       bond_fraction = segment_bond_UID / (pow(2.0, ceil( log2( num_segment_bonds + 1 ) )));
   }

   //float chain_fraction(pow(2.0, -chain_UID)); // chain_UID is never zero; they start at 1 and every residue has one
   chain_fraction = chain_UID / (pow(2.0, ceil( log2( num_chains + 1 ) )));
   
   pos_w = (float)rigid_domain_UID + bond_fraction;
   meta_w = (float)resSeq + chain_fraction;
}
#endif

double Residue::distance(const Vector3f p, const float bounding_value)
{
    float Xab(position.x - p.x);
    float Yab(position.y - p.y);
    float Zab(position.z - p.z);

    if (bounding_value)
    {
        Xab = Xab - bounding_value * round(Xab / bounding_value);
        Yab = Yab - bounding_value * round(Yab / bounding_value);
        Zab = Zab - bounding_value * round(Zab / bounding_value);
    }

    return sqrtf(Xab * Xab + Yab * Yab + Zab * Zab);
}


double Residue::distance(const Residue& r, const float bounding_value)
{
    return distance(r.position, bounding_value);
}

// TODO: if we never use this, we should remove it
void Residue::set_rotation(Quaternion q, Vector3f origin)
{
    Vector3f rel_pos = position - origin;
    rel_pos = q.rotateVector(rel_pos);
    // TODO: is adding component-wise an important optimisation?
    position = origin + rel_pos;
}

void Residue::set_rotation_about_center(Quaternion q, Vector3f center)
{
    // use cached relative position
    relativePosition = q.rotateVector(relativePosition);
    // TODO: is adding component-wise an important optimisation?
    position = center + relativePosition;
}

bool Residue::chain_neighbour_of(Residue &other)
{
    return chain_UID == other.chain_UID && abs(resSeq - other.resSeq) < 4;
}

bool Residue::same_rigid_domain_as(Residue &other)
{
    return rigid_domain_UID && rigid_domain_UID == other.rigid_domain_UID;
}

bool Residue::bonded_to(Residue &other)
{
    return segment_bond_UID && segment_bond_UID == other.segment_bond_UID;
}
