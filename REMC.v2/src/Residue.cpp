#include "Residue.h"

Residue::Residue()
{
    aminoAcidIndex = int(PADDER_IDENTIFIER);
#if FLEXIBLE_LINKS
    pseudo_angle = 0.0;
    e_angle = 0.0;
    update_e_angle = true;
#endif
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
    //sasa = r.sasa;
    isCrowder = r.isCrowder;
    chainId = r.chainId;
    resSeq = r.resSeq;
#if FLEXIBLE_LINKS
    pseudo_angle = r.pseudo_angle;
    e_angle = r.e_angle;
    update_e_angle = r.update_e_angle;
#endif
}

double Residue::distance(const Residue& r, const float bounding_value)
{
    float Xab(position.x - r.position.x);
    float Yab(position.y - r.position.y);
    float Zab(position.z - r.position.z);

    Xab = Xab - bounding_value * round(Xab / bounding_value);
    Yab = Yab - bounding_value * round(Yab / bounding_value);
    Zab = Zab - bounding_value * round(Zab / bounding_value);

    return sqrtf(Xab * Xab + Yab * Yab + Zab * Zab);
}

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