#include "Residue.h"

Residue::Residue()
{
    aminoAcidIndex = int(PADDER_IDENTIFIER);
#ifdef FLEXIBLE_LINKS
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
#ifdef FLEXIBLE_LINKS
    pseudo_angle = r.pseudo_angle;
    e_angle = r.e_angle;
    update_e_angle = r.update_e_angle;
#endif
}
