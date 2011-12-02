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


#ifdef FLEXIBLE_LINKS

float Residue::LJ_component(const Residue & rj, const double r, const AminoAcids & AminoAcidsData)
{
    float Eij(lambda * (AminoAcidsData.LJpotentials[aminoAcidIndex][rj.aminoAcidIndex] - e0));
    // sigmaij is the average atomic radius determined by the van der waal radius in kim2008
    float sigmaij(0.5f * (vanderWaalRadius + rj.vanderWaalRadius));
    //float r0 = powf(2.0f,(1.0f/6.0f));
    //float sigT = sigmaij / r ;
    double LJtmp(powf(sigmaij / r, 6.0f)); //sigT*sigT*sigT*sigT*sigT*sigT;
    double LJ(-4.0f * Eij * LJtmp * (LJtmp - 1.0f));

    if (Eij > 0.0f && r < sigmaij * r0_constant)  // attractive pairs
    {
        LJ = -LJ + 2.0f * Eij;
    }

    return LJ;
}

float Residue::DH_component(const Residue & rj, const double r)
{
    return electrostaticCharge * rj.electrostaticCharge * expf(-r / Xi) / r;
}

#endif
