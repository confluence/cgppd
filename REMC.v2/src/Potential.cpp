#include "Potential.h"

using namespace std;

double calculate_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData)
{
    double Eij(lambda * (AminoAcidsData.LJpotentials[ri.aminoAcidIndex][rj.aminoAcidIndex] - e0));
    // sigmaij is the average atomic radius determined by the van der waal radius in kim2008
    double sigmaij(0.5f * (ri.vanderWaalRadius + rj.vanderWaalRadius));
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

double calculate_DH(Residue &ri, Residue &rj, const double r)
{
    double DH (ri.electrostaticCharge * rj.electrostaticCharge * expf(-r / Xi) / r);
    return DH;
}

#if FLEXIBLE_LINKS
double calculate_bond(Residue &ri, Link &l, Residue &rj, const float bounding_value)
{
    if (l.update_e_bond)
    {
        // eqn 9: kim2008
        double r(rj.distance(ri, bounding_value));
        LOG(DEBUG_LOCAL_TRANSLATE, "\nri = (%f, %f, %f) rj = (%f, %f, %f) r = %f\n", ri.position.x, ri.position.y, ri.position.z, rj.position.x, rj.position.y, rj.position.z, r);
        l.pseudo_bond = r;
        l.e_bond = (r - R0) * (r - R0); // in angstroms
        l.update_e_bond = false;
    }

    return l.e_bond;
}

double calculate_angle(Residue &rh, Residue &ri, Residue &rj)
{
    if (ri.update_e_angle)
    {
        Vector3f ab = rh.position - ri.position;
        Vector3f cb = rj.position - ri.position;
        double theta = ab.angle(cb);
        ri.pseudo_angle = theta;
        // eqn 10: kim2008
        ri.e_angle = exp(-GammaAngle * (KAlpha * (theta - ThetaAlpha) * (theta - ThetaAlpha) + EpsilonAlpha)) +
                    exp(-GammaAngle * (KBeta *(theta - ThetaBeta) *(theta - ThetaBeta)));

        ri.update_e_angle = false;
    }

    return ri.e_angle;
}

double calculate_torsion(Residue &rh, Residue &ri, Link &l, Residue &rj, Residue &rk, TorsionalLookupMatrix &torsions)
{
    if (l.update_e_torsion)
    {
        // TODO: is this the most efficient way?
        Vector3f b1 = ri.position - rh.position;
        Vector3f b2 = rj.position - ri.position;
        Vector3f b3 = rk.position - rj.position;
        Vector3f b2xb3 = b2.cross(b3);
        double phi = atan2((b2.magnitude() * b1.dot(b2xb3)), b1.cross(b2).dot(b2xb3));

        l.pseudo_torsion = phi;
        // eqn 11: kim2008
        int r1 = ri.aminoAcidIndex;
        int r2 = rj.aminoAcidIndex;
        l.e_torsion = (1 + cos(phi - torsions.getSigma(r1, r2, 1))) * torsions.getV(r1, r2, 1) +
            (1 + cos(2 * phi - torsions.getSigma(r1, r2, 2))) * torsions.getV(r1, r2, 2) +
            (1 + cos(3 * phi - torsions.getSigma(r1, r2, 3))) * torsions.getV(r1, r2, 3) +
            (1 + cos(4 * phi - torsions.getSigma(r1, r2, 4))) * torsions.getV(r1, r2, 4);

        l.update_e_torsion = false;
    }

    return l.e_torsion;
}
#endif // FLEXIBLE_LINKS

Potential::Potential()
{
    LJ = 0.0f;
    DH = 0.0f;

#if FLEXIBLE_LINKS
    LJ_subtotal = 0.0f;
    DH_subtotal = 0.0f;
    bond = 0.0f;
    /* NB: angle potential is a sum of logs -- we multiply the components and take the log at the end. */
    angle = 1.0f;
    torsion = 0.0f;
#endif // FLEXIBLE_LINKS

#if COMPENSATE_KERNEL_SUM
    c_lj = 0.0f;
    c_dh = 0.0f;
#if FLEXIBLE_LINKS
    c_lj_subtotal = 0.0f;
    c_dh_subtotal = 0.0f;
    c_bond = 0.0f;
    c_torsion = 0.0f;
#endif // FLEXIBLE_LINKS
#endif // COMPENSATE_KERNEL_SUM
}

Potential::~Potential()
{
}

#if COMPENSATE_KERNEL_SUM
void Potential::kahan_sum(double &sum, const double i, double &c)
{
    double y(i - c);
    double t(sum + y);
    k.c = (t - sum) - y;
    sum = t;
}
#endif // COMPENSATE_KERNEL_SUM

void Potential::increment_LJ(const double LJ)
{
#if COMPENSATE_KERNEL_SUM
    kahan_sum(this->LJ, LJ, c_lj);
#else // if not COMPENSATE_KERNEL_SUM
    this->LJ += LJ;
#endif // COMPENSATE_KERNEL_SUM
}

void Potential::increment_DH(const double DH)
{
#if COMPENSATE_KERNEL_SUM
    kahan_sum(this->DH, DH, c_dh);
#else // if not COMPENSATE_KERNEL_SUM
    this->DH += DH;
#endif // COMPENSATE_KERNEL_SUM
}

#if FLEXIBLE_LINKS
void Potential::reset_LJ_subtotal()
{
    LJ_subtotal = 0.0f;
#if COMPENSATE_KERNEL_SUM
    c_lj_subtotal = 0.0f;
#endif // COMPENSATE_KERNEL_SUM
}

void Potential::reset_DH_subtotal()
{
    DH_subtotal = 0.0f;
#if COMPENSATE_KERNEL_SUM
    c_dh_subtotal = 0.0f;
#endif // COMPENSATE_KERNEL_SUM
}

void Potential::increment_LJ_subtotal(const double LJ)
{
#if COMPENSATE_KERNEL_SUM
    kahan_sum(this->LJ_subtotal, LJ, c_lj_subtotal);
#else // if not COMPENSATE_KERNEL_SUM
    this->LJ_subtotal += LJ;
#endif // COMPENSATE_KERNEL_SUM
}

void Potential::increment_DH_subtotal(const double DH)
{
#if COMPENSATE_KERNEL_SUM
    kahan_sum(this->DH_subtotal, DH, c_dh_subtotal);
#else // if not COMPENSATE_KERNEL_SUM
    this->DH_subtotal += DH;
#endif // COMPENSATE_KERNEL_SUM
}

void Potential::increment_bond(const double bond)
{
#if COMPENSATE_KERNEL_SUM
    kahan_sum(this->bond, bond, c_bond);
#else // if not COMPENSATE_KERNEL_SUM
    this->bond += bond;
#endif // COMPENSATE_KERNEL_SUM
}

void Potential::increment_angle(const double angle)
{
    this->angle *= angle;
}

void Potential::increment_torsion(const double torsion)
{
#if COMPENSATE_KERNEL_SUM
    kahan_sum(this->torsion, torsion, c_torsion);
#else // if not COMPENSATE_KERNEL_SUM
    this->torsion += torsion;
#endif // COMPENSATE_KERNEL_SUM
}
#endif // FLEXIBLE_LINKS

void Potential::increment(const Potential p)
{
    increment_LJ(p.LJ);
    increment_DH(p.DH);

#if FLEXIBLE_LINKS
    increment_bond(p.bond);
    increment_angle(p.angle);
    increment_torsion(p.torsion);
#endif // FLEXIBLE_LINKS
}

double Potential::total_LJ()
{
    return LJ * LJ_CONVERSION_FACTOR * KBTConversionFactor;
}

double Potential::total_DH()
{
    return DH * DH_constant_component * KBTConversionFactor;
}

#if FLEXIBLE_LINKS
double Potential::total_bond()
{
    LOG(DEBUG_LOCAL_TRANSLATE, "\nbond: %f K_spring: %d KBTConversionFactor: %f total bond: %f\n", bond, K_spring, KBTConversionFactor, bond * 0.5 * K_spring * KBTConversionFactor);
    return bond * 0.5 * K_spring * KBTConversionFactor;
}

double Potential::total_angle()
{
    return log(angle) * -GammaAngleReciprocal * KBTConversionFactor;
}

double Potential::total_torsion()
{
    return torsion * KBTConversionFactor;
}
#endif // FLEXIBLE_LINKS

double Potential::total()
{
#if FLEXIBLE_LINKS
    return (DH * DH_constant_component + LJ * LJ_CONVERSION_FACTOR + bond * 0.5 * K_spring + log(angle) * -GammaAngleReciprocal + torsion) * KBTConversionFactor;
#else // if not FLEXIBLE_LINKS
    return (LJ * LJ_CONVERSION_FACTOR + DH * DH_constant_component) * KBTConversionFactor;
#endif // FLEXIBLE_LINKS
}

void Potential::print_log()
{
#if FLEXIBLE_LINKS
    LOG(INFO_POTENTIAL, "\nDH: %f\tLJ: %f\tbond: %f\tangle: %f\ttorsion: %f\tTOTAL: %f\n", total_LJ(), total_DH(), total_bond(), total_angle(), total_torsion(), total());
#else // if not FLEXIBLE_LINKS
    LOG(INFO_POTENTIAL, "\nDH: %f\tLJ: %f\tTOTAL: %f\n", total_LJ(), total_DH(), total());
#endif // FLEXIBLE_LINKS
}