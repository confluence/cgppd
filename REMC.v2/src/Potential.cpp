#include "Potential.h"

using namespace std;

#if COMPENSATE_KERNEL_SUM
KahanTuple::KahanTuple()
{
    double dsum(0.0f);
    sum = &dsum;
    c = 0.0f;
}

KahanTuple::KahanTuple(double &dsum)
{
    sum = &dsum;
    c = 0.0f;
}
#endif // COMPENSATE_KERNEL_SUM

Potential::Potential()
{
#if !COMPENSATE_KERNEL_SUM
    LJ = 0.0f;
    DH = 0.0f;
#endif // not COMPENSATE_KERNEL_SUM

#if FLEXIBLE_LINKS
    /* NB: angle potential is a sum of logs -- we multiply the components and take the log at the end. */
    angle = 1.0f;

#if !COMPENSATE_KERNEL_SUM
    bond = 0.0f;
    torsion = 0.0f;
#endif // not COMPENSATE_KERNEL_SUM
#endif // FLEXIBLE_LINKS
}

Potential::~Potential()
{
}

#if COMPENSATE_KERNEL_SUM
double Potential::get(const KahanTuple &sum)
{
    return *sum.sum;
}
#else // if not COMPENSATE_KERNEL_SUM
double Potential::get(const double &sum)
{
    return sum;
}
#endif // COMPENSATE_KERNEL_SUM

void Potential::increment_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData)
{
    increment_LJ(ri, rj, r, AminoAcidsData, LJ);
}

void Potential::increment_DH(Residue &ri, Residue &rj, const double r)
{
    increment_DH(ri, rj, r, DH);
}

#if COMPENSATE_KERNEL_SUM
void Potential::increment_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData, KahanTuple &sum)
#else // if not COMPENSATE_KERNEL_SUM
void Potential::increment_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData, double &sum)
#endif // COMPENSATE_KERNEL_SUM
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

    increment_LJ(sum, LJ);
}

#if COMPENSATE_KERNEL_SUM
void Potential::increment_DH(Residue &ri, Residue &rj, const double r, KahanTuple &sum)
#else // if not COMPENSATE_KERNEL_SUM
void Potential::increment_DH(Residue &ri, Residue &rj, const double r, double &sum)
#endif // COMPENSATE_KERNEL_SUM
{
    double DH (ri.electrostaticCharge * rj.electrostaticCharge * expf(-r / Xi) / r);
    increment_DH(sum, DH);
}

#if FLEXIBLE_LINKS
void Potential::increment_bond(Residue &ri, Link &l, Residue &rj, const float bounding_value)
{
    if (l.update_e_bond)
    {
        // eqn 9: kim2008
        double r(rj.distance(ri, bounding_value));
        l.pseudo_bond = r;
        l.e_bond = (r - R0) * (r - R0); // in angstroms
        l.update_e_bond = false;
    }

    increment_bond(l.e_bond);
}

void Potential::increment_angle(Residue &rh, Residue &ri, Residue &rj)
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

    increment_angle(ri.e_angle);
}

void Potential::increment_torsion(Residue &rh, Residue &ri, Link &l, Residue &rj, Residue &rk, TorsionalLookupMatrix &torsions)
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

    increment_torsion(l.e_torsion);
}
#endif

#if COMPENSATE_KERNEL_SUM
void Potential::kahan_sum(KahanTuple &k, const double i)
{
    double y(i - k.c);
    double t(*k.sum + y);
    k.c = (t - *k.sum) - y;
    *k.sum = t;
}

void Potential::sum(KahanTuple &sum, const double i)
{
    kahan_sum(sum, i);
}
#else // if not COMPENSATE_KERNEL_SUM
void Potential::sum(double &sum, const double i)
{
    sum += i;
}
#endif // COMPENSATE_KERNEL_SUM

void Potential::increment_LJ(const double LJ)
{
    increment_LJ(this->LJ, LJ);
}

void Potential::increment_DH(const double DH)
{
    increment_DH(this->DH, DH);
}

#if COMPENSATE_KERNEL_SUM
void Potential::increment_LJ(KahanTuple &sum, const double LJ)
#else // if not COMPENSATE_KERNEL_SUM
void Potential::increment_LJ(double &sum, const double LJ)
#endif // COMPENSATE_KERNEL_SUM
{
    this->sum(sum, LJ);
}

#if COMPENSATE_KERNEL_SUM
void Potential::increment_DH(KahanTuple &sum, const double DH)
#else // if not COMPENSATE_KERNEL_SUM
void Potential::increment_DH(double &sum, const double DH)
#endif // COMPENSATE_KERNEL_SUM
{
    this->sum(sum, DH);
}

#if FLEXIBLE_LINKS
void Potential::increment_bond(const double bond)
{
    sum(this->bond, bond);
}

void Potential::increment_angle(const double angle)
{
    this->angle *= angle;
}

void Potential::increment_torsion(const double torsion)
{
    sum(this->torsion, torsion);
}
#endif

void Potential::increment(const Potential p)
{
    sum(LJ, get(p.LJ));
    sum(DH, get(p.DH));

#if FLEXIBLE_LINKS
    sum(bond, get(p.bond));
    angle *= p.angle;
    sum(torsion, get(p.torsion));
#endif
}

double Potential::total_LJ()
{
    return get(LJ) * LJ_CONVERSION_FACTOR * KBTConversionFactor;
}

double Potential::total_DH()
{
    return get(DH) * DH_constant_component * KBTConversionFactor;
}

#if FLEXIBLE_LINKS
double Potential::total_bond()
{
    return get(bond) * 0.5 * K_spring * KBTConversionFactor;
}

double Potential::total_angle()
{
    return log(angle) * -GammaAngleReciprocal * KBTConversionFactor;
}

double Potential::total_torsion()
{
    return get(torsion) * KBTConversionFactor;
}
#endif

double Potential::total()
{
#if FLEXIBLE_LINKS
    return (get(DH) * DH_constant_component + get(LJ) * LJ_CONVERSION_FACTOR + get(bond) * 0.5 * K_spring + log(angle) * -GammaAngleReciprocal + get(torsion)) * KBTConversionFactor;
#else
    return (get(LJ) * LJ_CONVERSION_FACTOR + get(DH) * DH_constant_component) * KBTConversionFactor;
#endif
}