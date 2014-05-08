#include "Potential.h"

using namespace std;

double calculate_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData)
{
    /* This calculation is described in kim2008, in equations 2, 3, 5 and 6
     *
     * Eij is the LJ interaction strength between residues i and j.
     * Sigmaij is the average atomic radius of residues i and j.
     */

    // Eij = LJ_lambda(eij - e0)
    double Eij(LJ_lambda * (AminoAcidsData.LJpotentials[ri.aminoAcidIndex][rj.aminoAcidIndex] - e0));

    // sigmaij = (sigmai + sigmaj) / 2
    double sigmaij(0.5f * (ri.vanderWaalRadius + rj.vanderWaalRadius));
    double LJtmp(pow(sigmaij / r, 6.0f)); // (sigmaij/r)^6

    // if attractive interactions (Eij < 0), or repulsive interactions (Eij > 0) and r >= r0ij:
    // uij(r) = -4Eij( (sigmaij/r)^12 - (sigmaij/r)^6 )
    double LJ(-4.0f * Eij * LJtmp * (LJtmp - 1.0f));

    // if repulsive interactions (Eij > 0) and r < r0ij:
    // uij(r) = 4Eij( (sigmaij/r)^12 - (sigmaij/r)^6 ) + 2Eij
    if (Eij > 0.0f && r < sigmaij * r0_constant)
    {
        LJ = -LJ + 2.0f * Eij;
    }

    // NOTE: in kim2008, uij(r) for attractive interactions is written as 4|Eij|( (sigmaij/r)^12 - (sigmaij/r)^6 )
    // Because Eij is negative, 4|Eij| = -4Eij

    LOG(DEBUG && (fabs(LJ) > 1.0), "!!! ");
    LOG(DEBUG, "M%d R%d / M%d R%d: r = %f LJ = %f\n", ri.moleculeId, ri.resSeq, rj.moleculeId, rj.resSeq, r, LJ);

    return LJ;
}

double calculate_DH(Residue &ri, Residue &rj, const double r)
{
    // kim2008, equation 4
    double DH (ri.electrostaticCharge * rj.electrostaticCharge * expf(-r / Xi) / r);
    return DH;
}

#if FLEXIBLE_LINKS
double calculate_bond(Residue * residues, Bond &bond, const float bounding_value)
{
    if (bond.update_potential)
    {
        Residue & ri = residues[bond.i];
        Residue & rj = residues[bond.j];
        // eqn 9: kim2008
        double r(rj.distance(ri, bounding_value));
        bond.length = r;
        bond.potential = (r - R0) * (r - R0); // in angstroms
        bond.update_potential = false;
    }

    return bond.potential;
}

double calculate_angle(Residue * residues, Angle &angle)
{

    if (angle.update_potential)
    {
        Residue & ri = residues[angle.i];
        Residue & rj = residues[angle.j];
        Residue & rk = residues[angle.k];

        Vector3f ab = ri.position - rj.position;
        Vector3f cb = rk.position - rj.position;
        double theta = ab.angle(cb);
        angle.theta = theta;
        // eqn 10: kim2008
        angle.potential = exp(-GammaAngle * (KAlpha * (theta - ThetaAlpha) * (theta - ThetaAlpha) + EpsilonAlpha)) +
                    exp(-GammaAngle * (KBeta *(theta - ThetaBeta) *(theta - ThetaBeta)));

        angle.update_potential = false;
    }

    return angle.potential;
}

double calculate_torsion(Residue * residues, Torsion &torsion, TorsionalLookupMatrix &torsion_data)
{
    if (torsion.update_potential)
    {
        Residue & ri = residues[torsion.i];
        Residue & rj = residues[torsion.j];
        Residue & rk = residues[torsion.k];
        Residue & rl = residues[torsion.l];

        // TODO: is this the most efficient way?
        Vector3f b1 = rj.position - ri.position;
        Vector3f b2 = rk.position - rj.position;
        Vector3f b3 = rl.position - rk.position;
        Vector3f b2xb3 = b2.cross(b3);
        double phi = atan2((b2.magnitude() * b1.dot(b2xb3)), b1.cross(b2).dot(b2xb3));

        torsion.phi = phi;
        // eqn 11: kim2008
        int r1 = rj.aminoAcidIndex;
        int r2 = rk.aminoAcidIndex;


        torsion.potential = (1 + cos(phi - torsion_data.getSigma(r1, r2, 1))) * torsion_data.getV(r1, r2, 1) +
            (1 + cos(2 * phi - torsion_data.getSigma(r1, r2, 2))) * torsion_data.getV(r1, r2, 2) +
            (1 + cos(3 * phi - torsion_data.getSigma(r1, r2, 3))) * torsion_data.getV(r1, r2, 3) +
            (1 + cos(4 * phi - torsion_data.getSigma(r1, r2, 4))) * torsion_data.getV(r1, r2, 4);

        torsion.update_potential = false;
    }

    return torsion.potential;
}
#endif // FLEXIBLE_LINKS

Potential::Potential() : 
#if FLEXIBLE_LINKS
/* NB: angle potential is a sum of logs -- we multiply the components and take the log at the end. */
bond(0.0f), angle(1.0f), torsion(0.0f),
#endif // FLEXIBLE_LINKS
#if COMPENSATE_KERNEL_SUM
c_lj(0.0f), c_dh(0.0f),
#if FLEXIBLE_LINKS
c_bond(0.0f), c_torsion(0.0f),
#endif // FLEXIBLE_LINKS
#endif // COMPENSATE_KERNEL_SUM
LJ(0.0f), DH(0.0f)
{
}

Potential::Potential(const double LJ, const double DH, const double bond, const double angle, const double torsion) : Potential()
{
    // Reverse constant factors to store raw totals
    this->LJ = LJ/RT_to_kcalmol;
    this->DH = DH/DH_CONVERSION_FACTOR;
#if FLEXIBLE_LINKS
    this->bond = bond/(0.5 * K_spring);
    this->angle = exp(angle/-GammaAngleReciprocal);
    this->torsion = torsion;
#endif // FLEXIBLE_LINKS
}

Potential::~Potential()
{
}

#if COMPENSATE_KERNEL_SUM
void Potential::kahan_sum(double &sum, const double i, double &c)
{
    double y(i - c);
    double t(sum + y);
    c = (t - sum) - y;
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

void Potential::increment(const Potential p) // TODO: replace this with an addition operator
{
    increment_LJ(p.LJ);
    increment_DH(p.DH);

#if FLEXIBLE_LINKS
    increment_bond(p.bond);
    increment_angle(p.angle);
    increment_torsion(p.torsion);
#endif // FLEXIBLE_LINKS
}

double Potential::total_LJ () const
{
    return LJ * RT_to_kcalmol;
}

double Potential::total_DH () const
{
    return DH * DH_CONVERSION_FACTOR;
}

#if FLEXIBLE_LINKS
double Potential::total_bond () const
{
    return bond * 0.5 * K_spring;
}

double Potential::total_angle () const
{
    return log(angle) * -GammaAngleReciprocal;
}

double Potential::total_torsion () const
{
    return torsion;
}
#endif // FLEXIBLE_LINKS

double Potential::total () const
{
#if FLEXIBLE_LINKS
    return DH * DH_CONVERSION_FACTOR + LJ * RT_to_kcalmol + bond * 0.5 * K_spring + log(angle) * -GammaAngleReciprocal + torsion;
#else // if not FLEXIBLE_LINKS
    return DH * DH_CONVERSION_FACTOR + LJ * RT_to_kcalmol;
#endif // FLEXIBLE_LINKS
}

bool Potential::almost_equal(const Potential p, const float eps)
{
#if FLEXIBLE_LINKS
    return (fabs(LJ - p.LJ) < eps && fabs(DH - p.DH) < eps && fabs(bond - p.bond) < eps && fabs(angle - p.angle) < eps && fabs(torsion - p.torsion) < eps);
#else // if not FLEXIBLE_LINKS
    return (fabs(LJ - p.LJ) < eps && fabs(DH - p.DH) < eps);
#endif // FLEXIBLE_LINKS
}

ostream& operator<<(ostream& os, const Potential p)
{
#if FLEXIBLE_LINKS
    os << fixed << setprecision(6) << "Potential(total=" << p.total() << ", LJ=" << p.total_LJ() << ", DH=" << p.total_DH() << ", bond=" << p.total_bond() << ", angle=" << p.total_angle() << ", torsion=" << p.total_torsion() << ")";
#else // if not FLEXIBLE_LINKS
    os << fixed << setprecision(6) << "Potential(total=" << p.total() << ", LJ=" << p.total_LJ() << ", DH=" << p.total_DH() << ")";
#endif // FLEXIBLE_LINKS
    return os;
};
