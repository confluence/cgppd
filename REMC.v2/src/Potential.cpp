#include "Potential.h"

using namespace std;

Potential::Potential()
{
    LJ = 0.0f;
    DH = 0.0f;

    c_lj =  0.0f;
    c_dh =  0.0f;

#if FLEXIBLE_LINKS
    bond = 0.0f;
    /* NB: angle potential is a sum of logs -- we multiply the components and take the log at the end. */
    angle = 1.0f;
    torsion = 0.0f;

    c_bond =  0.0f;
    c_torsion =  0.0f;
#endif
}

Potential::~Potential()
{
}

void Potential::kahan_sum(double &sum, const double i, double &c)
{
    double y(i - c);
    double t(sum + y);
    c = (t - sum) - y;
    sum = t;
}

void Potential::sum(double &sum, const double i, double &c)
{
#if COMPENSATE_KERNEL_SUM
    kahan_sum(sum, i, c);
#else
    sum += i;
#endif
}
void Potential::increment_LJ(const double LJ)
{
    sum(this->LJ, LJ, c_lj);
}

void Potential::increment_DH(const double DH)
{
    sum(this->DH, DH, c_dh);
}

#if FLEXIBLE_LINKS
void Potential::increment_bond(const double bond)
{
    sum(this->bond, bond, c_bond);
}

void Potential::increment_angle(const double angle)
{
    this->angle *= angle;
}

void Potential::increment_torsion(const double torsion)
{
    sum(this->torsion, torsion, c_torsion);
}
#endif

void Potential::increment(const Potential p)
{
    sum(LJ, p.LJ, c_lj);
    sum(DH, p.DH, c_dh);

#if FLEXIBLE_LINKS
    sum(bond, p.bond, c_bond);
    angle *= p.angle;
    sum(torsion, p.torsion, c_torsion);
#endif
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
//     cout << "bond " << (bond) << endl;
//     cout << "bond * 0.5 " << (bond * 0.5) << endl;
//     cout << "bond * 0.5 * K_spring " << (bond * 0.5 * K_spring) << endl;
//     cout << "bond * 0.5 * K_spring * KBTConversionFactor " << (bond * 0.5 * K_spring * KBTConversionFactor) << endl;
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
#endif

double Potential::total()
{
#if FLEXIBLE_LINKS
    return (DH * DH_constant_component + LJ * LJ_CONVERSION_FACTOR + bond * 0.5 * K_spring + log(angle) * -GammaAngleReciprocal + torsion) * KBTConversionFactor;
#else
    return (LJ * LJ_CONVERSION_FACTOR + DH * DH_constant_component) * KBTConversionFactor;
#endif
}