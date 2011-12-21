#ifndef POTENTIAL_H_
#define POTENTIAL_H_
#include "definitions.h"
#include "TorsionalLookupMatrix.h"
#include "Residue.h"
#include "Link.h"
#include <cmath>
#include <iostream>

#if COMPENSATE_KERNEL_SUM
struct KahanTuple
{
    double * sum;
    double c;

    KahanTuple();
    KahanTuple(&sum);
    ~KahanTuple();
};
#endif // COMPENSATE_KERNEL_SUM

class Potential
{
private:
#if COMPENSATE_KERNEL_SUM
    void kahan_sum(KahanTuple &sum, const double i);
    void sum(KahanTuple &sum, const double i);
#else // if not COMPENSATE_KERNEL_SUM
    void sum(double &sum, const double i);
#endif // COMPENSATE_KERNEL_SUM

public:

#if COMPENSATE_KERNEL_SUM
    KahanTuple LJ;
    KahanTuple DH;
#else // if not COMPENSATE_KERNEL_SUM
    double LJ;
    double DH;
#endif // COMPENSATE_KERNEL_SUM

#if FLEXIBLE_LINKS
    double angle;
#if COMPENSATE_KERNEL_SUM
    KahanTuple bond;
    KahanTuple torsion;
#else // if not COMPENSATE_KERNEL_SUM
    double bond;
    double torsion;
#endif // COMPENSATE_KERNEL_SUM
#endif // FLEXIBLE_LINKS

    Potential();
    ~Potential();

#if COMPENSATE_KERNEL_SUM
    double get(const KahanTuple &sum);
#else // if not COMPENSATE_KERNEL_SUM
    double get(const double &sum);
#endif // COMPENSATE_KERNEL_SUM

    /* Increment components individually after calculating */
    void increment_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData);
    void increment_DH(Residue &ri, Residue &rj, const double r);
#if COMPENSATE_KERNEL_SUM
    void increment_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData, KahanTuple &sum);
    void increment_DH(Residue &ri, Residue &rj, const double r, KahanTuple &sum);
#else // if not COMPENSATE_KERNEL_SUM
    void increment_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData, double &sum);
    void increment_DH(Residue &ri, Residue &rj, const double r, double &sum);
#endif // COMPENSATE_KERNEL_SUM
#if FLEXIBLE_LINKS
    void increment_bond(Residue &ri, Link &l, Residue &rj, const float bounding_value);
    void increment_angle(Residue &rh, Residue &ri, Residue &rj);
    void increment_torsion(Residue &rh, Residue &ri, Link &l, Residue &rj, Residue &rk, TorsionalLookupMatrix &torsions);
#endif // FLEXIBLE_LINKS

    /* Increment components individually*/
    void increment_LJ(const double LJ);
    void increment_DH(const double DH);
#if COMPENSATE_KERNEL_SUM
    void increment_LJ(KahanTuple &sum, const double LJ);
    void increment_DH(KahanTuple &sum, const double DH);
#else // if not COMPENSATE_KERNEL_SUM
    void increment_LJ(double &sum, const double LJ);
    void increment_DH(double &sum, const double DH);
#endif // COMPENSATE_KERNEL_SUM
#if FLEXIBLE_LINKS
    void increment_bond(const double bond);
    void increment_angle(const double angle);
    void increment_torsion(const double torsion);
#endif // FLEXIBLE_LINKS

    /* Increment all components */
    void increment(const Potential p);

    /* Component totals (for unit tests) */
    double total_LJ();
    double total_DH();
#if FLEXIBLE_LINKS
    double total_bond();
    double total_angle();
    double total_torsion();
#endif // FLEXIBLE_LINKS

    /* Output final total */
    double total();
};

#endif /*POTENTIAL_H_*/
