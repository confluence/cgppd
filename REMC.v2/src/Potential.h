#ifndef POTENTIAL_H_
#define POTENTIAL_H_
#include "definitions.h"
#include "TorsionalLookupMatrix.h"
#include "Residue.h"
#include "Link.h"
#include <cmath>
#include <iostream>

class Potential
{
private:
    double c_lj;
    double c_dh;
#if FLEXIBLE_LINKS
    double c_bond;
    double c_torsion;
#endif

    void kahan_sum(double &sum, const double i, double &c);
    void sum(double &sum, const double i, double &c);

public:
    double LJ;
    double DH;
#if FLEXIBLE_LINKS
    double bond;
    double angle;
    double torsion;
#endif

    Potential();
    ~Potential();

    /* Increment components individually after calculating */
    void increment_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData);
    void increment_DH(Residue &ri, Residue &rj, const double r);
    void increment_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData, double &sum, double &c);
    void increment_DH(Residue &ri, Residue &rj, const double r, double &sum, double &c);
#if FLEXIBLE_LINKS
    void increment_bond(Residue &ri, Link &l, Residue &rj, const float bounding_value);
    void increment_angle(Residue &rh, Residue &ri, Residue &rj);
    void increment_torsion(Residue &rh, Residue &ri, Link &l, Residue &rj, Residue &rk, TorsionalLookupMatrix &torsions);
#endif

    /* Increment components individually*/
    void increment_LJ(const double LJ);
    void increment_DH(const double DH);
    void increment_LJ(double &sum, const double LJ, double &c);
    void increment_DH(double &sum, const double DH, double &c);
#if FLEXIBLE_LINKS
    void increment_bond(const double bond);
    void increment_angle(const double angle);
    void increment_torsion(const double torsion);
#endif

    /* Increment all components */
    void increment(const Potential p);

    /* Component totals (for unit tests) */
    double total_LJ();
    double total_DH();
#if FLEXIBLE_LINKS
    double total_bond();
    double total_angle();
    double total_torsion();
#endif

    /* Output final total */
    double total();
};

#endif /*POTENTIAL_H_*/
