#ifndef POTENTIAL_H_
#define POTENTIAL_H_
#include "definitions.h"
#include "constants.h"
#include "TorsionalLookupMatrix.h"
#include "Residue.h"
#include "Geometry.h"
#include <cmath>
#include <iostream>
#include <iomanip>

/* Potential component calculations */
double calculate_LJ(Residue &ri, Residue &rj, const double r, const AminoAcids &AminoAcidsData);
double calculate_DH(Residue &ri, Residue &rj, const double r);
double calculate_bond(Residue * residues, Bond &bond, const float bounding_value);
double calculate_angle(Residue * residues, Angle &angle);
double calculate_torsion(Residue * residues, Torsion &torsion, TorsionalLookupMatrix &torsion_data);

class Potential
{
private:
#if COMPENSATE_KERNEL_SUM
    void kahan_sum(double &sum, const double i, double &c);

    double c_lj;
    double c_dh;
#if FLEXIBLE_LINKS
    double c_bond;
    double c_torsion;
#endif // FLEXIBLE_LINKS
#endif // COMPENSATE_KERNEL_SUM

public:
    double LJ;
    double DH;

#if FLEXIBLE_LINKS
    double bond;
    double angle;
    double torsion;
#endif // FLEXIBLE_LINKS

    Potential();
    Potential(const double LJ, const double DH, const double bond=0.0f, const double angle=0.0f, const double torsion=0.0f);
    ~Potential();
    
    bool operator == (const Potential p) {
#if FLEXIBLE_LINKS
        return (LJ == p.LJ) && (DH == p.DH) && (bond == p.bond) && (angle == p.angle) && (torsion = p.torsion);
#else
        return (LJ == p.LJ) && (DH == p.DH);
#endif // FLEXIBLE_LINKS
    };
    
    bool almost_equal(const Potential p, const float eps);

    friend ostream& operator<<(ostream& os, const Potential p);

    /* Increment components individually*/
    void increment_LJ(const double LJ);
    void increment_DH(const double DH);

#if FLEXIBLE_LINKS
    void increment_bond(const double bond);
    void increment_angle(const double angle);
    void increment_torsion(const double torsion);
#endif // FLEXIBLE_LINKS

    /* Increment all components */
    void increment(const Potential p);

    /* Component totals (for unit tests) */
    double total_LJ () const;
    double total_DH () const;
#if FLEXIBLE_LINKS
    double total_bond () const;
    double total_angle () const;
    double total_torsion () const;
#endif // FLEXIBLE_LINKS

    /* Output final total */
    double total () const;
};

#endif /*POTENTIAL_H_*/
