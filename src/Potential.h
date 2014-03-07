#ifndef POTENTIAL_H_
#define POTENTIAL_H_
#include "definitions.h"
#include "constants.h"
#include "TorsionalLookupMatrix.h"
#include "Residue.h"
#include "Geometry.h"
#include <cmath>
#include <iostream>

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
    double c_lj_subtotal;
    double c_dh_subtotal;
    double c_bond;
    double c_torsion;
#endif // FLEXIBLE_LINKS
#endif // COMPENSATE_KERNEL_SUM

public:
    double LJ;
    double DH;

#if FLEXIBLE_LINKS
    double LJ_subtotal;
    double DH_subtotal;

    double bond;
    double angle;
    double torsion;
#endif // FLEXIBLE_LINKS

    Potential();
    ~Potential();

    /* Increment components individually*/
    void increment_LJ(const double LJ);
    void increment_DH(const double DH);

#if FLEXIBLE_LINKS
    void reset_LJ_subtotal();
    void reset_DH_subtotal();
    void increment_LJ_subtotal(const double LJ);
    void increment_DH_subtotal(const double DH);

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

    /* Log */
    void to_string(char * destination);
    void print_log(const bool level, const char * prefix);
};

#endif /*POTENTIAL_H_*/
