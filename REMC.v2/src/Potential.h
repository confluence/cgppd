#ifndef POTENTIAL_H_
#define POTENTIAL_H_
#include "definitions.h"
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

    /* Increment components individually*/
    void increment_LJ(const double LJ);
    void increment_DH(const double DH);
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
