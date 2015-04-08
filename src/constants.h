#ifndef CONSTANTS_H_
    #define CONSTANTS_H_

    // LJ and DH constants

    // These conversions were incorrect in CGPPD v1. They should now be fixed.

    #define E_charge    1.602176487e-19f    // elemental charge (in coulombs)
    #define Ang         1e-10f              // one angstrom (in metres)
    #define D_water     80.0f               // dielectric constant for water
    #define Rgas        8.3144621f          // universal gas constant (in J / mol K)
    #define T_room      298.0f              // room temperature (in Kelvin)
    #define kcal        4184.0f             // one kilocalorie (in joules)
    #define N_A         6.02214129e23f      // Avogadro's constant (in 1 / mol)

    #define r0_constant         pow(2.0f, 1.0f/6.0f) // constant factor in r0_ij

    #define LJPDSOURCE          "data/pair_potentials_MJ1996" // contact potentials from MJ1996 model; see [K&H 2008]
    #define e0                  -2.27f          // RT units: experimental fitted value for MJ1996 model, see [K&H 2008]
    #define LJ_lambda              0.159f          // scalar: experimental fitted value for MJ1996 model, see [K&H 2008]

    #define Xi                      10.0f // Debye screening length (in angstroms)
    #define eps0                    8.85418782e-12f  // permittivity of free space (in farads / metre)

    #define RT_to_kcalmol       (T_room * Rgas / kcal) // conversion from RT units to kcal/mol
    #define DH_CONVERSION_FACTOR    (pow(E_charge, 2.0f) * N_A / (4.0f * float(M_PIl) * eps0 * kcal * Ang * D_water)) // as CCELEC in CHARMM

    // TODO: probably move these to the definitions file; they're implementation details. Maybe move some theoretical stuff from there to here.
    #define AA_COUNT        20
    #define AA_COUNTf       float(AA_COUNT)
    #define LJArraySize     (sizeof(float)*AA_COUNT*AA_COUNT)

    // pseudobond constants

    #define K_spring    378             // 378 kcal / (mol.A^2)
    #define R0          3.81f           // reference bond length in angstroms

    // pseudo angle constants

    #define GammaAngle      0.1f        // mol/kcal
    #define GammaAngleReciprocal    10.0f
    #define EpsilonAlpha    4.3f        // kcal/mol
    #define ThetaAlpha      1.6f        // radians
    #define ThetaBeta       2.27f       // radians
    #define KAlpha          106.4f      // kcal/(mol rad)^2
    #define KBeta           26.3f       // kcal/(mol rad)^2

    // membrane interaction constants

    #define Zreference  20.0f
    #define ephi0       -30.0f          //mV check units
    #define invDebye    1e9f            // inverse debye length water = 1/(10 angtrom)
    #define kappa       invDebye

    // these are used for linker validation
    #ifdef EnableLJOff
        #define LJ_OFF 1
    #else
        #define LJ_OFF 0
    #endif

    #ifdef EnableLJRepulsive
        #define LJ_REPULSIVE 1
    #else
        #define LJ_REPULSIVE 0
    #endif

    #if LJ_REPULSIVE
        #undef e0
        #define e0  0.0001f // k_b T
        #define ASSUME_POLYMER_FOLDING_TEST 1 // shortcut the non-bonded potential criteria
    #endif

    #if LJ_OFF
        //#undef LJ_lambda
        //#define LJ_lambda 0.0f // not efficient, but takes care of both GPU and CPU calculation
        // That was really dumb. We're going to do this properly so that we don't waste days on unnecessary calculations.
        #define ASSUME_POLYMER_FOLDING_TEST 1 // shortcut the non-bonded potential criteria
    #endif

#endif
