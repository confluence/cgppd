#ifndef CONSTANTS_H_
    #define CONSTANTS_H_

    // global constants

    #define K_b         1.3806504240e-23f   // boltzmann constant  J/K
    #define E_charge    1.602176487e-19f    // elemental charge in coulombs
    #define Angstrom    1e-10f              // one angstrom
    #define D_water     80.0f               // dialectric constant for water
    #define N_A         6.02214179e23f      // avogadros number
    #define Rgas        8.314472f           // universal gas constant
    #define PI          float(M_PIl)
    #define KBTConversionFactor  (1.0f/(294.0f*Rgas/4184.0f))  // use 294K for verification (should this be 298?)
    #define K_bTtoKcalmol  KBTConversionFactor
    #define KcalmoltoK_bt  (1.0f/K_bTtoKcalmol)
    #define BTU_to_J             0.948f
    #define Xi 10.0f
    #define LJ_CONVERSION_FACTOR 0.3507221006079f  // 1/(KBTConversionFactor)^2
    #define DH_CONVERSION_FACTOR BTU_to_J //10e7f //2.0f
    #define DH_constant_component (DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f)
    #define CF KBTConversionFactor
    #define r0_constant 1.122462048309372981433533049679f // constant factor of r0_ij in LJ calculation

    //#define LJPDSOURCE    "data/pair_potentials_BT"
    //#define e0            0.0372f     // K_b T : experimental fitted value, see [K & B 2008], BT model
    //#define lambda        1.700f      // scalar  : experimental fitted value, see [K & B 2008]
    //#define LJPDSOURCE    "data/pair_potentials_MJ1999"
    //#define e0            0.0328f     // K_b T : experimental fitted value, see [K & B 2008], MJ1999 model
    //#define lambda        1.480f      // scalar  : experimental fitted value, see [K & B 2008]
    #define LJPDSOURCE      "data/pair_potentials_MJ1996"
    #define e0              -2.27f          // K_b T : experimental fitted value, see [K & B 2008], MJ1996 model
    #define lambda          0.159f          // scalar  : experimental fitted value, see [K & B 2008]
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
    #define LJ_OFF true
    #define LJ_REPULSIVE false

    #if LJ_REPULSIVE
        #undef e0
        #define e0  0.0001f // k_b T
    #endif

    #if LJ_OFF
        #undef lambda
        #define lambda 0.0f // not efficient, but takes care of both GPU and CPU calculation
    #endif

#endif
