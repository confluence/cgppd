#ifndef TORSIONALLOOKUPMATRIX_H_
#define TORSIONALLOOKUPMATRIX_H_

#include "AminoAcids.h"
#include "definitions.h"

#define ELPP_NO_DEFAULT_LOG_FILE 1
#include "easylogging++.h"

using namespace std;

class TorsionalLookupMatrix
{
public:
    double V[20][20][4];
    double Sigma[20][20][4];

    TorsionalLookupMatrix();
    ~TorsionalLookupMatrix();

    bool loadData(const char* filename, AminoAcids &aminoAcids);
    double getV(int acidi, int acidj, int n);
    double getSigma(int acidi, int acidj, int n);
};

#endif /*TORSIONALLOOKUPMATRIX_H_*/
