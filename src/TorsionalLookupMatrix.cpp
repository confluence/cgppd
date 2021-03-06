#include "TorsionalLookupMatrix.h"

TorsionalLookupMatrix::TorsionalLookupMatrix()
{
}

TorsionalLookupMatrix::~TorsionalLookupMatrix()
{
}

double TorsionalLookupMatrix::getV(int ai, int aj, int n)
{
    return V[ai][aj][n - 1];
}

double TorsionalLookupMatrix::getSigma(int ai, int aj, int n)
{
    return Sigma[ai][aj][n - 1];
}

bool TorsionalLookupMatrix::loadData(const char* filename, AminoAcids &aminoAcids)
{
    ifstream input(filename);
    if (!input.good())
    {
        throw "Failed to open torsional data file.";
    }

    char line[256];

    char acid0[5] = {0};
    char acid1[5] = {0};
    double v;
    double sigma;
    int n;

    while (!input.eof())
    {
        input.getline(line, 255);
        // ignore blank lines and comments
        if (strcmp(line, "") != 0 && strncmp(line, "#", 1) != 0)
        {
            sscanf(line,"%3s %3s %lf %d %lf", acid0, acid1, &v, &n, &sigma);
            int a0 = aminoAcids.getAminoAcidIndex(acid0);
            int a1 = aminoAcids.getAminoAcidIndex(acid1);
            V[a0][a1][n - 1] = v;
            Sigma[a0][a1][n - 1] = sigma/57.29577951308232; // convert to radians
        }
    }
    input.close();
    return true;
}
