#ifndef AMINOACID_H_
#define AMINOACID_H_

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

class AminoAcid
{
public:
    AminoAcid();
    virtual ~AminoAcid();
    char name[20];
    char sname[4];
    char shortSymbol;
    float vanderWaalRadius;
    float electrostaticCharge;
    int index;
    char *getSNAME();
};
#endif /*AMINOACID_H_*/
