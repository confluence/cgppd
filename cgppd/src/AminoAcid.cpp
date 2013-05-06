#include "AminoAcid.h"

AminoAcid::AminoAcid()
{
    shortSymbol = ' ';
    vanderWaalRadius = 0.0f;
    electrostaticCharge = 0.0f;
    index = -1;
    strcpy(name, "");
    strcpy(sname, "");
}


AminoAcid::~AminoAcid()
{

}

char * AminoAcid::getSNAME()
{
    sname[0] = toupper(sname[0]);
    sname[1] = toupper(sname[1]);
    sname[2] = toupper(sname[2]);
    return sname;
}
