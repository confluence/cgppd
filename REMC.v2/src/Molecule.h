#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <fstream>
#include <cstring>
#include <vector>
#include "definitions.h"
#include "AminoAcids.h"
#include "vector3f.h"
#include "Quaternion.h"
#include "TorsionalLookupMatrix.h"
#include "Residue.h"
#include "Link.h"
#include "Segment.h"

using namespace std;

struct Molecule_E
{
    double LJ;
    double DH;
    double bond;
    double angle;
    double torsion;
};

class Molecule
{
public:
    // TODO: clean this up
    Molecule();
    ~Molecule();
    Molecule(const Molecule& m);
    Molecule operator = (const Molecule& m) {
        return Molecule(m);
    };
    void deleteResidues(); // TODO: actually use this
    void init(AminoAcid residues[], char *description, char *data);
    bool initFromPDB(const char* pdbfilename);
    void copy(const Molecule& m); // TODO: eliminate?
    void saveBeforeStateChange(const Molecule* m);
    void undoStateChange(const Molecule* m);
    void reserveResidueSpace(int size);
    char *print(); // TODO: eliminate?
    void saveAsPDB(const char *filename);
    Vector3f calculateCenter(); // TODO: actually use this
    void setPosition(Vector3f v);
    bool translate(Vector3f v);
    void setRotation(Quaternion q);
    bool rotateQ(const Vector3double Raxis, const double angle);
    void setMoleculeRoleIdentifier(float moleculeRoleIdentifier);
    bool amIACrowder; // TODO: REMOVE

    float getVolume();// TODO: REMOVE
    float calculateVolume();
    float calculateSASA();

    AminoAcids AminoAcidsData;
    Residue *Residues;
    int residueCount;
#ifdef FLEXIBLE_LINKS
    int linkCount;
    int segmentCount;
//     int linkerCount;
    Link *Links;
    Segment *Segments; // Boundaries of all domains and linkers
//     Segment **Linkers; // Boundaries of all flexible linkers (pointers into Segments array)
    TorsionalLookupMatrix torsions;
    Molecule_E E();

    double LJ; // cached inter-segment LJ component of this molecule
    double DH; // cached inter-segment DH component of this molecule
    bool update_E; // LJ and DH values need to be updated
#endif
    size_t chainCount;
    Vector3f position;
    Quaternion rotation;
    Vector3f center;
    int index;
    float volume;
    float moleculeRoleIdentifier;

    float translationalStep;// TODO: REMOVE
    float rotationalStep;// TODO: REMOVE
    char *filename;
    bool hasFilename;
    Vector3f xAxis;
    Vector3f yAxis;
    Vector3f zAxis;
    Vector3f translation;
};

#endif /*MOLECULE_H_*/



