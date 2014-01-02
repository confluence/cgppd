#ifndef SEGMENT_H_
#define SEGMENT_H_
#include <gsl/gsl_rng.h>

class Segment
{
public:
    int start; // index of starting residue
    int end; // index of ending residue
    int size; // length of this segment, including both end residues

    bool flexible; // this segment is a flexible linker

    double LJ; // cached internal LJ component of this segment if it is a linker
    double DH; // cached internal DH component of this segment if it is a linker
    bool update_LJ_and_DH; // LJ and DH values need to be updated

    Segment();
    uint random_residue_index(gsl_rng * r);
    uint random_residue_index_middle(gsl_rng * r);
};

#endif
