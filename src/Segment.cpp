#include "Segment.h"

Segment::Segment()
{
    start = 0;
    end = 0;
    size = 0;

    flexible = false;

    LJ = 0.0f;
    DH = 0.0f;
    update_LJ_and_DH = true;
}

uint Segment::random_residue_index(gsl_rng * r)
{
    return (int) gsl_rng_uniform_int(r, size) + start;
}

uint Segment::random_residue_index_middle(gsl_rng * r)
{

    if (size > 2)
    {
        return (int) gsl_rng_uniform_int(r, size - 2) + start + 1;
    }
    else
    {
        return -1;
    }
}