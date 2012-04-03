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