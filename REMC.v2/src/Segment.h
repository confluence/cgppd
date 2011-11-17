#ifndef SEGMENT_H_
#define SEGMENT_H_

class Segment
{
public:
    int start; // index of starting residue
    int end; // index of ending residue
    bool flexible; // this segment is a flexible linker

    Segment();
};

#endif
