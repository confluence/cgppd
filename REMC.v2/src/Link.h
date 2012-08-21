#ifndef LINK_H_
#define LINK_H_

class Link
{
public:
    Link();
    virtual ~Link();
    Link(const Link& l);
    Link operator = (Link l);

    /* These values are only cached so that they can be checked by unit tests. */
    float pseudo_bond;
    float pseudo_torsion;

    // Cached energy potential values

    /* Bond-stretching potential for this link, before multiplication by
     * (0.5 * K_spring), which is performed on the final sum over all links. */
    float e_bond;
    /* Torsion-angle potential for this link. */
    float e_torsion;

    // These are set to 1 if local MC moves necessitate an update of the cached values
    bool update_e_bond;
    bool update_e_torsion;

    float flexible; // true if this link is part of a flexible linker rather than a rigid domain
    bool dummy; // not actually a link; the space between the last residue of a chain and the first residue of the next
};

#endif /*LINK_H_*/
