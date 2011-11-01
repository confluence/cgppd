#ifndef LINK_H_
#define LINK_H_

class Link
{
public:
    Link();
    virtual ~Link();
    Link(const Link& l);
    Link operator = (Link l);

    //private:
    // Cached energy potential values

    /* Bond-stretching potential for this link, before multiplication by
     * (0.5 * K_spring), which is performed on the final sum over all links. */
    float e_bond;
    /* Angle potential for the previous residue, before natural log and division
     * by gamma angle, which are performed on the final product over all links.
     * TODO: consider moving this to the residue object. */
    float e_angle;
    /* Torsion-angle potential for this link. */
    float e_torsion;

    // These are set to 1 if local MC moves necessitate an update of the cached values
    bool update_e_bond;
    bool update_e_angle;
    bool update_e_torsion;

    float flexible; // 0 (inflexible) or 1 (flexible): determines fixed tertiary structures
    bool terminal; // true if a new chain starts from this link;
};

#endif /*LINK_H_*/
