#include "Link.h"

Link::Link()
{

    pseudo_bond = 0.0;
    pseudo_angle = 0.0;
    pseudo_torsion = 0.0;

    e_bond = 0.0;
    e_angle = 0.0;
    e_torsion = 0.0;

    update_e_bond = true;
    update_e_angle = true;
    update_e_torsion = true;

    flexible = false;
    dummy = false;
}

Link::~Link()
{
}

Link::Link(const Link & l)
{
    pseudo_bond = l.pseudo_bond;
    pseudo_angle = l.pseudo_angle;
    pseudo_torsion = l.pseudo_torsion;

    e_bond = l.e_bond;
    e_angle = l.e_angle;
    e_torsion = l.e_torsion;

    update_e_bond = l.update_e_bond;
    update_e_angle = l.update_e_angle;
    update_e_torsion = l.update_e_torsion;

    flexible = l.flexible;
    dummy = l.dummy;
}

Link Link::operator =(Link l)
{
    return Link(l);
}
