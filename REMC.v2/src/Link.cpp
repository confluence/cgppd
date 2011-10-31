#include "Link.h"

Link::Link()
{
    Angle = 0.0;
    BondLength = 0.0;
    TorsionAngle = 0.0;

    e_bond = 0.0;
    e_angle = 0.0;
    e_torsion = 0.0;

    update_e_bond = true;
    update_e_angle = true;
    update_e_torsion = true;

    flexible = false;
    terminal = false;
}

Link::~Link()
{
}

Link::Link(const Link & l)
{
    Angle = l.Angle;
    BondLength = l.BondLength;
    TorsionAngle = l.TorsionAngle;

    e_bond = l.e_bond;
    e_angle = l.e_angle;
    e_torsion = l.e_torsion;

    update_e_bond = l.update_e_bond;
    update_e_angle = l.update_e_angle;
    update_e_torsion = l.update_e_torsion;

    flexible = l.flexible;
    terminal = l.terminal;
}

Link Link::operator =(Link l)
{
    return Link(l);
}
