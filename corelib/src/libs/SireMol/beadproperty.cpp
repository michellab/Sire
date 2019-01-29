/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#include "beadproperty.hpp"

using namespace SireMol;
using namespace SireStream;

namespace SireMol
{
    template class BeadProperty<QString>;
    template class BeadProperty<qint64>;
    template class BeadProperty<double>;
    template class BeadProperty<QVariant>;
}

static const RegisterMetaType<BeadStringProperty> r_beadstring;
static const RegisterMetaType<BeadIntProperty> r_beadint;
static const RegisterMetaType<BeadFloatProperty> r_beadfloat;
static const RegisterMetaType<BeadVariantProperty> r_beadvariant;


///////
/////// Implementation of BeadProp
///////

static const RegisterMetaType<BeadProp> r_beadprop(MAGIC_ONLY,
                                                 "SireMol::BeadProp");
                                                   
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const BeadProp &beadprop)
{
    writeHeader(ds, r_beadprop, 1);
    
    SharedDataStream sds(ds);
    
    sds << beadprop.bdng << static_cast<const MolViewProperty&>(beadprop);
         
    return ds;
}

/** Extract from a binary datastream */QDataStream &operator>>(QDataStream &ds, BeadProp &beadprop)
{
    VersionID v = readHeader(ds, r_beadprop);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        ds >> beadprop.bdng >> static_cast<MolViewProperty&>(beadprop);
    }
    else
        throw version_error(v, "1", r_beadprop, CODELOC);
        
    return ds;
}

BeadProp::BeadProp() : MolViewProperty()
{}

BeadProp::BeadProp(const Beading &beading) : MolViewProperty(), bdng(beading)
{}

BeadProp::BeadProp(const BeadProp &other) : MolViewProperty(other), bdng(other.bdng)
{}

BeadProp::~BeadProp()
{}

BeadProp& BeadProp::operator=(const BeadProp &other)
{
    bdng = other.bdng;
    MolViewProperty::operator=(other);
    return *this;
}

bool BeadProp::operator==(const BeadProp &other) const
{
    return bdng.read().equals(other.bdng.read()) and MolViewProperty::operator==(other);
}

bool BeadProp::operator!=(const BeadProp &other) const
{
    return not BeadProp::operator==(other);
}

/** Return the beading property used to define the beads in the molecule */
const Beading& BeadProp::beading() const
{
    return bdng.read();
}

/** Set the beading property used to define the beads in the molecule */
void BeadProp::setBeading(const Beading &beading)
{
    bdng = beading;
}

/** Internal function used to get the number of beads when the 
    contained beading function is used on the passed molecule info */
int BeadProp::getNBeads(const MoleculeInfoData &molinfo) const
{
    return bdng.read().nBeads(molinfo);
}
