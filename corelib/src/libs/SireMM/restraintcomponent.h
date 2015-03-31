/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREMM_RESTRAINTCOMPONENT_H
#define SIREMM_RESTRAINTCOMPONENT_H

#include "SireFF/ffcomponent.h"
#include "SireFF/ff.h"

namespace SireMM
{
class RestraintComponent;
}

QDataStream& operator<<(QDataStream&, const SireMM::RestraintComponent&);
QDataStream& operator>>(QDataStream&, SireMM::RestraintComponent&);

namespace SireMM
{

using SireFF::FF;
using SireFF::FFName;

typedef SireFF::ComponentEnergy<RestraintComponent> RestraintEnergy;

/** This class represents a restraint component of a forcefield */
class SIREMM_EXPORT RestraintComponent : public SireFF::FFComponent
{
public:
    RestraintComponent(const FFName &ffname = FFName());
    RestraintComponent(const FFName &ffname, const QString &suffix);
    
    RestraintComponent(const SireCAS::Symbol &symbol);
    
    RestraintComponent(const RestraintComponent &other);
    
    ~RestraintComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return RestraintComponent::typeName();
    }
    
    RestraintComponent* clone() const
    {
        return new RestraintComponent(*this);
    }
    
    const RestraintComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const RestraintEnergy &nrg) const;
    void changeEnergy(FF &ff, const RestraintEnergy &nrg) const;
    
    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

}

Q_DECLARE_METATYPE( SireMM::RestraintComponent )

SIRE_EXPOSE_CLASS( SireMM::RestraintComponent )

#endif
