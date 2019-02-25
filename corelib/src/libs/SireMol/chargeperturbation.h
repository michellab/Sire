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

#ifndef SIREMOL_CHARGEPERTURBATION_H
#define SIREMOL_CHARGEPERTURBATION_H

#include "perturbation.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ChargePerturbation;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ChargePerturbation&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ChargePerturbation&);

namespace SireMol
{

/** This perturbation is used to scale charges from one value
    to another as a function of lambda
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ChargePerturbation
         : public SireBase::ConcreteProperty<ChargePerturbation,Perturbation>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ChargePerturbation&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ChargePerturbation&);

public: 
    ChargePerturbation();
    ChargePerturbation(const PropertyMap &map);
    
    ChargePerturbation(const Expression &mapping_function,
                       const PropertyMap &map = PropertyMap());
    
    ChargePerturbation(const ChargePerturbation &other);
    
    ~ChargePerturbation();
    
    static const char* typeName();
    
    ChargePerturbation& operator=(const ChargePerturbation &other);
    
    bool operator==(const ChargePerturbation &other) const;
    bool operator!=(const ChargePerturbation &other) const;

    QString toString() const;

    QSet<QString> requiredProperties() const;
    
    bool wouldChange(const Molecule &molecule, const Values &values) const;

protected:
    void perturbMolecule(MolEditor &molecule, const Values &values) const;
};

}

Q_DECLARE_METATYPE( SireMol::ChargePerturbation )

SIRE_EXPOSE_CLASS( SireMol::ChargePerturbation )

SIRE_END_HEADER

#endif
