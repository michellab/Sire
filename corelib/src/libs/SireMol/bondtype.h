/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREMOL_BONDTYPE_H
#define SIREMOL_BONDTYPE_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class BondType;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::BondType&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::BondType&);

namespace SireMol
{

/** This class represents a bond type (e.g. single, double etc.)

    @author Christopher Woods
*/
class SIREMOL_EXPORT BondType
    : public SireBase::ConcreteProperty<BondType,SireBase::Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const BondType&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, BondType&);

public:
    BondType();

    BondType(const QString &s);
    BondType(int sdf_value);

    BondType(const BondType &other);

    ~BondType();

    BondType& operator=(const BondType &other);

    bool operator==(const BondType &other) const;
    bool operator!=(const BondType &other) const;

    static const char* typeName();

    static BondType singleBond();
    static BondType doubleBond();
    static BondType tripleBond();
    static BondType aromaticBond();
    static BondType undefinedBond();

    QString toString() const;

    int value() const;
    int sdfValue() const;

    bool isDefined() const;
    bool isSingle() const;
    bool isDouble() const;
    bool isTriple() const;
    bool isAromatic() const;

private:
    /** The bond type. We use an integer in SDF format */
    qint32 bond_type;
};

}

Q_DECLARE_METATYPE( SireMol::BondType )

SIRE_EXPOSE_CLASS( SireMol::BondType )

SIRE_END_HEADER

#endif
