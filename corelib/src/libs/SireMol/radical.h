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

#ifndef SIREMOL_RADICAL_H
#define SIREMOL_RADICAL_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Radical;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Radical&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Radical&);

namespace SireMol
{

/** This class provides information about the radical type of an
    atom (e.g. not a radical, singlet, doublet etc)

    @author Christopher Woods
*/
class SIREMOL_EXPORT Radical
    : public SireBase::ConcreteProperty<Radical,SireBase::Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Radical&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, Radical&);

public:
    Radical();

    Radical(const QString &s);
    Radical(int value);

    Radical(const Radical &other);

    ~Radical();

    Radical& operator=(const Radical &other);

    bool operator==(const Radical &other) const;
    bool operator!=(const Radical &other) const;

    static const char* typeName();

    static Radical singlet();
    static Radical doublet();
    static Radical triplet();
    static Radical undefined();

    QString toString() const;

    int value() const;
    int sdfValue() const;

    bool isDefined() const;
    bool isSinglet() const;
    bool isDoublet() const;
    bool isTriplet() const;

private:
    /** The radical type. */
    qint32 radical_type;
};

}

Q_DECLARE_METATYPE( SireMol::Radical )

SIRE_EXPOSE_CLASS( SireMol::Radical )

SIRE_END_HEADER

#endif
