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

#ifndef SIREMOL_STEREOSCOPY_H
#define SIREMOL_STEREOSCOPY_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Stereoscopy;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Stereoscopy&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Stereoscopy&);

namespace SireMol
{

/** This class represents a bond type (e.g. single, double etc.)

    @author Christopher Woods
*/
class SIREMOL_EXPORT Stereoscopy
    : public SireBase::ConcreteProperty<Stereoscopy,SireBase::Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Stereoscopy&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, Stereoscopy&);

public:
    Stereoscopy();

    Stereoscopy(const QString &s);
    Stereoscopy(int sdf_value);

    Stereoscopy(const Stereoscopy &other);

    ~Stereoscopy();

    Stereoscopy& operator=(const Stereoscopy &other);

    bool operator==(const Stereoscopy &other) const;
    bool operator!=(const Stereoscopy &other) const;

    static const char* typeName();

    static Stereoscopy up();
    static Stereoscopy down();
    static Stereoscopy notStereo();
    static Stereoscopy undefined();

    QString toString() const;

    int value() const;
    int sdfValue() const;

    bool isDefined() const;
    bool isUp() const;
    bool isDown() const;
    bool isNotStereo() const;

private:
    /** The stereo type. We use an integer in SDF format */
    qint32 stereo_type;
};

}

Q_DECLARE_METATYPE( SireMol::Stereoscopy )

SIRE_EXPOSE_CLASS( SireMol::Stereoscopy )

SIRE_END_HEADER

#endif
