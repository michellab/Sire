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

#ifndef SIREIO_TINKER_H
#define SIREIO_TINKER_H

#include "iobase.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class Tinker;
}

QDataStream& operator<<(QDataStream&, const SireIO::Tinker&);
QDataStream& operator>>(QDataStream&, SireIO::Tinker&);

namespace SireIO
{

/** This class holds all of the sources and default values of the
    properties and parameters used by the Tinker reader/writer
    
    @author Christopher Woods
*/
class SIREIO_EXPORT TinkerParameters : public IOParametersBase
{
public:
    TinkerParameters();
    ~TinkerParameters();
};

/** This class is used to read and write files in Tinker format

    @author Christopher Woods
*/
class SIREIO_EXPORT Tinker : public SireBase::ConcreteProperty<Tinker,IOBase>
{

friend QDataStream& ::operator<<(QDataStream&, const Tinker&);
friend QDataStream& ::operator>>(QDataStream&, Tinker&);

public:
    Tinker();
    Tinker(const Tinker &other);
    
    ~Tinker();
    
    static const char* typeName();
    
    Tinker& operator=(const Tinker &other);
    
    bool operator==(const Tinker &other) const;
    bool operator!=(const Tinker &other) const;
    
    static const TinkerParameters& parameters();
    
    void loadParameters(const QString &prmfile);
    
protected:
    MoleculeGroup readMols(const QByteArray &data,
                           const PropertyMap &map) const;
                           
    QByteArray writeMols(const MoleculeGroup &molgroup,
                         const PropertyMap &map) const;
                         
    QByteArray writeMols(const Molecules &molecules,
                         const PropertyMap &map) const;
};

}

Q_DECLARE_METATYPE( SireIO::Tinker )

SIRE_EXPOSE_CLASS( SireIO::TinkerParameters )
SIRE_EXPOSE_CLASS( SireIO::Tinker )

SIRE_END_HEADER

#endif
