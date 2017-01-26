/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#ifndef SIREIO_AMBER2_H
#define SIREIO_AMBER2_H

#include "SireBase/propertymap.h"
#include "SireBase/shareddatapointer.hpp"

#include "SireMol/residuecutting.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class Amber2;
}

QDataStream& operator<<(QDataStream&, const SireIO::Amber2&);
QDataStream& operator>>(QDataStream&, SireIO::Amber2&);

namespace SireSystem
{
class System;
}

namespace SireIO
{

/** This class is used to read and write AMBER molecule files.
    The class will aim to support the range of file formats
    used by Amber. This class is based on the original Amber
    class written by Julien Michel
    
    @author Christopher Woods
*/
class SIREIO_EXPORT Amber2
{

friend QDataStream& ::operator<<(QDataStream&, const SireIO::Amber2&);
friend QDataStream& ::operator>>(QDataStream&, SireIO::Amber2&);
  
public:
    Amber2();
    Amber2(const Amber2 &other);
    ~Amber2();
  
    Amber2& operator=(const Amber2 &other);
    
    bool operator==(const Amber2 &other) const;
    bool operator!=(const Amber2 &other) const;
  
    static const char* typeName();
  
    const char* what() const;
  
    void set14Factors(double coul_14, double lj_14);
    
    double coulomb14Factor() const;
    double lj14Factor() const;

    SireSystem::System readRst7Parm7(const QStringList &rstlines,
                                     const QStringList &prmlines,
                                     const SireMol::CuttingFunction &cutting_function
                                             = SireMol::ResidueCutting(),
                                     const SireBase::PropertyMap &map
                                             = SireBase::PropertyMap()) const;
  
    SireSystem::System readRst7Parm7(const QString &rstfile,
                                     const QString &prmfile,
                                     const SireMol::CuttingFunction &cutting_function
                                             = SireMol::ResidueCutting(),
                                     const SireBase::PropertyMap &map
                                             = SireBase::PropertyMap()) const;
  
    SireSystem::System readRst7Parm7(const QString &rstfile,
                                     const QString &prmfile,
                                     const SireBase::PropertyMap &map,
                                     const SireMol::CuttingFunction &cutting_function
                                             = SireMol::ResidueCutting()) const;
  
    SireSystem::System readRst7Parm7(const QStringList &rstlines,
                                     const QStringList &prmlines,
                                     const SireBase::PropertyMap &map,
                                     const SireMol::CuttingFunction &cutting_function
                                             = SireMol::ResidueCutting()) const;

    void writeRst7Parm7(const SireSystem::System &system,
                        const QString &rstfile, const QString &prmfile,
                        const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

private:
    double coul_14scl;
    double lj_14scl;
};

}

Q_DECLARE_METATYPE( SireIO::Amber2 )

SIRE_EXPOSE_CLASS( SireIO::Amber2 )

SIRE_END_HEADER

#endif

