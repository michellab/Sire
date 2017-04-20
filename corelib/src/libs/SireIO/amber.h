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

#ifndef SIREIO_AMBER_H
#define SIREIO_AMBER_H

#include "SireBase/propertymap.h"
#include "SireBase/shareddatapointer.hpp"

#include "SireMol/atomidx.h"
#include "SireMol/molviewproperty.h"

#include "iobase.h"
#include "SireVol/space.h"


SIRE_BEGIN_HEADER

namespace SireIO
{
class Amber;
}

QDataStream& operator<<(QDataStream&, const SireIO::Amber&);
QDataStream& operator>>(QDataStream&, SireIO::Amber&);

namespace SireMol
{
class Molecules;
class MoleculeGroup;
}

namespace SireMM
{
class TwoAtomFunctions;
class ThreeAtomFunctions;
class FourAtomFunctions;
class CLJNBPairs;
}

namespace SireVol
{
class Space;
}

namespace SireIO
{
using boost::tuple;

using SireMol::Molecules;
using SireVol::SpacePtr;
using SireVol::Space;

/** This class is used to read in an AMBER top file and crd file 
    
    @author Julien Michel
*/
class SIREIO_EXPORT Amber
{

friend QDataStream& ::operator<<(QDataStream&, const SireIO::Amber&);
friend QDataStream& ::operator>>(QDataStream&, SireIO::Amber&);
  
public:
    Amber();
    Amber(const Amber &other);
    ~Amber();
  
    Amber& operator=(const Amber &other);
    
    bool operator==(const Amber &other) const;
    bool operator!=(const Amber &other) const;
  
    static const char* typeName();
  
    const char* what() const;
  
    void set14Factors(double coul_14, double lj_14);
    
    double coulomb14Factor() const;
    double lj14Factor() const;
  
    tuple<MoleculeGroup,SpacePtr> readCrdTop(const QString &crdfile,
                                             const QString &topfile,
                                             QString flag_cutting="perresidue") const;

    void writeCrd(const MoleculeGroup &mols, const Space &space, const QString &crdfile,
                  const PropertyMap &map = PropertyMap()) const;

private:
    double coul_14scl;
    double lj_14scl;
};

}

Q_DECLARE_METATYPE( SireIO::Amber )

SIRE_EXPOSE_CLASS( SireIO::Amber )

SIRE_END_HEADER

#endif

