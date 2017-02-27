/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include <QMutex>

#include "cuttingfunction.h"
#include "residuecutting.h"

#include "molecule.h"
#include "mover.hpp"
#include "moleditor.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

//////////////
////////////// Implementation of CuttingFunction
//////////////

static const RegisterMetaType<CuttingFunction> r_cutfunc( MAGIC_ONLY,
                                                          "SireMol::CuttingFunction" );
                                                          
/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const CuttingFunction &cutfunc)
{
    writeHeader(ds, r_cutfunc, 1);
    
    ds << static_cast<const Property&>(cutfunc);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, CuttingFunction &cutfunc)
{
    VersionID v = readHeader(ds, r_cutfunc);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(cutfunc);
    }
    else
        throw version_error( v, "1", r_cutfunc, CODELOC );
        
    return ds;
}

/** Constructor */
CuttingFunction::CuttingFunction()
                : Property()
{}

/** Copy constructor */
CuttingFunction::CuttingFunction(const CuttingFunction &other)
                : Property(other)
{}

/** Destructor */
CuttingFunction::~CuttingFunction()
{}

/** Apply this function to a molecule */
Molecule CuttingFunction::operator()(const Molecule &molecule) const
{
    MolStructureEditor moleditor( molecule );
    
    moleditor = this->operator()(moleditor);
    
    return moleditor.commit();
}

const ResidueCutting& CuttingFunction::null()
{
    return *(create_shared_null<ResidueCutting>());
}
