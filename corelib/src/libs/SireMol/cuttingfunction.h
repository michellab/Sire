/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMOL_CUTTINGFUNCTION_H
#define SIREMOL_CUTTINGFUNCTION_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class CuttingFunction;
}

QDataStream& operator<<(QDataStream&, const SireMol::CuttingFunction&);
QDataStream& operator>>(QDataStream&, SireMol::CuttingFunction&);

namespace SireMol
{

class Molecule;
class MolStructureEditor;

class ResidueCutting;

/** This is the base class of all cutting functions. These are
    functions that divide a molecule up into CutGroups.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT CuttingFunction : public SireBase::Property
{
public:
    CuttingFunction();
    CuttingFunction(const CuttingFunction &other);
    
    virtual ~CuttingFunction();
    
    static const char* typeName()
    {
        return "SireMol::CuttingFunction";
    }

    virtual CuttingFunction* clone() const=0;
    
    virtual Molecule operator()(const Molecule &molecule) const;
    virtual MolStructureEditor operator()(MolStructureEditor &moleditor) const=0;

    static const ResidueCutting& null();
};

typedef SireBase::PropPtr<CuttingFunction> CutFuncPtr;

}

SIRE_EXPOSE_CLASS( SireMol::CuttingFunction )

SIRE_EXPOSE_PROPERTY( SireMol::CutFuncPtr, SireMol::CuttingFunction )

SIRE_END_HEADER

/// need to include ResidueCutting header as it is the null object
#include "residuecutting.h"

#endif
