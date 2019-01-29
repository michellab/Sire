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

#ifndef SIREMOL_RESIDUECUTTING_H
#define SIREMOL_RESIDUECUTTING_H

#include "cuttingfunction.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ResidueCutting;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ResidueCutting&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ResidueCutting&);

namespace SireMol
{

/** This is a cutting function that divides up a molecule into 
    CutGroups based on residue - each residue is placed into 
    a different CutGroup
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ResidueCutting 
        : public SireBase::ConcreteProperty<ResidueCutting,CuttingFunction>
{
public:
    ResidueCutting();
    
    ResidueCutting(const ResidueCutting &other);
    
    ~ResidueCutting();
    
    static const char* typeName();
    
    ResidueCutting& operator=(const ResidueCutting &other);
    
    bool operator==(const ResidueCutting &other) const;
    bool operator!=(const ResidueCutting &other) const;
    
    MolStructureEditor operator()(MolStructureEditor &moleditor) const;
}; 

}

Q_DECLARE_METATYPE( SireMol::ResidueCutting );

SIRE_EXPOSE_CLASS( SireMol::ResidueCutting )

SIRE_END_HEADER

#endif
