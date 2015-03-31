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

#ifndef SIREMOL_ATOMCUTTING_H
#define SIREMOL_ATOMCUTTING_H

#include "cuttingfunction.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class AtomCutting;
}

QDataStream& operator<<(QDataStream&, const SireMol::AtomCutting&);
QDataStream& operator>>(QDataStream&, SireMol::AtomCutting&);

namespace SireMol
{

/** This is a cutting function that divides up a molecule into 
    CutGroups based on atoms - each atom is placed into 
    a different CutGroup
    
    @author Gaetano Calabro'
*/
class SIREMOL_EXPORT AtomCutting 
        : public SireBase::ConcreteProperty<AtomCutting,CuttingFunction>
{
public:
    AtomCutting();
    
    AtomCutting(const AtomCutting &other);
    
    ~AtomCutting();
    
    static const char* typeName();
    
    AtomCutting& operator=(const AtomCutting &other);
    
    bool operator==(const AtomCutting &other) const;
    bool operator!=(const AtomCutting &other) const;
    
    MolStructureEditor operator()(MolStructureEditor &moleditor) const;
}; 

}

Q_DECLARE_METATYPE( SireMol::AtomCutting );

SIRE_EXPOSE_CLASS( SireMol::AtomCutting )

SIRE_END_HEADER

#endif
