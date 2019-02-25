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

#include "atomcutting.h"

#include "moleditor.h"

#include "mover.hpp"
#include "selector.hpp"

#include "molecule.h"
#include "residue.h"
#include "cutgroup.h"
#include "atom.h"

#include "atomeditor.h"
#include "cgeditor.h"
#include "reseditor.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<AtomCutting> r_rescut;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                       const AtomCutting &rescut)
{
    writeHeader(ds, r_rescut, 1);
    
    ds << static_cast<const CuttingFunction&>(rescut);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AtomCutting &rescut)
{
    VersionID v = readHeader(ds, r_rescut);
    
    if (v == 1)
    {
        ds >> static_cast<CuttingFunction&>(rescut);
    }
    else
        throw version_error(v, "1", r_rescut, CODELOC);
        
    return ds;
}

/** Constructor */
AtomCutting::AtomCutting()
               : ConcreteProperty<AtomCutting,CuttingFunction>()
{}

/** Copy constructor */
AtomCutting::AtomCutting(const AtomCutting &other)
               : ConcreteProperty<AtomCutting,CuttingFunction>(other)
{}

/** Destructor */
AtomCutting::~AtomCutting()
{}

/** Copy assignment operator */
AtomCutting& AtomCutting::operator=(const AtomCutting&)
{
    return *this;
}

/** Comparison operator */
bool AtomCutting::operator==(const AtomCutting&) const
{
    return true;
}

/** Comparison operator */
bool AtomCutting::operator!=(const AtomCutting&) const
{
    return false;
}

/** Apply this function - this creates one CutGroup per atom */
MolStructureEditor AtomCutting::operator()(MolStructureEditor &moleditor) const
{
    //remove the existing CutGroups
    moleditor.removeAllCutGroups();
    
    //now create one CutGroup for each atom, giving it the same
    //number as the atom index
    
    int k=0;
    
    for (ResIdx i(0); i<moleditor.nResidues(); ++i)
    {
        ResStructureEditor reseditor = moleditor.residue(i);
        
        for (int j=0; j<reseditor.nAtoms(); ++j)
        {
            
	  moleditor.add( CGName(QString::number(k)) );          
            
	  reseditor.atom(j).reparent( CGIdx(k) );
            
	  k++;
        }
        
        //k=reseditor.nAtoms() + 1;
        
           
     }
    
    return moleditor;
}

const char* AtomCutting::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomCutting>() );
}
