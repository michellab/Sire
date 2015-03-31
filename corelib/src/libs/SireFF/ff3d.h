/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREFF_FF3D_H
#define SIREFF_FF3D_H

#include "forcetable.h"
#include "energytable.h"
#include "fieldtable.h"
#include "potentialtable.h"
#include "probe.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Symbol;
}

namespace SireFF
{

/** This class provides the virtual interface for the 3D
    forcefields. These are forcefields that use 3D coordinates
    for the atoms, and therefore you can calculate 3D forces 
    on the atoms. This class provides a virtual interface,
    and should be multiply inherited with FF to be used.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FF3D
{
public:
    FF3D();
    FF3D(const FF3D &other);
    
    virtual ~FF3D();
    
    static const char* typeName()
    {
        return "SireFF::FF3D";
    }
    
    /** Calculate all of the energies of the 
        molecules in the forcetable 'forcetable' due to the
        molecules' interactions in this forcefield */
    virtual void energy(EnergyTable &energytable, double scale_energy=1)=0;

    /** Calculate all of the energies acting on all of the 
        molecules in the energytable 'energytable' due to the
        specified component of the molecules' interactions in 
        this forcefield */
    virtual void energy(EnergyTable &energytable, const SireCAS::Symbol &component,
                             double scale_energy=1)=0;

    /** Calculate all of the forces acting on all of the 
        molecules in the forcetable 'forcetable' due to the
        molecules' interactions in this forcefield */
    virtual void force(ForceTable &forcetable, double scale_force=1)=0;

    /** Calculate all of the forces acting on all of the 
        molecules in the forcetable 'forcetable' due to the
        specified component of the molecules' interactions in 
        this forcefield */
    virtual void force(ForceTable &forcetable, const SireCAS::Symbol &component,
                       double scale_force=1)=0;
                       
    /** Calculate the fields acting at all of the points
        in 'fieldtable' due to the molecules in this forcefield */
    virtual void field(FieldTable &fieldtable, double scale_field=1)=0;
    
    /** Calculate the fields acting at all of the points
        in 'fieldtable' due to the specified component 'component'
        from the molecules in this forcefield */
    virtual void field(FieldTable &fieldtable, const SireCAS::Symbol &component,
                       double scale_field=1)=0;
                       
    /** Calculate the fields acting at all of the points
        in 'fieldtable' due to the molecules in this forcefield */
    virtual void field(FieldTable &fieldtable, const Probe &probe,
                       double scale_field=1)=0;
    
    /** Calculate the fields acting at all of the points
        in 'fieldtable' due to the specified component 'component'
        from the molecules in this forcefield */
    virtual void field(FieldTable &fieldtable, const SireCAS::Symbol &component,
                       const Probe &probe, double scale_field=1)=0;
                       
    /** Calculate the potential acting at the points in 'potentialtable'
        due to the molecules in this forcefield */
    virtual void potential(PotentialTable &potentialtable, double scale_potential=1)=0;
                       
    /** Calculate the potential acting at the points in 'potentialtable'
        due to the molecules in this forcefield */
    virtual void potential(PotentialTable &potentialtable, const Probe &probe,
                           double scale_potential=1)=0;
    
    /** Calculate the potential acting at the points in 'potentialtable'
        due to the specified component 'component' from the molecules in 
        this table */
    virtual void potential(PotentialTable &potentialtable, 
                           const SireCAS::Symbol &component,
                           double scale_potential=1)=0;
    
    /** Calculate the potential acting at the points in 'potentialtable'
        due to the specified component 'component' from the molecules in 
        this table */
    virtual void potential(PotentialTable &potentialtable, 
                           const SireCAS::Symbol &component,
                           const Probe &probe,
                           double scale_potential=1)=0;
};

}

SIRE_EXPOSE_CLASS( SireFF::FF3D )

SIRE_END_HEADER

#endif
