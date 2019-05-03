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

#include "spacewrapper.h"
#include "system.h"
#include "delta.h"

#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/moleditor.h"

#include "SireVol/space.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireFF;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<SpaceWrapper> r_spacewrapper;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                          const SpaceWrapper &spacewrapper)
{
    writeHeader(ds, r_spacewrapper, 1);
    
    SharedDataStream sds(ds);
    
    sds << spacewrapper.wrap_point << spacewrapper.molgroup
        << spacewrapper.map
        << static_cast<const MoleculeConstraint&>(spacewrapper);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SpaceWrapper &spacewrapper)
{
    VersionID v = readHeader(ds, r_spacewrapper);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        spacewrapper = SpaceWrapper();
        
        sds >> spacewrapper.wrap_point >> spacewrapper.molgroup
            >> spacewrapper.map
            >> static_cast<MoleculeConstraint&>(spacewrapper);
            
        spacewrapper.space_property = spacewrapper.map["space"];
        spacewrapper.coords_property = spacewrapper.map["coordinates"];
    }
    else
        throw version_error(v, "1", r_spacewrapper, CODELOC);
        
    return ds;
}

/** Constructor */
SpaceWrapper::SpaceWrapper() : ConcreteProperty<SpaceWrapper,MoleculeConstraint>()
{}

/** Construct to wrap all of the molecules in the group 'molgroup' 
    into the same periodic box as the point 'point' using the
    supplied property map to find the space and coordinate properties */
SpaceWrapper::SpaceWrapper(const PointRef &point,
                           const MoleculeGroup &wrap_group,
                           const PropertyMap &wrap_map)
             : ConcreteProperty<SpaceWrapper,MoleculeConstraint>(),
               wrap_point(point), molgroup(wrap_group), map(wrap_map),
               space_property(wrap_map["space"]),
               coords_property(wrap_map["coordinates"])
{}

/** Copy constructor */
SpaceWrapper::SpaceWrapper(const SpaceWrapper &other)
             : ConcreteProperty<SpaceWrapper,MoleculeConstraint>(other),
               wrap_point(other.wrap_point), molgroup(other.molgroup),
               map(other.map), space_property(other.space_property),
               coords_property(other.coords_property),
               spce(other.spce), changed_mols(other.changed_mols)
{}

/** Destructor */
SpaceWrapper::~SpaceWrapper()
{}

/** Copy assignment operator */
SpaceWrapper& SpaceWrapper::operator=(const SpaceWrapper &other)
{
    if (this != &other)
    {
        wrap_point = other.wrap_point;
        molgroup = other.molgroup;
        map = other.map;
        space_property = other.space_property;
        coords_property = other.coords_property;
        spce = other.spce;
        changed_mols = other.changed_mols;
        MoleculeConstraint::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool SpaceWrapper::operator==(const SpaceWrapper &other) const
{
    return (this == &other) or
           (wrap_point == other.wrap_point and
            map == other.map and
            molgroup == other.molgroup and
            MoleculeConstraint::operator==(other));
}

/** Comparison operator */
bool SpaceWrapper::operator!=(const SpaceWrapper &other) const
{
    return not this->operator==(other);
}

const char* SpaceWrapper::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SpaceWrapper>() );
}

/** Return the point that defines the center of the periodic box
    (the molecules will be wrapped so that they are in the same
    box as this point) */
const Point& SpaceWrapper::point() const
{
    return wrap_point.read();
}

/** Return the molecule group containing the molecules being wrapped */
const MoleculeGroup& SpaceWrapper::moleculeGroup() const
{
    return molgroup.read();
}

/** Return the property map used to find the coordinates and 
    space properties */
const PropertyMap& SpaceWrapper::propertyMap() const
{
    return map;
}

/** Set the baseline system for the constraint - this is 
    used to pre-calculate everything for the system
    and to check if the constraint is satisfied */
void SpaceWrapper::setSystem(const System &system)
{
    if (Constraint::wasLastSystem(system) and Constraint::wasLastSubVersion(system))
        return;

    if (wrap_point.read().usesMoleculesIn(system))
        wrap_point.edit().update(system);

    if (not molgroup.isNull())
    {
        if (system.contains(molgroup.read().number()))
            molgroup = system[molgroup.read().number()];
        else
            molgroup.edit().update(system.molecules());
    }

    if (not system.containsProperty(space_property))
        spce = SpacePtr();
    else
    {
        const Property &new_space = system.property(space_property);

        if (not new_space.isA<Space>())
            throw SireError::incompatible_error( QObject::tr(
                    "You cannot use a SpaceWrapper constraint with a "
                    "system property (%1) that is not derived from "
                    "Space (%2).")
                        .arg(space_property.toString())
                        .arg(new_space.toString()), CODELOC );
    
        spce = new_space.asA<Space>();
    }

    if (molgroup.isNull() or molgroup.read().isEmpty() or 
        spce.isNull() or (not spce.read().isPeriodic()))
    {
        Constraint::setSatisfied(system, true);
        return;
    }

    const Vector &center_point = wrap_point.read().point();
    
    const Molecules &molecules = molgroup.read().molecules();
    
    const Space &space = spce.read();

    changed_mols = Molecules();
    
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        Molecule molecule = it->molecule();
    
        const AtomCoords &coords = molecule.property(coords_property)
                                           .asA<AtomCoords>();
        
        //translate the molecule as a single entity (we don't want
        //molecules splitting over two sides of the box)
        CoordGroupArray new_coords = space.getMinimumImage(coords.array(),
                                                           center_point, true);
        
        if (new_coords.constData() != coords.constData())
        {
            //the molecule has moved
            molecule = molecule.edit().setProperty(coords_property,
                                                   AtomCoords(new_coords))
                                      .commit();
                                      
            changed_mols.add(molecule);
        }
    }

    Constraint::setSatisfied(system, changed_mols.isEmpty());
}

/** Return whether or not the changes in the passed
    delta *may* have changed the system since the last
    subversion 'subversion' */
bool SpaceWrapper::mayChange(const Delta &delta, quint32 last_subversion) const
{
    if (molgroup.isNull())
        return false;

    else if (not changed_mols.isEmpty())
        return true;

    else if (delta.sinceChanged(space_property, last_subversion))
    {
        if (delta.deltaSystem().containsProperty(space_property))
        {
            const Property &new_space = delta.deltaSystem().property(space_property);
            
            if (not new_space.isA<Space>())
                throw SireError::incompatible_error( QObject::tr(
                        "You cannot use a SpaceWrapper constraint with a "
                        "system property (%1) that is not derived from "
                        "Space (%2).")
                            .arg(space_property.toString())
                            .arg(new_space.toString()), CODELOC );
                            
            if (spce.isNull())
                return true;

            else if (not spce.read().equals(new_space))
                return true;
        }
    }
    
    return delta.sinceChanged(molgroup.read().molecules(), last_subversion) or
           delta.sinceChanged(wrap_point.read(), last_subversion);
}

/** Fully apply this constraint on the passed delta - this returns
    whether or not this constraint affects the delta */
bool SpaceWrapper::fullApply(Delta &delta)
{
    this->setSystem(delta.deltaSystem());

    if (changed_mols.isEmpty())
        return false;
    else
    {
        bool changed = delta.update(changed_mols);
        
        if (delta.deltaSystem().contains(molgroup.read().number()))
            molgroup = delta.deltaSystem()[molgroup.read().number()];
        
        changed_mols = Molecules();
        return changed;
    }
}

/** Apply this constraint based on the delta, knowing that the 
    last application of this constraint was on this system, 
    at subversion number last_subversion */
bool SpaceWrapper::deltaApply(Delta &delta, quint32 last_subversion)
{
    const System &system = delta.deltaSystem();

    if (delta.sinceChanged(space_property, last_subversion))
    {
        if (system.containsProperty(space_property))
        {
            const Property &new_space = system.property(space_property);
        
            if (new_space.isA<Space>())
            {
                if (spce.isNull())
                    return this->fullApply(delta);
                
                else if (not spce.read().equals(new_space))
                    return this->fullApply(delta);
            }
            else
                throw SireError::incompatible_error( QObject::tr(
                        "You cannot use a SpaceWrapper constraint with a "
                        "system property (%1) that is not derived from "
                        "Space (%2).")
                            .arg(space_property.toString())
                            .arg(new_space.toString()), CODELOC );
        }
        else if (not spce.isNull())
        {
            return this->fullApply(delta);
        }
    }
    
    if (molgroup.isNull() or spce.isNull() or (not spce.read().isPeriodic()))
    {
        return false;
    }
    
    if (delta.sinceChanged(wrap_point.read(), last_subversion))
        return this->fullApply(delta);

    const Space &space = spce.read();
    const Vector &center_point = wrap_point.read().point();

    if (delta.hasMoleculeChangeSince(last_subversion))
    {
        QList<MolNum> changed_molnums = delta.changedMoleculesSince(
                                                    molgroup.read().molecules(),
                                                    last_subversion);

        if (not changed_molnums.isEmpty())
        {
            if (system.contains(molgroup.read().number()))
                molgroup = system[molgroup.read().number()];
            else
                molgroup.edit().update(system.molecules());

            //some of the molecules have changed - see if they need
            //further changes to maintain the constraint
            const Molecules &molecules = molgroup.read().molecules();

            changed_mols = Molecules();
        
            foreach (MolNum molnum, changed_molnums)
            {
                Molecule molecule = molecules[molnum].molecule();
    
                const AtomCoords &coords = molecule.property(coords_property)
                                                   .asA<AtomCoords>();
                                                 
                //translate the molecule as a single entity (we don't want
                //molecules splitting over two sides of the box)
                CoordGroupArray new_coords = space.getMinimumImage(coords.array(),
                                                                   center_point, true);
                                                           
                if (new_coords.constData() != coords.constData())
                {
                    //the molecule has moved
                    molecule = molecule.edit().setProperty(coords_property,
                                                           AtomCoords(new_coords))
                                              .commit();
                                      
                    changed_mols.add(molecule);
                }
            }
    
            if (not changed_mols.isEmpty())
            {
                bool changed = delta.update(changed_mols);
        
                if (delta.deltaSystem().contains(molgroup.read().number()))
                    molgroup = delta.deltaSystem()[molgroup.read().number()];
        
                changed_mols = Molecules();
                return changed;
            }
        }
    }

    return false;
}
