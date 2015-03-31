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

#include "restraintff.h"

#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/moleculedata.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/viewsofmol.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"
#include "SireCAS/errors.h"
#include "SireFF/errors.h"
#include "SireBase/errors.h"

using namespace SireMM;
using namespace SireFF;
using namespace SireMol;
using namespace SireCAS;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<RestraintFF> r_restraintff;

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const RestraintFF &restraintff)
{
    writeHeader(ds, r_restraintff, 1);
    
    SharedDataStream sds(ds);
    
    sds << restraintff.restraints_by_idx
        << restraintff.spce
        << static_cast<const G1FF&>(restraintff);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, RestraintFF &restraintff)
{
    VersionID v = readHeader(ds, r_restraintff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> restraintff.restraints_by_idx
            >> restraintff.spce
            >> static_cast<G1FF&>(restraintff);

        restraintff.reindexRestraints();
        restraintff._pvt_updateName();
    }
    else
        throw version_error( v, "1", r_restraintff, CODELOC );
        
    return ds;
}

/** Force recalculation of the restraint energy from scratch */
void RestraintFF::mustNowRecalculateFromScratch()
{
    if (not recalc_from_scratch)
    {
        old_restraints_by_idx.clear();
        recalc_from_scratch = true;
    }
}

/** Internal function used to rebuild the list of available
    properties in this forcefield */
void RestraintFF::rebuildProperties()
{
    props.clear();
    
    props.setProperty( "space", spce );

    foreach (Symbol symbol, user_values.symbols())
    {
        props.setProperty( symbol.toString(), 
                           VariantProperty(user_values[symbol]) );
    }
}

/** Internal function called to reindex the restraints */
void RestraintFF::reindexRestraints()
{
    const Restraint3DPtr *restraints_array = restraints_by_idx.constData();
    quint32 nrestraints = restraints_by_idx.count();
    
    restraints_by_molnum.clear();
    user_values = Values();
    builtin_symbols.clear();
    
    for (quint32 i=0; i<nrestraints; ++i)
    {
        const Restraint3D &restraint = restraints_array[i].read();
        
        Molecules mols = restraint.molecules();
        
        for (Molecules::const_iterator it = mols.constBegin();
             it != mols.constEnd();
             ++it)
        {
            restraints_by_molnum[it->number()].append(i);
        }
    
        user_values += restraint.userValues();
        builtin_symbols += restraint.builtinSymbols();
    }
    
    this->rebuildProperties();
    
    this->mustNowRecalculateFromScratch();
}

/** Constructor */
RestraintFF::RestraintFF() : ConcreteProperty<RestraintFF,G1FF>(true), FF3D(),
                             recalc_from_scratch(true)
{
    this->_pvt_updateName();
    this->rebuildProperties();
}

/** Construct, giving the forcefield the specified name */
RestraintFF::RestraintFF(const QString &name)
            : ConcreteProperty<RestraintFF,G1FF>(true), FF3D(),
              recalc_from_scratch(true)
{
    FF::setName(name);
    this->rebuildProperties();
}

/** Copy constructor */
RestraintFF::RestraintFF(const RestraintFF &other)
            : ConcreteProperty<RestraintFF,G1FF>(other), FF3D(other),
              ffcomponents(other.ffcomponents),
              restraints_by_idx(other.restraints_by_idx),
              old_restraints_by_idx(other.old_restraints_by_idx),
              restraints_by_molnum(other.restraints_by_molnum),
              spce(other.spce),
              user_values(other.user_values),
              builtin_symbols(other.builtin_symbols),
              props(other.props),
              recalc_from_scratch(other.recalc_from_scratch)
{}

/** Destructor */
RestraintFF::~RestraintFF()
{}

const char* RestraintFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RestraintFF>() );
}

/** Copy assignment operator */
RestraintFF& RestraintFF::operator=(const RestraintFF &other)
{
    if (this != &other)
    {
        ffcomponents = other.ffcomponents;
        restraints_by_idx = other.restraints_by_idx;
        old_restraints_by_idx = other.old_restraints_by_idx;
        restraints_by_molnum = other.restraints_by_molnum;
        user_values = other.user_values;
        builtin_symbols = other.builtin_symbols;
        spce = other.spce;
        props = other.props;
        recalc_from_scratch = other.recalc_from_scratch;
        
        G1FF::operator=(other);
        FF3D::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool RestraintFF::operator==(const RestraintFF &other) const
{
    return this == &other or
           (restraints_by_idx == other.restraints_by_idx and
            G1FF::operator==(other));
}

/** Comparison operator */
bool RestraintFF::operator!=(const RestraintFF &other) const
{
    return not RestraintFF::operator==(other);
}

/** Internal function used to get the forcefield components */
const RestraintComponent& RestraintFF::_pvt_components() const
{
    return ffcomponents;
}

/** Return the components of this forcefield */
const RestraintComponent& RestraintFF::components() const
{
    return ffcomponents;
}

/** Recalculate the restraint energy */
void RestraintFF::recalculateEnergy()
{
    if (old_restraints_by_idx.isEmpty())
    {
        //recalculate the energy from scratch
        const Restraint3DPtr *restraints_array = restraints_by_idx.constData();
        const int nrestraints = restraints_by_idx.count();
        
        RestraintEnergy nrg(0);
        
        #pragma omp parallel
        {
            const Restraint3DPtr *my_restraints_array = restraints_array;
            
            RestraintEnergy my_nrg(0);
            
            #pragma omp for schedule(dynamic)
            for (int i=0; i<nrestraints; ++i)
            {
                my_nrg += RestraintEnergy( my_restraints_array[i].read().energy() );
            }
            
            #pragma omp critical
            {
                nrg += my_nrg;
            }
        }
        
        this->components().setEnergy(*this, nrg);

        recalc_from_scratch = false;
    }
    else
    {
        //evaluate the change in energy
        RestraintEnergy delta_nrg(0);
        
        const Restraint3DPtr *restraints_array = restraints_by_idx.constData();
        const quint32 nrestraints = restraints_by_idx.count();
        
        for (QHash<quint32,Restraint3DPtr>::const_iterator 
                                        it = old_restraints_by_idx.constBegin();
             it != old_restraints_by_idx.constEnd();
             ++it)
        {
            BOOST_ASSERT( it.key() < nrestraints );
            
            delta_nrg += RestraintEnergy( restraints_array[it.key()].read().energy() );
            delta_nrg -= RestraintEnergy( it.value().read().energy() );
        }
        
        this->components().changeEnergy(*this, delta_nrg);
        old_restraints_by_idx.clear();

        recalc_from_scratch = false;
    }
}

/** Update the restraints with the new molecule data in 'moldata' */
void RestraintFF::updateRestraints(const MoleculeData &moldata)
{
    //get the list of restraints potentially affected by this change in molecule
    QList<quint32> restraint_idxs = restraints_by_molnum.value(moldata.number());
    
    if (restraint_idxs.isEmpty())
        return;
        
    else if (restraint_idxs.count() == 1)
    {
        if (old_restraints_by_idx.count() >= restraints_by_idx.count() / 2)
        {
            this->mustNowRecalculateFromScratch();
        }

        quint32 idx = *(restraint_idxs.constBegin());
        
        bool save_old_state = not ( recalc_from_scratch or 
                                    old_restraints_by_idx.contains(idx) );
        
        Restraint3DPtr old_restraint;
        
        if (save_old_state)
            old_restraint = restraints_by_idx.at(idx);
        
        restraints_by_idx[idx].edit().update(moldata);
        
        if (save_old_state)
            old_restraints_by_idx.insert(idx, old_restraint);
    }
    else
    {
        std::auto_ptr<RestraintFF> old_state( this->clone() );
        
        try
        {
            if (old_restraints_by_idx.count() >= restraints_by_idx.count() / 2)
            {
                this->mustNowRecalculateFromScratch();
            }
            
            foreach (quint32 idx, restraint_idxs)
            {
                bool save_old_state = not ( recalc_from_scratch or
                                            old_restraints_by_idx.contains(idx) );
                                            
                Restraint3DPtr old_restraint;
                
                if (save_old_state)
                    old_restraint = restraints_by_idx.at(idx);
                    
                restraints_by_idx[idx].edit().update(moldata);
                
                if (save_old_state)
                    old_restraints_by_idx.insert(idx, old_restraint);
            }
        }
        catch(...)
        {
            this->copy( *old_state );
            throw;
        }
    }
}
    
/** Update the restraints with the new molecule data in 'molecules' */
void RestraintFF::updateRestraints(const Molecules &molecules)
{
    //get the list of restraints potentially affected by this change in molecules
    QList<quint32> restraint_idxs;
    
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        restraint_idxs += restraints_by_molnum.value(it.key());
    }
    
    if (restraint_idxs.isEmpty())
        return;
        
    else if (restraint_idxs.count() == 1)
    {
        if (old_restraints_by_idx.count() >= restraints_by_idx.count() / 2)
        {
            this->mustNowRecalculateFromScratch();
        }

        quint32 idx = *(restraint_idxs.constBegin());
        
        bool save_old_state = not ( recalc_from_scratch or 
                                    old_restraints_by_idx.contains(idx) );
        
        Restraint3DPtr old_restraint;
        
        if (save_old_state)
            old_restraint = restraints_by_idx.at(idx);
        
        restraints_by_idx[idx].edit().update(molecules);
        
        if (save_old_state)
            old_restraints_by_idx.insert(idx, old_restraint);
    }
    else
    {
        std::auto_ptr<RestraintFF> old_state( this->clone() );
        
        try
        {
            if (old_restraints_by_idx.count() >= restraints_by_idx.count() / 2)
            {
                this->mustNowRecalculateFromScratch();
            }

            foreach (quint32 idx, restraint_idxs)
            {
                bool save_old_state = not ( recalc_from_scratch or
                                            old_restraints_by_idx.contains(idx) );
                                            
                Restraint3DPtr old_restraint;
                
                if (save_old_state)
                    old_restraint = restraints_by_idx.at(idx);
                    
                restraints_by_idx[idx].edit().update(molecules);
                
                if (save_old_state)
                    old_restraints_by_idx.insert(idx, old_restraint);
            }
        }
        catch(...)
        {
            this->copy( *old_state );
            throw;
        }
    }
}

/** Internal function used to update the name of the forcefield */
void RestraintFF::_pvt_updateName()
{
    ffcomponents = RestraintComponent(this->name());
    G1FF::_pvt_updateName();
}

/** Internal function called when a molecule is added - this forcefield
    doesn't care about molecules being added - only restraints */
void RestraintFF::_pvt_added(const PartialMolecule&, const PropertyMap&)
{}
  
/** Internal function called when a molecule is remove - this forcefield
    will remove the first restraint it finds that contains this exact
    molecule. It will then remove that restraint (which may then
    cause other molecules to be removed) */              
void RestraintFF::_pvt_removed(const PartialMolecule &mol)
{
    MolNum molnum = mol.number();
    
    foreach (quint32 idx, restraints_by_molnum.value(molnum))
    {
        if (restraints_by_idx.at(idx).read().molecules().contains(mol))
        {
            //this restraint contains this view - we need to remove
            //this restraint, which will remove all molecule views
            //involved with this restraint. *However* we need to ensure
            //that only one copy of this view is removed, so we need
            //to re-add it to this forcefield so that it is available
            //to be removed when this restraint is removed
            this->add(mol);
            this->removeRestraintAt(idx);
            return;
        }
    }
}

/** This internal function is called when the molecule 'mol' has
    changed - any restraints which depend on this molecule will
    need to be updated */
void RestraintFF::_pvt_changed(const SireMol::Molecule &mol, bool auto_update)
{
    this->updateRestraints( mol.data() );
}

/** This internal function is called when the molecules in 'mols'
    have changed - any restraints depending on these molecules
    will need to be updated */
void RestraintFF::_pvt_changed(const QList<SireMol::Molecule> &mols, bool auto_update)
{
    this->updateRestraints( Molecules(mols) );
}

/** This internal function is called when all molecules (and hence 
    all restraints) are removed */
void RestraintFF::_pvt_removedAll()
{
    restraints_by_idx.clear();
    old_restraints_by_idx.clear();
    restraints_by_molnum.clear();
    user_values = Values();
    builtin_symbols.clear();
    recalc_from_scratch = true;
}

/** This internal function is called to see if changing properties would
    change the forcefield - as we can't change properties of restraints,
    this cannot change the forcefield */
bool RestraintFF::_pvt_wouldChangeProperties(SireMol::MolNum, 
                                             const PropertyMap&) const
{
    return false;
}

/** Internal function called when a molecule is added - this forcefield
    doesn't care about molecules being added - only restraints */
void RestraintFF::_pvt_added(const ViewsOfMol&, const PropertyMap&)
{}

/** Internal function called to remove all of the views in 'mol'. 
    Note that this only removes the first restraints that contain
    each view, and removing the restraint will then remove *all*
    molecules involved in this restraint */
void RestraintFF::_pvt_removed(const ViewsOfMol &mol)
{
    for (int i=0; i<mol.nViews(); ++i)
    {
        this->_pvt_removed(mol[i]);
    }
}

/** Internal function called when all copies of the view of the molecule
    in 'mol' have been removed. This removes all restraints that
    involve this view, and also removes all molecules that are 
    in those restraints */
void RestraintFF::_pvt_removedAll(const PartialMolecule &mol)
{
    MolNum molnum = mol.number();
    
    while (true)
    {
        bool more_to_remove = false;
    
        foreach (quint32 idx, restraints_by_molnum.value(molnum))
        {
            if (restraints_by_idx.at(idx).read().molecules().contains(mol))
            {
                //this restraint contains this view - we need to remove
                //this restraint, which will remove all molecule views
                //involved with this restraint
                this->removeRestraintAt(idx);
                more_to_remove = true;
                break;
            }
        }
        
        if (not more_to_remove)
            break;
    }
}

/** Internal function called to remove all restraints that involve
    the views in 'mol' - this removes all of the restraints, and also
    all of the molecules used in those restraints */
void RestraintFF::_pvt_removedAll(const ViewsOfMol &mol)
{
    for (int i=0; i<mol.nViews(); ++i)
    {
        this->_pvt_removedAll( mol[i] );
    }
}

/** Return the space used by this forcefield */
const Space& RestraintFF::space() const
{
    return spce.read();
}

/** Set the space used by all of the restraints in this forcefield. 
    
    \throw SireVol::incompatible_space
*/
bool RestraintFF::setSpace(const Space &space)
{
    if (not space.equals(spce))
    {
        QVector<Restraint3DPtr> new_restraints = restraints_by_idx;
        
        for (int i=0; i<new_restraints.count(); ++i)
        {
            new_restraints[i].edit().setSpace(space);
        }
        
        restraints_by_idx = new_restraints;
        spce = space;
        
        props.setProperty("space", spce);
        
        return true;
    }
    else
        return false;
}

/** Set the value of the user symbol 'symbol' to the value 'value'.
    This will only work if there is a restraint in this forcefield
    that has this symbol. This returns whether or not this
    changes the forcefield. This raises an exception if you
    are trying to set the value of a built-in symbol
    
    \throw SireError::invalid_arg
*/
bool RestraintFF::setValue(const Symbol &symbol, double value)
{
    if (builtin_symbols.contains(symbol))
        throw SireError::invalid_arg( QObject::tr(
            "You cannot set the value of the symbol %1 to %2 as this "
            "is one of the built-in symbols (%3) of the restraints in "
            "this forcefield (%4).")
                .arg(symbol.toString()).arg(value)
                .arg( Sire::toString(builtin_symbols) )
                .arg(this->toString()), CODELOC );

    if (user_values.contains(symbol))
    {
        if (user_values[symbol] == value)
            return false;
            
        QVector<Restraint3DPtr> new_restraints = restraints_by_idx;
        
        for (int i=0; i<new_restraints.count(); ++i)
        {
            if (new_restraints.at(i).read().hasValue(symbol))
                new_restraints[i].edit().setValue(symbol, value);
        }
        
        user_values.set(symbol, value);
        
        restraints_by_idx = new_restraints;
        
        props.setProperty( symbol.toString(), VariantProperty(value) );
        
        return true;
    }
    else
        return false;
}

/** Return the value of the user-supplied symbol 'symbol' */
double RestraintFF::getValue(const Symbol &symbol) const
{
    if (not user_values.contains(symbol))
        throw SireCAS::missing_symbol( QObject::tr(
            "There is no user-supplied value of the symbol %1 in the "
            "forcefield %2.")
                .arg(symbol.toString(), this->toString()), CODELOC );

    return user_values[symbol];
}

/** Return whether or not there is a user-supplied value with the 
    symbol 'symbol' */
bool RestraintFF::hasValue(const Symbol &symbol) const
{
    return user_values.contains(symbol);
}

/** Set the property 'name' to the value 'property'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
bool RestraintFF::setProperty(const QString &name, const Property &property)
{
    if (name == QLatin1String("space"))
    {
        return this->setSpace( property.asA<Space>() );
    }
    else if (props.hasProperty(name))
    {
        return this->setValue( Symbol(name), property.asA<VariantProperty>()
                                                     .convertTo<double>() );
    }
    else
    {
        throw SireBase::missing_property( QObject::tr(
                "There is no property called \"%1\" in the forcefield %2, "
                "so it cannot be set to the value %3. Available properties "
                "are [ %4 ].")
                    .arg(name, this->toString(), property.toString())
                    .arg( Sire::toString(props.propertyKeys()) ), CODELOC );
                    
        return false;
    }
}

/** Return the property called 'name' 

    \throw SireBase::missing_property
*/
const Property& RestraintFF::property(const QString &name) const
{
    return props.property(name);
}

/** Return whether or not this forcefield contains a property called 'name' */
bool RestraintFF::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return all of the properties of this forcefield */
const Properties& RestraintFF::properties() const
{
    return props;
}

/** Return all of the user-supplied symbols for the restraints 
    in this forcefield */
Symbols RestraintFF::userSymbols() const
{
    return user_values.symbols();
}

/** Return all of the built-in symbols used by the restraints
    in this forcefield */
Symbols RestraintFF::builtinSymbols() const
{
    return builtin_symbols;
}

/** Return all of the symbols used in this forcefield - this includes
    both the user-supplied symbols and the built-in symbols */
Symbols RestraintFF::symbols() const
{
    return userSymbols() + builtinSymbols();
}

/** Return all of the values for the user-supplied values in the
    restraints in this forcefield */
Values RestraintFF::userValues() const
{
    return user_values;
}

/** Return a copy of this forcefield where all of the restraints
    have been differentiated with respect to 'symbol'. The returned
    forcefield will contain the same molecules in the same state
    as they are in this forcefield, and will be called
    d(forcefield_name)/d(symbol)
*/
RestraintFF RestraintFF::differentiate(const Symbol &symbol) const
{
    RestraintFF diff( *this );
    
    diff.restraints_by_idx.clear();
    
    for (int i=0; i<restraints_by_idx.count(); ++i)
    {
        Restraint3DPtr diff_restraint = restraints_by_idx.at(i).read()
                                                               .differentiate(symbol);
    
        if (not diff_restraint.isNull())
            diff.restraints_by_idx.append(diff_restraint);
    }
    
    diff.reindexRestraints();
    
    return diff;
}

/** Return whether or not this forcefield contains the restraint
    'restraint' */
bool RestraintFF::contains(const Restraint3D &restraint) const
{
    for (int i=0; i<restraints_by_idx.count(); ++i)
    {
        if (restraints_by_idx.at(i).read().equals(restraint))
            return true;
    }
    
    return false;
}

/** Add the passed restraint to this forcefield. This does nothing
    if this restraint is already part of this forcefield. This returns
    whether or not this changes this forcefield */
bool RestraintFF::add(const Restraint3D &restraint)
{
    //make sure that there are no current user values that clash
    //with the builtin symbols of this restraint
    const Symbols user_symbols = user_values.symbols();
    
    if (restraint.builtinSymbols().intersects(user_symbols))
        throw SireError::incompatible_error( QObject::tr(
            "Cannot add the restraint %1 to the forcefield %2 as there "
            "is a clash between the built-in symbols of this restraint (%3) "
            "and the user-supplied symbols of this forcefield (%4).")
                .arg(restraint.toString(), this->toString())
                .arg( Sire::toString(restraint.builtinSymbols()) )
                .arg( Sire::toString(user_symbols) ), CODELOC );

    Restraint3DPtr new_restraint( restraint );

    //update the restraint with the current space
    if (not new_restraint.read().space().equals(spce))
    {
        new_restraint.edit().setSpace(spce);
    }
    
    //update the restraint with the current symbol values
    foreach (Symbol symbol, user_symbols)
    {
        new_restraint.edit().setValue(symbol, user_values[symbol]);
    }
    
    //update the restraint with the current versions of the molecules
    if (not this->isEmpty())
    {
        new_restraint.edit().update( this->molecules() );
    }
    
    //now see if we have this restraint already?
    if (not this->contains(new_restraint))
    {
        //we don't - add the new restraint
        std::auto_ptr<RestraintFF> old_state( this->clone() );
        
        try
        {
            //add the restraint onto the list of existing restraints...
            restraints_by_idx.append( new_restraint );
            
            //now add the molecules in this restraint
            this->add( new_restraint.read().molecules() );
            
            //finally, reindex the restraints
            this->reindexRestraints();
            
            return true;
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
    }
    
    return false;
}

/** Return the array of all restraints */
QVector<Restraint3DPtr> RestraintFF::restraints() const
{
    return restraints_by_idx;
}

/** Return the number of restraints in this forcefield */
int RestraintFF::nRestraints() const
{
    return restraints_by_idx.count();
}

/** Return the ith restraint

    \throw SireError::index_error
*/
const Restraint3D& RestraintFF::restraintAt(int i) const
{
    return restraints_by_idx.at( Index(i).map(this->nRestraints()) );
}

/** Remove the ith restraint

    \throw SireError::invalid_index
*/
void RestraintFF::removeRestraintAt(int i)
{
    i = Index(i).map( this->nRestraints() );
    
    std::auto_ptr<RestraintFF> old_state( this->clone() );
    
    try
    {
        //get the molecules in this restraint
        Molecules mols = restraints_by_idx.at(i).read().molecules();
    
        //lose this restraint from the list
        restraints_by_idx.remove(i);
        
        //reindex the restraints before removing the molecules so as
        //to ensure that we don't trigger infinite recursion as 
        //removing more molecules causes the restraint to be removed again
        this->reindexRestraints();
        
        this->remove(mols);
    }
    catch(...)
    {
        this->copy( *old_state );
        throw;
    }
}

/** Remove the restraint 'restraint' from this forcefield. This does
    nothing if this restraint is not in this forcefield. This returns
    whether or not the restraint was removed. */
bool RestraintFF::remove(const Restraint3D &restraint)
{
    for (int i=0; i<restraints_by_idx.count(); ++i)
    {
        if (restraints_by_idx.at(i).read().equals(restraint))
        {
            this->removeRestraintAt(i);
            return true;
        }
    }
    
    return false;
}

void RestraintFF::energy(EnergyTable &energytable, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "RestraintFF does not yet support energy calculations!"), CODELOC );
}

void RestraintFF::energy(EnergyTable &energytable, const Symbol &symbol,
                           double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "RestraintFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the forces on the molecules in 'forcetable' caused
    by the restraints in this forcefield and add them onto the
    forcetable, optionally scaled by 'scale_force' */
void RestraintFF::force(ForceTable &forcetable, double scale_force)
{
    for (int i=0; i<restraints_by_idx.count(); ++i)
    {
        restraints_by_idx.at(i).read().force(forcetable, scale_force);
    }
}

/** Calculate the forces on the molecules in 'forcetable' caused by 
    the energy component 'symbol' in this forcefield, and add them
    onto the forcetable, optionally scaled by 'scale_force'
    
    \throw SireFF::missing_component
*/
void RestraintFF::force(ForceTable &forcetable, const Symbol &symbol, double scale_force)
{
    if (symbol == ffcomponents.total())
        this->force(forcetable, scale_force);
    else
        throw SireFF::missing_component( QObject::tr(
            "There is no forcefield component represented by %1 in the "
            "forcefield %2. Available components are %3.")
                .arg(symbol.toString(), this->toString(), 
                     ffcomponents.total().toString()), CODELOC );
}

RestraintFF* RestraintFF::clone() const
{
    return new RestraintFF(*this);
}
               
void RestraintFF::field(FieldTable &fieldtable, double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate the field of a restraint has not "
                "been written."), CODELOC );
}

void RestraintFF::field(FieldTable &fieldtable, const Symbol &component,
                        double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate the field of a restraint has not "
                "been written."), CODELOC );
}

void RestraintFF::potential(PotentialTable &potentialtable, double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate the field of a restraint has not "
                "been written."), CODELOC );
}

void RestraintFF::potential(PotentialTable &potentialtable, const Symbol &component,
                            double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate the field of a restraint has not "
                "been written."), CODELOC );
}

void RestraintFF::field(FieldTable &fieldtable, const Probe &probe, double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate the field of a restraint has not "
                "been written."), CODELOC );
}

void RestraintFF::field(FieldTable &fieldtable, const Symbol &component,
                        const Probe &probe, double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate the field of a restraint has not "
                "been written."), CODELOC );
}

void RestraintFF::potential(PotentialTable &potentialtable, const Probe &probe,
                            double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate the field of a restraint has not "
                "been written."), CODELOC );
}

void RestraintFF::potential(PotentialTable &potentialtable, const Symbol &component,
                            const Probe &probe, double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate the field of a restraint has not "
                "been written."), CODELOC );
}
