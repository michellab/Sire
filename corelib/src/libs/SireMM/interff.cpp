/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include "interff.h"
#include "cljshiftfunction.h"
#include "cljcalculator.h"

#include "SireBase/booleanproperty.h"
#include "SireBase/lengthproperty.h"

#include "SireError/errors.h"
#include "SireBase/errors.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/molresid.h"
#include "SireMol/residue.h"
#include "SireMol/atomselection.h"
#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QElapsedTimer>
#include <QDebug>

using namespace SireMM;
using namespace SireMol;
using namespace SireFF;
using namespace SireBase;
using namespace SireStream;

namespace SireMM
{
    namespace detail
    {
        class InterFFData : public QSharedData
        {
        public:
            InterFFData() : QSharedData(), fixed_only(false),
                            parallel_calc(true), repro_sum(false)
            {}
            
            InterFFData(const InterFFData &other)
                 : QSharedData(),
                   cljfuncs(other.cljfuncs),
                   fixed_atoms(other.fixed_atoms),
                   cljcomps(other.cljcomps),
                   props(other.props),
                   fixed_only(other.fixed_only),
                   parallel_calc(other.parallel_calc),
                   repro_sum(other.repro_sum)
            {}
            
            ~InterFFData()
            {}

            /** The function(s) used to calculate energies */
            QVector<CLJFunctionPtr> cljfuncs;
            
            /** All of the fixed atoms, duplicated for each 
                of the CLJFunctions */
            QVector<CLJGrid> fixed_atoms;

            /** The energy components available for this forcefield */
            MultiCLJComponent cljcomps;
            
            /** All of the properties in this forcefield */
            Properties props;
            
            /** Whether or not to only calculate the energy with
                the fixed atoms */
            bool fixed_only;
            
            /** Whether or not to calculate energies in parallel */
            bool parallel_calc;
            
            /** Whether or not to sum energies using a reproducible sum */
            bool repro_sum;
        };
    }
}

static RegisterMetaType<InterFF> r_interff;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const InterFF &interff)
{
    writeHeader(ds, r_interff, 3);
    
    SharedDataStream sds(ds);
    
    sds << interff.cljgroup << interff.needs_accepting
        << interff.d->cljfuncs
        << interff.d->cljcomps
        << interff.d->fixed_atoms
        << interff.d->fixed_only << interff.d->parallel_calc
        << interff.d->repro_sum
        << static_cast<const G1FF&>(interff);

    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, InterFF &interff)
{
    VersionID v = readHeader(ds, r_interff);
    
    if (v == 3)
    {
        SharedDataStream sds(ds);
        
        sds >> interff.cljgroup >> interff.needs_accepting
            >> interff.d->cljfuncs >> interff.d->cljcomps
            >> interff.d->fixed_atoms
            >> interff.d->fixed_only >> interff.d->parallel_calc
            >> interff.d->repro_sum
            >> static_cast<G1FF&>(interff);
        
        interff.rebuildProps();
        interff._pvt_updateName();
    }
    else
        throw version_error(v, "3", r_interff, CODELOC);
    
    return ds;
}

/** Constructor */
InterFF::InterFF()
        : ConcreteProperty<InterFF,G1FF>(),
          cljgroup(CLJExtractor::EXTRACT_BY_CUTGROUP), needs_accepting(false)
{
    d = new detail::InterFFData();
    this->_pvt_updateName();
    this->setCLJFunction( CLJShiftFunction::defaultShiftFunction() );
}

/** Construct, specifying the name of the forcefield */
InterFF::InterFF(const QString &name)
        : ConcreteProperty<InterFF, G1FF>(),
          cljgroup(CLJExtractor::EXTRACT_BY_CUTGROUP), needs_accepting(false)
{
    d = new detail::InterFFData();
    G1FF::setName(name);
    this->setCLJFunction( CLJShiftFunction::defaultShiftFunction() );
}

/** Copy constructor */
InterFF::InterFF(const InterFF &other)
        : ConcreteProperty<InterFF,G1FF>(other),
          cljgroup(other.cljgroup), d(other.d),
          needs_accepting(other.needs_accepting)
{}

/** Destructor */
InterFF::~InterFF()
{}

/** Function used to set the CLJFunction used to calculate the energy */
void InterFF::setCLJFunction(const CLJFunction &func)
{
    if (d.constData()->cljfuncs.isEmpty())
    {
        d->cljfuncs.append(func);
        d->fixed_atoms.append( CLJGrid() );
        d->fixed_atoms[0].setCLJFunction(func);
        rebuildProps();
        this->mustNowRecalculateFromScratch();
    }
    else if (not d.constData()->cljfuncs.at(0).read().equals(func))
    {
        d->fixed_atoms[0].setCLJFunction(func);
        d->cljfuncs[0] = func;
        rebuildProps();
        this->mustNowRecalculateFromScratch();
    }
}

/** Return the function used to calculate the energy */
const CLJFunction& InterFF::cljFunction() const
{
    if (d.constData()->cljfuncs.isEmpty())
        throw SireError::program_bug( QObject::tr(
                "There should always be at least one CLJFunction in InterFF!"),
                    CODELOC );
    
    return d.constData()->cljfuncs.at(0).read();
}

/** Set the CLJFunction with key 'key' equal to 'cljfunc' */
void InterFF::setCLJFunction(QString key, const CLJFunction &cljfunc)
{
    if (key == "default")
    {
        this->setCLJFunction(cljfunc);
    }
    else
    {
        int idx = d->cljcomps.add(key);
        
        if (idx >= d->cljfuncs.count())
        {
            d->cljfuncs.resize( idx + 1 );
        }
        
        //duplicate the fixed atoms so that we have a set for this
        //CLJFunction
        while (d->fixed_atoms.count() < d->cljfuncs.count())
        {
            d->fixed_atoms.append( d->fixed_atoms.last() );
        }
        
        d->cljfuncs[idx] = cljfunc;
        d->fixed_atoms[idx].setCLJFunction(cljfunc);

        rebuildProps();

        this->mustNowRecalculateFromScratch();
    }
}

/** Remove the CLJ function with key 'key' - note that you cannot remove
    the 'default' CLJ function */
void InterFF::removeCLJFunctionAt(QString key)
{
    if (key != "default")
    {
        int idx = d->cljcomps.remove(key);
        
        if (idx > 0)
        {
            d->cljfuncs.removeAt(idx);
            d->fixed_atoms.removeAt(idx);
            
            rebuildProps();
            this->mustNowRecalculateFromScratch();
        }
    }
}

/** Function to remove all of the CLJFunctions (except for the "default" function) */
void InterFF::removeAllCLJFunctions()
{
    if (d->cljfuncs.count() < 2)
        return;

    d->cljcomps.removeAll();
    
    while (d->cljfuncs.count() > 1)
    {
        d->cljfuncs.removeLast();
    }
    
    while (d->fixed_atoms.count() > 1)
    {
        d->fixed_atoms.removeLast();
    }
    
    rebuildProps();
    this->mustNowRecalculateFromScratch();
}

/** Return the keys of all CLJFunctions added to this forcefield */
QStringList InterFF::cljFunctionKeys() const
{
    return d.constData()->cljcomps.keys();
}

/** Return the CLJFunction associated with the passed key */
const CLJFunction& InterFF::cljFunction(QString key) const
{
    return d.constData()->cljfuncs.at( d.constData()->cljcomps.indexOf(key) ).read();
}

/** Return the number of CLJ functions in this forcefield. There should always
    be at least one */
int InterFF::nCLJFunctions() const
{
    return d.constData()->cljfuncs.count();
}

/** Return the hash of all CLJFunctions in this forcefield, indexed by their key */
QHash<QString,CLJFunctionPtr> InterFF::cljFunctions() const
{
    QHash<QString,CLJFunctionPtr> funcs;
    
    foreach (QString key, d->cljcomps.keys())
    {
        funcs.insert( key, this->cljFunction(key) );
    }
    
    return funcs;
}

/** Internal function called when the name of the forcefield changes */
void InterFF::_pvt_updateName()
{
    d->cljcomps = d->cljcomps.rename(this->name());
    G1FF::_pvt_updateName();
}

const char* InterFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<InterFF>() );
}

const char* InterFF::what() const
{
    return InterFF::typeName();
}

/** Copy assignment operator */
InterFF& InterFF::operator=(const InterFF &other)
{
    if (this != &other)
    {
        cljgroup = other.cljgroup;
        needs_accepting = other.needs_accepting;
        d = other.d;
        G1FF::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool InterFF::operator==(const InterFF &other) const
{
    return (this == &other) or
           (G1FF::operator==(other) and d->fixed_atoms == other.d->fixed_atoms and
            d->cljfuncs == other.d->cljfuncs and d->cljcomps == other.d->cljcomps and
            cljgroup == other.cljgroup and needs_accepting == other.needs_accepting);
}

/** Comparison operator */
bool InterFF::operator!=(const InterFF &other) const
{
    return not operator==(other);
}

InterFF* InterFF::clone() const
{
    return new InterFF(*this);
}

/** Return the energy components of this forcefield */
const MultiCLJComponent& InterFF::components() const
{
    return d->cljcomps;
}

/** Internal function used to rebuild the properties object that 
    stores all of the properties of this forcefield */
void InterFF::rebuildProps()
{
    //collect all of the properties from all of the CLJFunctions - note that
    //the first 'default' CLJFunction has precedence on the value of
    //properties
    d->props = d.constData()->cljfuncs.at(0).read().properties();
    
    for (int i=1; i<d.constData()->cljfuncs.count(); ++i)
    {
        Properties p = d.constData()->cljfuncs.at(i).read().properties();
        
        foreach (QString key, p.propertyKeys())
        {
            if (not d->props.hasProperty(key))
            {
                d->props.setProperty(key, p.property(key));
            }
        }
    }
    
    d->props.setProperty("cljFunction", this->cljFunction());
    d->props.setProperty("useGrid", BooleanProperty(d->fixed_atoms[0].usesGrid()));
    d->props.setProperty("gridBuffer", LengthProperty(d->fixed_atoms[0].gridBuffer()));
    d->props.setProperty("gridSpacing", LengthProperty(d->fixed_atoms[0].gridSpacing()));
    d->props.setProperty("fixedOnly", BooleanProperty(d->fixed_only));
    d->props.setProperty("parallelCalculation", BooleanProperty(d->parallel_calc));
    d->props.setProperty("reproducibleCalculation", BooleanProperty(d->repro_sum));

    for (int i=0; i<d->fixed_atoms.count(); ++i)
    {
        d->fixed_atoms[i].setUseParallelCalculation(this->usesParallelCalculation());
        d->fixed_atoms[i].setUseReproducibleCalculation(this->usesReproducibleCalculation());
    }
}

/** Set the forcefield property called 'name' to the value 'property'. Note that
    this only affects the "default" CLJFunction. Additional functions must
    be configured before adding them to the forcefield */
bool InterFF::setProperty(const QString &name, const Property &property)
{
    if (name == "cljFunction")
    {
        if (not cljFunction().equals(property))
        {
            this->setCLJFunction(property.asA<CLJFunction>());
            return true;
        }
        else
            return false;
    }
    else if (name == "useGrid")
    {
        bool use_grid = property.asA<BooleanProperty>().value();
        
        if (use_grid != d.constData()->fixed_atoms[0].usesGrid())
        {
            for (int i=0; i<d->fixed_atoms.count(); ++i)
            {
                d->fixed_atoms[i].setUseGrid(use_grid);
            }
            
            d->props.setProperty("useGrid", BooleanProperty(use_grid));
            return true;
        }
        else
            return false;
    }
    else if (name == "gridBuffer")
    {
        Length buffer = property.asA<LengthProperty>().value();
        
        if (buffer != d.constData()->fixed_atoms[0].gridBuffer())
        {
            for (int i=0; i<d->fixed_atoms.count(); ++i)
            {
                d->fixed_atoms[i].setGridBuffer(buffer);
            }
            
            d->props.setProperty("gridBuffer", property);
            return true;
        }
        else
            return false;
    }
    else if (name == "gridSpacing")
    {
        Length spacing = property.asA<LengthProperty>().value();
        
        if (spacing != d.constData()->fixed_atoms[0].gridSpacing())
        {
            for (int i=0; i<d->fixed_atoms.count(); ++i)
            {
                d->fixed_atoms[i].setGridSpacing(spacing);
            }
            
            d->props.setProperty("gridSpacing", property);
            return true;
        }
        else
            return false;
    }
    else if (name == "fixedOnly")
    {
        bool fixed_only = property.asA<BooleanProperty>().value();
        
        if (fixed_only != d.constData()->fixed_only)
        {
            d->fixed_only = fixed_only;
            d->props.setProperty("fixedOnly", property);
            this->mustNowRecalculateFromScratch();
            return true;
        }
        else
            return false;
    }
    else if (name == "parallelCalculation")
    {
        bool parallel_calc = property.asA<BooleanProperty>().value();
        
        if (parallel_calc != d.constData()->parallel_calc)
        {
            d->parallel_calc = parallel_calc;
            
            for (int i=0; i<d->fixed_atoms.count(); ++i)
            {
                d->fixed_atoms[i].setUseParallelCalculation(parallel_calc);
            }
            
            d->props.setProperty("parallelCalculation", property);
            return true;
        }
        else
            return false;
    }
    else if (name == "reproducibleCalculation")
    {
        bool repro_sum = property.asA<BooleanProperty>().value();
        
        if (repro_sum != d.constData()->repro_sum)
        {
            d->repro_sum = repro_sum;
            
            for (int i=0; i<d->fixed_atoms.count(); ++i)
            {
                d->fixed_atoms[i].setUseReproducibleCalculation(repro_sum);
            }
            
            d->props.setProperty("reproducibleCalculation", property);
            return true;
        }
        else
            return false;
    }
    else
    {
        bool found_property = false;
        bool changed_property = false;
    
        for (int i=0; i<d.constData()->cljfuncs.count(); ++i)
        {
            PropertyPtr old_prop;
            bool this_func_has_property = false;

            try
            {
                PropertyPtr old_prop = d.constData()->cljfuncs.at(i).read().property(name);
                found_property = true;
                this_func_has_property = true;
            }
            catch(...)
            {}
            
            if (this_func_has_property and not property.equals(old_prop.read()))
            {
                //need to set the property
                CLJFunctionPtr new_func = d->cljfuncs[i].read().setProperty(name, property);
                d->cljfuncs[i] = new_func;
                d->fixed_atoms[i].setCLJFunction(new_func.read());
                changed_property = true;
            }
        }

        if (changed_property)
        {
            this->rebuildProps();
            this->mustNowRecalculateFromScratch();
            return true;
        }
        else if (not found_property)
        {
            throw SireBase::missing_property( QObject::tr(
                    "No property at the key '%1' in this forcefield. Available "
                    "properties are %2.").arg(name).arg(Sire::toString(this->propertyKeys())),
                        CODELOC );
        }
        
        return false;
    }
}

/** Return the value of the forcefield property with name 'name' */
const Property& InterFF::property(const QString &name) const
{
    return d->props.property(name);
}

/** Return whether or not this forcefield contains the property 'property' */
bool InterFF::containsProperty(const QString &name) const
{
    return d->props.hasProperty(name);
}

/** Return all of the properties of this function */
const Properties& InterFF::properties() const
{
    return d->props;
}

/** Add the passed atoms as fixed atoms to the forcefield */
void InterFF::addFixedAtoms(const CLJAtoms &atoms)
{
    d->fixed_atoms[0].addFixedAtoms(atoms);
    
    for (int i=1; i<d->fixed_atoms.count(); ++i)
    {
        d->fixed_atoms[i] = d->fixed_atoms[0];
        d->fixed_atoms[i].setCLJFunction( d->cljfuncs.at(i) );
    }
    
    this->mustNowRecalculateFromScratch();
}

/** Add the passed atoms as fixed atoms to the forcefield */
void InterFF::addFixedAtoms(const MoleculeView &molecule, const PropertyMap &map)
{
    this->addFixedAtoms( CLJAtoms(molecule,map) );
}

/** Add the passed molecules as fixed atoms to the forcefield */
void InterFF::addFixedAtoms(const Molecules &molecules, const PropertyMap &map)
{
    this->addFixedAtoms( CLJAtoms(molecules,map) );
}

/** Set the fixed atoms equal to 'atoms' */
void InterFF::setFixedAtoms(const CLJAtoms &atoms)
{
    d->fixed_atoms[0].setFixedAtoms(atoms);
    
    for (int i=1; i<d->fixed_atoms.count(); ++i)
    {
        d->fixed_atoms[i] = d->fixed_atoms[0];
        d->fixed_atoms[i].setCLJFunction( d->cljfuncs.at(i).read() );
    }
    
    this->mustNowRecalculateFromScratch();
}

/** Set the fixed atoms equal to 'molecule' */
void InterFF::setFixedAtoms(const MoleculeView &molecule, const PropertyMap &map)
{
    this->setFixedAtoms( CLJAtoms(molecule,map) );
}

/** Set the fixed atoms equal to 'molecules' */
void InterFF::setFixedAtoms(const Molecules &molecules, const PropertyMap &map)
{
    this->setFixedAtoms( CLJAtoms(molecules,map) );
}

/** Set whether or not the energy calculation is only between the mobile and 
    fixed atoms (i.e. the mobile-mobile interaction is ignored) */
void InterFF::setFixedOnly(bool on)
{
    if (d.constData()->fixed_only != on)
    {
        d->fixed_only = on;
        d->props.setProperty("fixedOnly", BooleanProperty(on));
        this->mustNowRecalculateFromScratch();
    }
}

/** Set whether or not a grid is used to optimise energy calculations with the fixed atoms */
void InterFF::setUseGrid(bool on)
{
    if (this->usesGrid() != on)
    {
        for (int i=0; i<d->fixed_atoms.count(); ++i)
        {
            d->fixed_atoms[i].setUseGrid(on);
        }
        
        if (this->usesGrid() == on)
        {
            this->mustNowRecalculateFromScratch();
            d->props.setProperty("useGrid", BooleanProperty(on));
        }
        else
        {
            qDebug() << "Switching on the grid failed as none of the CLJFunctions "
                        "in this forcefield support use of a grid.";
        }
    }
}

/** Turn on the use of the grid */
void InterFF::enableGrid()
{
    this->setUseGrid(true);
}

/** Turn off use of the grid */
void InterFF::disableGrid()
{
    this->setUseGrid(false);
}

/** Return whether or not the grid is used */
bool InterFF::usesGrid() const
{
    for (int i=0; i<d.constData()->fixed_atoms.count(); ++i)
    {
        if (d->fixed_atoms.at(i).usesGrid())
            return true;
    }

    return false;
}

/** Set whether or not to use a multicore parallel algorithm
    to calculate the energy */
void InterFF::setUseParallelCalculation(bool on)
{
    if (on != usesParallelCalculation())
    {
        d->parallel_calc = on;
        d->props.setProperty("parallelCalculation", BooleanProperty(on));
        
        for (int i=0; i<d->fixed_atoms.count(); ++i)
        {
            d->fixed_atoms[i].setUseParallelCalculation(on);
        }
    }
}

/** Turn on use of a multicore parallel calculation of the energy.
    This is on by default, and spreads the energy calculations over
    available cores */
void InterFF::enableParallelCalculation()
{
    this->setUseParallelCalculation(true);
}

/** Turn off use of a multicore parallel calculation of the energy.
    This may be quicker if you have few atoms in the forcefield,
    or if you are only planning on allocating one core per forcefield */
void InterFF::disableParallelCalculation()
{
    this->setUseParallelCalculation(false);
}

/** Return whether or not a parallel algorithm is used to calculate energies */
bool InterFF::usesParallelCalculation() const
{
    return d->parallel_calc;
}

/** Turn on an energy summing algorithm that guarantees the same energy
    regardless of whether a single core or multicore calculation is being
    performed (i.e. rounding errors in both cases will be identical) */
void InterFF::enableReproducibleCalculation()
{
    setUseReproducibleCalculation(true);
}

/** Turn off an energy summing algorithm that guarantees the same energy
    regardless of whether a single core or multicore calculation is being
    performed (i.e. rounding errors in both cases will not be identical) */
void InterFF::disableReproducibleCalculation()
{
    setUseReproducibleCalculation(false);
}

/** Switch on or off use of an energy summing algorithm that guarantees the 
    same energy regardless of whether a single core or multicore calculation 
    is being performed */
void InterFF::setUseReproducibleCalculation(bool on)
{
    if (on != d->repro_sum)
    {
        d->repro_sum = on;
        
        for (int i=0; i<d->fixed_atoms.count(); ++i)
        {
            d->fixed_atoms[i].setUseReproducibleCalculation(on);
        }
        
        d->props.setProperty("reproducibleCalculation", BooleanProperty(on));
    }
}

/** Return whether or not a reproducible energy summing algorithm is being
    used to accumulate the energies */
bool InterFF::usesReproducibleCalculation() const
{
    return d->repro_sum;
}

/** Return whether or not only the energy between the mobile and fixed
    atoms is being calculated */
bool InterFF::fixedOnly() const
{
    return d->fixed_only;
}

/** Set the buffer used when using a grid. This is the distance
    added around the maximum extent of the atoms when working out the
    dimension of the grid */
void InterFF::setGridBuffer(Length buffer)
{
    if (d.constData()->fixed_atoms[0].gridBuffer() != buffer)
    {
        for (int i=0; i<d->fixed_atoms.count(); ++i)
        {
            d->fixed_atoms[i].setGridBuffer(buffer);
        }
        
        d->props.setProperty("gridBuffer", LengthProperty(buffer));
        
        if (usesGrid())
            this->mustNowRecalculateFromScratch();
    }
}

/** Return the buffer used when working out the dimension of the grid */
Length InterFF::gridBuffer() const
{
    return d->fixed_atoms[0].gridBuffer();
}

/** Set the spacing between grid points */
void InterFF::setGridSpacing(Length spacing)
{
    if (d.constData()->fixed_atoms[0].gridSpacing() != spacing)
    {
        for (int i=0; i<d->fixed_atoms.count(); ++i)
        {
            d->fixed_atoms[i].setGridSpacing(spacing);
        }
        
        d->props.setProperty("gridSpacing", LengthProperty(spacing));
        
        if (usesGrid())
            this->mustNowRecalculateFromScratch();
    }
}

/** Return spacing between grid points */
Length InterFF::gridSpacing() const
{
    return d->fixed_atoms[0].gridSpacing();
}

/** Return the grid used to calculate the energy with fixed atoms. This will
    only be set after the first energy calculation that uses the grid */
GridInfo InterFF::grid() const
{
    return d->fixed_atoms[0].grid();
}

/** Internal function used to regrid the atoms */
void InterFF::regridAtoms()
{
    if (usesGrid() and not d.constData()->fixed_atoms.isEmpty())
    {
        if (cljgroup.needsAccepting())
        {
            cljgroup.accept();
            needs_accepting = false;
            
            this->mustNowRecalculateFromScratch();
        }
        
        for (int i=0; i<d->fixed_atoms.count(); ++i)
        {
            d->fixed_atoms[i].setGridDimensions( cljgroup.cljBoxes().atoms() );
        }
        
        this->setDirty();
    }
}

/** Signal that this forcefield must now be recalculated from scratch */
void InterFF::mustNowRecalculateFromScratch()
{
    cljgroup.mustRecalculateFromScratch();
    needs_accepting = false;

    this->setDirty();
}

/** Signal to completely do everything from scratch */
void InterFF::mustNowReallyRecalculateFromScratch()
{
    cljgroup.mustReallyRecalculateFromScratch();
    needs_accepting = false;
    
    this->regridAtoms();
    this->setDirty();
}

/** Recalculate the energy of this forcefield */
void InterFF::recalculateEnergy()
{
    if (cljgroup.recalculatingFromScratch())
    {
        //calculate the energy from first principles and regenerate the
        //grid if needed
        cljgroup.accept();
        needs_accepting = false;
        
        if (cljgroup.isEmpty())
        {
            //no atoms
            if (d.constData()->cljcomps.count() == 1)
            {
                d.constData()->cljcomps.setEnergy(*this, MultiCLJEnergy(0,0));
            }
            else
            {
                int n = d.constData()->cljcomps.count();
                
                d.constData()->cljcomps.setEnergy(*this,
                        MultiCLJEnergy( QVector<double>(n,0.0), QVector<double>(n,0.0) ) );
            }
            
            this->setClean();
            return;
        }
        
        //calculate the energy from scratch
        if (d.constData()->cljcomps.count() == 1)
        {
            tuple<double,double> nrgs(0,0);
            
            if (not d.constData()->fixed_only)
            {
                if (d.constData()->parallel_calc)
                {
                    CLJCalculator calc(d->repro_sum);
                    nrgs = calc.calculate(cljFunction(), cljgroup.cljBoxes());
                }
                else
                {
                    nrgs = cljFunction().calculate(cljgroup.cljBoxes());
                }
            }

            if (not d.constData()->fixed_atoms[0].isEmpty())
            {
                this->regridAtoms();
                
                tuple<double,double> grid_nrgs = d.constData()->fixed_atoms[0]
                                                                    .calculate(cljgroup.cljBoxes());
                nrgs.get<0>() += grid_nrgs.get<0>();
                nrgs.get<1>() += grid_nrgs.get<1>();
            }
            
            d.constData()->cljcomps.setEnergy(*this, MultiCLJEnergy(nrgs.get<0>(), nrgs.get<1>()));
        }
        else
        {
            tuple< QVector<double>,QVector<double> > nrgs;
            
            if (not d.constData()->fixed_only)
            {
                if (d.constData()->parallel_calc)
                {
                    CLJCalculator calc(d->repro_sum);
                    nrgs = calc.calculate(d.constData()->cljfuncs, cljgroup.cljBoxes());
                }
                else
                {
                    nrgs = CLJFunction::multiCalculate(d.constData()->cljfuncs,
                                                       cljgroup.cljBoxes());
                }
            }
            else
            {
                nrgs.get<0>() = QVector<double>(d.constData()->cljfuncs.count(), 0.0);
                nrgs.get<1>() = QVector<double>(d.constData()->cljfuncs.count(), 0.0);
            }
            
            if (not d.constData()->fixed_atoms[0].isEmpty()) // if 0 is empty, they are all empty
            {
                this->regridAtoms();
                
                for (int i=0; i<d.constData()->fixed_atoms.count(); ++i)
                {
                    tuple<double,double> grid_nrgs = d.constData()->fixed_atoms[i]
                                                                    .calculate(cljgroup.cljBoxes());

                    nrgs.get<0>()[i] += grid_nrgs.get<0>();
                    nrgs.get<1>()[i] += grid_nrgs.get<1>();
                }
            }
            
            d.constData()->cljcomps.setEnergy(*this, MultiCLJEnergy(nrgs.get<0>(), nrgs.get<1>()));
        }
        
        this->setClean();
    }
    else if (cljgroup.needsAccepting())
    {
        //we can calculate using just the change in energy
        if (d.constData()->cljfuncs.count() == 1)
        {
            tuple<double,double> delta_nrgs(0,0);

            CLJAtoms changed_atoms;
            
            if (cljgroup.isSingleIDChange() or d.constData()->fixed_only)
            {
                //the changed atoms don't interact with each other
                changed_atoms = cljgroup.changedAtoms();
            
                if (not d.constData()->fixed_only)
                {
                    //calculate the change in energy using the molecules in changed_atoms
                    if (d.constData()->parallel_calc)
                    {
                        CLJCalculator calc(d.constData()->repro_sum);
                        delta_nrgs = calc.calculate(cljFunction(),
                                                    changed_atoms, cljgroup.cljBoxes());
                    }
                    else
                    {
                        delta_nrgs = cljFunction().calculate(changed_atoms, cljgroup.cljBoxes());
                    }
                }
            }
            else
            {
                //the changed atoms interact with one another
                tuple<CLJAtoms,CLJAtoms,CLJAtoms> changes = cljgroup.mergeChanges();
                
                changed_atoms = changes.get<0>();
                const CLJAtoms &old_atoms = changes.get<1>();
                const CLJAtoms &new_atoms = changes.get<2>();
                
                if (d.constData()->parallel_calc)
                {
                    CLJCalculator calc(d.constData()->repro_sum);
                    delta_nrgs = calc.calculate(cljFunction(),
                                                changed_atoms, cljgroup.cljBoxes());
                }
                else
                {
                    delta_nrgs = cljFunction().calculate(changed_atoms, cljgroup.cljBoxes());
                }

                //add on the change in interaction within changed_atoms
                tuple<double,double> old_changes = cljFunction().calculate(old_atoms);
                tuple<double,double> new_changes = cljFunction().calculate(new_atoms);
                
                delta_nrgs.get<0>() += new_changes.get<0>() - old_changes.get<0>();
                delta_nrgs.get<1>() += new_changes.get<1>() - old_changes.get<1>();
            }
            
            if (not d.constData()->fixed_atoms[0].isEmpty())
            {
                tuple<double,double> grid_deltas = d.constData()->fixed_atoms[0]
                                                                .calculate(changed_atoms);
                
                delta_nrgs.get<0>() += grid_deltas.get<0>();
                delta_nrgs.get<1>() += grid_deltas.get<1>();
            }
            
            d.constData()->cljcomps.changeEnergy(*this,
                                        MultiCLJEnergy(delta_nrgs.get<0>(), delta_nrgs.get<1>()));
        }
        else
        {
            tuple< QVector<double>,QVector<double> > delta_nrgs;

            CLJAtoms changed_atoms;
            
            if (cljgroup.isSingleIDChange() or d.constData()->fixed_only)
            {
                //the changed atoms don't interact with each other
                changed_atoms = cljgroup.changedAtoms();
            
                if (not d.constData()->fixed_only)
                {
                    //calculate the change in energy using the molecules in changed_atoms
                    if (d.constData()->parallel_calc)
                    {
                        CLJCalculator calc(d.constData()->repro_sum);
                        delta_nrgs = calc.calculate(d.constData()->cljfuncs,
                                                    changed_atoms, cljgroup.cljBoxes());
                    }
                    else
                    {
                        delta_nrgs = CLJFunction::multiCalculate(d.constData()->cljfuncs,
                                                    changed_atoms, cljgroup.cljBoxes());
                    }
                }
            }
            else
            {
                //the changed atoms interact with one another
                tuple<CLJAtoms,CLJAtoms,CLJAtoms> changes = cljgroup.mergeChanges();
                
                changed_atoms = changes.get<0>();
                const CLJAtoms &old_atoms = changes.get<1>();
                const CLJAtoms &new_atoms = changes.get<2>();
                
                if (d.constData()->parallel_calc)
                {
                    CLJCalculator calc(d.constData()->repro_sum);
                    delta_nrgs = calc.calculate(d.constData()->cljfuncs,
                                                changed_atoms, cljgroup.cljBoxes());
                }
                else
                {
                    delta_nrgs = CLJFunction::multiCalculate(d.constData()->cljfuncs,
                                                             changed_atoms, cljgroup.cljBoxes());
                }

                //add on the change in interaction within changed_atoms
                tuple< QVector<double>,QVector<double> > old_changes
                            = CLJFunction::multiCalculate(d.constData()->cljfuncs,old_atoms);
                tuple< QVector<double>,QVector<double> > new_changes
                            = CLJFunction::multiCalculate(d.constData()->cljfuncs,new_atoms);
                
                for (int i=0; i<d.constData()->cljfuncs.count(); ++i)
                {
                    delta_nrgs.get<0>()[i] += new_changes.get<0>()[i] - old_changes.get<0>()[i];
                    delta_nrgs.get<1>()[i] += new_changes.get<1>()[i] - old_changes.get<1>()[i];
                }
            }
            
            if (not d.constData()->fixed_atoms[0].isEmpty()) // if 0 is empty, they are all empty
            {
                for (int i=0; i<d.constData()->fixed_atoms.count(); ++i)
                {
                    tuple<double,double> grid_deltas = d.constData()->fixed_atoms[i]
                                                                .calculate(changed_atoms);
                
                    //only add on grid energies for the default CLJ function
                    delta_nrgs.get<0>()[i] += grid_deltas.get<0>();
                    delta_nrgs.get<1>()[i] += grid_deltas.get<1>();
                }
            }
            
            d.constData()->cljcomps.changeEnergy(*this,
                                        MultiCLJEnergy(delta_nrgs.get<0>(), delta_nrgs.get<1>()));
        }

        //the CLJGroup needs to be accepted before we can change anything else
        needs_accepting = true;

        this->setClean();
    }
    else
    {
        //recalculate everything from scratch as this has been requested
        //calculate the energy from scratch
        cljgroup.accept();
        needs_accepting = false;

        if (d.constData()->cljfuncs.count() == 1)
        {
            tuple<double,double> nrgs(0,0);
            
            if (not d.constData()->fixed_only)
            {
                if (d.constData()->parallel_calc)
                {
                    CLJCalculator calc(d.constData()->repro_sum);
                    nrgs = calc.calculate(cljFunction(), cljgroup.cljBoxes());
                }
                else
                {
                    nrgs = cljFunction().calculate(cljgroup.cljBoxes());
                }
            }

            if (not d.constData()->fixed_atoms[0].isEmpty())
            {
                this->regridAtoms();
                tuple<double,double> grid_nrgs = d.constData()->fixed_atoms[0]
                                                        .calculate(cljgroup.cljBoxes());
                nrgs.get<0>() += grid_nrgs.get<0>();
                nrgs.get<1>() += grid_nrgs.get<1>();
            }
            
            d.constData()->cljcomps.setEnergy(*this, MultiCLJEnergy(nrgs.get<0>(), nrgs.get<1>()));
        }
        else
        {
            tuple< QVector<double>,QVector<double> > nrgs;
            
            if (not d.constData()->fixed_only)
            {
                if (d.constData()->parallel_calc)
                {
                    CLJCalculator calc(d.constData()->repro_sum);
                    nrgs = calc.calculate(d.constData()->cljfuncs, cljgroup.cljBoxes());
                }
                else
                {
                    nrgs = CLJFunction::multiCalculate(d.constData()->cljfuncs,
                                                       cljgroup.cljBoxes());
                }
            }
            else
            {
                nrgs.get<0>() = QVector<double>( d.constData()->cljfuncs.count(), 0.0 );
                nrgs.get<1>() = QVector<double>( d.constData()->cljfuncs.count(), 0.0 );
            }

            if (not d.constData()->fixed_atoms[0].isEmpty())
            {
                this->regridAtoms();
                
                for (int i=0; i<d.constData()->fixed_atoms.count(); ++i)
                {
                    tuple<double,double> grid_nrgs = d.constData()->fixed_atoms[i]
                                                                .calculate(cljgroup.cljBoxes());

                    //we only calculate the grid energy for the default CLJ function
                    nrgs.get<0>()[i] += grid_nrgs.get<0>();
                    nrgs.get<1>()[i] += grid_nrgs.get<1>();
                }
            }
            
            d.constData()->cljcomps.setEnergy(*this, MultiCLJEnergy(nrgs.get<0>(), nrgs.get<1>()));
        }
        
        this->setClean();
    }
}

/** Function called to add a molecule to this forcefield */
void InterFF::_pvt_added(const SireMol::PartialMolecule &mol, const SireBase::PropertyMap &map)
{
    if (needs_accepting)
    {
        cljgroup.accept();
        needs_accepting = false;
    }

    cljgroup.add(mol, map);
    setDirty();
}

/** Function called to remove a molecule from this forcefield */
void InterFF::_pvt_removed(const SireMol::PartialMolecule &mol)
{
    if (needs_accepting)
    {
        cljgroup.accept();
        needs_accepting = false;
    }

    cljgroup.remove(mol);
    setDirty();
}

/** Function called to indicate that the passed molecule has changed */
void InterFF::_pvt_changed(const Molecule &molecule, bool auto_update)
{
    if (needs_accepting)
    {
        cljgroup.accept();
        needs_accepting = false;
    }

    cljgroup.update(molecule);
    setDirty();
}

/** Function called to indicate that a list of molecules in this forcefield have changed */
void InterFF::_pvt_changed(const QList<SireMol::Molecule> &molecules, bool auto_update)
{
    if (needs_accepting)
    {
        cljgroup.accept();
        needs_accepting = false;
    }

    foreach (const Molecule &molecule, molecules)
    {
        cljgroup.update(molecule);
    }
    
    setDirty();
}

/** Function called to indicate that all molecules in this forcefield have been removed */
void InterFF::_pvt_removedAll()
{
    if (needs_accepting)
    {
        cljgroup.accept();
        needs_accepting = false;
    }

    cljgroup.removeAll();
    this->setDirty();
}

/** Function called to query whether or not a change in source properties would
    change the properties needed by this forcefield for the molecule with number 'molnum' */
bool InterFF::_pvt_wouldChangeProperties(SireMol::MolNum molnum,
                                         const SireBase::PropertyMap &map) const
{
    return cljgroup.mapForMolecule(molnum) != map;
}

/** Return whether or not this forcefield is using a temporary workspace that 
    needs to be accepted */
bool InterFF::needsAccepting() const
{
    return needs_accepting or G1FF::needsAccepting();
}

/** Tell the forcefield that the last move was accepted. This tells the
    forcefield to make permanent any temporary changes that were used a workspace
    to avoid memory allocation during a move */
void InterFF::accept()
{
    if (needs_accepting)
    {
        cljgroup.accept();
        needs_accepting = false;
    }
    
    G1FF::accept();
}
