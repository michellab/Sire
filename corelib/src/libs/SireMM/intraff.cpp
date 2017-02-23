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

#include "intraff.h"

#include "cljshiftfunction.h"
#include "cljcalculator.h"

#include "SireBase/booleanproperty.h"
#include "SireBase/lengthproperty.h"
#include "SireBase/refcountdata.h"

#include "SireError/errors.h"
#include "SireBase/errors.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/molresid.h"
#include "SireMol/residue.h"
#include "SireMol/atomselection.h"
#include "SireMol/selector.hpp"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QElapsedTimer>
#include <QDebug>

using namespace SireMM;
using namespace SireMol;
using namespace SireFF;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits;

QDataStream& operator<<(QDataStream &ds,
                        const SharedDataPointer<SireMM::detail::IntraFFMolData> &ptr)
{
    SharedDataStream sds(ds);
    sds << ptr;
    return ds;
}

QDataStream& operator>>(QDataStream &ds, SharedDataPointer<SireMM::detail::IntraFFMolData> &ptr)
{
    SharedDataStream sds(ds);
    sds >> ptr;
    return ds;
}

namespace SireMM
{
    namespace detail
    {
        /** This class holds (mostly) const data about the forcefield */
        class IntraFFData : public RefCountData
        {
        public:
            IntraFFData() : RefCountData(),
                            parallel_calc(true), repro_sum(false)
            {}
            
            IntraFFData(const IntraFFData &other)
                 : RefCountData(),
                   cljfuncs(other.cljfuncs),
                   cljcomps(other.cljcomps),
                   props(other.props),
                   parallel_calc(other.parallel_calc),
                   repro_sum(other.repro_sum)
            {}
            
            ~IntraFFData()
            {}

            /** The function(s) used to calculate energies */
            QVector<CLJFunctionPtr> cljfuncs;

            /** The energy components available for this forcefield */
            MultiCLJComponent cljcomps;
            
            /** All of the properties in this forcefield */
            Properties props;
            
            /** Whether or not to calculate energies in parallel */
            bool parallel_calc;
            
            /** Whether or not to sum energies using a reproducible sum */
            bool repro_sum;
        };
        
        /** This class holds all of the information about a particular 
            molecule in the forcefield */
        class IntraFFMolData : public RefCountData
        {
        public:
            /** The CLJGroup for this molecule */
            CLJGroup cljgroup;
            
            /** All of the CLJ functions that have been specialised for
                this molecule */
            QVector<CLJFunctionPtr> cljfuncs;
            
            /** The connectivity of the molecule */
            Connectivity cty;
            
            /** All of the energies for this molecule */
            MultiCLJEnergy nrg;
            
            /** The property name of the connectivity property */
            PropertyName connectivity_property;
            
            /** The current version of the connectivity property */
            quint64 connectivity_version;
            
            /** Whether or not the energy needs to be calculated */
            bool needs_energy_calc;
            
            /** Whether or not this group needs accepting */
            bool needs_accepting;

            IntraFFMolData() : RefCountData(), connectivity_version(0),
                               needs_energy_calc(false), needs_accepting(false)
            {}
            
            IntraFFMolData(const MoleculeView &molview,
                           const PropertyMap &map,
                           const QVector<CLJFunctionPtr> &funcs)
                    : RefCountData(), needs_energy_calc(true), needs_accepting(false)
            {
                cljgroup = CLJGroup(CLJAtoms::USE_ATOMIDX, CLJExtractor::EXTRACT_BY_CUTGROUP);
                cljgroup.setBoxLength( 7.5*angstrom );
                cljgroup.add(molview, map);
                connectivity_property = map["connectivity"];
                cty = molview.data().property(connectivity_property).asA<Connectivity>();
                connectivity_version = molview.data().version(connectivity_property);
                
                setCLJFunctions(funcs);
            }
            
            IntraFFMolData(const IntraFFMolData &other)
                : RefCountData(), cljgroup(other.cljgroup), cljfuncs(other.cljfuncs),
                  cty(other.cty), nrg(other.nrg),
                  connectivity_property(other.connectivity_property),
                  connectivity_version(other.connectivity_version),
                  needs_energy_calc(other.needs_energy_calc),
                  needs_accepting(other.needs_accepting)
            {}
            
            ~IntraFFMolData()
            {}
            
            void mustNowRecalculateFromScratch()
            {
                cljgroup.mustRecalculateFromScratch();
                needs_accepting = false;
                needs_energy_calc = true;
                nrg = MultiCLJEnergy(0,0);
            }
            
            void mustReallyRecalculateFromScratch()
            {
                cljgroup.mustReallyRecalculateFromScratch();
                needs_accepting = false;
                needs_energy_calc = true;
                nrg = MultiCLJEnergy(0,0);
            }
            
            void checkForChangeInConnectivity(const MoleculeView &molview)
            {
                if (molview.data().version(connectivity_property) != connectivity_version)
                {
                    //the connectivity may have changed
                    Connectivity new_cty = molview.data().property(connectivity_property)
                                                         .asA<Connectivity>();
                    
                    connectivity_version = molview.data().version(connectivity_property);
                    
                    if (cty != new_cty)
                    {
                        QVector<CLJFunctionPtr> newfuncs = cljfuncs;
                        
                        for (int i=0; i<newfuncs.count(); ++i)
                        {
                            newfuncs[i].edit().asA<CLJIntraFunction>().setConnectivity(new_cty);
                        }
                        
                        cljfuncs = newfuncs;
                        cty = new_cty;
                        mustReallyRecalculateFromScratch();
                    }
                }
            }
            
            void setCLJFunctions(const QVector<CLJFunctionPtr> &funcs)
            {
                QVector<CLJFunctionPtr> newfuncs = funcs;
                
                for (int i=0; i<funcs.count(); ++i)
                {
                    newfuncs[i].edit().asA<CLJIntraFunction>().setConnectivity(cty);
                }
                
                cljfuncs = newfuncs;
                mustNowRecalculateFromScratch();
            }
        };
    }
}

QDataStream& operator<<(QDataStream &ds, const SireMM::detail::IntraFFMolData &intraff)
{
    quint32 version = 1;
    
    ds << version;
    
    SharedDataStream sds(ds);
    
    sds << intraff.cljgroup << intraff.cljfuncs << intraff.cty
        << intraff.nrg << intraff.connectivity_property
        << intraff.connectivity_version << intraff.needs_energy_calc
        << intraff.needs_accepting;
    
    return ds;
}

QDataStream& operator>>(QDataStream &ds, SireMM::detail::IntraFFMolData &intraff)
{
    quint32 version;
    
    ds >> version;
    
    if (version == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> intraff.cljgroup >> intraff.cljfuncs >> intraff.cty
            >> intraff.nrg >> intraff.connectivity_property
            >> intraff.connectivity_version >> intraff.needs_energy_calc
            >> intraff.needs_accepting;
    }
    else
        throw version_error( QObject::tr(
                "Wrong version of SireMM::detail::IntraFFData. Version %1, while "
                "only version 1 is supported.").arg(version), CODELOC );
    
    return ds;
}

static const RegisterMetaType<IntraFF> r_intraff;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const IntraFF &intraff)
{
    writeHeader(ds, r_intraff, 1);
    
    SharedDataStream sds(ds);
    
    sds << intraff.moldata
        << intraff.needs_accepting
        << intraff.d->cljfuncs
        << intraff.d->cljcomps
        << intraff.d->parallel_calc
        << intraff.d->repro_sum
        << static_cast<const G1FF&>(intraff);
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, IntraFF &intraff)
{
    VersionID v = readHeader(ds, r_intraff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> intraff.moldata
            >> intraff.needs_accepting
            >> intraff.d->cljfuncs
            >> intraff.d->cljcomps
            >> intraff.d->parallel_calc
            >> intraff.d->repro_sum
            >> static_cast<G1FF&>(intraff);
        
        intraff.rebuildProps();
    }
    else
        throw version_error(v, "1", r_intraff, CODELOC);
    
    return ds;
}


/** Constructor */
IntraFF::IntraFF()
        : ConcreteProperty<IntraFF,G1FF>(), needs_accepting(false)
{
    d = new detail::IntraFFData();
    this->_pvt_updateName();
    this->setCLJFunction( CLJIntraShiftFunction::defaultShiftFunction()
                                    .read().asA<CLJIntraFunction>() );
}

/** Construct, specifying the name of the forcefield */
IntraFF::IntraFF(const QString &name)
        : ConcreteProperty<IntraFF, G1FF>(),
          needs_accepting(false)
{
    d = new detail::IntraFFData();
    G1FF::setName(name);
    this->setCLJFunction( CLJIntraShiftFunction::defaultShiftFunction()
                                    .read().asA<CLJIntraFunction>() );
}

/** Copy constructor */
IntraFF::IntraFF(const IntraFF &other)
        : ConcreteProperty<IntraFF,G1FF>(other),
          moldata(other.moldata), d(other.d),
          needs_accepting(other.needs_accepting)
{}

/** Destructor */
IntraFF::~IntraFF()
{}

/** Function used to set the CLJIntraFunction used to calculate 
    the intramolecular energy */
void IntraFF::setCLJFunction(const CLJIntraFunction &func)
{
    if (d.constData()->cljfuncs.isEmpty())
    {
        d->cljfuncs.append(func);
    }
    else if (not d.constData()->cljfuncs.at(0).read().equals(func))
    {
        d->cljfuncs[0] = func;
    }

    if (not moldata.isEmpty())
    {
        for (MolData::iterator it = moldata.begin();
             it != moldata.end();
             ++it)
        {
            it.value()->setCLJFunctions(d->cljfuncs);
        }
    }

    rebuildProps();
    this->mustNowRecalculateFromScratch();
}

/** Return the function used to calculate the energy */
const CLJIntraFunction& IntraFF::cljFunction() const
{
    if (d.constData()->cljfuncs.isEmpty())
        throw SireError::program_bug( QObject::tr(
                "There should always be at least one CLJFunction in InterFF!"),
                    CODELOC );
    
    return d.constData()->cljfuncs.at(0).read().asA<CLJIntraFunction>();
}

/** Set the CLJFunction with key 'key' equal to 'cljfunc' */
void IntraFF::setCLJFunction(QString key, const CLJIntraFunction &cljfunc)
{
    if (key == "all")
    {
        for (int i=0; i<d->cljfuncs.count(); ++i)
        {
            d->cljfuncs[i] = cljfunc;
        }

        if (not moldata.isEmpty())
        {
            for (MolData::iterator it = moldata.begin();
                 it != moldata.end();
                 ++it)
            {
                it.value()->setCLJFunctions(d->cljfuncs);
            }
        }
        
        rebuildProps();
        
        this->mustNowRecalculateFromScratch();
    }
    else if (key == "default")
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
        
        d->cljfuncs[idx] = cljfunc;

        if (not moldata.isEmpty())
        {
            for (MolData::iterator it = moldata.begin();
                 it != moldata.end();
                 ++it)
            {
                it.value()->setCLJFunctions(d->cljfuncs);
            }
        }

        rebuildProps();
        this->mustNowRecalculateFromScratch();
    }
}

/** Remove the CLJ function with key 'key' - note that you cannot remove
    the 'default' CLJ function */
void IntraFF::removeCLJFunctionAt(QString key)
{
    if (key != "default")
    {
        int idx = d->cljcomps.remove(key);
        
        if (idx > 0)
        {
            d->cljfuncs.removeAt(idx);

            if (not moldata.isEmpty())
            {
                for (MolData::iterator it = moldata.begin();
                     it != moldata.end();
                     ++it)
                {
                    it.value()->setCLJFunctions(d->cljfuncs);
                }
            }
            
            rebuildProps();
            this->mustNowRecalculateFromScratch();
        }
    }
}

/** Function to remove all of the CLJFunctions (except for the "default" function) */
void IntraFF::removeAllCLJFunctions()
{
    if (d->cljfuncs.count() < 2)
        return;

    d->cljcomps.removeAll();
    
    while (d->cljfuncs.count() > 1)
    {
        d->cljfuncs.removeLast();
    }

    if (not moldata.isEmpty())
    {
        for (MolData::iterator it = moldata.begin();
             it != moldata.end();
             ++it)
        {
            it.value()->setCLJFunctions(d->cljfuncs);
        }
    }

    rebuildProps();
    this->mustNowRecalculateFromScratch();
}

/** Return the keys of all CLJFunctions added to this forcefield */
QStringList IntraFF::cljFunctionKeys() const
{
    return d.constData()->cljcomps.keys();
}

/** Return the CLJFunction associated with the passed key */
const CLJIntraFunction& IntraFF::cljFunction(QString key) const
{
    return d.constData()->cljfuncs.at( d.constData()->cljcomps.indexOf(key) )
                                  .read().asA<CLJIntraFunction>();
}

/** Return the number of CLJ functions in this forcefield. There should always
    be at least one */
int IntraFF::nCLJFunctions() const
{
    return d.constData()->cljfuncs.count();
}

/** Return the hash of all CLJFunctions in this forcefield, indexed by their key */
QHash<QString,CLJFunctionPtr> IntraFF::cljFunctions() const
{
    QHash<QString,CLJFunctionPtr> funcs;
    
    foreach (QString key, d->cljcomps.keys())
    {
        funcs.insert( key, this->cljFunction(key) );
    }
    
    return funcs;
}

/** Internal function called when the name of the forcefield changes */
void IntraFF::_pvt_updateName()
{
    d->cljcomps = d->cljcomps.rename(this->name());
    G1FF::_pvt_updateName();
}

const char* IntraFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<IntraFF>() );
}

const char* IntraFF::what() const
{
    return IntraFF::typeName();
}

/** Copy assignment operator */
IntraFF& IntraFF::operator=(const IntraFF &other)
{
    if (this != &other)
    {
        moldata = other.moldata;
        needs_accepting = other.needs_accepting;
        d = other.d;
        G1FF::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool IntraFF::operator==(const IntraFF &other) const
{
    return (this == &other) or
           (G1FF::operator==(other) and
            d->cljfuncs == other.d->cljfuncs and d->cljcomps == other.d->cljcomps and
            moldata == other.moldata and needs_accepting == other.needs_accepting);
}

/** Comparison operator */
bool IntraFF::operator!=(const IntraFF &other) const
{
    return not operator==(other);
}

IntraFF* IntraFF::clone() const
{
    return new IntraFF(*this);
}

/** Return the energy components of this forcefield */
const MultiCLJComponent& IntraFF::components() const
{
    return d->cljcomps;
}

/** Internal function used to rebuild the properties object that 
    stores all of the properties of this forcefield */
void IntraFF::rebuildProps()
{
    //collect all of the properties from all of the CLJFunctions - note that
    //the first 'default' CLJFunction has precedence on the value of
    //properties
    d->props = d.constData()->cljfuncs.at(0).read().properties();
    
    for (QString key : d.constData()->cljcomps.keys())
    {
        const int idx = d.constData()->cljcomps.indexOf(key);
    
        d->props.setProperty( QString("cljFunction[%1]").arg(key),
                              d.constData()->cljfuncs.at(idx) );
    
        Properties p = d.constData()->cljfuncs.at(idx).read().properties();
        
        foreach (QString propkey, p.propertyKeys())
        {
            d->props.setProperty( QString("%1[%2]").arg(propkey).arg(key),
                                  p.property(propkey) );
        }
    }
    
    d->props.setProperty("cljFunction", this->cljFunction());
    d->props.setProperty("parallelCalculation", BooleanProperty(d->parallel_calc));
    d->props.setProperty("reproducibleCalculation", BooleanProperty(d->repro_sum));
}

/** Set the forcefield property called 'name' to the value 'property'. Note that
    this only affects the "default" CLJFunction. Additional functions must
    be configured before adding them to the forcefield */
bool IntraFF::setProperty(const QString &name, const Property &property)
{
    if (name == "cljFunction")
    {
        if (not cljFunction().equals(property))
        {
            this->setCLJFunction(property.asA<CLJIntraFunction>());
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
            d->props.setProperty("reproducibleCalculation", property);
            return true;
        }
        else
            return false;
    }
    else
    {
        //see if the property is in the default CLJFunction
        if (cljFunction().containsProperty(name))
        {
            CLJFunctionPtr newfunc = cljFunction().setProperty(name, property);
            
            if ( not cljFunction().equals( newfunc.read() ) )
            {
                this->setCLJFunction(newfunc.read().asA<CLJIntraFunction>());
                return true;
            }
            else
                return false;
        }

        //ok, it doesn't. Now look to see if the property has an index, meaning
        //that it may belong to one of the extra CLJ functions
        auto subscr = getSubscriptedProperty(name);
        QString cljname = subscr.get<0>();
        QString cljkey = subscr.get<1>();
    
        if (cljname == "cljFunction")
        {
            if (cljkey == "all")
            {
                this->setCLJFunction("all", property.asA<CLJIntraFunction>());
            }
            else if (not cljFunction(cljkey).equals(property))
            {
                this->setCLJFunction(cljkey, property.asA<CLJIntraFunction>());
                return true;
            }
            else
                return false;
        }
        else if (cljkey.length() > 0)
        {
            bool found_property = false;
            bool changed_property = false;

            if (cljkey == "all")
            {
                //set this property in all cljfunctions (if possible)
                for (int i=0; i<d.constData()->cljfuncs.count(); ++i)
                {
                    PropertyPtr old_prop;
                    bool this_func_has_property = false;

                    try
                    {
                        PropertyPtr old_prop = d.constData()->cljfuncs.at(i)
                                                                      .read().property(cljname);
                        found_property = true;
                        this_func_has_property = true;
                    }
                    catch(...)
                    {}
            
                    if (this_func_has_property and not property.equals(old_prop.read()))
                    {
                        //need to set the property
                        CLJFunctionPtr new_func = d->cljfuncs[i].read()
                                                                .setProperty(cljname, property);
                        d->cljfuncs[i] = new_func;
                        changed_property = true;
                    }
                }
            }
            else if (cljkey == "default")
            {
                //set this property in the default cljfunction
                if (cljFunction().containsProperty(cljname))
                {
                    found_property = true;
                
                    CLJFunctionPtr newfunc = cljFunction().setProperty(cljname, property);
            
                    if ( not cljFunction().equals( newfunc.read() ) )
                    {
                        this->setCLJFunction(newfunc.read().asA<CLJIntraFunction>());
                        return true;
                    }
                    else
                        return false;
                }
            }
            else
            {
                if (d.constData()->cljcomps.hasKey(cljkey))
                {
                    found_property = true;

                    if (cljFunction(cljkey).containsProperty(cljname))
                    {
                        CLJFunctionPtr newfunc = cljFunction(cljkey).setProperty(cljname, property);
                        
                        if ( not cljFunction(cljkey).equals( newfunc.read() ) )
                        {
                            this->setCLJFunction(cljkey, newfunc.read().asA<CLJIntraFunction>());
                            return true;
                        }
                        else
                            return false;
                    }
                }
            }

            if (changed_property)
            {
                if (not moldata.isEmpty())
                {
                    for (MolData::iterator it = moldata.begin();
                         it != moldata.end();
                         ++it)
                    {
                        it.value()->setCLJFunctions(d->cljfuncs);
                    }
                }
        
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
        }
        else
            throw SireBase::missing_property( QObject::tr(
                    "No property at the key '%1' in this forcefield. Available "
                    "properties are %2.").arg(name).arg(Sire::toString(this->propertyKeys())),
                        CODELOC );
        
        return false;
    }
}

/** Return the value of the forcefield property with name 'name' */
const Property& IntraFF::property(const QString &name) const
{
    return d->props.property(name);
}

/** Return whether or not this forcefield contains the property 'property' */
bool IntraFF::containsProperty(const QString &name) const
{
    return d->props.hasProperty(name);
}

/** Return all of the properties of this function */
const Properties& IntraFF::properties() const
{
    return d->props;
}

/** Set whether or not to use a multicore parallel algorithm
    to calculate the energy */
void IntraFF::setUseParallelCalculation(bool on)
{
    if (on != usesParallelCalculation())
    {
        d->parallel_calc = on;
        d->props.setProperty("parallelCalculation", BooleanProperty(on));
    }
}

/** Turn on use of a multicore parallel calculation of the energy.
    This is on by default, and spreads the energy calculations over
    available cores */
void IntraFF::enableParallelCalculation()
{
    this->setUseParallelCalculation(true);
}

/** Turn off use of a multicore parallel calculation of the energy.
    This may be quicker if you have few atoms in the forcefield,
    or if you are only planning on allocating one core per forcefield */
void IntraFF::disableParallelCalculation()
{
    this->setUseParallelCalculation(false);
}

/** Return whether or not a parallel algorithm is used to calculate energies */
bool IntraFF::usesParallelCalculation() const
{
    return d->parallel_calc;
}

/** Turn on an energy summing algorithm that guarantees the same energy
    regardless of whether a single core or multicore calculation is being
    performed (i.e. rounding errors in both cases will be identical) */
void IntraFF::enableReproducibleCalculation()
{
    setUseReproducibleCalculation(true);
}

/** Turn off an energy summing algorithm that guarantees the same energy
    regardless of whether a single core or multicore calculation is being
    performed (i.e. rounding errors in both cases will not be identical) */
void IntraFF::disableReproducibleCalculation()
{
    setUseReproducibleCalculation(false);
}

/** Switch on or off use of an energy summing algorithm that guarantees the 
    same energy regardless of whether a single core or multicore calculation 
    is being performed */
void IntraFF::setUseReproducibleCalculation(bool on)
{
    if (on != d->repro_sum)
    {
        d->repro_sum = on;
        d->props.setProperty("reproducibleCalculation", BooleanProperty(on));
    }
}

/** Return whether or not a reproducible energy summing algorithm is being
    used to accumulate the energies */
bool IntraFF::usesReproducibleCalculation() const
{
    return d->repro_sum;
}

/** Signal that this forcefield must now be recalculated from scratch */
void IntraFF::mustNowRecalculateFromScratch()
{
    if (not moldata.isEmpty())
    {
        for (MolData::iterator it = moldata.begin();
             it != moldata.end();
             ++it)
        {
            it.value()->mustNowRecalculateFromScratch();
        }
    }

    needs_accepting = false;

    this->setDirty();
}

/** Signal to completely do everything from scratch */
void IntraFF::mustNowReallyRecalculateFromScratch()
{
    if (not moldata.isEmpty())
    {
        for (MolData::iterator it = moldata.begin();
             it != moldata.end();
             ++it)
        {
            it.value()->mustReallyRecalculateFromScratch();
        }
    }

    needs_accepting = false;

    this->setDirty();
}

/** Calculate the intramolecular energy of a particular molecule */
MultiCLJEnergy IntraFF::calcEnergy(detail::IntraFFMolData &mol) const
{
    if (mol.needs_energy_calc)
    {
        if (mol.cljgroup.recalculatingFromScratch())
        {
            if (mol.cljgroup.needsAccepting() or mol.needs_accepting)
            {
                mol.cljgroup.accept();
                mol.needs_accepting = false;
            }
        
            if (mol.cljfuncs.count() == 1)
            {
                tuple<double,double> nrg;
            
                if (usesParallelCalculation())
                {
                    CLJCalculator cljcalc( usesReproducibleCalculation() );
                    
                    nrg = cljcalc.calculate( mol.cljfuncs.at(0).read(),
                                             mol.cljgroup.cljBoxes() );
                }
                else
                {
                    nrg = mol.cljfuncs.at(0).read().calculate( mol.cljgroup.cljBoxes() );
                }
                
                mol.nrg = MultiCLJEnergy( nrg.get<0>(), nrg.get<1>() );
                mol.needs_energy_calc = false;
                return mol.nrg;
            }
            else
            {
                tuple< QVector<double>,QVector<double> > nrg;
                
                if (usesParallelCalculation())
                {
                    CLJCalculator cljcalc( usesReproducibleCalculation() );
                    
                    nrg = cljcalc.calculate( mol.cljfuncs, mol.cljgroup.cljBoxes() );
                }
                else
                {
                    nrg = CLJFunction::multiCalculate( mol.cljfuncs,
                                                       mol.cljgroup.cljBoxes() );
                }
                
                mol.nrg = MultiCLJEnergy( nrg.get<0>(), nrg.get<1>() );
                mol.needs_energy_calc = false;
                return mol.nrg;
            }
        }
        else if (mol.cljgroup.needsAccepting())
        {
            //we can calculate the energy using a delta
            
            //get the changed atoms  tuple<changedAtoms(),oldAtoms(),newAtoms()>
            tuple<CLJAtoms,CLJAtoms,CLJAtoms> change = mol.cljgroup.mergeChanges();
            
            if (mol.cljfuncs.count() == 1)
            {
                tuple<double,double> delta, oldnrg, newnrg;

                if (usesParallelCalculation())
                {
                    CLJCalculator cljcalc(usesReproducibleCalculation());
                    
                    //get the change in energy from the changed atoms with the
                    //atoms that have not changed
                    delta = cljcalc.calculate( mol.cljfuncs.at(0).read(),
                                               change.get<0>(), mol.cljgroup.cljBoxes() );
                    
                    //get the energy of the old atoms interacting with each other
                    oldnrg = cljcalc.calculate( mol.cljfuncs.at(0).read(),
                                                change.get<1>() );
                    
                    //get the energy of the new atoms interacting with each other
                    newnrg = cljcalc.calculate( mol.cljfuncs.at(0).read(),
                                                change.get<2>() );
                }
                else
                {
                    //get the change in energy from the changed atoms with the
                    //atoms that have not changed
                    delta = mol.cljfuncs.at(0).read().calculate(
                                               change.get<0>(), mol.cljgroup.cljBoxes() );
                    
                    //get the energy of the old atoms interacting with each other
                    oldnrg = mol.cljfuncs.at(0).read().calculate( change.get<1>() );
                    
                    //get the energy of the new atoms interacting with each other
                    newnrg = mol.cljfuncs.at(0).read().calculate( change.get<2>() );
                }
                
                mol.nrg += MultiCLJEnergy(delta.get<0>(),delta.get<1>()) +
                           MultiCLJEnergy(newnrg.get<0>(),newnrg.get<1>()) -
                           MultiCLJEnergy(oldnrg.get<0>(),oldnrg.get<1>());
                
                mol.needs_energy_calc = false;
                mol.needs_accepting = true;
                
                return mol.nrg;
            }
            else  // delta calculation with more than one clj function
            {
                tuple< QVector<double>,QVector<double> > delta, oldnrg, newnrg;
                
                if (usesParallelCalculation())
                {
                    CLJCalculator cljcalc(usesReproducibleCalculation());
                    
                    //get the change in energy from the changed atoms with the
                    //atoms that have not changed
                    delta = cljcalc.calculate( mol.cljfuncs,
                                               change.get<0>(), mol.cljgroup.cljBoxes() );
                    
                    //get the energy of the old atoms interacting with each other
                    oldnrg = cljcalc.calculate( mol.cljfuncs,
                                                change.get<1>() );
                    
                    //get the energy of the new atoms interacting with each other
                    newnrg = cljcalc.calculate( mol.cljfuncs,
                                                change.get<2>() );
                }
                else
                {
                    //get the change in energy from the changed atoms with the
                    //atoms that have not changed
                    delta = CLJFunction::multiCalculate( mol.cljfuncs,
                                                         change.get<0>(),
                                                         mol.cljgroup.cljBoxes() );
                    
                    //get the energy of the old atoms interacting with each other
                    oldnrg = CLJFunction::multiCalculate( mol.cljfuncs, change.get<1>() );
                    
                    //get the energy of the new atoms interacting with each other
                    newnrg = CLJFunction::multiCalculate( mol.cljfuncs, change.get<2>() );
                }
                
                mol.nrg += MultiCLJEnergy( delta.get<0>(), delta.get<1>() ) +
                           MultiCLJEnergy( newnrg.get<0>(), newnrg.get<1>() ) -
                           MultiCLJEnergy( oldnrg.get<0>(), oldnrg.get<1>() );
                
                mol.needs_energy_calc = false;
                mol.needs_accepting = true;
                
                return mol.nrg;
            }
        }
        else
        {
            //something weird - recalculate from scratch
            mol.mustNowRecalculateFromScratch();
            return calcEnergy(mol);
        }
    }
    else
    {
        //something weird - recalculate from scratch
        mol.mustNowRecalculateFromScratch();
        return calcEnergy(mol);
    }
}

/** Recalculate the energy of this forcefield */
void IntraFF::recalculateEnergy()
{
    MultiCLJEnergy cljenergy;
    
    if (d.constData()->cljfuncs.count() == 1)
    {
        cljenergy = MultiCLJEnergy(0,0);
    }
    else
    {
        QVector<double> zero( d.constData()->cljfuncs.count(), 0.0 );
        cljenergy = MultiCLJEnergy(zero, zero);
    }
    
    for (MolData::iterator it = moldata.begin();
         it != moldata.end();
         ++it)
    {
        if (it.value().constData()->needs_energy_calc)
        {
            cljenergy += this->calcEnergy(*(it.value()));
            needs_accepting = needs_accepting or it.value().constData()->needs_accepting;
        }
        else
        {
            cljenergy += it.value().constData()->nrg;
        }
    }

    d.constData()->cljcomps.setEnergy(*this, cljenergy);
    setClean();
}

/** Function called to add a molecule to this forcefield */
void IntraFF::_pvt_added(const SireMol::PartialMolecule &mol, const SireBase::PropertyMap &map)
{
    MolData::iterator it = moldata.find(mol.number());
    
    if (it == moldata.end())
    {
        //need to add this as a new molecule
        moldata.insert( mol.number(), SharedDataPointer<detail::IntraFFMolData>(
                                new detail::IntraFFMolData(mol, map, d.constData()->cljfuncs) ) );
    }
    else
    {
        if (it.value()->needs_accepting)
        {
            it.value()->cljgroup.accept();
            it.value()->needs_accepting = false;
        }
        
        it.value()->cljgroup.add(mol, map);
        it.value()->checkForChangeInConnectivity(mol);
        it.value()->needs_energy_calc = true;
    }

    setDirty();
}

/** Function called to remove a molecule from this forcefield */
void IntraFF::_pvt_removed(const SireMol::PartialMolecule &mol)
{
    MolData::iterator it = moldata.find(mol.number());
    
    if (it != moldata.end())
    {
        if (it.value()->needs_accepting)
        {
            it.value()->cljgroup.accept();
            it.value()->needs_accepting = false;
        }
        
        it.value()->cljgroup.remove(mol);
        
        if (it.value()->cljgroup.isEmpty())
        {
            moldata.remove(mol.number());
        }
        else
        {
            it.value()->checkForChangeInConnectivity(mol);
            it.value()->needs_energy_calc = true;
        }

        setDirty();
    }
}

/** Function called to indicate that the passed molecule has changed */
void IntraFF::_pvt_changed(const Molecule &molecule, bool auto_update)
{
    MolData::iterator it = moldata.find(molecule.number());
    
    if (it != moldata.end())
    {
        if (it.value()->needs_accepting)
        {
            it.value()->cljgroup.accept();
            it.value()->needs_accepting = false;
        }

        it.value()->cljgroup.update(molecule);
        it.value()->checkForChangeInConnectivity(molecule);
        it.value()->needs_energy_calc = true;

        setDirty();
    }
}

/** Function called to indicate that a list of molecules in this forcefield have changed */
void IntraFF::_pvt_changed(const QList<SireMol::Molecule> &molecules, bool auto_update)
{
    foreach (const Molecule &molecule, molecules)
    {
        this->_pvt_changed(molecule, auto_update);
    }
}

/** Function called to indicate that all molecules in this forcefield have been removed */
void IntraFF::_pvt_removedAll()
{
    moldata.clear();
    this->setDirty();
}

/** Function called to query whether or not a change in source properties would
    change the properties needed by this forcefield for the molecule with number 'molnum' */
bool IntraFF::_pvt_wouldChangeProperties(SireMol::MolNum molnum,
                                         const SireBase::PropertyMap &map) const
{
    MolData::const_iterator it = moldata.constFind(molnum);
    
    if (it != moldata.constEnd())
    {
        return it.value()->cljgroup.mapForMolecule(molnum) != map;
    }
    else
    {
        return false;
    }
}

/** Return whether or not this forcefield is using a temporary workspace that 
    needs to be accepted */
bool IntraFF::needsAccepting() const
{
    return needs_accepting or G1FF::needsAccepting();
}

/** Tell the forcefield that the last move was accepted. This tells the
    forcefield to make permanent any temporary changes that were used a workspace
    to avoid memory allocation during a move */
void IntraFF::accept()
{
    if (needs_accepting)
    {
        for (MolData::iterator it = moldata.begin();
             it != moldata.end();
             ++it)
        {
            if (it.value()->needs_accepting)
            {
                it.value()->cljgroup.accept();
                it.value()->needs_accepting = false;
            }
        }
        
        needs_accepting = false;
    }
    
    G1FF::accept();
}
