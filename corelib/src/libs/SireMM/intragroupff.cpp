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

#include "intragroupff.h"

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
                        const SharedDataPointer<SireMM::detail::IntraGroupFFMolData> &ptr)
{
    SharedDataStream sds(ds);
    sds << ptr;
    return ds;
}

QDataStream& operator>>(QDataStream &ds,
                        SharedDataPointer<SireMM::detail::IntraGroupFFMolData> &ptr)
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
        class IntraGroupFFData : public RefCountData
        {
        public:
            IntraGroupFFData() : RefCountData(),
                                 parallel_calc(true), repro_sum(false)
            {}
            
            IntraGroupFFData(const IntraGroupFFData &other)
                 : RefCountData(),
                   cljfuncs(other.cljfuncs),
                   cljcomps(other.cljcomps),
                   props(other.props),
                   parallel_calc(other.parallel_calc),
                   repro_sum(other.repro_sum)
            {}
            
            ~IntraGroupFFData()
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
        class IntraGroupFFMolData : public RefCountData
        {
        public:
            /** The CLJGroup for this molecule in group 0 */
            CLJGroup cljgroup0;
            
            /** The CLJGroup for this molecule in group 1 */
            CLJGroup cljgroup1;
            
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

            IntraGroupFFMolData() : RefCountData(), connectivity_version(0),
                                    needs_energy_calc(false), needs_accepting(false)
            {}
            
            IntraGroupFFMolData(quint32 group_id,
                                const MoleculeView &molview,
                                const PropertyMap &map,
                                const QVector<CLJFunctionPtr> &funcs)
                    : RefCountData(), needs_energy_calc(true), needs_accepting(false)
            {
                cljgroup0 = CLJGroup(CLJAtoms::USE_ATOMIDX, CLJExtractor::EXTRACT_BY_CUTGROUP);
                cljgroup1 = CLJGroup(CLJAtoms::USE_ATOMIDX, CLJExtractor::EXTRACT_BY_CUTGROUP);
                
                cljgroup0.setBoxLength( 7.5*angstrom );
                cljgroup1.setBoxLength( 7.5*angstrom );
                
                if (group_id == 0)
                    cljgroup0.add(molview, map);
                else
                    cljgroup1.add(molview, map);

                connectivity_property = map["connectivity"];
                
                cty = molview.data().property(connectivity_property).asA<Connectivity>();
                connectivity_version = molview.data().version(connectivity_property);
                
                setCLJFunctions(funcs);
            }
            
            IntraGroupFFMolData(const IntraGroupFFMolData &other)
                : RefCountData(), cljgroup0(other.cljgroup0),
                  cljgroup1(other.cljgroup1), cljfuncs(other.cljfuncs),
                  cty(other.cty), nrg(other.nrg),
                  connectivity_property(other.connectivity_property),
                  connectivity_version(other.connectivity_version),
                  needs_energy_calc(other.needs_energy_calc),
                  needs_accepting(other.needs_accepting)
            {}
            
            ~IntraGroupFFMolData()
            {}
            
            void assertMoleculeCompatible(const MoleculeView &molview,
                                          const CLJGroup &cljgroup) const
            {
                if (cljgroup.isEmpty())
                    return;
                
                Molecules mols = cljgroup.molecules();
                
                if (mols.nMolecules() != 1 or mols.nViews() != 1)
                {
                    throw SireError::program_bug( QObject::tr(
                            "Strange - the IntraGroupFF cljgroup should only contain "
                            "a single molecule with a single view. Instead, we see "
                            "nMolecules() = %1 and nViews() = %2.")
                                .arg(mols.nMolecules())
                                .arg(mols.nViews()), CODELOC );
                }
                
                const MoleculeData &moldata = mols.constBegin().value().data();
                
                if (moldata.number() != molview.data().number() or
                    moldata.version() != molview.data().version())
                {
                    throw SireError::incompatible_error( QObject::tr(
                            "You are trying to add a second part of the molecule "
                            "that is not the *same* molecule as that added from the "
                            "first part. The first part added was molecule %1, version %2, "
                            "while the second part is molecule %3, version %4.")
                                .arg(molview.data().number().toString())
                                .arg(molview.data().version())
                                .arg(moldata.number().toString())
                                .arg(moldata.version()), CODELOC );
                }
            }
            
            void assertConnectivityCompatible(const PropertyMap &map) const
            {
                PropertyName new_connectivity = map["connectivity"];
                
                if (new_connectivity != connectivity_property)
                {
                    throw SireError::incompatible_error( QObject::tr(
                            "You cannot add the second part of the molecule with a "
                            "different connectivity property! The first part was added "
                            "with a 'connectivity' property equal to '%1', while you are "
                            "now adding the second part of the molecule with a "
                            "'connectivity' property equal to '%2'.")
                                .arg(connectivity_property.toString())
                                .arg(new_connectivity.toString()), CODELOC );
                }
            }
            
            void add(quint32 group_id,
                     const MoleculeView &molview,
                     const PropertyMap &map)
            {
                if (group_id == 0)
                {
                    if (cljgroup0.isEmpty())
                    {
                        //check that this molecule is compatible with that in cljgroup1
                        assertMoleculeCompatible(molview, cljgroup1);
                        assertConnectivityCompatible(map);
                    }
                    
                    cljgroup0.add(molview, map);
                }
                else
                {
                    if (cljgroup1.isEmpty())
                    {
                        //check that this molecule is compatible with that in cljgroup0
                        assertMoleculeCompatible(molview, cljgroup0);
                        assertConnectivityCompatible(map);
                    }

                    cljgroup1.add(molview, map);
                }
            }
            
            void mustNowRecalculateFromScratch()
            {
                cljgroup0.mustRecalculateFromScratch();
                cljgroup1.mustRecalculateFromScratch();
                needs_accepting = false;
                needs_energy_calc = true;
                nrg = MultiCLJEnergy(0,0);
            }
            
            void mustReallyRecalculateFromScratch()
            {
                cljgroup0.mustReallyRecalculateFromScratch();
                cljgroup1.mustReallyRecalculateFromScratch();
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
                        mustNowRecalculateFromScratch();
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

QDataStream& operator<<(QDataStream &ds, const SireMM::detail::IntraGroupFFMolData &intraff)
{
    quint32 version = 1;
    
    ds << version;
    
    SharedDataStream sds(ds);
    
    sds << intraff.cljgroup0 << intraff.cljgroup1 << intraff.cljfuncs << intraff.cty
        << intraff.nrg << intraff.connectivity_property
        << intraff.connectivity_version << intraff.needs_energy_calc
        << intraff.needs_accepting;
    
    return ds;
}

QDataStream& operator>>(QDataStream &ds, SireMM::detail::IntraGroupFFMolData &intraff)
{
    quint32 version;
    
    ds >> version;
    
    if (version == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> intraff.cljgroup0 >> intraff.cljgroup1 >> intraff.cljfuncs >> intraff.cty
            >> intraff.nrg >> intraff.connectivity_property
            >> intraff.connectivity_version >> intraff.needs_energy_calc
            >> intraff.needs_accepting;
    }
    else
        throw version_error( QObject::tr(
                "Wrong version of SireMM::detail::IntraGroupFFData. Version %1, while "
                "only version 1 is supported.").arg(version), CODELOC );
    
    return ds;
}

static const RegisterMetaType<IntraGroupFF> r_intragroupff;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const IntraGroupFF &intraff)
{
    writeHeader(ds, r_intragroupff, 1);
    
    SharedDataStream sds(ds);
    
    sds << intraff.moldata
        << intraff.needs_accepting
        << intraff.d->cljfuncs
        << intraff.d->cljcomps
        << intraff.d->parallel_calc
        << intraff.d->repro_sum
        << static_cast<const G2FF&>(intraff);
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, IntraGroupFF &intraff)
{
    VersionID v = readHeader(ds, r_intragroupff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> intraff.moldata
            >> intraff.needs_accepting
            >> intraff.d->cljfuncs
            >> intraff.d->cljcomps
            >> intraff.d->parallel_calc
            >> intraff.d->repro_sum
            >> static_cast<G2FF&>(intraff);
    }
    else
        throw version_error(v, "1", r_intragroupff, CODELOC);
    
    return ds;
}


/** Constructor */
IntraGroupFF::IntraGroupFF()
             : ConcreteProperty<IntraGroupFF,G2FF>(), needs_accepting(false)
{
    d = new detail::IntraGroupFFData();
    this->_pvt_updateName();
    this->setCLJFunction( CLJIntraShiftFunction::defaultShiftFunction()
                                    .read().asA<CLJIntraFunction>() );
}

/** Construct, specifying the name of the forcefield */
IntraGroupFF::IntraGroupFF(const QString &name)
             : ConcreteProperty<IntraGroupFF, G2FF>(),
               needs_accepting(false)
{
    d = new detail::IntraGroupFFData();
    G2FF::setName(name);
    this->setCLJFunction( CLJIntraShiftFunction::defaultShiftFunction()
                                    .read().asA<CLJIntraFunction>() );
}

/** Copy constructor */
IntraGroupFF::IntraGroupFF(const IntraGroupFF &other)
             : ConcreteProperty<IntraGroupFF,G2FF>(other),
               moldata(other.moldata), d(other.d),
               needs_accepting(other.needs_accepting)
{}

/** Destructor */
IntraGroupFF::~IntraGroupFF()
{}

/** Function used to set the CLJIntraFunction used to calculate 
    the intramolecular energy */
void IntraGroupFF::setCLJFunction(const CLJIntraFunction &func)
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
const CLJIntraFunction& IntraGroupFF::cljFunction() const
{
    if (d.constData()->cljfuncs.isEmpty())
        throw SireError::program_bug( QObject::tr(
                "There should always be at least one CLJFunction in InterGroupFF!"),
                    CODELOC );
    
    return d.constData()->cljfuncs.at(0).read().asA<CLJIntraFunction>();
}

/** Set the CLJFunction with key 'key' equal to 'cljfunc' */
void IntraGroupFF::setCLJFunction(QString key, const CLJIntraFunction &cljfunc)
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
void IntraGroupFF::removeCLJFunctionAt(QString key)
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
void IntraGroupFF::removeAllCLJFunctions()
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
QStringList IntraGroupFF::cljFunctionKeys() const
{
    return d.constData()->cljcomps.keys();
}

/** Return the CLJFunction associated with the passed key */
const CLJIntraFunction& IntraGroupFF::cljFunction(QString key) const
{
    return d.constData()->cljfuncs.at( d.constData()->cljcomps.indexOf(key) )
                                  .read().asA<CLJIntraFunction>();
}

/** Return the number of CLJ functions in this forcefield. There should always
    be at least one */
int IntraGroupFF::nCLJFunctions() const
{
    return d.constData()->cljfuncs.count();
}

/** Return the hash of all CLJFunctions in this forcefield, indexed by their key */
QHash<QString,CLJFunctionPtr> IntraGroupFF::cljFunctions() const
{
    QHash<QString,CLJFunctionPtr> funcs;
    
    foreach (QString key, d->cljcomps.keys())
    {
        funcs.insert( key, this->cljFunction(key) );
    }
    
    return funcs;
}

/** Internal function called when the name of the forcefield changes */
void IntraGroupFF::_pvt_updateName()
{
    d->cljcomps = d->cljcomps.rename(this->name());
    G2FF::_pvt_updateName();
}

const char* IntraGroupFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<IntraGroupFF>() );
}

const char* IntraGroupFF::what() const
{
    return IntraGroupFF::typeName();
}

/** Copy assignment operator */
IntraGroupFF& IntraGroupFF::operator=(const IntraGroupFF &other)
{
    if (this != &other)
    {
        moldata = other.moldata;
        needs_accepting = other.needs_accepting;
        d = other.d;
        G2FF::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool IntraGroupFF::operator==(const IntraGroupFF &other) const
{
    return (this == &other) or
           (G2FF::operator==(other) and
            d->cljfuncs == other.d->cljfuncs and d->cljcomps == other.d->cljcomps and
            moldata == other.moldata and needs_accepting == other.needs_accepting);
}

/** Comparison operator */
bool IntraGroupFF::operator!=(const IntraGroupFF &other) const
{
    return not operator==(other);
}

IntraGroupFF* IntraGroupFF::clone() const
{
    return new IntraGroupFF(*this);
}

/** Return the energy components of this forcefield */
const MultiCLJComponent& IntraGroupFF::components() const
{
    return d->cljcomps;
}

/** Internal function used to rebuild the properties object that 
    stores all of the properties of this forcefield */
void IntraGroupFF::rebuildProps()
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
bool IntraGroupFF::setProperty(const QString &name, const Property &property)
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
const Property& IntraGroupFF::property(const QString &name) const
{
    return d->props.property(name);
}

/** Return whether or not this forcefield contains the property 'property' */
bool IntraGroupFF::containsProperty(const QString &name) const
{
    return d->props.hasProperty(name);
}

/** Return all of the properties of this function */
const Properties& IntraGroupFF::properties() const
{
    return d->props;
}

/** Set whether or not to use a multicore parallel algorithm
    to calculate the energy */
void IntraGroupFF::setUseParallelCalculation(bool on)
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
void IntraGroupFF::enableParallelCalculation()
{
    this->setUseParallelCalculation(true);
}

/** Turn off use of a multicore parallel calculation of the energy.
    This may be quicker if you have few atoms in the forcefield,
    or if you are only planning on allocating one core per forcefield */
void IntraGroupFF::disableParallelCalculation()
{
    this->setUseParallelCalculation(false);
}

/** Return whether or not a parallel algorithm is used to calculate energies */
bool IntraGroupFF::usesParallelCalculation() const
{
    return d->parallel_calc;
}

/** Turn on an energy summing algorithm that guarantees the same energy
    regardless of whether a single core or multicore calculation is being
    performed (i.e. rounding errors in both cases will be identical) */
void IntraGroupFF::enableReproducibleCalculation()
{
    setUseReproducibleCalculation(true);
}

/** Turn off an energy summing algorithm that guarantees the same energy
    regardless of whether a single core or multicore calculation is being
    performed (i.e. rounding errors in both cases will not be identical) */
void IntraGroupFF::disableReproducibleCalculation()
{
    setUseReproducibleCalculation(false);
}

/** Switch on or off use of an energy summing algorithm that guarantees the 
    same energy regardless of whether a single core or multicore calculation 
    is being performed */
void IntraGroupFF::setUseReproducibleCalculation(bool on)
{
    if (on != d->repro_sum)
    {
        d->repro_sum = on;
        d->props.setProperty("reproducibleCalculation", BooleanProperty(on));
    }
}

/** Return whether or not a reproducible energy summing algorithm is being
    used to accumulate the energies */
bool IntraGroupFF::usesReproducibleCalculation() const
{
    return d->repro_sum;
}

/** Signal that this forcefield must now be recalculated from scratch */
void IntraGroupFF::mustNowRecalculateFromScratch()
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
void IntraGroupFF::mustNowReallyRecalculateFromScratch()
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
MultiCLJEnergy IntraGroupFF::calcEnergy(detail::IntraGroupFFMolData &mol) const
{
    if (mol.needs_energy_calc)
    {
        if (mol.cljgroup0.recalculatingFromScratch() or mol.cljgroup1.recalculatingFromScratch())
        {
            if (mol.needs_accepting)
            {
                mol.cljgroup0.accept();
                mol.cljgroup1.accept();
                mol.needs_accepting = false;
            }
            
            if (mol.cljgroup0.needsAccepting())
            {
                mol.cljgroup0.accept();
            }

            if (mol.cljgroup1.needsAccepting())
            {
                mol.cljgroup1.accept();
            }
        
            if (mol.cljfuncs.count() == 1)
            {
                tuple<double,double> nrg;
            
                if (usesParallelCalculation())
                {
                    CLJCalculator cljcalc( usesReproducibleCalculation() );
                    
                    nrg = cljcalc.calculate( mol.cljfuncs.at(0).read(),
                                             mol.cljgroup0.cljBoxes(),
                                             mol.cljgroup1.cljBoxes() );
                }
                else
                {
                    nrg = mol.cljfuncs.at(0).read().calculate( mol.cljgroup0.cljBoxes(),
                                                               mol.cljgroup1.cljBoxes() );
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
                    
                    nrg = cljcalc.calculate( mol.cljfuncs, mol.cljgroup0.cljBoxes(),
                                                           mol.cljgroup1.cljBoxes() );
                }
                else
                {
                    nrg = CLJFunction::multiCalculate( mol.cljfuncs,
                                                       mol.cljgroup0.cljBoxes(),
                                                       mol.cljgroup1.cljBoxes() );
                }
                
                mol.nrg = MultiCLJEnergy( nrg.get<0>(), nrg.get<1>() );
                mol.needs_energy_calc = false;
                return mol.nrg;
            }
        }
        else if (mol.cljgroup0.needsAccepting() and mol.cljgroup1.needsAccepting())
        {
            //they both need accepting - just recalculate from scratch for the moment
            //until I have time or need to write a more optimised calculation
            mol.mustNowRecalculateFromScratch();
            return this->calcEnergy(mol);
        }
        else if (mol.cljgroup0.needsAccepting())
        {
            //we can calculate the energy using a delta
            CLJAtoms changed_atoms = mol.cljgroup0.changedAtoms();
            
            if (mol.cljfuncs.count() == 1)
            {
                tuple<double,double> delta;

                if (usesParallelCalculation())
                {
                    CLJCalculator cljcalc(usesReproducibleCalculation());
                    
                    delta = cljcalc.calculate( mol.cljfuncs.at(0).read(),
                                               changed_atoms, mol.cljgroup1.cljBoxes() );
                }
                else
                {
                    delta = mol.cljfuncs.at(0).read().calculate(
                                               changed_atoms, mol.cljgroup1.cljBoxes() );
                }
                
                mol.nrg += MultiCLJEnergy(delta.get<0>(),delta.get<1>());
                
                mol.needs_energy_calc = false;
                mol.needs_accepting = true;
                
                return mol.nrg;
            }
            else  // delta calculation with more than one clj function
            {
                tuple< QVector<double>,QVector<double> > delta;
                
                if (usesParallelCalculation())
                {
                    CLJCalculator cljcalc(usesReproducibleCalculation());
                    
                    delta = cljcalc.calculate( mol.cljfuncs,
                                               changed_atoms, mol.cljgroup1.cljBoxes() );
                }
                else
                {
                    delta = CLJFunction::multiCalculate( mol.cljfuncs,
                                                         changed_atoms,
                                                         mol.cljgroup1.cljBoxes() );
                }
                
                mol.nrg += MultiCLJEnergy( delta.get<0>(), delta.get<1>() );
                
                mol.needs_energy_calc = false;
                mol.needs_accepting = true;
                
                return mol.nrg;
            }
        }
        else if (mol.cljgroup1.needsAccepting())
        {
            //we can calculate the energy using a delta
            CLJAtoms changed_atoms = mol.cljgroup1.changedAtoms();
            
            if (mol.cljfuncs.count() == 1)
            {
                tuple<double,double> delta;

                if (usesParallelCalculation())
                {
                    CLJCalculator cljcalc(usesReproducibleCalculation());
                    
                    delta = cljcalc.calculate( mol.cljfuncs.at(0).read(),
                                               changed_atoms, mol.cljgroup0.cljBoxes() );
                }
                else
                {
                    delta = mol.cljfuncs.at(0).read().calculate(
                                               changed_atoms, mol.cljgroup0.cljBoxes() );
                }
                
                mol.nrg += MultiCLJEnergy(delta.get<0>(),delta.get<1>());
                
                mol.needs_energy_calc = false;
                mol.needs_accepting = true;
                
                return mol.nrg;
            }
            else  // delta calculation with more than one clj function
            {
                tuple< QVector<double>,QVector<double> > delta;
                
                if (usesParallelCalculation())
                {
                    CLJCalculator cljcalc(usesReproducibleCalculation());
                    
                    delta = cljcalc.calculate( mol.cljfuncs,
                                               changed_atoms, mol.cljgroup0.cljBoxes() );
                }
                else
                {
                    delta = CLJFunction::multiCalculate( mol.cljfuncs,
                                                         changed_atoms,
                                                         mol.cljgroup0.cljBoxes() );
                }
                
                mol.nrg += MultiCLJEnergy( delta.get<0>(), delta.get<1>() );
                
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
void IntraGroupFF::recalculateEnergy()
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
void IntraGroupFF::_pvt_added(quint32 group_id,
                              const SireMol::PartialMolecule &mol,
                              const SireBase::PropertyMap &map)
{
    MolData::iterator it = moldata.find(mol.number());
    
    if (it == moldata.end())
    {
        //need to add this as a new molecule
        moldata.insert( mol.number(), SharedDataPointer<detail::IntraGroupFFMolData>(
                                new detail::IntraGroupFFMolData(group_id, mol,
                                                                map, d.constData()->cljfuncs) ) );
    }
    else
    {
        if (it.value()->needs_accepting)
        {
            it.value()->cljgroup0.accept();
            it.value()->cljgroup1.accept();
            it.value()->needs_accepting = false;
        }
        
        it.value()->add(group_id, mol, map);
        it.value()->checkForChangeInConnectivity(mol);
        it.value()->needs_energy_calc = true;
    }

    setDirty();
}

/** Function called to remove a molecule from this forcefield */
void IntraGroupFF::_pvt_removed(quint32 group_id, const SireMol::PartialMolecule &mol)
{
    MolData::iterator it = moldata.find(mol.number());
    
    if (it != moldata.end())
    {
        if (it.value()->needs_accepting)
        {
            it.value()->cljgroup0.accept();
            it.value()->cljgroup1.accept();
            it.value()->needs_accepting = false;
        }
        
        if (group_id == 0)
            it.value()->cljgroup0.remove(mol);
        else
            it.value()->cljgroup1.remove(mol);
        
        if (it.value()->cljgroup0.isEmpty() and it.value()->cljgroup1.isEmpty())
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
void IntraGroupFF::_pvt_changed(quint32 group_id,
                                const Molecule &molecule, bool auto_update)
{
    MolData::iterator it = moldata.find(molecule.number());
    
    if (it != moldata.end())
    {
        if (it.value()->needs_accepting)
        {
            it.value()->cljgroup0.accept();
            it.value()->cljgroup1.accept();
            it.value()->needs_accepting = false;
        }

        if (group_id == 0)
            it.value()->cljgroup0.update(molecule);
        else
            it.value()->cljgroup1.update(molecule);
        
        it.value()->checkForChangeInConnectivity(molecule);
        it.value()->needs_energy_calc = true;

        setDirty();
    }
}

/** Function called to indicate that a list of molecules in this forcefield have changed */
void IntraGroupFF::_pvt_changed(quint32 group_id,
                                const QList<SireMol::Molecule> &molecules, bool auto_update)
{
    foreach (const Molecule &molecule, molecules)
    {
        this->_pvt_changed(group_id, molecule, auto_update);
    }
}

/** Function called to indicate that all molecules in this forcefield have been removed */
void IntraGroupFF::_pvt_removedAll(quint32 group_id)
{
    QSet<MolNum> mols_to_remove;

    for (MolData::iterator it = moldata.begin();
         it != moldata.end();
         ++it)
    {
        if (group_id == 0)
        {
            it.value()->cljgroup0.removeAll();
            if (it.value()->cljgroup1.isEmpty())
            {
                mols_to_remove.insert(it.key());
            }
        }
        else
        {
            it.value()->cljgroup1.removeAll();
            if (it.value()->cljgroup0.isEmpty())
            {
                mols_to_remove.insert(it.key());
            }
        }
    }

    foreach (MolNum molnum, mols_to_remove)
    {
        moldata.remove(molnum);
    }

    this->setDirty();
}

/** Function called to query whether or not a change in source properties would
    change the properties needed by this forcefield for the molecule with number 'molnum' */
bool IntraGroupFF::_pvt_wouldChangeProperties(quint32 group_id,
                                              SireMol::MolNum molnum,
                                              const SireBase::PropertyMap &map) const
{
    MolData::const_iterator it = moldata.constFind(molnum);
    
    if (it != moldata.constEnd())
    {
        if (group_id == 0)
            return it.value()->cljgroup0.mapForMolecule(molnum) != map;
        else
            return it.value()->cljgroup1.mapForMolecule(molnum) != map;
    }
    else
    {
        return false;
    }
}

/** Return whether or not this forcefield is using a temporary workspace that 
    needs to be accepted */
bool IntraGroupFF::needsAccepting() const
{
    return needs_accepting or G2FF::needsAccepting();
}

/** Tell the forcefield that the last move was accepted. This tells the
    forcefield to make permanent any temporary changes that were used a workspace
    to avoid memory allocation during a move */
void IntraGroupFF::accept()
{
    if (needs_accepting)
    {
        for (MolData::iterator it = moldata.begin();
             it != moldata.end();
             ++it)
        {
            if (it.value()->needs_accepting)
            {
                it.value()->cljgroup0.accept();
                it.value()->cljgroup1.accept();
                it.value()->needs_accepting = false;
            }
        }
        
        needs_accepting = false;
    }
    
    G2FF::accept();
}
