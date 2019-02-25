/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#include "titrator.h"

#include "SireError/errors.h"
#include "SireBase/errors.h"

#include "SireID/index.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/molidx.h"
#include "SireMol/atomcharges.h"
#include "SireMol/evaluator.h"
#include "SireMol/editor.hpp"
#include "SireMol/mover.hpp"

#include "SireMaths/rangenerator.h"

#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

using namespace SireMove;
using namespace SireMol;
using namespace SireFF;
using namespace SireID;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Titrator> r_titrator;

QDataStream &operator<<(QDataStream &ds, const Titrator &titrator)
{
    writeHeader(ds, r_titrator, 1);
    
    SharedDataStream sds(ds);
    
    sds << titrator.mgname << titrator.mgnum << titrator.mgversion
        << titrator.chgs << titrator.desired_chgs
        << titrator.neutral_template << titrator.negative_template
        << titrator.positive_template << titrator.neutral_properties
        << titrator.negative_properties << titrator.positive_properties
        << titrator.neutral_map << titrator.negative_map
        << titrator.positive_map << titrator.propmap
        << titrator.pos_charge << titrator.neg_charge
        << static_cast<const Property&>(titrator);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, Titrator &titrator)
{
    VersionID v = readHeader(ds, r_titrator);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> titrator.mgname >> titrator.mgnum >> titrator.mgversion
            >> titrator.chgs >> titrator.desired_chgs
            >> titrator.neutral_template >> titrator.negative_template
            >> titrator.positive_template >> titrator.neutral_properties
            >> titrator.negative_properties >> titrator.positive_properties
            >> titrator.neutral_map >> titrator.negative_map
            >> titrator.positive_map >> titrator.propmap
            >> titrator.pos_charge >> titrator.neg_charge
            >> static_cast<Property&>(titrator);
    }
    else
        throw version_error(v, "1", r_titrator, CODELOC);
    
    return ds;
}

/** Constructor */
Titrator::Titrator()
         : ConcreteProperty<Titrator,Property>(), pos_charge(0), neg_charge(0)
{}

/** Construct to titrate the molecules in the passed molecule group */
Titrator::Titrator(const MoleculeGroup &group)
         : ConcreteProperty<Titrator,Property>(), mgname(group.name()),
           pos_charge(0), neg_charge(0)
{}

/** Copy constructor */
Titrator::Titrator(const Titrator &other)
         : ConcreteProperty<Titrator,Property>(other),
           mgname(other.mgname), mgnum(other.mgnum),
           mgversion(other.mgversion), chgs(other.chgs),
           desired_chgs(other.desired_chgs), neutral_template(other.neutral_template),
           negative_template(other.negative_template),
           positive_template(other.positive_template),
           neutral_properties(other.neutral_properties),
           negative_properties(other.negative_properties),
           positive_properties(other.positive_properties),
           neutral_map(other.neutral_map), negative_map(other.negative_map),
           positive_map(other.positive_map), propmap(other.propmap),
           pos_charge(other.pos_charge), neg_charge(other.neg_charge)
{}

/** Destructor */
Titrator::~Titrator()
{}

const char* Titrator::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Titrator>() );
}

const char* Titrator::what() const
{
    return Titrator::typeName();
}

/** Copy assignment operator */
Titrator& Titrator::operator=(const Titrator &other)
{
    if (this != &other)
    {
        mgname = other.mgname;
        mgnum = other.mgnum;
        mgversion = other.mgversion;
        chgs = other.chgs;
        desired_chgs = other.desired_chgs;
        neutral_template = other.neutral_template;
        negative_template = other.negative_template;
        positive_template = other.positive_template;
        neutral_properties = other.neutral_properties;
        negative_properties = other.negative_properties;
        positive_properties = other.positive_properties;
        neutral_map = other.neutral_map;
        negative_map = other.negative_map;
        positive_map = other.positive_map;
        propmap = other.propmap;
        pos_charge = other.pos_charge;
        neg_charge = other.neg_charge;
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool Titrator::operator==(const Titrator &other) const
{
    if (this == &other)
        return true;

    return mgname == other.mgname and
           mgnum == other.mgnum and
           mgversion == other.mgversion and
           chgs == other.chgs and
           desired_chgs == other.desired_chgs and
           neutral_template.equals(other.neutral_template) and
           negative_template.equals(other.negative_template) and
           positive_template.equals(other.positive_template) and
           neutral_properties == other.neutral_properties and
           negative_properties == other.negative_properties and
           positive_properties == other.positive_properties and
           neutral_map == other.neutral_map and
           negative_map == other.negative_map and
           positive_map == other.positive_map and
           propmap == other.propmap and
           pos_charge == other.pos_charge and
           neg_charge == other.neg_charge and
           Property::operator==(other);
}

/** Comparison operator */
bool Titrator::operator!=(const Titrator &other) const
{
    return not Titrator::operator==(other);
}

/** Clear the current state */
void Titrator::clearState()
{
    mgnum = MGNum();
    mgversion = Version();
    
    chgs.clear();
    desired_chgs.clear();
}

/** Set the molecule group containing the molecules whose charge
    state will be changed */
void Titrator::setMoleculeGroup(const MoleculeGroup &group)
{
    if (group.name() != mgname)
    {
        mgname = group.name();
        clearState();
    }
}

/** Set the template for the positive ion state. The collection of properties
    listed in "properties" are copied from "positive_ion" (once mapped by
    "map") to the molecule to put it into the "positive_ion" state */
void Titrator::setPositiveTemplate(const Molecule &positive_ion,
                                   const QStringList &properties,
                                   const PropertyMap &map)
{
    foreach (QString property, properties)
    {
        if (not positive_ion.hasProperty(map[property]))
        {
            throw SireBase::missing_property( QObject::tr(
                    "You have asked to use the property %1 (mapped from %2) from the "
                    "passed positive ion template, but the only available properties are %3.")
                        .arg(map[property].toString(), property)
                        .arg(Sire::toString(positive_ion.propertyKeys())), CODELOC );
        }
    }

    //check also that there is a charge property of type "AtomCharges"
    const Property &p = positive_ion.property( map["charge"] );
    
    if (not p.isA<AtomCharges>())
    {
        throw SireBase::missing_property( QObject::tr(
                "The positive ion template is missing the AtomCharges charge property..."),
                    CODELOC );
    }
    
    double chg = positive_ion.evaluate().charge(map);
    
    if (chg <= 0)
        throw SireError::invalid_arg( QObject::tr(
                "The charge on the positive ion template is not positive! %1")
                    .arg(chg), CODELOC );
    
    pos_charge = qint32( chg+0.5 );

    positive_template = positive_ion;
    positive_properties = properties;
    positive_map = map;
    
    clearState();
}

/** Set the template for the negative ion state. The collection of properties
    listed in "properties" are copied from "negative_ion" (once mapped by
    "map") to the molecule to put it into the "negative_ion" state */
void Titrator::setNegativeTemplate(const Molecule &negative_ion,
                                   const QStringList &properties,
                                   const PropertyMap &map)
{
    foreach (QString property, properties)
    {
        if (not negative_ion.hasProperty(map[property]))
        {
            throw SireBase::missing_property( QObject::tr(
                    "You have asked to use the property %1 (mapped from %2) from the "
                    "passed negative ion template, but the only available properties are %3.")
                        .arg(map[property].toString(), property)
                        .arg(Sire::toString(negative_ion.propertyKeys())), CODELOC );
        }
    }

    //check also that there is a charge property of type "AtomCharges"
    const Property &p = negative_ion.property( map["charge"] );
    
    if (not p.isA<AtomCharges>())
    {
        throw SireBase::missing_property( QObject::tr(
                "The negative ion template is missing the AtomCharges charge property..."),
                    CODELOC );
    }

    double chg = negative_ion.evaluate().charge(map);
    
    if (chg >= 0)
        throw SireError::invalid_arg( QObject::tr(
                "The charge on the negative ion template is not negative! %1")
                    .arg(chg), CODELOC );
    
    neg_charge = qint32( chg-0.5 );

    negative_template = negative_ion;
    negative_properties = properties;
    negative_map = map;
    
    clearState();
}

/** Set the template for the neutral state. The collection of properties
    listed in "properties" are copied from "neutral_mol" (once mapped by
    "map") to the molecule to put it into the "neutral" state */
void Titrator::setNeutralTemplate(const Molecule &neutral_mol,
                                  const QStringList &properties,
                                  const PropertyMap &map)
{
    foreach (QString property, properties)
    {
        if (not neutral_mol.hasProperty(map[property]))
        {
            throw SireBase::missing_property( QObject::tr(
                    "You have asked to use the property %1 (mapped from %2) from the "
                    "passed neutral molecule template, but the only available properties are %3.")
                        .arg(map[property].toString(), property)
                        .arg(Sire::toString(neutral_mol.propertyKeys())), CODELOC );
        }
    }

    //check also that there is a charge property of type "AtomCharges"
    const Property &p = neutral_mol.property( map["charge"] );
    
    if (not p.isA<AtomCharges>())
    {
        throw SireBase::missing_property( QObject::tr(
                "The neutral molecule template is missing the AtomCharges charge property..."),
                    CODELOC );
    }

    double chg = neutral_mol.evaluate().charge(map);
    
    if (not SireMaths::isZero(chg))
        throw SireError::invalid_arg( QObject::tr(
                "The charge on the neutral molecule template is not zero! %1")
                    .arg(chg), CODELOC );

    neutral_template = neutral_mol;
    neutral_properties = properties;
    neutral_map = map;
    
    clearState();
}

/** Set the template for the positive ion state, where all properties except
    map["coordinates"] will be copied to the molecule */
void Titrator::setPositiveTemplate(const Molecule &positive_ion, const PropertyMap &map)
{
    this->setPositiveTemplate(positive_ion, QStringList(), map);
}

/** Set the template for the negative ion state, where all properties except
    map["coordinates"] will be copied to the molecule */
void Titrator::setNegativeTemplate(const Molecule &negative_ion, const PropertyMap &map)
{
    this->setNegativeTemplate(negative_ion, QStringList(), map);
}

/** Set the template for the neutral state, where all properties except
    map["coordinates"] will be copied to the molecule */
void Titrator::setNeutralTemplate(const Molecule &neutral_mol, const PropertyMap &map)
{
    this->setNeutralTemplate(neutral_mol, QStringList(), map);
}

/** Return the number of ions */
int Titrator::nIons() const
{
    int nions = 0;
    
    if (desired_chgs.isEmpty())
    {
        foreach (int charge, chgs)
        {
            if (charge != 0)
                nions += 1;
        }
    }
    else
    {
        foreach (int charge, desired_chgs)
        {
            if (charge != 0)
                nions += 1;
        }
    }
        
    return nions;
}

/** Return the number of neutral molecules */
int Titrator::nNeutrals() const
{
    int nmols = 0;

    if (desired_chgs.isEmpty())
    {
        foreach (int charge, chgs)
        {
            if (charge == 0)
                nmols += 1;
        }
    }
    else
    {
        foreach (int charge, desired_chgs)
        {
            if (charge == 0)
                nmols += 1;
        }
    }
    
    return nmols;
}

/** Return the number of positive ions */
int Titrator::nPositiveIons() const
{
    int nions = 0;
    
    if (desired_chgs.isEmpty())
    {
        foreach (int charge, chgs)
        {
            if (charge > 0)
                nions += 1;
        }
    }
    else
    {
        foreach (int charge, desired_chgs)
        {
            if (charge > 0)
                nions += 1;
        }
    }
    
    return nions;
}

/** Return the number of negative ions */
int Titrator::nNegativeIons() const
{
    int nions = 0;
    
    if (desired_chgs.isEmpty())
    {
        foreach (int charge, chgs)
        {
            if (charge < 0)
                nions += 1;
        }
    }
    else
    {
        foreach (int charge, desired_chgs)
        {
            if (charge < 0)
                nions += 1;
        }
    }
    
    return nions;
}

/** Return the index of the molecule at ion_index 'ion_index' */
int Titrator::getIonIndex(int ion_index) const
{
    ion_index = Index(ion_index).map( this->nIons() );
    
    if (desired_chgs.isEmpty())
    {
        for (int i=0; i<chgs.count(); ++i)
        {
            if (chgs[i] != 0)
            {
                if (ion_index == 0)
                    return i;
                
                ion_index -= 1;
            }
        }
    }
    else
    {
        for (int i=0; i<desired_chgs.count(); ++i)
        {
            if (desired_chgs[i] != 0)
            {
                if (ion_index == 0)
                    return i;
                
                ion_index -= 1;
            }
        }
    }

    throw SireError::program_bug( QObject::tr(
                    "It should not be possible to not find an ion!"), CODELOC );
    
    return 0;
}

/** Return the index of the molecule at positive ion_index 'ion_index' */
int Titrator::getPositiveIonIndex(int ion_index) const
{
    ion_index = Index(ion_index).map( this->nPositiveIons() );
    
    if (desired_chgs.isEmpty())
    {
        for (int i=0; i<chgs.count(); ++i)
        {
            if (chgs[i] > 0)
            {
                if (ion_index == 0)
                    return i;
                
                ion_index -= 1;
            }
        }
    }
    else
    {
        for (int i=0; i<desired_chgs.count(); ++i)
        {
            if (desired_chgs[i] > 0)
            {
                if (ion_index == 0)
                    return i;
                
                ion_index -= 1;
            }
        }
    }

    throw SireError::program_bug( QObject::tr(
                    "It should not be possible to not find a positive ion!"), CODELOC );
    
    return 0;
}

/** Return the index of the molecule at negative ion_index 'ion_index' */
int Titrator::getNegativeIonIndex(int ion_index) const
{
    ion_index = Index(ion_index).map( this->nNegativeIons() );
    
    if (desired_chgs.isEmpty())
    {
        for (int i=0; i<chgs.count(); ++i)
        {
            if (chgs[i] < 0)
            {
                if (ion_index == 0)
                    return i;
                
                ion_index -= 1;
            }
        }
    }
    else
    {
        for (int i=0; i<desired_chgs.count(); ++i)
        {
            if (desired_chgs[i] < 0)
            {
                if (ion_index == 0)
                    return i;
                
                ion_index -= 1;
            }
        }
    }

    throw SireError::program_bug( QObject::tr(
                    "It should not be possible to not find a negative ion!"), CODELOC );
    
    return 0;
}

/** Return the index of the molecule at neutral 'ion_index' */
int Titrator::getNeutralIndex(int neutral_index) const
{
    neutral_index = Index(neutral_index).map( this->nNeutrals() );
    
    if (desired_chgs.isEmpty())
    {
        for (int i=0; i<chgs.count(); ++i)
        {
            if (chgs[i] == 0)
            {
                if (neutral_index == 0)
                    return i;
                
                neutral_index -= 1;
            }
        }
    }
    else
    {
        for (int i=0; i<desired_chgs.count(); ++i)
        {
            if (desired_chgs[i] == 0)
            {
                if (neutral_index == 0)
                    return i;
                
                neutral_index -= 1;
            }
        }
    }

    throw SireError::program_bug( QObject::tr(
                    "It should not be possible to not find a neutral molecule!"), CODELOC );
    
    return 0;
}

/** Return the charge of the 'ith' molecule */
int Titrator::getCharge(int i) const
{
    i = Index(i).map(chgs.count());
    
    if (desired_chgs.isEmpty())
        return chgs[i];
    else
        return desired_chgs[i];
}

/** Swap the charge of the ith and jth molecules */
void Titrator::swapCharge(int i, int j)
{
    i = Index(i).map(chgs.count());
    j = Index(j).map(chgs.count());
    
    if (i != j)
    {
        if (desired_chgs.isEmpty())
        {
            desired_chgs = chgs;
        }
        
        int old_chg = desired_chgs[i];
        desired_chgs[i] = desired_chgs[j];
        desired_chgs[j] = old_chg;
    }
}

Molecule copyIntoMol(Molecule molecule, const Molecule &templ,
                     QStringList properties, const PropertyMap &templ_map,
                     const PropertyMap &map)
{
    if (properties.isEmpty())
    {
        properties = templ.propertyKeys();
        
        if (properties.contains("coordinates"))
            properties.removeAt( properties.indexOf("coordinates") );
    }

    if (not properties.contains("charge"))
        properties.append("charge");

    MolEditor editor = molecule.edit();
    
    foreach (QString property, properties)
    {
        editor.setProperty(map[property], templ.property(templ_map[property]));
    }

    return editor.commit();
}

/** Randomise all of the charges, ensuring there are npositive positive charges
    and nnegative negative charges */
void Titrator::randomiseCharge(int npositive, int nnegative)
{
    if (chgs.isEmpty())
        return;

    const RanGenerator &rand = RanGenerator::global();
    
    desired_chgs = QVector<qint32>(chgs.count(), 0);
    
    while (npositive + nnegative > desired_chgs.count())
    {
        npositive -= 1;
        
        if (npositive + nnegative > desired_chgs.count())
            nnegative -= 1;
    }
    
    while (npositive > 0)
    {
        int idx = rand.randInt(desired_chgs.count());
        
        if (desired_chgs[idx] == 0)
        {
            desired_chgs[idx] = pos_charge;
            npositive -= 1;
        }
    }
    
    while (nnegative > 0)
    {
        int idx = rand.randInt(desired_chgs.count());
        
        if (desired_chgs[idx] == 0)
        {
            desired_chgs[idx] = neg_charge;
            nnegative -= 1;
        }
    }
}

/** Randomise all of the charges - this ensures that there are 2*ncharges charges
    in the system (ncharges positive, and ncharges negative, with the rest neutral) */
void Titrator::randomiseCharge(int ncharges)
{
    this->randomiseCharge(ncharges, ncharges);
}

/** Apply the set of charges to the passed system */
double Titrator::applyTo(System &system)
{
    MoleculeGroup group = system.group(mgname);
    
    if (group.number() != mgnum or group.version().majorVersion() != mgversion.majorVersion())
        //if there is a change in group, or a change in the molecules in a group,
        //then we need to reinitialise this object
        clearState();
    
    PropertyName chgkey = propmap["charge"];
    
    AtomCharges positive_chgs = positive_template.property( positive_map["charge"] )
                                                 .asA<AtomCharges>();
    
    AtomCharges negative_chgs = negative_template.property( negative_map["charge"] )
                                                 .asA<AtomCharges>();
    
    AtomCharges neutral_chgs = neutral_template.property( neutral_map["charge"] )
                                               .asA<AtomCharges>();
    
    if (chgs.isEmpty())
    {
        desired_chgs.clear();
        
        //go through all of the molecules and work out what charge they have. We do this
        //by checking the "charge" property of the molecule and seeing which of the three
        //templates it agrees with. If it agrees with none of them, then we will force the
        //molecule to have a neutral charge and will assign the neutral molecule
        int nmols = group.nMolecules();
        
        chgs = QVector<qint32>(nmols, 0);
        
        Molecules mols;
        
        for (int i=0; i<nmols; ++i)
        {
            Molecule molecule = group[MolIdx(i)].molecule();
            
            AtomCharges molchgs;
            
            try
            {
                molchgs = molecule.property(chgkey).asA<AtomCharges>();
            }
            catch(...)
            {}
            
            if (molchgs.array() == positive_chgs.array())
            {
                //this is a positively charged molecule - copy across all of the
                //other positive ion properties
                molecule = copyIntoMol(molecule, positive_template, positive_properties,
                                       positive_map, propmap);
            
                chgs[i] = pos_charge;
            }
            else if (molchgs.array() == negative_chgs.array())
            {
                //this is a negatively charged molecule - copy across all of the
                //other negative ion properties
                molecule = copyIntoMol(molecule, negative_template, negative_properties,
                                       negative_map, propmap);
                
                chgs[i] = neg_charge;
            }
            else
            {
                //this is a neutral molecule - copy across all of the neutral molecule
                //properties
                molecule = copyIntoMol(molecule, neutral_template, neutral_properties,
                                       neutral_map, propmap);
                
                chgs[i] = 0;
            }
            
            mols.add(molecule);
        }
        
        group.update(mols);
    }
    else
    {
        if (not desired_chgs.isEmpty())
        {
            Molecules mols;
        
            for (int i=0; i<chgs.count(); ++i)
            {
                if (desired_chgs.at(i) != chgs.at(i))
                {
                    //we need to change the charge of the molecule
                    qint32 newchg = desired_chgs.at(i);

                    Molecule molecule = group[MolIdx(i)].molecule();
                    
                    if (newchg > 0)
                    {
                        //switch the molecule to the positive ion
                        molecule = copyIntoMol(molecule, positive_template, positive_properties,
                                               positive_map, propmap);
                    }
                    else if (newchg < 0)
                    {
                        //switch the molecule to the negative ion
                        molecule = copyIntoMol(molecule, negative_template, negative_properties,
                                               negative_map, propmap);
                    }
                    else
                    {
                        //switch the molecule to the neutral state
                        molecule = copyIntoMol(molecule, neutral_template, neutral_properties,
                                               neutral_map, propmap);
                    }
                    
                    mols.add(molecule);
                }
            }
            
            group.update(mols);
            
            chgs = desired_chgs;
            desired_chgs.clear();
        }
    }

    mgnum = group.number();
    mgversion = group.version();
    
    system.update(group);
    
    PropertyName tiprop = propmap["titrator"];
    
    if (tiprop.hasSource())
        system.setProperty(tiprop.source(), *this );
    
    return 0;
}
