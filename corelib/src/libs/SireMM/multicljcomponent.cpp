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

#include "multicljcomponent.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include <QRegExp>

using namespace SireMM;
using namespace SireFF;
using namespace SireCAS;
using namespace SireStream;

namespace SireMM
{
    boost::tuple<QString,QString> getSubscriptedProperty(QString name)
    {
        QRegExp key_regexp("\\[(.*)\\]$");
 
        int idx = key_regexp.indexIn(name);
        if (idx != -1)
        {
            QString cljkey = key_regexp.cap(1).trimmed();
            QString cljprop = name;
            cljprop.truncate(idx);
            cljprop = cljprop.trimmed();
            
            return boost::tuple<QString,QString>(cljprop,cljkey);
        }
        else
            return boost::tuple<QString,QString>(name, QString::null);
    }

}

//////////
////////// Implementation of MultiCLJEnergy
//////////

QDataStream &operator<<(QDataStream &ds, const MultiCLJEnergy &nrg)
{
    quint32 version = 1;
    
    ds << version;
    
    SharedDataStream sds(ds);
    
    sds << nrg.cnrgs << nrg.ljnrgs << nrg.coulomb() << nrg.lj();
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, MultiCLJEnergy &nrg)
{
    quint32 version;
    
    ds >> version;
    
    if (version == 1)
    {
        SharedDataStream sds(ds);
        
        double cnrg, ljnrg;
        QVector<double> cnrgs, ljnrgs;
        
        sds >> cnrgs >> ljnrgs >> cnrg >> ljnrg;
        
        if (cnrgs.isEmpty())
        {
            nrg = MultiCLJEnergy(cnrg, ljnrg);
        }
        else
        {
            nrg = MultiCLJEnergy(cnrgs, ljnrgs);
        }
    }
    else
        throw version_error( QObject::tr(
                "Unsupported version of MultiCLJEnergy being loaded (%1). "
                "Only version 1 is supported.").arg(version), CODELOC );
    
    return ds;
}

void MultiCLJEnergy::assertValidCoulombIndex(quint32 i) const
{
    if (i >= cnrgs.count())
        throw SireError::invalid_index( QObject::tr(
            "There is no CoulombEnergy at index '%1'. Number of energies equals %2.")
                .arg(i).arg(cnrgs.count()), CODELOC );
}

void MultiCLJEnergy::assertValidLJIndex(quint32 i) const
{
    if (i >= ljnrgs.count())
        throw SireError::invalid_index( QObject::tr(
            "There is no LJEnergy at index '%1'. Number of energies equals %2.")
                .arg(i).arg(ljnrgs.count()), CODELOC );
}

//////////
////////// Implementation of MultiCLJComponent
//////////

static const RegisterMetaType<MultiCLJComponent> r_cljcomp;

QDataStream &operator<<(QDataStream &ds, const MultiCLJComponent &cljcomp)
{
    writeHeader(ds, r_cljcomp, 1);
    
    SharedDataStream sds(ds);
    
    sds << cljcomp.comps << cljcomp.key_to_idx
        << static_cast<const FFComponent&>(cljcomp);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, MultiCLJComponent &cljcomp)
{
    VersionID v = readHeader(ds, r_cljcomp);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> cljcomp.comps >> cljcomp.key_to_idx
            >> static_cast<FFComponent&>(cljcomp);
    }
    else
        throw version_error(v, "1", r_cljcomp, CODELOC);
    
    return ds;
}

/** Construct with just a single, default, CLJComponent for the 
    forcefield with the passed name */
MultiCLJComponent::MultiCLJComponent(const FFName &name)
                  : FFComponent(name, QLatin1String("CLJ"))
{
    comps.append( CLJComponent(name) );
    key_to_idx.insert("default", 0);
}

/** Copy constructor */
MultiCLJComponent::MultiCLJComponent(const MultiCLJComponent &other)
                  : FFComponent(other), comps(other.comps), key_to_idx(other.key_to_idx)
{}

/** Destructor */
MultiCLJComponent::~MultiCLJComponent()
{}

/** Return a copy of this MultiCLJComponent that has been renamed for
    the passed forcefield */
MultiCLJComponent MultiCLJComponent::rename(const FFName &name) const
{
    MultiCLJComponent ret(name);
    
    ret.comps = QVector<CLJComponent>(comps.count());
    ret.key_to_idx = key_to_idx;
    
    for (QHash<QString,quint32>::const_iterator it = key_to_idx.constBegin();
         it != key_to_idx.constEnd();
         ++it)
    {
        ret.comps[it.value()] = CLJComponent(name, it.key());
    }
    
    return ret;
}

/** Return a string representation of these components */
QString MultiCLJComponent::toString() const
{
    return QObject::tr("MultiCLJComponent( %1 )").arg( Sire::toString(symbols()) );
}

/** Copy assignment operator */
MultiCLJComponent& MultiCLJComponent::operator=(const MultiCLJComponent &other)
{
    if (this != &other)
    {
        FFComponent::operator=(other);
        comps = other.comps;
        key_to_idx = other.key_to_idx;
    }
    
    return *this;
}

/** Comparison operator */
bool MultiCLJComponent::operator==(const MultiCLJComponent &other) const
{
    return FFComponent::operator==(other) and comps == other.comps
              and key_to_idx == other.key_to_idx;
}

/** Comparison operator */
bool MultiCLJComponent::operator!=(const MultiCLJComponent &other) const
{
    return not operator==(other);
}

/** Return the default coulomb component */
const CoulombComponent& MultiCLJComponent::coulomb() const
{
    return comps.at(0).coulomb();
}

/** Return the coulomb component for the key 'key' */
const CoulombComponent& MultiCLJComponent::coulomb(QString key) const
{
    assertValidKey(key);
    return comps.at( key_to_idx.value(key) ).coulomb();
}

/** Return the default LJ component */
const LJComponent& MultiCLJComponent::lj() const
{
    return comps.at(0).lj();
}

/** Return the LJ component for the key 'key' */
const LJComponent& MultiCLJComponent::lj(QString key) const
{
    assertValidKey(key);
    return comps.at( key_to_idx.value(key) ).lj();
}

/** Return the default total component */
const CLJComponent& MultiCLJComponent::total() const
{
    return comps.at(0).total();
}

/** Return the total component for the key 'key' */
const CLJComponent& MultiCLJComponent::total(QString key) const
{
    assertValidKey(key);
    return comps.at( key_to_idx.value(key) ).total();
}

const char* MultiCLJComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MultiCLJComponent>() );
}

const char* MultiCLJComponent::what() const
{
    return MultiCLJComponent::typeName();
}

MultiCLJComponent* MultiCLJComponent::clone() const
{
    return new MultiCLJComponent(*this);
}

void MultiCLJComponent::assertValidKey(QString key) const
{
    if (not key_to_idx.contains(key))
        throw SireError::invalid_key( QObject::tr(
                "There is no CLJ component with key '%1'. Available keys are %2.")
                    .arg(key).arg(Sire::toString(keys())), CODELOC );
}

/** Set the energy in the forcefield from the passed 'cljnrg' object */
void MultiCLJComponent::setEnergy(FF &ff, const MultiCLJEnergy &value) const
{
    if (comps.count() == 1)
    {
        FFComponent::setEnergy(ff, total(), value.total());
        FFComponent::setEnergy(ff, coulomb(), value.coulomb());
        FFComponent::setEnergy(ff, lj(), value.lj());
    }
    else
    {
        for (int i=0; i<comps.count(); ++i)
        {
            const CLJComponent &comp = comps.constData()[i];
        
            FFComponent::setEnergy(ff, comp.total(), value.total(i));
            FFComponent::setEnergy(ff, comp.coulomb(), value.coulomb(i));
            FFComponent::setEnergy(ff, comp.lj(), value.lj(i));
        }
    }
}

/** Change the energy in the forcefield from the passed 'cljnrg' object */
void MultiCLJComponent::changeEnergy(FF &ff, const MultiCLJEnergy &delta) const
{
    if (comps.count() == 1)
    {
        FFComponent::changeEnergy(ff, total(), delta.total());
        FFComponent::changeEnergy(ff, coulomb(), delta.coulomb());
        FFComponent::changeEnergy(ff, lj(), delta.lj());
    }
    else
    {
        for (int i=0; i<comps.count(); ++i)
        {
            const CLJComponent &comp = comps.constData()[i];
        
            FFComponent::changeEnergy(ff, comp.total(), delta.total(i));
            FFComponent::changeEnergy(ff, comp.coulomb(), delta.coulomb(i));
            FFComponent::changeEnergy(ff, comp.lj(), delta.lj(i));
        }
    }
}

/** Return all of the symbols used by these components */
SireCAS::Symbols MultiCLJComponent::symbols() const
{
    Symbols symbs;
    
    for (int i=0; i<comps.count(); ++i)
    {
        symbs.add( comps.at(i).symbols() );
    }
    
    return symbs;
}

/** Add new CLJ components with the passed key. This returns the index
    of the new component */
int MultiCLJComponent::add(QString key)
{
    if (not key_to_idx.contains(key))
    {
        comps.append( CLJComponent(this->forceFieldName(),key) );
        key_to_idx.insert(key, comps.count()-1);
        return comps.count()-1;
    }
    else
        return key_to_idx.value(key);
}

/** Remove the component with key 'key'. This returns the index of the component
    that was removed, or -1 if there was no such component */
int MultiCLJComponent::remove(QString key)
{
    if (key == "default")
        return -1;

    else if (key_to_idx.contains(key))
    {
        quint32 idx = key_to_idx.value(key);
        comps.removeAt(idx);
        key_to_idx.remove(key);
        
        QMutableHashIterator<QString, quint32> it(key_to_idx);
        
        while (it.hasNext())
        {
            it.next();
            
            if (it.value() > idx)
            {
                it.value() -= 1;
            }
        }
        
        return idx;
    }
    else
        return -1;
}

/** Remove all keys from the index (except 'default') */
void MultiCLJComponent::removeAll()
{
    while (comps.count() > 1)
    {
        comps.removeLast();
    }
    
    key_to_idx.clear();
    key_to_idx.insert( "default", 0 );
}

/** Return the index of the component with key 'key' */
int MultiCLJComponent::indexOf(QString key) const
{
    if (not key_to_idx.contains(key))
    {
        throw SireError::invalid_key( QObject::tr(
                "There is no CLJComponent with key '%1'. Available keys are %2.")
                    .arg(key).arg(Sire::toString(key_to_idx.keys())), CODELOC );
    }
    
    return key_to_idx.value(key);
}

/** Return whether or not there is a function with the associated key */
bool MultiCLJComponent::hasKey(QString key) const
{
    return key_to_idx.contains(key);
}

/** Return all of the keys for the different CLJComponents */
QStringList MultiCLJComponent::keys() const
{
    return key_to_idx.keys();
}

/** Return the number of CLJComponents / number of keys */
int MultiCLJComponent::count() const
{
    return comps.count();
}

/** Return the number of CLJComponents / number of keys */
int MultiCLJComponent::size() const
{
    return count();
}

/** Return the number of CLJComponents / number of keys */
int MultiCLJComponent::nKeys() const
{
    return count();
}
