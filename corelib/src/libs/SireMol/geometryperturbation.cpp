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

#include "geometryperturbation.h"

#include "molecule.h"
#include "moleditor.h"
#include "mover.hpp"
#include "core.h"

#include "SireCAS/values.h"
#include "SireCAS/identities.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireCAS;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

///////////
/////////// Implementation of GeometryPerturbation
///////////

static const RegisterMetaType<GeometryPerturbation> r_geompert( MAGIC_ONLY,
                                                     GeometryPerturbation::typeName() );

QDataStream &operator<<(QDataStream &ds,
                                       const GeometryPerturbation &geompert)
{
    writeHeader(ds, r_geompert, 1);

    ds << static_cast<const Perturbation&>(geompert);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, GeometryPerturbation &geompert)
{
    VersionID v = readHeader(ds, r_geompert);

    if (v == 1)
    {
        ds >> static_cast<Perturbation&>(geompert);
    }
    else
        throw version_error(v, "1", r_geompert, CODELOC);

    return ds;
}

/** Constructor */
GeometryPerturbation::GeometryPerturbation(const PropertyMap &map) : Perturbation(map)
{}

/** Constructor */
GeometryPerturbation::GeometryPerturbation(const Expression &mapping_function,
                                           const PropertyMap &map)
                     : Perturbation(mapping_function, map)
{}

/** Copy constructor */
GeometryPerturbation::GeometryPerturbation(const GeometryPerturbation &other)
                     : Perturbation(other)
{}

/** Destructor */
GeometryPerturbation::~GeometryPerturbation()
{}

/** Copy assignment operator */
GeometryPerturbation& GeometryPerturbation::operator=(const GeometryPerturbation &other)
{
    Perturbation::operator=(other);
    return *this;
}

/** Comparison operator */
bool GeometryPerturbation::operator==(const GeometryPerturbation &other) const
{
    return Perturbation::operator==(other);
}

/** Comparison operator */
bool GeometryPerturbation::operator!=(const GeometryPerturbation &other) const
{
    return Perturbation::operator!=(other);
}

/** Return the properties required or affected by this perturbation */
QSet<QString> GeometryPerturbation::requiredProperties() const
{
    QSet<QString> props;

    PropertyName coords_property = propertyMap()["coordinates"];

    if (coords_property.hasSource())
        props.insert( coords_property.source() );

    return props;
}

void GeometryPerturbation::perturbMolecule(MolEditor &molecule,
                                           const Values &values) const
{
    Mover<Molecule> molmover = molecule.move();

    this->perturbMolecule(molmover, values);

    molecule = molmover.commit().edit();
}

const char* GeometryPerturbation::typeName()
{
    return "SireMol::GeometryPerturbation";
}

Q_GLOBAL_STATIC( QMutex, globalMutex )
Q_GLOBAL_STATIC( SharedPolyPointer<GeometryPerturbation>, perturbationPtr );

const NullGeometryPerturbation& GeometryPerturbation::null()
{
    SharedPolyPointer<GeometryPerturbation> *ptr = perturbationPtr();

    if (ptr->constData() == 0)
    {
        QMutexLocker lkr( globalMutex() );

        if (ptr->constData() == 0)
            *ptr = static_cast<GeometryPerturbation*>(new NullGeometryPerturbation());
    }

    return ptr->constData()->asA<NullGeometryPerturbation>();
}

//////////
////////// Implementation of NullGeometryPerturbation
//////////

static const RegisterMetaType<NullGeometryPerturbation> r_nullpert;

QDataStream &operator<<(QDataStream &ds,
                                       const NullGeometryPerturbation &nullpert)
{
    writeHeader(ds, r_nullpert, 1);

    ds << static_cast<const GeometryPerturbation&>(nullpert);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                       NullGeometryPerturbation &nullpert)
{
    VersionID v = readHeader(ds, r_nullpert);

    if (v == 1)
    {
        ds >> static_cast<GeometryPerturbation&>(nullpert);
    }
    else
        throw version_error(v, "1", r_nullpert, CODELOC);

    return ds;
}

NullGeometryPerturbation::NullGeometryPerturbation()
        : ConcreteProperty<NullGeometryPerturbation,GeometryPerturbation>()
{}

NullGeometryPerturbation::NullGeometryPerturbation(const NullGeometryPerturbation &other)
        : ConcreteProperty<NullGeometryPerturbation,GeometryPerturbation>(other)
{}

NullGeometryPerturbation::~NullGeometryPerturbation()
{}

const char* NullGeometryPerturbation::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullGeometryPerturbation>() );
}

NullGeometryPerturbation&
NullGeometryPerturbation::operator=(const NullGeometryPerturbation &other)
{
    GeometryPerturbation::operator=(other);
    return *this;
}

bool NullGeometryPerturbation::operator==(const NullGeometryPerturbation &other) const
{
    return GeometryPerturbation::operator==(other);
}

bool NullGeometryPerturbation::operator!=(const NullGeometryPerturbation &other) const
{
    return GeometryPerturbation::operator!=(other);
}

QSet<Symbol> NullGeometryPerturbation::requiredSymbols() const
{
    return QSet<Symbol>();
}

QSet<QString> NullGeometryPerturbation::requiredProperties() const
{
    return QSet<QString>();
}

bool NullGeometryPerturbation::wouldChange(const Molecule&, const Values&) const
{
    return false;
}

void NullGeometryPerturbation::perturbMolecule(MolEditor&, const Values&) const
{}

void NullGeometryPerturbation::perturbMolecule(Mover<Molecule>&, const Values&) const
{}

///////////
/////////// Implementation of GeometryPerturbations
///////////

static const RegisterMetaType<GeometryPerturbations> r_geomperts;

QDataStream &operator<<(QDataStream &ds,
                                       const GeometryPerturbations &geomperts)
{
    writeHeader(ds, r_geomperts, 1);

    SharedDataStream sds(ds);

    sds << geomperts.perts << static_cast<const GeometryPerturbation&>(geomperts);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                       GeometryPerturbations &geomperts)
{
    VersionID v = readHeader(ds, r_geomperts);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> geomperts.perts >> static_cast<GeometryPerturbation&>(geomperts);
    }
    else
        throw version_error(v, "1", r_geomperts, CODELOC);

    return ds;
}

/** Constructor */
GeometryPerturbations::GeometryPerturbations()
    : ConcreteProperty<GeometryPerturbations,GeometryPerturbation>(Expression())
{}

/** Construct to hold just a single perturbation */
GeometryPerturbations::GeometryPerturbations(const GeometryPerturbation &perturbation)
    : ConcreteProperty<GeometryPerturbations,GeometryPerturbation>(Expression())
{
    if (perturbation.isA<GeometryPerturbations>())
        perts = perturbation.asA<GeometryPerturbations>().perts;
    else
        perts.append( perturbation );
}

/** Construct to hold all of the passed perturbations */
GeometryPerturbations::GeometryPerturbations(const QList<GeomPertPtr> &perturbations)
    : ConcreteProperty<GeometryPerturbations,GeometryPerturbation>(Expression())
{
    for (QList<GeomPertPtr>::const_iterator it = perturbations.constBegin();
         it != perturbations.constEnd();
         ++it)
    {
        if ( (*it)->isA<GeometryPerturbations>() )
        {
            perts += (*it)->asA<GeometryPerturbations>().perts;
        }
        else
        {
            perts += *it;
        }
    }
}

/** Copy constructor */
GeometryPerturbations::GeometryPerturbations(const GeometryPerturbations &other)
       : ConcreteProperty<GeometryPerturbations,GeometryPerturbation>(other),
         perts(other.perts)
{}

/** Destructor */
GeometryPerturbations::~GeometryPerturbations()
{}

const char* GeometryPerturbations::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GeometryPerturbations>() );
}

/** Copy assignment operator */
GeometryPerturbations&
GeometryPerturbations::operator=(const GeometryPerturbations &other)
{
    if (this != &other)
    {
        perts = other.perts;
        GeometryPerturbation::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool GeometryPerturbations::operator==(const GeometryPerturbations &other) const
{
    return perts == other.perts and GeometryPerturbation::operator==(other);
}

/** Comparison operator */
bool GeometryPerturbations::operator!=(const GeometryPerturbations &other) const
{
    return not GeometryPerturbations::operator==(other);
}

QString GeometryPerturbations::toString() const
{
    if (perts.isEmpty())
        return QObject::tr("GeometryPerturbations::null");

    QStringList lines;

    lines.append( QObject::tr("GeometryPerturbations:") );

    foreach (GeomPertPtr pert, perts)
    {
        lines.append( QString("  %1").arg(pert->toString()) );
    }

    return lines.join("\n");
}

/** Return the geometry perturbations in this collection */
QList<GeomPertPtr> GeometryPerturbations::perturbations() const
{
    return perts;
}

/** Return a re-created version of this set of perturbations where all child
    perturbations are changed to use the passed mapping function */
PerturbationPtr GeometryPerturbations::recreate(const Expression &mapping_function) const
{
    QList<GeomPertPtr> new_perts;

    for (QList<GeomPertPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        new_perts.append(
            it->read().recreate(mapping_function)->asA<GeometryPerturbation>() );
    }

    GeometryPerturbations ret(*this);
    ret.perts = new_perts;

    return ret;
}

/** Return a re-created version of this set of perturbations where all child
    perturbations are changed to use the passed property map */
PerturbationPtr GeometryPerturbations::recreate(const PropertyMap &map) const
{
    QList<GeomPertPtr> new_perts;

    for (QList<GeomPertPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        new_perts.append(
            it->read().recreate(map)->asA<GeometryPerturbation>() );
    }

    GeometryPerturbations ret(*this);
    ret.perts = new_perts;

    return ret;
}

/** Return a re-created version of this set of perturbations where all child
    perturbations are changed to use the passed mapping function and property map */
PerturbationPtr GeometryPerturbations::recreate(const Expression &mapping_function,
                                                const PropertyMap &map) const
{
    QList<GeomPertPtr> new_perts;

    for (QList<GeomPertPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        new_perts.append(
            it->read().recreate(mapping_function,map)->asA<GeometryPerturbation>() );
    }

    GeometryPerturbations ret(*this);
    ret.perts = new_perts;

    return ret;
}

/** Substitute the identities in 'identities' in all of the mapping functions
    used by this perturbation. This is useful if, for example, you want to
    switch from using 'lambda' to control the perturbation to using 'alpha', e.g.

    alpha_perturbations = lambda_perturbations.substitute( lam == Expression(alpha) );
*/
PerturbationPtr GeometryPerturbations::substitute(const Identities &identities) const
{
    QList<GeomPertPtr> new_perts;

    for (QList<GeomPertPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        new_perts.append(
            it->read().substitute(identities)->asA<GeometryPerturbation>() );
    }

    GeometryPerturbations ret(*this);
    ret.perts = new_perts;

    return ret;
}

PerturbationPtr GeometryPerturbations::substitute(const SireCAS::Symbol &old_symbol,
                                                  const SireCAS::Symbol &new_symbol) const
{
    return Perturbation::substitute(old_symbol, new_symbol);
}

/** Return the list of all child perturbations (and children of children) */
QList<PerturbationPtr> GeometryPerturbations::children() const
{
    QList<PerturbationPtr> kids;

    for (QList<GeomPertPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        kids += it->read().children();
    }

    return kids;
}

/** Return all of the symbols that need to be input to these perturbations */
QSet<Symbol> GeometryPerturbations::requiredSymbols() const
{
    QSet<Symbol> syms;

    for (QList<GeomPertPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        syms += it->read().requiredSymbols();
    }

    return syms;
}

/** Return all of the properties that are needed or affected by
    these perturbations */
QSet<QString> GeometryPerturbations::requiredProperties() const
{
    QSet<QString> props;

    for (QList<GeomPertPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        props += it->read().requiredProperties();
    }

    return props;
}

/** Return whether or not these perturbations with the passed values would
    change the molecule 'molecule' */
bool GeometryPerturbations::wouldChange(const Molecule &molecule,
                                        const Values &values) const
{
    try
    {
        for (QList<GeomPertPtr>::const_iterator it = perts.constBegin();
             it != perts.constEnd();
             ++it)
        {
            if (it->read().wouldChange(molecule,values))
                return true;
        }

        return false;
    }
    catch(...)
    {
        //if an error occured, then the molecule won't be changed
        return false;
    }
}

void GeometryPerturbations::perturbMolecule(Mover<Molecule> &molecule,
                                            const SireCAS::Values &values) const
{
    for (QList<GeomPertPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        it->read().perturbMolecule(molecule, values);
    }
}

///////////
/////////// Implementation of BondPerturbation
///////////

static const RegisterMetaType<BondPerturbation> r_bondpert;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const BondPerturbation &bondpert)
{
    writeHeader(ds, r_bondpert, 1);

    SharedDataStream sds(ds);

    sds << bondpert.bondid << bondpert.start_size.to(angstrom)
                           << bondpert.end_size.to(angstrom)
                           << static_cast<const GeometryPerturbation&>(bondpert);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, BondPerturbation &bondpert)
{
    VersionID v = readHeader(ds, r_bondpert);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        double start_size, end_size;

        sds >> bondpert.bondid >> start_size >> end_size
            >> static_cast<GeometryPerturbation&>(bondpert);

        bondpert.start_size = start_size * angstrom;
        bondpert.end_size = end_size * angstrom;
    }
    else
        throw version_error(v, "1", r_bondpert, CODELOC);

    return ds;
}

/** Constructor */
BondPerturbation::BondPerturbation()
                 : ConcreteProperty<BondPerturbation,GeometryPerturbation>(),
                   start_size(0), end_size(0)
{}

/** Construct to perturb the bond 'bond' from 'start' to 'end' */
BondPerturbation::BondPerturbation(const BondID &bond,
                                   const Length &start, const Length &end,
                                   const PropertyMap &map)
                 : ConcreteProperty<BondPerturbation,GeometryPerturbation>(map),
                   bondid(bond), start_size(start), end_size(end)
{}

/** Construct to perturb the bond 'bond' from 'start' to 'end'
    using the passed mapping function */
BondPerturbation::BondPerturbation(const BondID &bond,
                                   const Length &start, const Length &end,
                                   const Expression &mapping_function,
                                   const PropertyMap &map)
     : ConcreteProperty<BondPerturbation,GeometryPerturbation>(mapping_function, map),
       bondid(bond), start_size(start), end_size(end)
{}

/** Construct to perturb the bond between atoms 'atom0' and 'atom1'
    from 'start' to 'end' */
BondPerturbation::BondPerturbation(const AtomID &atom0, const AtomID &atom1,
                                   const Length &start, const Length &end,
                                   const PropertyMap &map)
                 : ConcreteProperty<BondPerturbation,GeometryPerturbation>(map),
                   bondid(atom0,atom1), start_size(start), end_size(end)
{}

/** Construct to perturb the bond between atoms 'atom0' and 'atom1'
    from 'start' to 'end' using the passed mapping function */
BondPerturbation::BondPerturbation(const AtomID &atom0, const AtomID &atom1,
                                   const Length &start, const Length &end,
                                   const Expression &mapping_function,
                                   const PropertyMap &map)
      : ConcreteProperty<BondPerturbation,GeometryPerturbation>(mapping_function, map),
        bondid(atom0,atom1), start_size(start), end_size(end)
{}

/** Copy constructor */
BondPerturbation::BondPerturbation(const BondPerturbation &other)
                 : ConcreteProperty<BondPerturbation,GeometryPerturbation>(other),
                   bondid(other.bondid), start_size(other.start_size),
                   end_size(other.end_size)
{}

/** Destructor */
BondPerturbation::~BondPerturbation()
{}

const char* BondPerturbation::typeName()
{
    return "SireMol::BondPerturbation";
}

/** Copy assignment operator */
BondPerturbation& BondPerturbation::operator=(const BondPerturbation &other)
{
    if (this != &other)
    {
        bondid = other.bondid;
        start_size = other.start_size;
        end_size = other.end_size;

        GeometryPerturbation::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool BondPerturbation::operator==(const BondPerturbation &other) const
{
    return bondid == other.bondid and start_size == other.start_size and
           end_size == other.end_size and GeometryPerturbation::operator==(other);
}

/** Comparison operator */
bool BondPerturbation::operator!=(const BondPerturbation &other) const
{
    return not BondPerturbation::operator==(other);
}

QString BondPerturbation::toString() const
{
    return QObject::tr("BondPerturbation( %1 from %2 A to %3 A ")
                .arg(bondid.toString())
                .arg(start_size.to(angstrom))
                .arg(end_size.to(angstrom));
}

/** Return the ID that identifies that bond that will be perturbed */
const BondID& BondPerturbation::bond() const
{
    return bondid;
}

/** Return the start length of the bond */
const SireUnits::Dimension::Length& BondPerturbation::start() const
{
    return start_size;
}

/** Return the end length of the bond */
const SireUnits::Dimension::Length& BondPerturbation::end() const
{
    return end_size;
}

/** Return whether or not this perturbation with the passed values would
    change the molecule 'molecule' */
bool BondPerturbation::wouldChange(const Molecule &molecule, const Values &values) const
{
    try
    {
        Values new_vals = values + ( symbols().initial() == start_size.value() ) +
                                   ( symbols().final() == end_size.value() );

        Length new_length = Length( mappingFunction().evaluate(new_vals) );

        Length old_length( bondid.length(molecule, propertyMap()) );

        return std::abs(new_length - old_length) > 0.000001;
    }
    catch(...)
    {
        return false;
    }
}

/** Apply this perturbation

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
void BondPerturbation::perturbMolecule(Mover<Molecule> &molecule,
                                       const Values &values) const
{
    //calculate the desired value of the bond
    Values new_vals = values + ( symbols().initial() == start_size.value() ) +
                               ( symbols().final() == end_size.value() );

    Length new_length = Length( mappingFunction().evaluate(new_vals) );

    Length old_length( bondid.length(molecule, propertyMap()) );

    if (std::abs(new_length - old_length) > 0.000001)
        molecule.set(bondid, new_length, propertyMap()).commit();
}

///////////
/////////// Implementation of AnglePerturbation
///////////

static const RegisterMetaType<AnglePerturbation> r_anglepert;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const AnglePerturbation &anglepert)
{
    writeHeader(ds, r_anglepert, 1);

    SharedDataStream sds(ds);

    sds << anglepert.angleid << anglepert.start_size.to(degrees)
                             << anglepert.end_size.to(degrees)
                             << static_cast<const GeometryPerturbation&>(anglepert);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AnglePerturbation &anglepert)
{
    VersionID v = readHeader(ds, r_anglepert);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        double start_size, end_size;

        sds >> anglepert.angleid >> start_size >> end_size
            >> static_cast<GeometryPerturbation&>(anglepert);

        anglepert.start_size = start_size * degrees;
        anglepert.end_size = end_size * degrees;
    }
    else
        throw version_error(v, "1", r_anglepert, CODELOC);

    return ds;
}

/** Constructor */
AnglePerturbation::AnglePerturbation()
                  : ConcreteProperty<AnglePerturbation,GeometryPerturbation>(),
                    start_size(0), end_size(0)
{}

/** Construct to perturb the angle 'angle' from 'start' to 'end' */
AnglePerturbation::AnglePerturbation(const AngleID &angle,
                                     const Angle &start, const Angle &end,
                                     const PropertyMap &map)
                  : ConcreteProperty<AnglePerturbation,GeometryPerturbation>(map),
                    angleid(angle), start_size(start), end_size(end)
{}

/** Construct to perturb the angle 'angle' from 'start' to 'end'
    using the passed mapping function */
AnglePerturbation::AnglePerturbation(const AngleID &angle,
                                     const Angle &start, const Angle &end,
                                     const Expression &mapping_function,
                                     const PropertyMap &map)
     : ConcreteProperty<AnglePerturbation,GeometryPerturbation>(mapping_function, map),
       angleid(angle), start_size(start), end_size(end)
{}

/** Construct to perturb the angle between atoms 'atom0', 'atom1' and 'atom2'
    from 'start' to 'end' */
AnglePerturbation::AnglePerturbation(const AtomID &atom0, const AtomID &atom1,
                                     const AtomID &atom2,
                                     const Angle &start, const Angle &end,
                                     const PropertyMap &map)
                 : ConcreteProperty<AnglePerturbation,GeometryPerturbation>(map),
                   angleid(atom0,atom1,atom2), start_size(start), end_size(end)
{}

/** Construct to perturb the angle between atoms 'atom0', 'atom1' and 'atom2'
    from 'start' to 'end' using the passed mapping function */
AnglePerturbation::AnglePerturbation(const AtomID &atom0, const AtomID &atom1,
                                     const AtomID &atom2,
                                     const Angle &start, const Angle &end,
                                     const Expression &mapping_function,
                                     const PropertyMap &map)
      : ConcreteProperty<AnglePerturbation,GeometryPerturbation>(mapping_function, map),
        angleid(atom0,atom1,atom2), start_size(start), end_size(end)
{}

/** Copy constructor */
AnglePerturbation::AnglePerturbation(const AnglePerturbation &other)
                 : ConcreteProperty<AnglePerturbation,GeometryPerturbation>(other),
                   angleid(other.angleid), start_size(other.start_size),
                   end_size(other.end_size)
{}

/** Destructor */
AnglePerturbation::~AnglePerturbation()
{}

const char* AnglePerturbation::typeName()
{
    return "SireMol::AnglePerturbation";
}

QString AnglePerturbation::toString() const
{
    return QObject::tr("AnglePerturbation( %1 from %2° to %3° ")
                .arg(angleid.toString())
                .arg(start_size.to(degrees))
                .arg(end_size.to(degrees));
}

/** Copy assignment operator */
AnglePerturbation& AnglePerturbation::operator=(const AnglePerturbation &other)
{
    if (this != &other)
    {
        angleid = other.angleid;
        start_size = other.start_size;
        end_size = other.end_size;

        GeometryPerturbation::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool AnglePerturbation::operator==(const AnglePerturbation &other) const
{
    return angleid == other.angleid and start_size == other.start_size and
           end_size == other.end_size and GeometryPerturbation::operator==(other);
}

/** Comparison operator */
bool AnglePerturbation::operator!=(const AnglePerturbation &other) const
{
    return not AnglePerturbation::operator==(other);
}

/** Return the ID that identifies that angle that will be perturbed */
const AngleID& AnglePerturbation::angle() const
{
    return angleid;
}

/** Return the start length of the angle */
const SireUnits::Dimension::Angle& AnglePerturbation::start() const
{
    return start_size;
}

/** Return the end length of the angle */
const SireUnits::Dimension::Angle& AnglePerturbation::end() const
{
    return end_size;
}

/** Return whether or not this perturbation with the passed values would
    change the molecule 'molecule' */
bool AnglePerturbation::wouldChange(const Molecule &molecule, const Values &values) const
{
    try
    {
        Values new_vals = values + ( symbols().initial() == start_size.value() ) +
                                   ( symbols().final() == end_size.value() );

        Angle new_size = Angle( mappingFunction().evaluate(new_vals) );

        Angle old_size( angleid.size(molecule, propertyMap()) );

        return std::abs(new_size.value() - old_size.value()) > 0.0001;
    }
    catch(...)
    {
        return false;
    }
}

/** Apply this perturbation

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
void AnglePerturbation::perturbMolecule(Mover<Molecule> &molecule,
                                       const Values &values) const
{
    //calculate the desired value of the angle
    Values new_vals = values + ( symbols().initial() == start_size.value() ) +
                               ( symbols().final() == end_size.value() );

    Angle old_size( angleid.size(molecule, propertyMap()) );
    Angle new_angle = Angle( mappingFunction().evaluate(new_vals) );

    if (std::abs(new_angle.value() - old_size.value()) > 0.0001)
        molecule.set(angleid, new_angle, propertyMap()).commit();
}

///////////
/////////// Implementation of DihedralPerturbation
///////////

static const RegisterMetaType<DihedralPerturbation> r_dihedralpert;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const DihedralPerturbation &dihedralpert)
{
    writeHeader(ds, r_dihedralpert, 1);

    SharedDataStream sds(ds);

    sds << dihedralpert.dihedralid << dihedralpert.start_size.to(degrees)
                        << dihedralpert.end_size.to(degrees)
                        << static_cast<const GeometryPerturbation&>(dihedralpert);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, DihedralPerturbation &dihedralpert)
{
    VersionID v = readHeader(ds, r_dihedralpert);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        double start_size, end_size;

        sds >> dihedralpert.dihedralid >> start_size >> end_size
            >> static_cast<GeometryPerturbation&>(dihedralpert);

        dihedralpert.start_size = start_size * degrees;
        dihedralpert.end_size = end_size * degrees;
    }
    else
        throw version_error(v, "1", r_dihedralpert, CODELOC);

    return ds;
}

/** Constructor */
DihedralPerturbation::DihedralPerturbation()
                     : ConcreteProperty<DihedralPerturbation,GeometryPerturbation>(),
                       start_size(0), end_size(0)
{}

/** Construct to perturb the dihedral 'dihedral' from 'start' to 'end' */
DihedralPerturbation::DihedralPerturbation(const DihedralID &dihedral,
                                     const Angle &start, const Angle &end,
                                     const PropertyMap &map)
                     : ConcreteProperty<DihedralPerturbation,GeometryPerturbation>(map),
                       dihedralid(dihedral), start_size(start), end_size(end)
{}

/** Construct to perturb the dihedral 'dihedral' from 'start' to 'end'
    using the passed mapping function */
DihedralPerturbation::DihedralPerturbation(const DihedralID &dihedral,
                                     const Angle &start, const Angle &end,
                                     const Expression &mapping_function,
                                     const PropertyMap &map)
     : ConcreteProperty<DihedralPerturbation,GeometryPerturbation>(mapping_function, map),
       dihedralid(dihedral), start_size(start), end_size(end)
{}

/** Construct to perturb the dihedral between atoms 'atom0', 'atom1', 'atom2' and 'atom3'
    from 'start' to 'end' */
DihedralPerturbation::DihedralPerturbation(const AtomID &atom0, const AtomID &atom1,
                                     const AtomID &atom2, const AtomID &atom3,
                                     const Angle &start, const Angle &end,
                                     const PropertyMap &map)
                 : ConcreteProperty<DihedralPerturbation,GeometryPerturbation>(map),
                   dihedralid(atom0,atom1,atom2,atom3), start_size(start), end_size(end)
{}

/** Construct to perturb the dihedral between atoms 'atom0', 'atom1', 'atom2' and 'atom3'
    from 'start' to 'end' using the passed mapping function */
DihedralPerturbation::DihedralPerturbation(const AtomID &atom0, const AtomID &atom1,
                                     const AtomID &atom2, const AtomID &atom3,
                                     const Angle &start, const Angle &end,
                                     const Expression &mapping_function,
                                     const PropertyMap &map)
    : ConcreteProperty<DihedralPerturbation,GeometryPerturbation>(mapping_function, map),
      dihedralid(atom0,atom1,atom2,atom3), start_size(start), end_size(end)
{}

/** Copy constructor */
DihedralPerturbation::DihedralPerturbation(const DihedralPerturbation &other)
                 : ConcreteProperty<DihedralPerturbation,GeometryPerturbation>(other),
                   dihedralid(other.dihedralid), start_size(other.start_size),
                   end_size(other.end_size)
{}

/** Destructor */
DihedralPerturbation::~DihedralPerturbation()
{}

const char* DihedralPerturbation::typeName()
{
    return "SireMol::DihedralPerturbation";
}

QString DihedralPerturbation::toString() const
{
    return QObject::tr("DihedralPerturbation( %1 from %2° to %3° ")
                .arg(dihedralid.toString())
                .arg(start_size.to(degrees))
                .arg(end_size.to(degrees));
}

/** Copy assignment operator */
DihedralPerturbation& DihedralPerturbation::operator=(const DihedralPerturbation &other)
{
    if (this != &other)
    {
        dihedralid = other.dihedralid;
        start_size = other.start_size;
        end_size = other.end_size;

        GeometryPerturbation::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool DihedralPerturbation::operator==(const DihedralPerturbation &other) const
{
    return dihedralid == other.dihedralid and start_size == other.start_size and
           end_size == other.end_size and GeometryPerturbation::operator==(other);
}

/** Comparison operator */
bool DihedralPerturbation::operator!=(const DihedralPerturbation &other) const
{
    return not DihedralPerturbation::operator==(other);
}

/** Return the ID that identifies that dihedral that will be perturbed */
const DihedralID& DihedralPerturbation::dihedral() const
{
    return dihedralid;
}

/** Return the start length of the dihedral */
const SireUnits::Dimension::Angle& DihedralPerturbation::start() const
{
    return start_size;
}

/** Return the end length of the dihedral */
const SireUnits::Dimension::Angle& DihedralPerturbation::end() const
{
    return end_size;
}

/** Return whether or not this perturbation with the passed values would
    change the molecule 'molecule' */
bool DihedralPerturbation::wouldChange(const Molecule &molecule,
                                       const Values &values) const
{
    try
    {
        Values new_vals = values + ( symbols().initial() == start_size.value() ) +
                                   ( symbols().final() == end_size.value() );

        Angle new_size = Angle( mappingFunction().evaluate(new_vals) );

        Angle old_size( dihedralid.size(molecule, propertyMap()) );

        return std::abs(new_size.value() - old_size.value()) > 0.0001;
    }
    catch(...)
    {
        return false;
    }
}

/** Apply this perturbation

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
void DihedralPerturbation::perturbMolecule(Mover<Molecule> &molecule,
                                       const Values &values) const
{
    //calculate the desired value of the dihedral
    Values new_vals = values + ( symbols().initial() == start_size.value() ) +
                               ( symbols().final() == end_size.value() );

    Angle new_dihedral = Angle( mappingFunction().evaluate(new_vals) );

    molecule.set(dihedralid, new_dihedral, propertyMap()).commit();
}
