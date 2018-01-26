/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Christopher Woods
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

#include "mmdetail.h"

#include "SireBase/propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireFF;
using namespace SireBase;
using namespace SireStream;
using namespace SireMM;

static const RegisterMetaType<MMDetail> r_mm;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const MMDetail &mm)
{
    writeHeader(ds, r_mm, 1);
    
    ds << static_cast<const FFDetail&>(mm);
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, MMDetail &mm)
{
    VersionID v = readHeader(ds, r_mm);
    
    if (v == 1)
    {
        ds >> static_cast<FFDetail&>(mm);
    }
    else
        throw version_error(v, "1", r_mm, CODELOC);
    
    return ds;
}

/** Null constructor */
MMDetail::MMDetail() : ConcreteProperty<MMDetail,FFDetail>()
{}

/** Construct passing in all of the required values */
MMDetail::MMDetail(QString name, QString combining_rules,
                  double scale14elec, double scale14vdw,
                  QString elecstyle, QString vdwstyle,
                  QString bondstyle, QString anglestyle,
                  QString dihedralstyle)
         : ConcreteProperty<MMDetail,FFDetail>(name)
{
    //convert all of the passed data into a canonical form
    if (combining_rules.startsWith("arith", Qt::CaseInsensitive))
    {
        setProperty("combining_rules", wrap("arithmetic"));
    }
    else if (combining_rules.startsWith("geo", Qt::CaseInsensitive))
    {
        setProperty("combining_rules", wrap("geometric"));
    }
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot understand the required combining rules from '%1'. These should "
            "be either 'arithmetic' or 'geometric'.").arg(combining_rules), CODELOC );
    
    setProperty("scale14elec", wrap(scale14elec));
    setProperty("scale14vdw", wrap(scale14vdw));
    
    if (elecstyle.startsWith("coul", Qt::CaseInsensitive))
    {
        setProperty("elecstyle", wrap("coulomb"));
    }
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot understand the required electrostatic model from '%1'. This should "
            "be 'coulomb'").arg(elecstyle), CODELOC );
    
    if (vdwstyle.startsWith("lj", Qt::CaseInsensitive))
    {
        setProperty("vdwstyle", wrap("lj"));
    }
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot understand the required vdw model from '%1'. This should "
            "be 'lj'").arg(vdwstyle), CODELOC );
    
    if (bondstyle.startsWith("harm", Qt::CaseInsensitive))
    {
        setProperty("bondstyle", wrap("harmonic"));
    }
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot understand the required bond model from '%1'. This should "
            "be 'harmonic'").arg(bondstyle), CODELOC );
    
    if (anglestyle.startsWith("harm", Qt::CaseInsensitive))
    {
        setProperty("anglestyle", wrap("harmonic"));
    }
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot understand the required angle model from '%1'. This should "
            "be 'harmonic'").arg(anglestyle), CODELOC );
    
    if (dihedralstyle.startsWith("cos", Qt::CaseInsensitive))
    {
        setProperty("dihedralstyle", wrap("cosine"));
    }
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot understand the required dihedral model from '%1'. This should "
            "be 'cosine'").arg(dihedralstyle), CODELOC );
    
    this->operator=( FFDetail::registerForceField(*this).read().asA<MMDetail>() );
}

/** Copy constructor */
MMDetail::MMDetail(const MMDetail &other) : ConcreteProperty<MMDetail,FFDetail>(other)
{}

/** Destructor */
MMDetail::~MMDetail()
{}

/** Copy assignment operator */
MMDetail& MMDetail::operator=(const MMDetail &other)
{
    FFDetail::operator=(other);
    return *this;
}

/** Comparison operator */
bool MMDetail::operator==(const MMDetail &other) const
{
    return FFDetail::operator==(other);
}

/** Comparison operator */
bool MMDetail::operator!=(const MMDetail &other) const
{
    return FFDetail::operator!=(other);
}

MMDetail* MMDetail::clone() const
{
    return new MMDetail(*this);
}

const char* MMDetail::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MMDetail>() );
}

const char* MMDetail::what() const
{
    return MMDetail::typeName();
}

/** Return a string form of this forcefield */
QString MMDetail::toString() const
{
    if (this->isNull())
        return QObject::tr("MMDetail::null");

    return QObject::tr("MM ForceField{ %1,\n"
                       "               combining_rules = %2,\n"
                       "               1-4 scaling = %3, %4,\n"
                       "               nonbonded = %5, %6,\n"
                       "               bond = %7, angle = %8,\n"
                       "               dihedral = %9 }")
                .arg(name()).arg(combiningRules())
                .arg(electrostatic14ScaleFactor())
                .arg(vdw14ScaleFactor())
                .arg(electrostaticStyle())
                .arg(vdwStyle())
                .arg(bondStyle())
                .arg(angleStyle())
                .arg(dihedralStyle());
}

/** Return the combining rules for this forcefield */
QString MMDetail::combiningRules() const
{
    return property("combining_rules").asAString();
}

/** Return whether or not this forcefield uses arithmetic combining rules */
bool MMDetail::usesArithmeticCombiningRules() const
{
    return combiningRules() == "arithmetic";
}

/** Return whether or not this forcefield uses arithmetic combining rules */
bool MMDetail::usesGeometricCombiningRules() const
{
    return combiningRules() == "geometric";
}

/** Return the electrostatic 1-4 scale factor */
double MMDetail::electrostatic14ScaleFactor() const
{
    return property("scale14elec").asADouble();
}

/** Return the vdw 1-4 scale factor */
double MMDetail::vdw14ScaleFactor() const
{
    return property("scale14vdw").asADouble();
}

/** Return the electrostatic model used by this forcefield */
QString MMDetail::electrostaticStyle() const
{
    return property("elecstyle").asAString();
}

/** Return whether or not this forcefield uses coulomb charges */
bool MMDetail::usesCoulombCharges() const
{
    return electrostaticStyle() == "coulomb";
}

/** Return the vdw model used by this forcefield */
QString MMDetail::vdwStyle() const
{
    return property("vdwstyle").asAString();
}

/** Return whether or not this forcefield uses the Lennard Jones vdw model */
bool MMDetail::usesLJTerm() const
{
    return vdwStyle() == "lj";
}

/** Return the bond model used by this forcefield */
QString MMDetail::bondStyle() const
{
    return property("bondstyle").asAString();
}

/** Return whether or not this forcefield uses harmonic bonds */
bool MMDetail::usesHarmonicBonds() const
{
    return bondStyle() == "harmonic";
}

/** Return the angle model used by this forcefield */
QString MMDetail::angleStyle() const
{
    return property("anglestyle").asAString();
}

/** Return whether or not this forcefield uses harmonic angles */
bool MMDetail::usesHarmonicAngles() const
{
    return angleStyle() == "harmonic";
}

/** Return the dihedral model uses by this forcefield */
QString MMDetail::dihedralStyle() const
{
    return property("dihedralstyle").asAString();
}

/** Return whether or not this forcefield uses cosine-series dihedrals */
bool MMDetail::usesCosineDihedrals() const
{
    return dihedralStyle() == "cosine";
}
