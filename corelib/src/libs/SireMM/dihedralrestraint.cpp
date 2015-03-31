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

#include "dihedralrestraint.h"

#include "SireFF/forcetable.h"

#include "SireCAS/symbols.h"
#include "SireCAS/values.h"
#include "SireCAS/conditional.h"
#include "SireCAS/power.h"

#include "SireID/index.h"

#include "SireUnits/units.h"
#include "SireUnits/angle.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireCAS/errors.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireFF;
using namespace SireID;
using namespace SireBase;
using namespace SireCAS;
using namespace SireMaths;
using namespace SireStream;
using namespace SireUnits;
using namespace SireUnits::Dimension;

////////////
//////////// Implementation of DihedralRestraint
////////////

static const RegisterMetaType<DihedralRestraint> r_dihrest;

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, 
                                      const DihedralRestraint &dihrest)
{
    writeHeader(ds, r_dihrest, 1);
    
    SharedDataStream sds(ds);
    
    sds << dihrest.p[0] << dihrest.p[1] << dihrest.p[2] << dihrest.p[3]
        << dihrest.force_expression
        << static_cast<const ExpressionRestraint3D&>(dihrest);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, DihedralRestraint &dihrest)
{
    VersionID v = readHeader(ds, r_dihrest);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> dihrest.p[0] >> dihrest.p[1] >> dihrest.p[2] >> dihrest.p[3]
            >> dihrest.force_expression
            >> static_cast<ExpressionRestraint3D&>(dihrest);

        dihrest.intra_molecule_points = Point::areIntraMoleculePoints(dihrest.p[0],
                                                                   dihrest.p[1]) and
                                        Point::areIntraMoleculePoints(dihrest.p[0],
                                                                   dihrest.p[2]) and
                                        Point::areIntraMoleculePoints(dihrest.p[0],
                                                                   dihrest.p[3]);
    }
    else
        throw version_error( v, "1", r_dihrest, CODELOC );
        
    return ds;
}

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getPhiSymbol, ("phi") );

/** Return the symbol that represents the dihedral angle between the points (phi) */
const Symbol& DihedralRestraint::phi()
{
    return *(getPhiSymbol());
}

/** Constructor */
DihedralRestraint::DihedralRestraint()
                  : ConcreteProperty<DihedralRestraint,ExpressionRestraint3D>()
{}

void DihedralRestraint::calculatePhi()
{
    if (this->restraintFunction().isFunction(phi()))
    {
        SireUnits::Dimension::Angle angle;
        
        if (intra_molecule_points)
            //we don't use the space when calculating intra-molecular angles
            angle = Vector::dihedral( p[0].read().point(), p[1].read().point(),
                                      p[2].read().point(), p[3].read().point() );
        else
            angle = this->space().calcDihedral( p[0].read().point(),
                                                p[1].read().point(),
                                                p[2].read().point(),
                                                p[3].read().point() );
                                                  
        ExpressionRestraint3D::_pvt_setValue( phi(), angle );
    }
}

/** Construct a restraint that acts on the angle within the 
    three points 'point0', 'point1' and 'point2' (theta == a(012)),
    restraining the angle within these points using the expression 
    'restraint' */
DihedralRestraint::DihedralRestraint(const PointRef &point0, const PointRef &point1,
                                     const PointRef &point2, const PointRef &point3,
                                     const Expression &restraint)
                  : ConcreteProperty<DihedralRestraint,ExpressionRestraint3D>(restraint)
{
    p[0] = point0;
    p[1] = point1;
    p[2] = point2;
    p[3] = point3;

    force_expression = this->restraintFunction().differentiate(phi());
        
    if (force_expression.isConstant())
        force_expression = force_expression.evaluate(Values());
    
    intra_molecule_points = Point::areIntraMoleculePoints(p[0], p[1]) and
                            Point::areIntraMoleculePoints(p[0], p[2]) and
                            Point::areIntraMoleculePoints(p[0], p[3]);

    this->calculatePhi();
}

/** Construct a restraint that acts on the angle within the 
    three points 'point0', 'point1' and 'point2' (theta == a(012)),
    restraining the angle within these points using the expression 
    'restraint' */
DihedralRestraint::DihedralRestraint(const PointRef &point0, const PointRef &point1,
                                     const PointRef &point2, const PointRef &point3,
                                     const Expression &restraint, const Values &values)
     : ConcreteProperty<DihedralRestraint,ExpressionRestraint3D>(restraint, values)
{
    p[0] = point0;
    p[1] = point1;
    p[2] = point2;
    p[3] = point3;

    force_expression = this->restraintFunction().differentiate(phi());
        
    if (force_expression.isConstant())
        force_expression = force_expression.evaluate(Values());
    
    intra_molecule_points = Point::areIntraMoleculePoints(p[0], p[1]) and
                            Point::areIntraMoleculePoints(p[0], p[2]) and
                            Point::areIntraMoleculePoints(p[0], p[3]);

    this->calculatePhi();
}

/** Internal constructor used to construct a restraint using the specified
    points, energy expression and force expression */
DihedralRestraint::DihedralRestraint(const PointRef &point0, const PointRef &point1,
                                     const PointRef &point2, const PointRef &point3,
                                     const Expression &nrg_restraint,
                                     const Expression &force_restraint)
           : ConcreteProperty<DihedralRestraint,ExpressionRestraint3D>(nrg_restraint),
             force_expression(force_restraint)
{
    p[0] = point0;
    p[1] = point1;
    p[2] = point2;
    p[3] = point3;

    if (force_expression.isConstant())
    {
        force_expression = force_expression.evaluate(Values());
    }
    else
    {
        if (not this->restraintFunction().symbols().contains(force_expression.symbols()))
            throw SireError::incompatible_error( QObject::tr(
                "You cannot use a force function which uses more symbols "
                "(%1) than the energy function (%2).")
                    .arg( Sire::toString(force_expression.symbols()),
                          Sire::toString(restraintFunction().symbols()) ), CODELOC );
    }
    
    intra_molecule_points = Point::areIntraMoleculePoints(p[0], p[1]) and
                            Point::areIntraMoleculePoints(p[0], p[2]) and
                            Point::areIntraMoleculePoints(p[0], p[3]);

    this->calculatePhi();
}

/** Copy constructor */
DihedralRestraint::DihedralRestraint(const DihedralRestraint &other)
                  : ConcreteProperty<DihedralRestraint,ExpressionRestraint3D>(other),
                    force_expression(other.force_expression),
                    intra_molecule_points(other.intra_molecule_points)
{
    for (int i=0; i<this->nPoints(); ++i)
    {
        p[i] = other.p[i];
    }
}

/** Destructor */
DihedralRestraint::~DihedralRestraint()
{}

/** Copy assignment operator */
DihedralRestraint& DihedralRestraint::operator=(const DihedralRestraint &other)
{
    if (this != &other)
    {
        ExpressionRestraint3D::operator=(other);

        for (int i=0; i<this->nPoints(); ++i)
        {
            p[i] = other.p[i];
        }
        
        force_expression = other.force_expression;
        intra_molecule_points = other.intra_molecule_points;
    }
    
    return *this;
}

/** Comparison operator */
bool DihedralRestraint::operator==(const DihedralRestraint &other) const
{
    return this == &other or
           ( ExpressionRestraint3D::operator==(other) and
             p[0] == other.p[0] and p[1] == other.p[1] and
             p[2] == other.p[2] and p[3] == other.p[3] and
             force_expression == other.force_expression);
}

/** Comparison operator */
bool DihedralRestraint::operator!=(const DihedralRestraint &other) const
{
    return not DihedralRestraint::operator==(other);
}

const char* DihedralRestraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DihedralRestraint>() );
}

/** This restraint involves four points */
int DihedralRestraint::nPoints() const
{
    return 4;
}

/** Return the ith point */
const Point& DihedralRestraint::point(int i) const
{
    i = Index(i).map( this->nPoints() );

    return p[i].read();
}

/** Return the first point */
const Point& DihedralRestraint::point0() const
{
    return p[0].read();
}

/** Return the second point */
const Point& DihedralRestraint::point1() const
{
    return p[1].read();
}

/** Return the third point */
const Point& DihedralRestraint::point2() const
{
    return p[2].read();
}

/** Return the fourth point */
const Point& DihedralRestraint::point3() const
{
    return p[3].read();
}

/** Return the built-in symbols for this restraint */
Symbols DihedralRestraint::builtinSymbols() const
{
    if (this->restraintFunction().isFunction(phi()))
        return phi();
    else
        return Symbols();
}

/** Return the values of the built-in symbols of this restraint */
Values DihedralRestraint::builtinValues() const
{
    if (this->restraintFunction().isFunction(phi()))
        return phi() == this->values()[phi()];
    else
        return Values();
}

/** Return the differential of this restraint with respect to
    the symbol 'symbol'
    
    \throw SireCAS::unavailable_differential
*/
RestraintPtr DihedralRestraint::differentiate(const Symbol &symbol) const
{
    if (this->restraintFunction().isFunction(symbol))
        return DihedralRestraint( p[0], p[1], p[2], p[3],
                                  restraintFunction().differentiate(symbol),
                                  this->values() );
    else
        return NullRestraint();
}

/** Set the space used to evaluate the energy of this restraint

    \throw SireVol::incompatible_space
*/
void DihedralRestraint::setSpace(const Space &new_space)
{
    if (not this->space().equals(new_space))
    {
        DihedralRestraint old_state(*this);
        
        try
        {
            for (int i=0; i<this->nPoints(); ++i)
            {
                p[i].edit().setSpace(new_space);
            }
            
            Restraint3D::setSpace(new_space);
            
            this->calculatePhi();
        }
        catch(...)
        {
            DihedralRestraint::operator=(old_state);
            throw;
        }
    }
}

/** Return the function used to calculate the restraint force */
const Expression& DihedralRestraint::differentialRestraintFunction() const
{
    return force_expression;
}

/** Calculate the force acting on the molecule in the forcetable 'forcetable' 
    caused by this restraint, and add it on to the forcetable scaled by 
    'scale_force' */
void DihedralRestraint::force(MolForceTable &forcetable, double scale_force) const
{
    bool in_p0 = p[0].read().contains(forcetable.molNum());
    bool in_p1 = p[1].read().contains(forcetable.molNum());
    bool in_p2 = p[2].read().contains(forcetable.molNum());
    bool in_p3 = p[3].read().contains(forcetable.molNum());
    
    if (not (in_p0 or in_p1 or in_p2 or in_p3))
        //this molecule is not affected by the restraint
        return;
        
    throw SireError::incomplete_code( QObject::tr(
            "Haven't yet written the code to calculate forces caused "
            "by a dihedral restraint."), CODELOC );
}

/** Calculate the force acting on the molecules in the forcetable 'forcetable' 
    caused by this restraint, and add it on to the forcetable scaled by 
    'scale_force' */
void DihedralRestraint::force(ForceTable &forcetable, double scale_force) const
{
    bool in_p0 = p[0].read().usesMoleculesIn(forcetable);
    bool in_p1 = p[1].read().usesMoleculesIn(forcetable);
    bool in_p2 = p[2].read().usesMoleculesIn(forcetable);
    bool in_p3 = p[3].read().usesMoleculesIn(forcetable);
    
    if (not (in_p0 or in_p1 or in_p2 or in_p3))
        //this molecule is not affected by the restraint
        return;

    throw SireError::incomplete_code( QObject::tr(
            "Haven't yet written the code to calculate forces caused "
            "by a dihedral restraint."), CODELOC );
}

/** Update the points of this restraint using new molecule data from 'moldata'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void DihedralRestraint::update(const MoleculeData &moldata)
{
    if (this->contains(moldata.number()))
    {
        DihedralRestraint old_state(*this);
        
        try
        {
            for (int i=0; i<this->nPoints(); ++i)
            {
                p[i].edit().update(moldata);
            }
            
            this->calculatePhi();
        }
        catch(...)
        {
            DihedralRestraint::operator=(old_state);
            throw;
        }
    }
}
            
/** Update the points of this restraint using new molecule data from 'molecules'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void DihedralRestraint::update(const Molecules &molecules)
{
    if (this->usesMoleculesIn(molecules))
    {
        DihedralRestraint old_state(*this);
        
        try
        {
            for (int i=0; i<this->nPoints(); ++i)
            {
                p[i].edit().update(molecules);
            }
            
            this->calculatePhi();
        }
        catch(...)
        {
            DihedralRestraint::operator=(old_state);
            throw;
        }
    }
}

/** Return the molecules used in this restraint */
Molecules DihedralRestraint::molecules() const
{
    Molecules mols;
    
    for (int i=0; i<this->nPoints(); ++i)
    {
        mols += p[i].read().molecules();
    }
    
    return mols;
}

/** Return whether or not this restraint affects the molecule
    with number 'molnum' */
bool DihedralRestraint::contains(MolNum molnum) const
{
    return p[0].read().contains(molnum) or p[1].read().contains(molnum) or
           p[2].read().contains(molnum) or p[3].read().contains(molnum);
}

/** Return whether or not this restraint affects the molecule
    with ID 'molid' */
bool DihedralRestraint::contains(const MolID &molid) const
{
    return p[0].read().contains(molid) or p[1].read().contains(molid) or
           p[2].read().contains(molid) or p[3].read().contains(molid);
}
    
/** Return whether or not this restraint involves any of the molecules
    that are in the forcetable 'forcetable' */
bool DihedralRestraint::usesMoleculesIn(const ForceTable &forcetable) const
{
    return p[0].read().usesMoleculesIn(forcetable) or 
           p[1].read().usesMoleculesIn(forcetable) or
           p[2].read().usesMoleculesIn(forcetable) or
           p[3].read().usesMoleculesIn(forcetable);
}

/** Return whether or not this restraint involves any of the molecules
    in 'molecules' */
bool DihedralRestraint::usesMoleculesIn(const Molecules &molecules) const
{
    return p[0].read().usesMoleculesIn(molecules) or
           p[1].read().usesMoleculesIn(molecules) or
           p[2].read().usesMoleculesIn(molecules) or
           p[3].read().usesMoleculesIn(molecules);
}

static Expression harmonicFunction(double force_constant)
{
    if (SireMaths::isZero(force_constant))
        return 0;
    else
        return force_constant * pow(DihedralRestraint::phi(), 2);
}

static Expression diffHarmonicFunction(double force_constant)
{
    if (SireMaths::isZero(force_constant))
        return 0;
    else
        return (2*force_constant) * DihedralRestraint::phi();
}

/** Return a distance restraint that applies a harmonic potential between 
    the points 'point0' and 'point1' using a force constant 'force_constant' */
DihedralRestraint DihedralRestraint::harmonic(
                                        const PointRef &point0,
                                        const PointRef &point1,
                                        const PointRef &point2,
                                        const PointRef &point3,
                                        const HarmonicAngleForceConstant &force_constant)
{
    return DihedralRestraint(point0, point1, point2, point3,
                             ::harmonicFunction(force_constant), 
                             ::diffHarmonicFunction(force_constant));
}

static Expression halfHarmonicFunction(double force_constant, double angle)
{
    if ( SireMaths::isZero(force_constant) )
        return 0;
        
    else if ( angle <= 0 )
        //this is just a harmonic function
        return ::harmonicFunction(force_constant);

    else
    {
        const Symbol &phi = DihedralRestraint::phi();
        return Conditional( 
                GreaterThan(phi, angle), force_constant * pow(phi-angle, 2), 0 );
    }
}

static Expression diffHalfHarmonicFunction(double force_constant, double angle)
{
    if ( SireMaths::isZero(force_constant) )
        return 0;
    
    else if (angle <= 0)
        //this is just a harmonic function
        return ::diffHarmonicFunction(force_constant);
    
    else
    {
        const Symbol &phi = DihedralRestraint::phi();
        return Conditional( GreaterThan(phi, angle), 
                                (2*force_constant) * (phi-angle), 0 );
    }
}

/** Return a distance restraint that applied a half-harmonic potential 
    between the points 'point0' and 'point1' above a distance 'distance'
    using a force constant 'force_constant' */
DihedralRestraint DihedralRestraint::halfHarmonic(
                                        const PointRef &point0,
                                        const PointRef &point1,
                                        const PointRef &point2,
                                        const PointRef &point3,
                                        const Angle &angle,
                                        const HarmonicAngleForceConstant &force_constant)
{
    double ang = angle.to(radians);

    return DihedralRestraint(point0, point1, point2, point3,
                          ::halfHarmonicFunction(force_constant, ang),
                          ::diffHalfHarmonicFunction(force_constant, ang));
}
