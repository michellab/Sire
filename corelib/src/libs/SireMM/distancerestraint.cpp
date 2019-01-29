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

#include "distancerestraint.h"

#include "SireFF/forcetable.h"

#include "SireCAS/symbols.h"
#include "SireCAS/values.h"
#include "SireCAS/conditional.h"
#include "SireCAS/power.h"

#include "SireID/index.h"

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
using namespace SireUnits::Dimension;

////////////
//////////// Implementation of DistanceRestraint
////////////

static const RegisterMetaType<DistanceRestraint> r_distrest;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                      const DistanceRestraint &distrest)
{
    writeHeader(ds, r_distrest, 1);
    
    SharedDataStream sds(ds);
    
    sds << distrest.p[0] << distrest.p[1]
        << distrest.force_expression
        << static_cast<const ExpressionRestraint3D&>(distrest);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, DistanceRestraint &distrest)
{
    VersionID v = readHeader(ds, r_distrest);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> distrest.p[0] >> distrest.p[1]
            >> distrest.force_expression
            >> static_cast<ExpressionRestraint3D&>(distrest);

        distrest.intra_molecule_points = Point::areIntraMoleculePoints(distrest.p[0],
                                                                       distrest.p[1]);
    }
    else
        throw version_error( v, "1", r_distrest, CODELOC );
        
    return ds;
}

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getRSymbol, ("r") );

/** Return the symbol that represents the distance between the
    two points ("r") */
const Symbol& DistanceRestraint::r()
{
    return *(getRSymbol());
}

/** Constructor */
DistanceRestraint::DistanceRestraint()
                  : ConcreteProperty<DistanceRestraint,ExpressionRestraint3D>()
{}

/** Construct a restraint that acts between the two points 'point0' and 'point1',
    restraining the distance between these points using the expression 
    'restraint' */
DistanceRestraint::DistanceRestraint(const PointRef &point0, const PointRef &point1,
                                     const Expression &restraint)
         : ConcreteProperty<DistanceRestraint,ExpressionRestraint3D>(restraint)
{
    p[0] = point0;
    p[1] = point1;

    force_expression = this->restraintFunction().differentiate(r());
        
    if (force_expression.isConstant())
        force_expression = force_expression.evaluate(Values());
    
    intra_molecule_points = Point::areIntraMoleculePoints(p[0], p[1]);
    
    this->calculateR();
}

/** Construct a restraint that acts between the two points 'point0' and 'point1',
    restraining the distance between these points using the expression 
    'restraint', with supplied values for this expression in 'values'.
    Note that any extra values in 'values' that aren't in the expression
    'restraint' are ignored */
DistanceRestraint::DistanceRestraint(const PointRef &point0, const PointRef &point1,
                                     const Expression &restraint, const Values &values)
         : ConcreteProperty<DistanceRestraint,ExpressionRestraint3D>(restraint, values)
{
    p[0] = point0;
    p[1] = point1;

    force_expression = this->restraintFunction().differentiate(r());
        
    if (force_expression.isConstant())
        force_expression = force_expression.evaluate(Values());
    
    intra_molecule_points = Point::areIntraMoleculePoints(p[0], p[1]);
    
    this->calculateR();
}

/** Internal constructor used to construct a restraint between the two 
    points 'point0' and 'point1' that restrains the distance between
    the two points using the function 'nrg_restraint' and
    calculates the force using the function 'force_restraint'. It 
    is a good idea to ensure that 'force_restraint' really is the
    differential of 'nrg_restraint' with respect to r() */
DistanceRestraint::DistanceRestraint(const PointRef &point0, const PointRef &point1,
                                     const Expression &nrg_restraint,
                                     const Expression &force_restraint)
      : ConcreteProperty<DistanceRestraint,ExpressionRestraint3D>(nrg_restraint),
        force_expression(force_restraint)
{
    p[0] = point0;
    p[1] = point1;

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
    
    intra_molecule_points = Point::areIntraMoleculePoints(p[0], p[1]);
    
    this->calculateR();
}

/** Copy constructor */
DistanceRestraint::DistanceRestraint(const DistanceRestraint &other)
      : ConcreteProperty<DistanceRestraint,ExpressionRestraint3D>(other),
        force_expression(other.force_expression),
        intra_molecule_points(other.intra_molecule_points)
{
    for (int i=0; i<2; ++i)
    {
        p[i] = other.p[i];
    }
}

/** Destructor */
DistanceRestraint::~DistanceRestraint()
{}

/** Copy assignment operator */
DistanceRestraint& DistanceRestraint::operator=(const DistanceRestraint &other)
{
    if (this != &other)
    {
        ExpressionRestraint3D::operator=(other);

        for (int i=0; i<2; ++i)
        {
            p[i] = other.p[i];
        }
        
        force_expression = other.force_expression;
        intra_molecule_points = other.intra_molecule_points;
    }
    
    return *this;
}

/** Comparison operator */
bool DistanceRestraint::operator==(const DistanceRestraint &other) const
{
    return this == &other or
           ( ExpressionRestraint3D::operator==(other) and
             p[0] == other.p[0] and p[1] == other.p[1] and
             force_expression == other.force_expression);
}

/** Comparison operator */
bool DistanceRestraint::operator!=(const DistanceRestraint &other) const
{
    return not DistanceRestraint::operator==(other);
}

const char* DistanceRestraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DistanceRestraint>() );
}

/** This restraint involves two points */
int DistanceRestraint::nPoints() const
{
    return 2;
}

/** Return the ith point */
const Point& DistanceRestraint::point(int i) const
{
    i = Index(i).map( this->nPoints() );

    return p[i].read();
}

/** Return the first point */
const Point& DistanceRestraint::point0() const
{
    return p[0].read();
}

/** Return the second point */
const Point& DistanceRestraint::point1() const
{
    return p[1].read();
}

/** Calculate the distance between the two points of this restraint */
void DistanceRestraint::calculateR()
{
    if (this->restraintFunction().isFunction(r()))
    {
        double distance;

        if (intra_molecule_points)
            //we don't use the space when calculating intra-molecular distances
            distance = Vector::distance( p[0].read().point(), p[1].read().point() );
        else
            distance = this->space().calcDist( p[0].read().point(),
                                               p[1].read().point() );

        ExpressionRestraint3D::_pvt_setValue(r(), distance);
    }
}

/** Set the space used to evaluate the energy of this restraint

    \throw SireVol::incompatible_space
*/
void DistanceRestraint::setSpace(const Space &new_space)
{
    if (not this->space().equals(new_space))
    {
        PointPtr old_p0 = p[0];
    
        try
        {
            p[0].edit().setSpace(new_space);
            p[1].edit().setSpace(new_space);
            
            ExpressionRestraint3D::setSpace(new_space);
            
            this->calculateR();
        }
        catch(...)
        {
            p[0] = old_p0;
            throw;
        }
    }
}

/** Return the built-in symbols of this restraint */
Symbols DistanceRestraint::builtinSymbols() const
{
    if (this->restraintFunction().isFunction(r()))
        return r();
    else
        return Symbols();
}

/** Return the built-in values of this restraint */
Values DistanceRestraint::builtinValues() const
{
    if (this->restraintFunction().isFunction(r()))
        return (r() == this->values()[r()]);
    else
        return Values();
}

/** Return the restraint that is the differential of this restraint
    with respect to the symbol 'symbol'
    
    \throw SireCAS::unavailable_differential
*/
RestraintPtr DistanceRestraint::differentiate(const Symbol &symbol) const
{
    if (this->restraintFunction().isFunction(symbol))
    {
        return DistanceRestraint( p[0], p[1],
                                  this->restraintFunction().differentiate(symbol),
                                  this->values() );
    }
    else
        return NullRestraint();
}

/** Return the function used to calculate the restraint force */
const Expression& DistanceRestraint::differentialRestraintFunction() const
{
    return force_expression;
}

template<class T>
static void addForce(const Point &p0, const Point &p1, const Space &space,
                     bool intra_molecule_points,
                     bool in_p0, bool in_p1, const double force,
                     T &forcetable)
{
    if ( SireMaths::isZero(force) )
        return;
    
    Vector delta;

    if (intra_molecule_points)
        delta = p1.point() - p0.point();
    else
    {
        const Vector &point0 = p0.point();
        Vector point1 = p1.point();
    
        point1 = space.getMinimumImage(point1, point0);
    
        delta = point1 - point0;
    }

    const double distance = delta.length();
    
    if ( SireMaths::isZero(distance) )
        return;

    const Vector f = (force / distance) * delta;
    
    if (not f.isZero())
    {
        if (in_p0)
            p0.addForce(forcetable, -f);

        if (in_p1)
            p1.addForce(forcetable, f);
    }
}

/** Calculate the force acting on the molecule in the forcetable 'forcetable' 
    caused by this restraint, and add it on to the forcetable scaled by 
    'scale_force' */
void DistanceRestraint::force(MolForceTable &forcetable, double scale_force) const
{
    bool in_p0 = p[0].read().contains(forcetable.molNum());
    bool in_p1 = p[1].read().contains(forcetable.molNum());
    
    if (in_p0 or in_p1)
    {
        const double force = scale_force * force_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[0], p[1], space(), intra_molecule_points,
                 in_p0, in_p1, force, forcetable);
    }
}

/** Calculate the force acting on the molecules in the forcetable 'forcetable' 
    caused by this restraint, and add it on to the forcetable scaled by 
    'scale_force' */
void DistanceRestraint::force(ForceTable &forcetable, double scale_force) const
{
    bool in_p0 = p[0].read().usesMoleculesIn(forcetable);
    bool in_p1 = p[1].read().usesMoleculesIn(forcetable);
    
    if (in_p0 or in_p1)
    {
        const double force = scale_force * force_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[0], p[1], space(), intra_molecule_points,
                 in_p0, in_p1, force, forcetable);
    }
}

/** Update the points of this restraint using new molecule data from 'moldata'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void DistanceRestraint::update(const MoleculeData &moldata)
{
    bool in_p0 = p[0].read().contains(moldata.number());
    bool in_p1 = p[1].read().contains(moldata.number());
    
    if (in_p0 and in_p1)
    {
        PointPtr old_p0 = p[1];
        
        try
        {
            p[0].edit().update(moldata);
            p[1].edit().update(moldata);
            
            this->calculateR();
        }
        catch(...)
        {
            p[0] = old_p0;
            throw;
        }
    }
    else if (in_p0)
    {
        p[0].edit().update(moldata);
        this->calculateR();
    }
    else if (in_p1)
    {
        p[1].edit().update(moldata);
        this->calculateR();
    }
}
            
/** Update the points of this restraint using new molecule data from 'molecules'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void DistanceRestraint::update(const Molecules &molecules)
{
    bool in_p0 = p[0].read().usesMoleculesIn(molecules);
    bool in_p1 = p[1].read().usesMoleculesIn(molecules);
    
    if (in_p0 and in_p1)
    {
        PointPtr old_p0 = p[0];
        
        try
        {
            p[0].edit().update(molecules);
            p[1].edit().update(molecules);
            this->calculateR();
        }
        catch(...)
        {
            p[0] = old_p0;
            throw;
        }
    }
    else if (in_p0)
    {
        p[0].edit().update(molecules);
        this->calculateR();
    }
    else if (in_p1)
    {
        p[1].edit().update(molecules);
        this->calculateR();
    }
}

/** Return the molecules used in this restraint */
Molecules DistanceRestraint::molecules() const
{
    int n0 = p[0].read().nMolecules();
    int n1 = p[1].read().nMolecules();
    
    if (n0 != 0 and n1 != 0)
    {
        Molecules mols = p[0].read().molecules();
        mols += p[1].read().molecules();
        
        return mols;
    }
    else if (n0 != 0)
    {
        return p[0].read().molecules();
    }
    else if (n1 != 0)
    {
        return p[1].read().molecules();
    }
    else
        return Molecules();
}

/** Return whether or not this restraint affects the molecule
    with number 'molnum' */
bool DistanceRestraint::contains(MolNum molnum) const
{
    return p[0].read().contains(molnum) or p[1].read().contains(molnum);
}

/** Return whether or not this restraint affects the molecule
    with ID 'molid' */
bool DistanceRestraint::contains(const MolID &molid) const
{
    return p[0].read().contains(molid) or p[1].read().contains(molid);
}
    
/** Return whether or not this restraint involves any of the molecules
    that are in the forcetable 'forcetable' */
bool DistanceRestraint::usesMoleculesIn(const ForceTable &forcetable) const
{
    return p[0].read().usesMoleculesIn(forcetable) or 
           p[1].read().usesMoleculesIn(forcetable);
}

/** Return whether or not this restraint involves any of the molecules
    in 'molecules' */
bool DistanceRestraint::usesMoleculesIn(const Molecules &molecules) const
{
    return p[0].read().usesMoleculesIn(molecules) or
           p[1].read().usesMoleculesIn(molecules);
}

static Expression harmonicFunction(double force_constant)
{
    if (SireMaths::isZero(force_constant))
        return 0;
    else
        return force_constant * pow(DistanceRestraint::r(), 2);
}

static Expression diffHarmonicFunction(double force_constant)
{
    if (SireMaths::isZero(force_constant))
        return 0;
    else
        return (2*force_constant) * DistanceRestraint::r();
}

/** Return a distance restraint that applies a harmonic potential between 
    the points 'point0' and 'point1' using a force constant 'force_constant' */
DistanceRestraint DistanceRestraint::harmonic(
                                  const PointRef &point0,
                                  const PointRef &point1,
                                  const HarmonicDistanceForceConstant &force_constant)
{
    return DistanceRestraint(point0, point1,
                             ::harmonicFunction(force_constant), 
                             ::diffHarmonicFunction(force_constant));
}

static Expression halfHarmonicFunction(double force_constant, double length)
{
    if ( SireMaths::isZero(force_constant) )
        return 0;
        
    else if ( length <= 0 )
        //this is just a harmonic function
        return ::harmonicFunction(force_constant);

    else
    {
        const Symbol &r = DistanceRestraint::r();
        return Conditional( 
                GreaterThan(r, length), force_constant * pow(r-length, 2), 0 );
    }
}

static Expression diffHalfHarmonicFunction(double force_constant, double distance)
{
    if ( SireMaths::isZero(force_constant) )
        return 0;
    
    else if (distance <= 0)
        //this is just a harmonic function
        return ::diffHarmonicFunction(force_constant);
    
    else
    {
        const Symbol &r = DistanceRestraint::r();
        return Conditional( GreaterThan(r, distance), 
                                (2*force_constant) * (r-distance), 0 );
    }
}

/** Return a distance restraint that applied a half-harmonic potential 
    between the points 'point0' and 'point1' above a distance 'distance'
    using a force constant 'force_constant' */
DistanceRestraint DistanceRestraint::halfHarmonic(
                                      const PointRef &point0,
                                      const PointRef &point1,
                                      const Length &distance,
                                      const HarmonicDistanceForceConstant &force_constant)
{
    return DistanceRestraint(point0, point1,
                             ::halfHarmonicFunction(force_constant, distance),
                             ::diffHalfHarmonicFunction(force_constant, distance));
}

////////////
//////////// Implementation of DoubleDistanceRestraint
////////////

static const RegisterMetaType<DoubleDistanceRestraint> r_doubledistrest;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                      const DoubleDistanceRestraint &doubledistrest)
{
    writeHeader(ds, r_doubledistrest, 1);
    
    SharedDataStream sds(ds);
    
    sds << doubledistrest.p[0] << doubledistrest.p[1]
        << doubledistrest.p[2] << doubledistrest.p[3]
        << doubledistrest.force01_expression
        << doubledistrest.force23_expression
        << static_cast<const ExpressionRestraint3D&>(doubledistrest);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, 
                                      DoubleDistanceRestraint &doubledistrest)
{
    VersionID v = readHeader(ds, r_doubledistrest);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> doubledistrest.p[0] >> doubledistrest.p[1]
            >> doubledistrest.p[2] >> doubledistrest.p[3]
            >> doubledistrest.force01_expression
            >> doubledistrest.force23_expression
            >> static_cast<ExpressionRestraint3D&>(doubledistrest);

        doubledistrest.intra_molecule_points01 = Point::areIntraMoleculePoints(
                                                                    doubledistrest.p[0],
                                                                    doubledistrest.p[1]);
                                                                    
        doubledistrest.intra_molecule_points23 = Point::areIntraMoleculePoints(
                                                                    doubledistrest.p[2],
                                                                    doubledistrest.p[3]);
    }
    else
        throw version_error( v, "1", r_doubledistrest, CODELOC );
        
    return ds;
}

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR01Symbol, ("r01") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR23Symbol, ("r23") );

/** Return the symbol that represents the distance between the
    points 0 and 1 ("r01") */
const Symbol& DoubleDistanceRestraint::r01()
{
    return *(getR01Symbol());
}

/** Return the symbol that represents the distance between the
    points 2 and 3 ("r23") */
const Symbol& DoubleDistanceRestraint::r23()
{
    return *(getR23Symbol());
}

/** Constructor */
DoubleDistanceRestraint::DoubleDistanceRestraint()
          : ConcreteProperty<DoubleDistanceRestraint,ExpressionRestraint3D>()
{}

void DoubleDistanceRestraint::calculateR()
{
    if (this->restraintFunction().isFunction(r01()))
    {
        double d01;
        
        if (intra_molecule_points01)
            //we don't use the space when calculating intra-molecular distances
            d01 = Vector::distance( p[0].read().point(), p[1].read().point() );
        else
            d01 = this->space().calcDist( p[0].read().point(),
                                          p[1].read().point() );
        
        ExpressionRestraint3D::_pvt_setValue( r01(), d01 );
    }
    
    if (this->restraintFunction().isFunction(r23()))
    {
        double d23;
        
        if (intra_molecule_points23)
            d23 = Vector::distance( p[2].read().point(), p[3].read().point() );
        else
            d23 = this->space().calcDist( p[2].read().point(),
                                          p[3].read().point() );

        ExpressionRestraint3D::_pvt_setValue( r23(), d23 );
    }
}

/** Construct a restraint that acts on the two distances defined
    using the passed four points, using the expression 'restraint' */
DoubleDistanceRestraint::DoubleDistanceRestraint(const PointRef &point0, 
                                                 const PointRef &point1,
                                                 const PointRef &point2,
                                                 const PointRef &point3,
                                                 const Expression &restraint)
  : ConcreteProperty<DoubleDistanceRestraint,ExpressionRestraint3D>(restraint)
{
    p[0] = point0;
    p[1] = point1;
    p[2] = point2;
    p[3] = point3;

    force01_expression = this->restraintFunction().differentiate(r01());
    force23_expression = this->restraintFunction().differentiate(r23());
        
    if (force01_expression.isConstant())
        force01_expression = force01_expression.evaluate(Values());
        
    if (force23_expression.isConstant())
        force23_expression = force23_expression.evaluate(Values());
    
    intra_molecule_points01 = Point::areIntraMoleculePoints(p[0], p[1]);
    intra_molecule_points23 = Point::areIntraMoleculePoints(p[2], p[3]);

    this->calculateR();
}

/** Construct a restraint that acts on the two distances defined
    using the passed four points, using the expression 'restraint'
    and supplied user values in 'values' */
DoubleDistanceRestraint::DoubleDistanceRestraint(const PointRef &point0, 
                                                 const PointRef &point1,
                                                 const PointRef &point2,
                                                 const PointRef &point3,
                                                 const Expression &restraint,
                                                 const Values &values)
  : ConcreteProperty<DoubleDistanceRestraint,ExpressionRestraint3D>(restraint, values)
{
    p[0] = point0;
    p[1] = point1;
    p[2] = point2;
    p[3] = point3;

    force01_expression = this->restraintFunction().differentiate(r01());
    force23_expression = this->restraintFunction().differentiate(r23());
        
    if (force01_expression.isConstant())
        force01_expression = force01_expression.evaluate(Values());
        
    if (force23_expression.isConstant())
        force23_expression = force23_expression.evaluate(Values());
    
    intra_molecule_points01 = Point::areIntraMoleculePoints(p[0], p[1]);
    intra_molecule_points23 = Point::areIntraMoleculePoints(p[2], p[3]);

    this->calculateR();
}

/** Copy constructor */
DoubleDistanceRestraint::DoubleDistanceRestraint(const DoubleDistanceRestraint &other)
          : ConcreteProperty<DoubleDistanceRestraint,ExpressionRestraint3D>(other),
            force01_expression(other.force01_expression),
            force23_expression(other.force23_expression),
            intra_molecule_points01(other.intra_molecule_points01),
            intra_molecule_points23(other.intra_molecule_points23)
{
    for (int i=0; i<this->nPoints(); ++i)
    {
        p[i] = other.p[i];
    }
}

/** Destructor */
DoubleDistanceRestraint::~DoubleDistanceRestraint()
{}

/** Copy assignment operator */
DoubleDistanceRestraint& DoubleDistanceRestraint::operator=(
                                                    const DoubleDistanceRestraint &other)
{
    if (this != &other)
    {
        ExpressionRestraint3D::operator=(other);

        for (int i=0; i<this->nPoints(); ++i)
        {
            p[i] = other.p[i];
        }
        
        force01_expression = other.force01_expression;
        force23_expression = other.force23_expression;
        intra_molecule_points01 = other.intra_molecule_points01;
        intra_molecule_points23 = other.intra_molecule_points23;
    }
    
    return *this;
}

/** Comparison operator */
bool DoubleDistanceRestraint::operator==(const DoubleDistanceRestraint &other) const
{
    return this == &other or
           ( ExpressionRestraint3D::operator==(other) and
             p[0] == other.p[0] and p[1] == other.p[1] and
             p[2] == other.p[2] and p[3] == other.p[3] and
             force01_expression == other.force01_expression and
             force23_expression == other.force23_expression);
}

/** Comparison operator */
bool DoubleDistanceRestraint::operator!=(const DoubleDistanceRestraint &other) const
{
    return not DoubleDistanceRestraint::operator==(other);
}

const char* DoubleDistanceRestraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DoubleDistanceRestraint>() );
}

/** This restraint involves four points */
int DoubleDistanceRestraint::nPoints() const
{
    return 4;
}

/** Return the ith point */
const Point& DoubleDistanceRestraint::point(int i) const
{
    i = Index(i).map( this->nPoints() );

    return p[i].read();
}

/** Return the first point */
const Point& DoubleDistanceRestraint::point0() const
{
    return p[0].read();
}

/** Return the second point */
const Point& DoubleDistanceRestraint::point1() const
{
    return p[1].read();
}

/** Return the third point */
const Point& DoubleDistanceRestraint::point2() const
{
    return p[2].read();
}

/** Return the fourth point */
const Point& DoubleDistanceRestraint::point3() const
{
    return p[3].read();
}

/** Return the built in symbols for this restraint */
Symbols DoubleDistanceRestraint::builtinSymbols() const
{
    Symbols symbols;

    if (this->restraintFunction().isFunction(r01()))
        symbols += r01();
    
    if (this->restraintFunction().isFunction(r23()))
        symbols += r23();
        
    return symbols;
}

/** Return the built in values for this restraint */
Values DoubleDistanceRestraint::builtinValues() const
{
    Values vals;
    
    if (this->restraintFunction().isFunction(r01()))
        vals.set( r01(), this->values()[r01()] );
        
    if (this->restraintFunction().isFunction(r23()))
        vals.set( r23(), this->values()[r23()] );
        
    return vals;
}

/** Return the differential of this restraint with respect to the 
    symbol 'symbol'
    
    \throw SireCAS::unavailable_differential
*/
RestraintPtr DoubleDistanceRestraint::differentiate(const Symbol &symbol) const
{
    if (this->restraintFunction().isFunction(symbol))
    {
        return DoubleDistanceRestraint( p[0], p[1], p[2], p[3],
                                        restraintFunction().differentiate(symbol),
                                        this->values() );
    }
    else
        return NullRestraint();
}

/** Set the space used to evaluate the energy of this restraint

    \throw SireVol::incompatible_space
*/
void DoubleDistanceRestraint::setSpace(const Space &new_space)
{
    if (not this->space().equals(new_space))
    {
        DoubleDistanceRestraint old_state(*this);
    
        try
        {
            for (int i=0; i<this->nPoints(); ++i)
            {
                p[i].edit().setSpace(new_space);
            }
            
            ExpressionRestraint3D::setSpace(new_space);
            
            this->calculateR();
        }
        catch(...)
        {
            DoubleDistanceRestraint::operator=(old_state);
            throw;
        }
    }
}

/** Return the function used to calculate the restraint force along the 
    distance r01 */
const Expression& DoubleDistanceRestraint::differentialRestraintFunction01() const
{
    return force01_expression;
}

/** Return the function used to calculate the restraint force along the 
    distance r23 */
const Expression& DoubleDistanceRestraint::differentialRestraintFunction23() const
{
    return force23_expression;
}

/** Calculate the force acting on the molecule in the forcetable 'forcetable' 
    caused by this restraint, and add it on to the forcetable scaled by 
    'scale_force' */
void DoubleDistanceRestraint::force(MolForceTable &forcetable, double scale_force) const
{
    bool in_p0 = p[0].read().contains(forcetable.molNum());
    bool in_p1 = p[1].read().contains(forcetable.molNum());
    bool in_p2 = p[2].read().contains(forcetable.molNum());
    bool in_p3 = p[3].read().contains(forcetable.molNum());

    if (in_p0 or in_p1)
    {
        const double force = scale_force * force01_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[0], p[1], space(), intra_molecule_points01,
                 in_p0, in_p1, force, forcetable);
    }
    
    if (in_p2 or in_p3)
    {
        const double force = scale_force * force23_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[2], p[3], space(), intra_molecule_points23,
                 in_p2, in_p3, force, forcetable);
    }
}

/** Calculate the force acting on the molecules in the forcetable 'forcetable' 
    caused by this restraint, and add it on to the forcetable scaled by 
    'scale_force' */
void DoubleDistanceRestraint::force(ForceTable &forcetable, double scale_force) const
{
    bool in_p0 = p[0].read().usesMoleculesIn(forcetable);
    bool in_p1 = p[1].read().usesMoleculesIn(forcetable);
    bool in_p2 = p[2].read().usesMoleculesIn(forcetable);
    bool in_p3 = p[3].read().usesMoleculesIn(forcetable);

    if (in_p0 or in_p1)
    {
        const double force = scale_force * force01_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[0], p[1], space(), intra_molecule_points01,
                 in_p0, in_p1, force, forcetable);
    }
    
    if (in_p2 or in_p3)
    {
        const double force = scale_force * force23_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[2], p[3], space(), intra_molecule_points23,
                 in_p2, in_p3, force, forcetable);
    }
}

/** Update the points of this restraint using new molecule data from 'moldata'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void DoubleDistanceRestraint::update(const MoleculeData &moldata)
{
    if (this->contains(moldata.number()))
    {
        DoubleDistanceRestraint old_state(*this);
        
        try
        {
            for (int i=0; i<this->nPoints(); ++i)
            {
                p[i].edit().update(moldata);
            }
            
            this->calculateR();
        }
        catch(...)
        {
            DoubleDistanceRestraint::operator=(old_state);
            throw;
        }
    }
}
            
/** Update the points of this restraint using new molecule data from 'molecules'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void DoubleDistanceRestraint::update(const Molecules &molecules)
{
    if (this->usesMoleculesIn(molecules))
    {
        DoubleDistanceRestraint old_state(*this);
        
        try
        {
            for (int i=0; i<this->nPoints(); ++i)
            {
                p[i].edit().update(molecules);
            }
            
            this->calculateR();
        }
        catch(...)
        {
            DoubleDistanceRestraint::operator=(old_state);
            throw;
        }
    }
}

/** Return the molecules used in this restraint */
Molecules DoubleDistanceRestraint::molecules() const
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
bool DoubleDistanceRestraint::contains(MolNum molnum) const
{
    return p[0].read().contains(molnum) or p[1].read().contains(molnum) or
           p[2].read().contains(molnum) or p[2].read().contains(molnum);
}

/** Return whether or not this restraint affects the molecule
    with ID 'molid' */
bool DoubleDistanceRestraint::contains(const MolID &molid) const
{
    return p[0].read().contains(molid) or p[1].read().contains(molid) or
           p[2].read().contains(molid) or p[3].read().contains(molid);
}
    
/** Return whether or not this restraint involves any of the molecules
    that are in the forcetable 'forcetable' */
bool DoubleDistanceRestraint::usesMoleculesIn(const ForceTable &forcetable) const
{
    return p[0].read().usesMoleculesIn(forcetable) or 
           p[1].read().usesMoleculesIn(forcetable) or
           p[2].read().usesMoleculesIn(forcetable) or
           p[3].read().usesMoleculesIn(forcetable);
}

/** Return whether or not this restraint involves any of the molecules
    in 'molecules' */
bool DoubleDistanceRestraint::usesMoleculesIn(const Molecules &molecules) const
{
    return p[0].read().usesMoleculesIn(molecules) or
           p[1].read().usesMoleculesIn(molecules) or
           p[2].read().usesMoleculesIn(molecules) or
           p[3].read().usesMoleculesIn(molecules);
}

////////////
//////////// Implementation of TripleDistanceRestraint
////////////

static const RegisterMetaType<TripleDistanceRestraint> r_tripledistrest;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                      const TripleDistanceRestraint &tripledistrest)
{
    writeHeader(ds, r_tripledistrest, 1);
    
    SharedDataStream sds(ds);
    
    sds << tripledistrest.p[0] << tripledistrest.p[1]
        << tripledistrest.p[2] << tripledistrest.p[3]
        << tripledistrest.p[4] << tripledistrest.p[5]
        << tripledistrest.force01_expression
        << tripledistrest.force23_expression
        << tripledistrest.force45_expression
        << static_cast<const ExpressionRestraint3D&>(tripledistrest);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, 
                                      TripleDistanceRestraint &tripledistrest)
{
    VersionID v = readHeader(ds, r_tripledistrest);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> tripledistrest.p[0] >> tripledistrest.p[1]
            >> tripledistrest.p[2] >> tripledistrest.p[3]
            >> tripledistrest.p[4] >> tripledistrest.p[5]
            >> tripledistrest.force01_expression
            >> tripledistrest.force23_expression
            >> tripledistrest.force45_expression
            >> static_cast<ExpressionRestraint3D&>(tripledistrest);

        tripledistrest.intra_molecule_points01 = Point::areIntraMoleculePoints(
                                                                    tripledistrest.p[0],
                                                                    tripledistrest.p[1]);
                                                                    
        tripledistrest.intra_molecule_points23 = Point::areIntraMoleculePoints(
                                                                    tripledistrest.p[2],
                                                                    tripledistrest.p[3]);
                                                                    
        tripledistrest.intra_molecule_points45 = Point::areIntraMoleculePoints(
                                                                    tripledistrest.p[4],
                                                                    tripledistrest.p[5]);
    }
    else
        throw version_error( v, "1", r_tripledistrest, CODELOC );
        
    return ds;
}

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR45Symbol, ("r45") );

/** Return the symbol that represents the distance between the
    points 0 and 1 ("r01") */
const Symbol& TripleDistanceRestraint::r01()
{
    return *(getR01Symbol());
}

/** Return the symbol that represents the distance between the
    points 2 and 3 ("r23") */
const Symbol& TripleDistanceRestraint::r23()
{
    return *(getR23Symbol());
}

/** Return the symbol that represents the distance between the
    points 4 and 5 ("r45") */
const Symbol& TripleDistanceRestraint::r45()
{
    return *(getR45Symbol());
}

/** Constructor */
TripleDistanceRestraint::TripleDistanceRestraint()
    : ConcreteProperty<TripleDistanceRestraint,ExpressionRestraint3D>()
{}

void TripleDistanceRestraint::calculateR()
{
    if (this->restraintFunction().isFunction(r01()))
    {
        double d01;
        
        if (intra_molecule_points01)
            //we don't use the space when calculating intra-molecular distances
            d01 = Vector::distance( p[0].read().point(), p[1].read().point() );
        else
            d01 = this->space().calcDist( p[0].read().point(),
                                          p[1].read().point() );
        
        ExpressionRestraint3D::_pvt_setValue( r01(), d01 );
    }
    
    if (this->restraintFunction().isFunction(r23()))
    {
        double d23;
        
        if (intra_molecule_points23)
            d23 = Vector::distance( p[2].read().point(), p[3].read().point() );
        else
            d23 = this->space().calcDist( p[2].read().point(),
                                          p[3].read().point() );

        ExpressionRestraint3D::_pvt_setValue( r23(), d23 );
    }
    
    if (this->restraintFunction().isFunction(r45()))
    {
        double d45;
        
        if (intra_molecule_points45)
            d45 = Vector::distance( p[4].read().point(), p[5].read().point() );
        else
            d45 = this->space().calcDist( p[4].read().point(),
                                          p[5].read().point() );

        ExpressionRestraint3D::_pvt_setValue( r45(), d45 );
    }
}

/** Construct a restraint that acts on the three distances defined
    using the passed six points, using the expression 'restraint' */
TripleDistanceRestraint::TripleDistanceRestraint(const PointRef &point0, 
                                                 const PointRef &point1,
                                                 const PointRef &point2,
                                                 const PointRef &point3,
                                                 const PointRef &point4,
                                                 const PointRef &point5,
                                                 const Expression &restraint)
    : ConcreteProperty<TripleDistanceRestraint,ExpressionRestraint3D>(restraint)
{
    p[0] = point0;
    p[1] = point1;
    p[2] = point2;
    p[3] = point3;
    p[4] = point4;
    p[5] = point5;

    force01_expression = this->restraintFunction().differentiate(r01());
    force23_expression = this->restraintFunction().differentiate(r23());
    force45_expression = this->restraintFunction().differentiate(r45());
        
    if (force01_expression.isConstant())
        force01_expression = force01_expression.evaluate(Values());
        
    if (force23_expression.isConstant())
        force23_expression = force23_expression.evaluate(Values());

    if (force45_expression.isConstant())
        force45_expression = force45_expression.evaluate(Values());
    
    intra_molecule_points01 = Point::areIntraMoleculePoints(p[0], p[1]);
    intra_molecule_points23 = Point::areIntraMoleculePoints(p[2], p[3]);
    intra_molecule_points45 = Point::areIntraMoleculePoints(p[4], p[5]);

    this->calculateR();
}

/** Construct a restraint that acts on the three distances defined
    using the passed six points, using the expression 'restraint' */
TripleDistanceRestraint::TripleDistanceRestraint(const PointRef &point0, 
                                                 const PointRef &point1,
                                                 const PointRef &point2,
                                                 const PointRef &point3,
                                                 const PointRef &point4,
                                                 const PointRef &point5,
                                                 const Expression &restraint,
                                                 const Values &values)
    : ConcreteProperty<TripleDistanceRestraint,ExpressionRestraint3D>(restraint, values)
{
    p[0] = point0;
    p[1] = point1;
    p[2] = point2;
    p[3] = point3;
    p[4] = point4;
    p[5] = point5;

    force01_expression = this->restraintFunction().differentiate(r01());
    force23_expression = this->restraintFunction().differentiate(r23());
    force45_expression = this->restraintFunction().differentiate(r45());
        
    if (force01_expression.isConstant())
        force01_expression = force01_expression.evaluate(Values());
        
    if (force23_expression.isConstant())
        force23_expression = force23_expression.evaluate(Values());

    if (force45_expression.isConstant())
        force45_expression = force45_expression.evaluate(Values());
    
    intra_molecule_points01 = Point::areIntraMoleculePoints(p[0], p[1]);
    intra_molecule_points23 = Point::areIntraMoleculePoints(p[2], p[3]);
    intra_molecule_points45 = Point::areIntraMoleculePoints(p[4], p[5]);

    this->calculateR();
}

/** Copy constructor */
TripleDistanceRestraint::TripleDistanceRestraint(const TripleDistanceRestraint &other)
    : ConcreteProperty<TripleDistanceRestraint,ExpressionRestraint3D>(other),
      force01_expression(other.force01_expression),
      force23_expression(other.force23_expression),
      force45_expression(other.force45_expression),
      intra_molecule_points01(other.intra_molecule_points01),
      intra_molecule_points23(other.intra_molecule_points23),
      intra_molecule_points45(other.intra_molecule_points45)
{
    for (int i=0; i<this->nPoints(); ++i)
    {
        p[i] = other.p[i];
    }
}

/** Destructor */
TripleDistanceRestraint::~TripleDistanceRestraint()
{}

/** Copy assignment operator */
TripleDistanceRestraint& TripleDistanceRestraint::operator=(
                                                const TripleDistanceRestraint &other)
{
    if (this != &other)
    {
        ExpressionRestraint3D::operator=(other);

        for (int i=0; i<this->nPoints(); ++i)
        {
            p[i] = other.p[i];
        }
        
        force01_expression = other.force01_expression;
        force23_expression = other.force23_expression;
        force45_expression = other.force45_expression;
        intra_molecule_points01 = other.intra_molecule_points01;
        intra_molecule_points23 = other.intra_molecule_points23;
        intra_molecule_points45 = other.intra_molecule_points45;
    }
    
    return *this;
}

/** Comparison operator */
bool TripleDistanceRestraint::operator==(const TripleDistanceRestraint &other) const
{
    return this == &other or
           ( ExpressionRestraint3D::operator==(other) and
             p[0] == other.p[0] and p[1] == other.p[1] and
             p[2] == other.p[2] and p[3] == other.p[3] and
             p[4] == other.p[4] and p[5] == other.p[5] and
             force01_expression == other.force01_expression and
             force23_expression == other.force23_expression and
             force45_expression == other.force45_expression);
}

/** Comparison operator */
bool TripleDistanceRestraint::operator!=(const TripleDistanceRestraint &other) const
{
    return not TripleDistanceRestraint::operator==(other);
}

const char* TripleDistanceRestraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TripleDistanceRestraint>() );
}

/** This restraint involves six points */
int TripleDistanceRestraint::nPoints() const
{
    return 6;
}

/** Return the ith point */
const Point& TripleDistanceRestraint::point(int i) const
{
    i = Index(i).map( this->nPoints() );

    return p[i].read();
}

/** Return the first point */
const Point& TripleDistanceRestraint::point0() const
{
    return p[0].read();
}

/** Return the second point */
const Point& TripleDistanceRestraint::point1() const
{
    return p[1].read();
}

/** Return the third point */
const Point& TripleDistanceRestraint::point2() const
{
    return p[2].read();
}

/** Return the fourth point */
const Point& TripleDistanceRestraint::point3() const
{
    return p[3].read();
}

/** Return the fifth point */
const Point& TripleDistanceRestraint::point4() const
{
    return p[4].read();
}

/** Return the sixth point */
const Point& TripleDistanceRestraint::point5() const
{
    return p[5].read();
}

/** Return the built-in symbols of this restraint */
Symbols TripleDistanceRestraint::builtinSymbols() const
{
    Symbols symbols;
    
    if (this->restraintFunction().isFunction(r01()))
        symbols += r01();
        
    if (this->restraintFunction().isFunction(r23()))
        symbols += r23();

    if (this->restraintFunction().isFunction(r45()))
        symbols += r45();

    return symbols;
}

/** Return the values of the built-in symbols of this restraint */
Values TripleDistanceRestraint::builtinValues() const
{
    Values vals;
    
    if (this->restraintFunction().isFunction(r01()))
    {
        vals.set( r01(), this->values()[r01()] );
    }
    
    if (this->restraintFunction().isFunction(r23()))
    {
        vals.set( r23(), this->values()[r23()] );
    }
    
    if (this->restraintFunction().isFunction(r45()))
    {
        vals.set( r45(), this->values()[r45()] );
    }

    return vals;
}

/** Return the differential of this restraint with respect to 
    the symbol 'symbol' */
RestraintPtr TripleDistanceRestraint::differentiate(const Symbol &symbol) const
{
    if (this->restraintFunction().isFunction(symbol))
    {
        return TripleDistanceRestraint( p[0], p[1], p[2], p[3], p[4], p[5],
                                        restraintFunction().differentiate(symbol),
                                        this->values() );
    }
    else
        return NullRestraint();
}

/** Set the space used to evaluate the energy of this restraint

    \throw SireVol::incompatible_space
*/
void TripleDistanceRestraint::setSpace(const Space &new_space)
{
    if (not this->space().equals(new_space))
    {
        TripleDistanceRestraint old_state(*this);
    
        try
        {
            for (int i=0; i<this->nPoints(); ++i)
            {
                p[i].edit().setSpace(new_space);
            }
            
            ExpressionRestraint3D::setSpace(new_space);
            
            this->calculateR();
        }
        catch(...)
        {
            TripleDistanceRestraint::operator=(old_state);
            throw;
        }
    }
}

/** Return the function used to calculate the restraint force along the 
    distance r01 */
const Expression& TripleDistanceRestraint::differentialRestraintFunction01() const
{
    return force01_expression;
}

/** Return the function used to calculate the restraint force along the 
    distance r23 */
const Expression& TripleDistanceRestraint::differentialRestraintFunction23() const
{
    return force23_expression;
}

/** Return the function used to calculate the restraint force along the 
    distance r45 */
const Expression& TripleDistanceRestraint::differentialRestraintFunction45() const
{
    return force45_expression;
}

/** Calculate the force acting on the molecule in the forcetable 'forcetable' 
    caused by this restraint, and add it on to the forcetable scaled by 
    'scale_force' */
void TripleDistanceRestraint::force(MolForceTable &forcetable, double scale_force) const
{
    bool in_p0 = p[0].read().contains(forcetable.molNum());
    bool in_p1 = p[1].read().contains(forcetable.molNum());
    bool in_p2 = p[2].read().contains(forcetable.molNum());
    bool in_p3 = p[3].read().contains(forcetable.molNum());
    bool in_p4 = p[4].read().contains(forcetable.molNum());
    bool in_p5 = p[5].read().contains(forcetable.molNum());

    if (in_p0 or in_p1)
    {
        const double force = scale_force * force01_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[0], p[1], space(), intra_molecule_points01,
                 in_p0, in_p1, force, forcetable);
    }
    
    if (in_p2 or in_p3)
    {
        const double force = scale_force * force23_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[2], p[3], space(), intra_molecule_points23,
                 in_p2, in_p3, force, forcetable);
    }
    
    if (in_p4 or in_p5)
    {
        const double force = scale_force * force45_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[4], p[5], space(), intra_molecule_points45,
                 in_p4, in_p5, force, forcetable);
    }
}

/** Calculate the force acting on the molecules in the forcetable 'forcetable' 
    caused by this restraint, and add it on to the forcetable scaled by 
    'scale_force' */
void TripleDistanceRestraint::force(ForceTable &forcetable, double scale_force) const
{
    bool in_p0 = p[0].read().usesMoleculesIn(forcetable);
    bool in_p1 = p[1].read().usesMoleculesIn(forcetable);
    bool in_p2 = p[2].read().usesMoleculesIn(forcetable);
    bool in_p3 = p[3].read().usesMoleculesIn(forcetable);
    bool in_p4 = p[4].read().usesMoleculesIn(forcetable);
    bool in_p5 = p[5].read().usesMoleculesIn(forcetable);

    if (in_p0 or in_p1)
    {
        const double force = scale_force * force01_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[0], p[1], space(), intra_molecule_points01,
                 in_p0, in_p1, force, forcetable);
    }
    
    if (in_p2 or in_p3)
    {
        const double force = scale_force * force23_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[2], p[3], space(), intra_molecule_points23,
                 in_p2, in_p3, force, forcetable);
    }
    
    if (in_p4 or in_p5)
    {
        const double force = scale_force * force45_expression.evaluate( 
                                                ExpressionRestraint3D::values() );

        addForce(p[4], p[5], space(), intra_molecule_points45,
                 in_p4, in_p5, force, forcetable);
    }
}

/** Update the points of this restraint using new molecule data from 'moldata'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void TripleDistanceRestraint::update(const MoleculeData &moldata)
{
    if (this->contains(moldata.number()))
    {
        TripleDistanceRestraint old_state(*this);
        
        try
        {
            for (int i=0; i<this->nPoints(); ++i)
            {
                p[i].edit().update(moldata);
            }
            
            this->calculateR();
        }
        catch(...)
        {
            TripleDistanceRestraint::operator=(old_state);
            throw;
        }
    }
}
            
/** Update the points of this restraint using new molecule data from 'molecules'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void TripleDistanceRestraint::update(const Molecules &molecules)
{
    if (this->usesMoleculesIn(molecules))
    {
        TripleDistanceRestraint old_state(*this);
        
        try
        {
            for (int i=0; i<this->nPoints(); ++i)
            {
                p[i].edit().update(molecules);
            }
            
            this->calculateR();
        }
        catch(...)
        {
            TripleDistanceRestraint::operator=(old_state);
            throw;
        }
    }
}

/** Return the molecules used in this restraint */
Molecules TripleDistanceRestraint::molecules() const
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
bool TripleDistanceRestraint::contains(MolNum molnum) const
{
    return p[0].read().contains(molnum) or p[1].read().contains(molnum) or
           p[2].read().contains(molnum) or p[2].read().contains(molnum) or
           p[4].read().contains(molnum) or p[5].read().contains(molnum);
}

/** Return whether or not this restraint affects the molecule
    with ID 'molid' */
bool TripleDistanceRestraint::contains(const MolID &molid) const
{
    return p[0].read().contains(molid) or p[1].read().contains(molid) or
           p[2].read().contains(molid) or p[3].read().contains(molid) or
           p[4].read().contains(molid) or p[5].read().contains(molid);
}
    
/** Return whether or not this restraint involves any of the molecules
    that are in the forcetable 'forcetable' */
bool TripleDistanceRestraint::usesMoleculesIn(const ForceTable &forcetable) const
{
    return p[0].read().usesMoleculesIn(forcetable) or 
           p[1].read().usesMoleculesIn(forcetable) or
           p[2].read().usesMoleculesIn(forcetable) or
           p[3].read().usesMoleculesIn(forcetable) or
           p[4].read().usesMoleculesIn(forcetable) or
           p[5].read().usesMoleculesIn(forcetable);
}

/** Return whether or not this restraint involves any of the molecules
    in 'molecules' */
bool TripleDistanceRestraint::usesMoleculesIn(const Molecules &molecules) const
{
    return p[0].read().usesMoleculesIn(molecules) or
           p[1].read().usesMoleculesIn(molecules) or
           p[2].read().usesMoleculesIn(molecules) or
           p[3].read().usesMoleculesIn(molecules) or
           p[4].read().usesMoleculesIn(molecules) or
           p[5].read().usesMoleculesIn(molecules);
}
