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

#ifndef SIREMM_ANGLERESTRAINT_H
#define SIREMM_ANGLERESTRAINT_H

#include "SireFF/point.h"

#include "restraint.h"

#include "SireCAS/expression.h"
#include "SireCAS/symbol.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class AngleRestraint;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::AngleRestraint&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::AngleRestraint&);

namespace SireMM
{

using SireFF::Point;
using SireFF::PointRef;

using SireCAS::Expression;
using SireCAS::Symbol;
using SireCAS::Symbols;

// typedef the unit of a harmonic force constant ( MolarEnergy / Angle^2 )
typedef SireUnits::Dimension::PhysUnit<1,2,-2,0,0,-1,-2> HarmonicAngleForceConstant;

/** This is a restraint that operates on the angle between 
    three SireMM::Point objects (e.g. three atoms in a molecule)
    
    @author Christopher Woods
*/
class SIREMM_EXPORT AngleRestraint 
            : public SireBase::ConcreteProperty<AngleRestraint,ExpressionRestraint3D>
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const AngleRestraint&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, AngleRestraint&);

public:
    AngleRestraint();
    
    AngleRestraint(const PointRef &point0, const PointRef &point1,
                   const PointRef &point2, const Expression &restraint);

    
    AngleRestraint(const PointRef &point0, const PointRef &point1,
                   const PointRef &point2, const Expression &restraint,
                   const Values &values);

    AngleRestraint(const AngleRestraint &other);
    
    ~AngleRestraint();
    
    AngleRestraint& operator=(const AngleRestraint &other);
    
    bool operator==(const AngleRestraint &other) const;
    bool operator!=(const AngleRestraint &other) const;
    
    static const char* typeName();
    
    const Point& point(int i) const;
    
    const Point& point0() const;
    const Point& point1() const;
    const Point& point2() const;
    
    int nPoints() const;
    
    static const Symbol& theta();

    Symbols builtinSymbols() const;
    Values builtinValues() const;
    
    RestraintPtr differentiate(const Symbol &symbol) const;

    void setSpace(const Space &space);

    const Expression& differentialRestraintFunction() const;

    void force(MolForceTable &forcetable, double scale_force=1) const;
    void force(ForceTable &forcetable, double scale_force=1) const;

    void update(const MoleculeData &moldata);
                
    void update(const Molecules &molecules);

    Molecules molecules() const;
    
    bool contains(MolNum molnum) const;
    bool contains(const MolID &molid) const;
    
    bool usesMoleculesIn(const ForceTable &forcetable) const;
    bool usesMoleculesIn(const Molecules &molecules) const;

    static AngleRestraint harmonic(const PointRef &point0, const PointRef &point1,
                                   const PointRef &point2,
                                   const HarmonicAngleForceConstant &force_constant);
                                       
    static AngleRestraint halfHarmonic(const PointRef &point0, const PointRef &point1,
                                       const PointRef &point2,
                                       const SireUnits::Dimension::Angle &angle,
                                       const HarmonicAngleForceConstant &force_constant);

protected:
    AngleRestraint(const PointRef &point0, const PointRef &point1,
                   const PointRef &point2,
                   const Expression &nrg_restraint,
                   const Expression &force_restraint);

private:
    void calculateTheta();

    /** The three points between which the restraint is calculated */
    SireFF::PointPtr p[3];
    
    /** The expression used to calculate the force */
    Expression force_expression;
    
    /** Whether or not all three points are within the same molecule */
    bool intra_molecule_points;
};

}

Q_DECLARE_METATYPE( SireMM::AngleRestraint )

SIRE_EXPOSE_CLASS( SireMM::AngleRestraint )

SIRE_END_HEADER

#endif
