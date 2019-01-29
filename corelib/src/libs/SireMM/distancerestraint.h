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

#ifndef SIREMM_DISTANCERESTRAINT_H
#define SIREMM_DISTANCERESTRAINT_H

#include "SireFF/point.h"

#include "restraint.h"

#include "SireCAS/expression.h"
#include "SireCAS/symbol.h"
#include "SireCAS/values.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class DistanceRestraint;
class DoubleDistanceRestraint;
class TripleDistanceRestraint;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::DistanceRestraint&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::DistanceRestraint&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::DoubleDistanceRestraint&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::DoubleDistanceRestraint&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::TripleDistanceRestraint&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::TripleDistanceRestraint&);

namespace SireMM
{

using SireFF::Point;
using SireFF::PointRef;

using SireCAS::Expression;
using SireCAS::Symbol;
using SireCAS::Symbols;

// typedef the unit of a harmonic force constant ( MolarEnergy / Length^2 )
typedef SireUnits::Dimension::PhysUnit<1,0,-2,0,0,-1,0> HarmonicDistanceForceConstant;

/** This is a restraint that operates on the distance between 
    two SireMM::Point objects (e.g. two atoms in a molecule,
    a point in space and the center of a molecule, the
    center of geometry of one molecule with the center of 
    mass of the solvent)
    
    @author Christopher Woods
*/
class SIREMM_EXPORT DistanceRestraint 
        : public SireBase::ConcreteProperty<DistanceRestraint,ExpressionRestraint3D>
{

friend QDataStream& ::operator<<(QDataStream&, const DistanceRestraint&);
friend QDataStream& ::operator>>(QDataStream&, DistanceRestraint&);

public:
    DistanceRestraint();
    
    DistanceRestraint(const PointRef &point0, const PointRef &point1,
                      const Expression &restraint);
    
    DistanceRestraint(const PointRef &point0, const PointRef &point1,
                      const Expression &restraint, const Values &values);

    DistanceRestraint(const DistanceRestraint &other);
    
    ~DistanceRestraint();
    
    DistanceRestraint& operator=(const DistanceRestraint &other);
    
    bool operator==(const DistanceRestraint &other) const;
    bool operator!=(const DistanceRestraint &other) const;
    
    static const char* typeName();
    
    const Point& point(int i) const;
    
    const Point& point0() const;
    const Point& point1() const;
    
    int nPoints() const;
    
    static const Symbol& r();

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

    static DistanceRestraint
                    harmonic(const PointRef &point0,
                             const PointRef &point1,
                             const HarmonicDistanceForceConstant &force_constant);
                                      
    static DistanceRestraint 
                    halfHarmonic(const PointRef &point0,
                                 const PointRef &point1,
                                 const SireUnits::Dimension::Length &distance,
                                 const HarmonicDistanceForceConstant &force_constant);

protected:
    DistanceRestraint(const PointRef &point0, const PointRef &point1,
                      const Expression &nrg_restraint,
                      const Expression &force_restraint);

private:
    void calculateR();

    /** The two points between which the restraint is calculated */
    SireFF::PointPtr p[2];
    
    /** The expression used to calculate the force */
    Expression force_expression;
    
    /** Whether or not these two points are both within the same molecule */
    bool intra_molecule_points;
};

/** This class provides a restraint that operates on a pair
    of inter-point distances (e.g. the difference between
    two bond lengths). You need to supply four points to this
    restraint which are used to calculate the two distances
    (r01 between points 0 and 1, and r23 between points 2 and 3)
    
    @author Christopher Woods
*/
class SIREMM_EXPORT DoubleDistanceRestraint 
    : public SireBase::ConcreteProperty<DoubleDistanceRestraint,ExpressionRestraint3D>
{

friend QDataStream& ::operator<<(QDataStream&, const DoubleDistanceRestraint&);
friend QDataStream& ::operator>>(QDataStream&, DoubleDistanceRestraint&);

public:
    DoubleDistanceRestraint();
    
    DoubleDistanceRestraint(const PointRef &point0, const PointRef &point1,
                            const PointRef &point2, const PointRef &point3,
                            const Expression &restraint);

    DoubleDistanceRestraint(const PointRef &point0, const PointRef &point1,
                            const PointRef &point2, const PointRef &point3,
                            const Expression &restraint, const Values &values);

    DoubleDistanceRestraint(const DoubleDistanceRestraint &other);
    
    ~DoubleDistanceRestraint();
    
    DoubleDistanceRestraint& operator=(const DoubleDistanceRestraint &other);
    
    bool operator==(const DoubleDistanceRestraint &other) const;
    bool operator!=(const DoubleDistanceRestraint &other) const;
    
    static const char* typeName();
    
    int nPoints() const;
    
    const Point& point(int i) const;
    
    const Point& point0() const;
    const Point& point1() const;
    const Point& point2() const;
    const Point& point3() const;
    
    static const Symbol& r01();
    static const Symbol& r23();

    Symbols builtinSymbols() const;
    Values builtinValues() const;

    RestraintPtr differentiate(const Symbol &symbol) const;

    void setSpace(const Space &space);

    const Expression& differentialRestraintFunction01() const;
    const Expression& differentialRestraintFunction23() const;

    void force(MolForceTable &forcetable, double scale_force=1) const;
    void force(ForceTable &forcetable, double scale_force=1) const;

    void update(const MoleculeData &moldata);
                
    void update(const Molecules &molecules);

    Molecules molecules() const;
    
    bool contains(MolNum molnum) const;
    bool contains(const MolID &molid) const;
    
    bool usesMoleculesIn(const ForceTable &forcetable) const;
    bool usesMoleculesIn(const Molecules &molecules) const;

private:
    void calculateR();

    /** The four points between which the restraint is calculated */
    SireFF::PointPtr p[4];
    
    /** The expression used to calculate the force between points 0 and 1 */
    Expression force01_expression;
    
    /** The expression used to calculate the force between points 2 and 3 */
    Expression force23_expression;
    
    /** Whether or not points 0 and 1 are both within the same molecule */
    bool intra_molecule_points01;
    
    /** Whether or not points 2 and 3 are both within the same molecule */    
    bool intra_molecule_points23;
};

/** This class provides a restraint that operates on a triple
    of inter-point distances (e.g. the differences between
    three bond lengths). You need to supply six points to this
    restraint which are used to calculate the three distances
    (r01 between points 0 and 1, r23 between points 2 and 3 
     and r45 between points 4 and 5)
    
    @author Christopher Woods
*/
class SIREMM_EXPORT TripleDistanceRestraint 
       : public SireBase::ConcreteProperty<TripleDistanceRestraint,ExpressionRestraint3D>
{

friend QDataStream& ::operator<<(QDataStream&, const TripleDistanceRestraint&);
friend QDataStream& ::operator>>(QDataStream&, TripleDistanceRestraint&);

public:
    TripleDistanceRestraint();
    
    TripleDistanceRestraint(const PointRef &point0, const PointRef &point1,
                            const PointRef &point2, const PointRef &point3,
                            const PointRef &point4, const PointRef &point5,
                            const Expression &restraint);
    
    TripleDistanceRestraint(const PointRef &point0, const PointRef &point1,
                            const PointRef &point2, const PointRef &point3,
                            const PointRef &point4, const PointRef &point5,
                            const Expression &restraint, const Values &values);

    TripleDistanceRestraint(const TripleDistanceRestraint &other);
    
    ~TripleDistanceRestraint();
    
    TripleDistanceRestraint& operator=(const TripleDistanceRestraint &other);
    
    bool operator==(const TripleDistanceRestraint &other) const;
    bool operator!=(const TripleDistanceRestraint &other) const;
    
    static const char* typeName();
    
    int nPoints() const;
    
    const Point& point(int i) const;
    
    const Point& point0() const;
    const Point& point1() const;
    const Point& point2() const;
    const Point& point3() const;
    const Point& point4() const;
    const Point& point5() const;
    
    static const Symbol& r01();
    static const Symbol& r23();
    static const Symbol& r45();

    Symbols builtinSymbols() const;
    Values builtinValues() const;

    RestraintPtr differentiate(const Symbol &symbol) const;

    void setSpace(const Space &space);

    const Expression& differentialRestraintFunction01() const;
    const Expression& differentialRestraintFunction23() const;
    const Expression& differentialRestraintFunction45() const;

    void force(MolForceTable &forcetable, double scale_force=1) const;
    void force(ForceTable &forcetable, double scale_force=1) const;

    void update(const MoleculeData &moldata);
    void update(const Molecules &molecules);

    Molecules molecules() const;
    
    bool contains(MolNum molnum) const;
    bool contains(const MolID &molid) const;
    
    bool usesMoleculesIn(const ForceTable &forcetable) const;
    bool usesMoleculesIn(const Molecules &molecules) const;

private:
    void calculateR();

    /** The six points between which the restraint is calculated */
    SireFF::PointPtr p[6];
    
    /** The expression used to calculate the force between points 0 and 1 */
    Expression force01_expression;
    
    /** The expression used to calculate the force between points 2 and 3 */
    Expression force23_expression;
    
    /** The expression used to calculate the force between points 4 and 5 */
    Expression force45_expression;
    
    /** Whether or not points 0 and 1 are both within the same molecule */
    bool intra_molecule_points01;
    
    /** Whether or not points 2 and 3 are both within the same molecule */    
    bool intra_molecule_points23;
    
    /** Whether or not points 4 and 5 are both within the same molecule */    
    bool intra_molecule_points45;
};

}

Q_DECLARE_METATYPE( SireMM::DistanceRestraint )
Q_DECLARE_METATYPE( SireMM::DoubleDistanceRestraint )
Q_DECLARE_METATYPE( SireMM::TripleDistanceRestraint )

Q_DECLARE_METATYPE( SireMM::HarmonicDistanceForceConstant )

SIRE_EXPOSE_CLASS( SireMM::DistanceRestraint )
SIRE_EXPOSE_CLASS( SireMM::DoubleDistanceRestraint )
SIRE_EXPOSE_CLASS( SireMM::TripleDistanceRestraint )

SIRE_END_HEADER

#endif
