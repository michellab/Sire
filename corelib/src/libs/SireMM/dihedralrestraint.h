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

#ifndef SIREMM_DIHEDRALRESTRAINT_H
#define SIREMM_DIHEDRALRESTRAINT_H

#include "restraint.h"
#include "anglerestraint.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class DihedralRestraint;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::DihedralRestraint&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::DihedralRestraint&);

namespace SireMM
{

/** This is a restraint that operates on the dihedral angle between 
    four SireMM::Point objects (e.g. four atoms in a molecule)
    
    @author Christopher Woods
*/
class SIREMM_EXPORT DihedralRestraint 
            : public SireBase::ConcreteProperty<DihedralRestraint,ExpressionRestraint3D>
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const DihedralRestraint&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, DihedralRestraint&);

public:
    DihedralRestraint();
    
    DihedralRestraint(const PointRef &point0, const PointRef &point1,
                      const PointRef &point2, const PointRef &point3,
                      const Expression &restraint);
    
    DihedralRestraint(const PointRef &point0, const PointRef &point1,
                      const PointRef &point2, const PointRef &point3,
                      const Expression &restraint, const Values &values);

    DihedralRestraint(const DihedralRestraint &other);
    
    ~DihedralRestraint();
    
    DihedralRestraint& operator=(const DihedralRestraint &other);
    
    bool operator==(const DihedralRestraint &other) const;
    bool operator!=(const DihedralRestraint &other) const;
    
    static const char* typeName();
    
    const Point& point(int i) const;
    
    const Point& point0() const;
    const Point& point1() const;
    const Point& point2() const;
    const Point& point3() const;
    
    int nPoints() const;
    
    static const Symbol& phi();

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

    static DihedralRestraint harmonic(
                                   const PointRef &point0, const PointRef &point1,
                                   const PointRef &point2, const PointRef &point3,
                                   const HarmonicAngleForceConstant &force_constant);
                                       
    static DihedralRestraint halfHarmonic(
                                       const PointRef &point0, const PointRef &point1,
                                       const PointRef &point2, const PointRef &point3,
                                       const SireUnits::Dimension::Angle &angle,
                                       const HarmonicAngleForceConstant &force_constant);

protected:
    DihedralRestraint(const PointRef &point0, const PointRef &point1,
                      const PointRef &point2, const PointRef &point3,
                      const Expression &nrg_restraint,
                      const Expression &force_restraint);

private:
    void calculatePhi();

    /** The four points between which the restraint is calculated */
    SireFF::PointPtr p[4];
    
    /** The expression used to calculate the force */
    Expression force_expression;
    
    /** Whether or not all four points are within the same molecule */
    bool intra_molecule_points;
};

}

Q_DECLARE_METATYPE( SireMM::DihedralRestraint )

SIRE_EXPOSE_CLASS( SireMM::DihedralRestraint )

SIRE_END_HEADER

#endif
