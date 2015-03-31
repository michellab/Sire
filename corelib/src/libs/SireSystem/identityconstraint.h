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

#ifndef SIRESYSTEM_IDENTITYCONSTRAINT_H
#define SIRESYSTEM_IDENTITYCONSTRAINT_H

#include <QVector>

#include "moleculeconstraint.h"

#include "SireBase/propertymap.h"

#include "SireMaths/vector.h"

#include "SireMol/moleculegroup.h"

#include "SireFF/point.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class IdentityConstraint;
}

QDataStream& operator<<(QDataStream&, const SireSystem::IdentityConstraint&);
QDataStream& operator>>(QDataStream&, SireSystem::IdentityConstraint&);

namespace SireSystem
{

using SireBase::PropertyMap;

using SireMaths::Vector;

using SireMol::MoleculeGroup;
using SireMol::Molecules;
using SireMol::MolID;
using SireMol::MolNum;

namespace detail
{
class IdentityConstraintPvt;
}

/** An identity constraint provides a method of constraining
    the identity of molecules based on *where* they are located.
    
    For example, it can be useful to be able to identify a 
    water in a binding pocket. However, in a normal simulation,
    we don't identify waters by location, but by their index 
    (e.g. this is the first water, this is the second etc.).
    
    This means that the identity of the water in the binding
    pocket can change, e.g. it can start with the fifth water
    in the pocket, but during the simulation the fifth water
    may diffuse out of the pocket, and the twentieth water
    would diffuse in its place. The identity of the water
    in the binding pocket will thus have changed from the 
    fifth water to the twentieth water.
    
    An identity constraint works by constantly monitoring 
    the locations of the waters, and so it can detect when 
    the fifth water is displaced by the twentieth water. When
    it detects that this has occured, the constraint swaps
    the coordinates of the fifth and twentieth waters,
    thereby ensuring that the fifth water *stays* in the pocket.
    This doesn't affect the energy or the statistics of the 
    system, as waters are indistinguishable (there are N!
    equivalent configurations of N waters - we only see them
    as N! different configurations as we identify each water,
    when really they are indistinguishable).
    
    The idea of constraining the identity of molecules
    was first presented by M. Tyka, R. Sessions and A. Clarke in 
   
    "Absolute Free-Energy Calculations of Liquids Using 
                a Harmonic Reference State" 
                
    J. Chem. Phys. B,  2007, 111, 9571-9580
    
            doi://10.1021/jp072357w
            
    They used the method to constrain the identity of *all*
    molecules, so that harmonic restraints can be applied
    to them all.
    
    This identity constraint is more general, and allows
    the identification of a subset of molecules to be
    constrained (e.g. just the waters in a binding pocket).
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT IdentityConstraint
               : public SireBase::ConcreteProperty<IdentityConstraint,MoleculeConstraint>
{

friend QDataStream& ::operator<<(QDataStream&, const IdentityConstraint&);
friend QDataStream& ::operator>>(QDataStream&, IdentityConstraint&);

public:
    IdentityConstraint();
    
    IdentityConstraint(const MoleculeGroup &molgroup,
                       const PropertyMap &map = PropertyMap());
    
    IdentityConstraint(const SireFF::PointRef &point,
                       const MoleculeGroup &molgroup,
                       const PropertyMap &map = PropertyMap());
                       
    IdentityConstraint(const QList<SireFF::PointPtr> &points,
                       const MoleculeGroup &molgroup,
                       const PropertyMap &map = PropertyMap());
    
    IdentityConstraint(const QVector<SireFF::PointPtr> &points,
                       const MoleculeGroup &molgroup,
                       const PropertyMap &map = PropertyMap());
    
    IdentityConstraint(const IdentityConstraint &other);
    
    ~IdentityConstraint();
    
    IdentityConstraint& operator=(const IdentityConstraint &other);
    
    bool operator==(const IdentityConstraint &other) const;
    bool operator!=(const IdentityConstraint &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    const MoleculeGroup& moleculeGroup() const;
    
    QVector<SireFF::PointPtr> points() const;

    const PropertyMap& propertyMap() const;
    
    void useManyPointsAlgorithm();
    void useFewPointsAlgorithm();
    void useSinglePointAlgorithm();

    static SireMol::MolGroupPtr constrain(const MoleculeGroup &molgroup,
                                          const SireFF::PointRef &point,
                                          const PropertyMap &map = PropertyMap());
                                   
    static SireMol::MolGroupPtr constrain(const MoleculeGroup &molgroup,
                                          const QVector<SireFF::PointPtr> &points,
                                          const PropertyMap &map = PropertyMap());

    static SireMol::MolGroupPtr constrain(const MoleculeGroup &molgroup,
                                          const QList<SireFF::PointPtr> &points,
                                          const PropertyMap &map = PropertyMap());

protected:
    void setSystem(const System &system);
    bool mayChange(const Delta &delta, quint32 last_subversion) const;

    bool fullApply(Delta &delta);
    bool deltaApply(Delta &delta, quint32 last_subversion);

private:
    /** The helper class that is used to implement the constraint */
    SireBase::SharedPolyPointer<detail::IdentityConstraintPvt> d;

    /** The space property used by the identity points */
    SireBase::PropertyName space_property;

    /** The molecules that need to change to maintain the constraint */
    Molecules changed_mols;
};

}

Q_DECLARE_METATYPE( SireSystem::IdentityConstraint )

SIRE_EXPOSE_CLASS( SireSystem::IdentityConstraint )

SIRE_END_HEADER

#endif
