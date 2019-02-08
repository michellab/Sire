/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIREMOL_EVALUATOR_H
#define SIREMOL_EVALUATOR_H

#include <QVector>
#include <QHash>

#include "moleculeview.h"
#include "atomselection.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Evaluator;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Evaluator&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Evaluator&);

namespace SireMaths
{
class AxisSet;
class Sphere;
class Vector;
}

namespace SireVol
{
class AABox;
}

namespace SireMol
{

class AtomMatcher;
class AtomSelection;
class BondID;
class AngleID;
class DihedralID;

using SireBase::Property;
using SireBase::PropertyMap;
using SireBase::PropertyName;

using SireMaths::AxisSet;
using SireMaths::Sphere;
using SireMaths::Vector;

using SireVol::AABox;

/** This class is used to add a nice API to the MoleculeView based classes to
    allow the evaluation of various properties of the molecule (without the
    need to clutter up the MoleculeView-based classes' APIs).

    e.g. can type mol.evaluate().center() rather than mol.center()

    @author Christopher Woods
*/
class SIREMOL_EXPORT Evaluator
            : public SireBase::ConcreteProperty<Evaluator,MoleculeView>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Evaluator&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, Evaluator&);

public:
    Evaluator();

    Evaluator(const MoleculeView &molecule);

    Evaluator(const MoleculeData &moldata);

    Evaluator(const MoleculeView &molecule,
              const AtomSelection &selected_atoms);

    Evaluator(const MoleculeData &moldata,
              const AtomSelection &selected_atoms);

    Evaluator(const Evaluator &other);

    ~Evaluator();

    Evaluator& operator=(const Evaluator &other);
    Evaluator& operator=(const MoleculeView &other);

    static const char* typeName();

    Evaluator* clone() const;

    QString toString() const;

    bool isEmpty() const;
    bool selectedAll() const;

    bool hasProperty(const PropertyName&) const
    {
        return false;
    }

    bool hasMetadata(const PropertyName&) const
    {
        return false;
    }

    bool hasMetadata(const PropertyName&, const PropertyName&) const
    {
        return false;
    }

    QStringList propertyKeys() const
    {
        return QStringList();
    }

    QStringList metadataKeys() const
    {
        return QStringList();
    }

    QStringList metadataKeys(const PropertyName&) const
    {
        return QStringList();
    }

    AtomSelection selection() const;

    SireUnits::Dimension::MolarMass mass(const PropertyMap &map = PropertyMap()) const;

    SireUnits::Dimension::Charge charge(const PropertyMap &map = PropertyMap()) const;

    Vector center(const PropertyMap &map = PropertyMap()) const;
    AABox aaBox(const PropertyMap &map = PropertyMap()) const;
    Sphere boundingSphere(const PropertyMap &map = PropertyMap()) const;

    Vector centroid(const PropertyMap &map = PropertyMap()) const;
    Vector centerOfGeometry(const PropertyMap &map = PropertyMap()) const;
    Vector centerOfMass(const PropertyMap &map = PropertyMap()) const;

    AxisSet principalAxes(const PropertyMap &map = PropertyMap()) const;
    AxisSet principalAxes(Vector &principal_moments,
                          const PropertyMap &map = PropertyMap()) const;

    AxisSet alignmentAxes(const MoleculeView &other,
                          const AtomMatcher &matcher,
                          const PropertyMap &map = PropertyMap()) const;

    AxisSet alignmentAxes(const MoleculeView &other,
                          const AtomMatcher &matcher,
                          const PropertyMap &map0,
                          const PropertyMap &map1) const;

    SireUnits::Dimension::Length measure(const AtomID &atom0,
                                         const AtomID &atom1,
                                         const PropertyMap &map=PropertyMap()) const;

    SireUnits::Dimension::Length measure(const BondID &bond,
                                         const PropertyMap &map=PropertyMap()) const;

    SireUnits::Dimension::Angle measure(const AtomID &atom0, const AtomID &atom1,
                                        const AtomID &atom2,
                                        const PropertyMap &map=PropertyMap()) const;

    SireUnits::Dimension::Angle measure(const AngleID &angle,
                                        const PropertyMap &map=PropertyMap()) const;

    SireUnits::Dimension::Angle measure(const AtomID &atom0, const AtomID &atom1,
                                        const AtomID &atom2, const AtomID &atom3,
                                        const PropertyMap &map=PropertyMap()) const;

    SireUnits::Dimension::Angle measure(const DihedralID &dihedral,
                                        const PropertyMap &map=PropertyMap()) const;

    SireUnits::Dimension::Length rmsd(const MoleculeView &other,
                                      const PropertyMap &map=PropertyMap()) const;

    SireUnits::Dimension::Length rmsd(const MoleculeView &other,
                                      const PropertyMap &map0,
                                      const PropertyMap &map1) const;

    SireUnits::Dimension::Length rmsd(const MoleculeView &other,
                                      const AtomMatcher &atommatcher,
                                      const PropertyMap &map=PropertyMap()) const;

    SireUnits::Dimension::Length rmsd(const MoleculeView &other,
                                      const AtomMatcher &atommatcher,
                                      const PropertyMap &map0,
                                      const PropertyMap &map1) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const PropertyMap &map=PropertyMap(),
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const AtomMatcher &atommatcher,
                                   const PropertyMap &map=PropertyMap(),
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const PropertyMap &map0,
                                   const PropertyMap &map1,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const AtomMatcher &atommatcher,
                                   const PropertyMap &map0,
                                   const PropertyMap &map1,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const SireUnits::Dimension::Time &timeout,
                                   const PropertyMap &map=PropertyMap(),
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const AtomMatcher &atommatcher,
                                   const SireUnits::Dimension::Time &timeout,
                                   const PropertyMap &map=PropertyMap(),
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const SireUnits::Dimension::Time &timeout,
                                   const PropertyMap &map0,
                                   const PropertyMap &map1,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const AtomMatcher &atommatcher,
                                   const SireUnits::Dimension::Time &timeout,
                                   const PropertyMap &map0,
                                   const PropertyMap &map1,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   bool match_light_atoms,
                                   const PropertyMap &map=PropertyMap(),
                                   int min_heavy_protons=6,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const AtomMatcher &atommatcher,
                                   bool match_light_atoms,
                                   const PropertyMap &map=PropertyMap(),
                                   int min_heavy_protons=6,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   bool match_light_atoms,
                                   const PropertyMap &map0,
                                   const PropertyMap &map1,
                                   int min_heavy_protons=6,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const AtomMatcher &atommatcher,
                                   bool match_light_atoms,
                                   const PropertyMap &map0,
                                   const PropertyMap &map1,
                                   int min_heavy_protons=6,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const SireUnits::Dimension::Time &timeout,
                                   bool match_light_atoms,
                                   const PropertyMap &map=PropertyMap(),
                                   int min_heavy_protons=6,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const AtomMatcher &atommatcher,
                                   const SireUnits::Dimension::Time &timeout,
                                   bool match_light_atoms,
                                   const PropertyMap &map=PropertyMap(),
                                   int min_heavy_protons=6,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const SireUnits::Dimension::Time &timeout,
                                   bool match_light_atoms,
                                   const PropertyMap &map0,
                                   const PropertyMap &map1,
                                   int min_heavy_protons=6,
                                   bool verbose=false) const;

    QHash<AtomIdx,AtomIdx> findMCS(const MoleculeView &other,
                                   const AtomMatcher &atommatcher,
                                   const SireUnits::Dimension::Time &timeout,
                                   bool match_light_atoms,
                                   const PropertyMap &map0,
                                   const PropertyMap &map1,
                                   int min_heavy_protons=6,
                                   bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const PropertyMap &map=PropertyMap(),
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const AtomMatcher &atommatcher,
                                                    const PropertyMap &map=PropertyMap(),
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const PropertyMap &map0,
                                                    const PropertyMap &map1,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const AtomMatcher &atommatcher,
                                                    const PropertyMap &map0,
                                                    const PropertyMap &map1,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const SireUnits::Dimension::Time &timeout,
                                                    const PropertyMap &map=PropertyMap(),
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const AtomMatcher &atommatcher,
                                                    const SireUnits::Dimension::Time &timeout,
                                                    const PropertyMap &map=PropertyMap(),
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const SireUnits::Dimension::Time &timeout,
                                                    const PropertyMap &map0,
                                                    const PropertyMap &map1,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const AtomMatcher &atommatcher,
                                                    const SireUnits::Dimension::Time &timeout,
                                                    const PropertyMap &map0,
                                                    const PropertyMap &map1,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    bool match_light_atoms,
                                                    const PropertyMap &map=PropertyMap(),
                                                    int min_heavy_protons=6,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const AtomMatcher &atommatcher,
                                                    bool match_light_atoms,
                                                    const PropertyMap &map=PropertyMap(),
                                                    int min_heavy_protons=6,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    bool match_light_atoms,
                                                    const PropertyMap &map0,
                                                    const PropertyMap &map1,
                                                    int min_heavy_protons=6,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const AtomMatcher &atommatcher,
                                                    bool match_light_atoms,
                                                    const PropertyMap &map0,
                                                    const PropertyMap &map1,
                                                    int min_heavy_protons=6,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const SireUnits::Dimension::Time &timeout,
                                                    bool match_light_atoms,
                                                    const PropertyMap &map=PropertyMap(),
                                                    int min_heavy_protons=6,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const AtomMatcher &atommatcher,
                                                    const SireUnits::Dimension::Time &timeout,
                                                    bool match_light_atoms,
                                                    const PropertyMap &map=PropertyMap(),
                                                    int min_heavy_protons=6,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const SireUnits::Dimension::Time &timeout,
                                                    bool match_light_atoms,
                                                    const PropertyMap &map0,
                                                    const PropertyMap &map1,
                                                    int min_heavy_protons=6,
                                                    bool verbose=false) const;

    QVector<QHash<AtomIdx,AtomIdx> > findMCSmatches(const MoleculeView &other,
                                                    const AtomMatcher &atommatcher,
                                                    const SireUnits::Dimension::Time &timeout,
                                                    bool match_light_atoms,
                                                    const PropertyMap &map0,
                                                    const PropertyMap &map1,
                                                    int min_heavy_protons=6,
                                                    bool verbose=false) const;

private:

    /** The atoms over which the properties will be
        evaluated */
    AtomSelection selected_atoms;
};

}

Q_DECLARE_METATYPE(SireMol::Evaluator)

SIRE_EXPOSE_CLASS( SireMol::Evaluator )

SIRE_END_HEADER

#endif
