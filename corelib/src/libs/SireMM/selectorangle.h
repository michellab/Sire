/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREMM_SELECTORANGLE_H
#define SIREMM_SELECTORANGLE_H

#include "angle.h"

#include "SireMol/selector.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{
class SelectorAngle;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SelectorAngle&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SelectorAngle&);

namespace SireMM
{

/** This provides a Selector<T>-style interface for multiple angles */
class SIREMM_EXPORT SelectorAngle :
    public SireBase::ConcreteProperty<SelectorAngle, SireMol::MoleculeView>
{

friend QDataStream& ::operator<<(QDataStream&, const SelectorAngle&);
friend QDataStream& ::operator>>(QDataStream&, SelectorAngle&);

public:
    typedef QList<Angle>::const_iterator const_iterator;

    SelectorAngle();
    SelectorAngle(const SireMol::MoleculeData &molecule,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorAngle(const MoleculeView &molecule,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorAngle(const SireMol::MoleculeView &molecule,
                  const QList<SireMol::AngleID> &angles);

    SelectorAngle(const SireMol::MoleculeData &molecule,
                  const SireMol::AtomID &atom,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorAngle(const SireMol::MoleculeData &molecule,
                  const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorAngle(const SireMol::MoleculeData &molecule,
                  const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                  const SireMol::AtomID &atom2,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorAngle(const SireMol::MoleculeView &molecule,
                  const SireMol::AtomID &atom,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorAngle(const SireMol::MoleculeView &molecule,
                  const SireMol::AngleID &angle,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorAngle(const SireMol::MoleculeView &molecule,
                  const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorAngle(const SireMol::MoleculeView &molecule,
                  const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                  const SireMol::AtomID &atom2,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorAngle(const SireMol::Selector<SireMol::Atom> &atoms,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorAngle(const SireMol::Selector<SireMol::Atom> &atoms0,
                  const SireMol::Selector<SireMol::Atom> &atoms1,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorAngle(const SireMol::Selector<SireMol::Atom> &atoms0,
                  const SireMol::Selector<SireMol::Atom> &atoms1,
                  const SireMol::Selector<SireMol::Atom> &atoms2,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorAngle(const SelectorAngle &other);

    virtual ~SelectorAngle();

    static const char* typeName();

    virtual const char* what() const
    {
        return SelectorAngle::typeName();
    }

    virtual SelectorAngle* clone() const
    {
        return new SelectorAngle(*this);
    }

    SelectorAngle& operator=(const SelectorAngle &angle);

    bool operator==(const SelectorAngle &other) const;
    bool operator!=(const SelectorAngle &other) const;

    QString toString() const;

    SireMol::MolViewPtr operator[](int i) const;
    SireMol::MolViewPtr operator[](const SireBase::Slice &slice) const;
    SireMol::MolViewPtr operator[](const QList<qint64> &idxs) const;
    SireMol::MolViewPtr operator[](const SireMol::AngleID &angle) const;

    Angle operator()(int i) const;
    SelectorAngle operator()(int i, int j) const;
    SelectorAngle operator()(const SireBase::Slice &slice) const;
    SelectorAngle operator()(const QList<qint64> &idxs) const;
    SelectorAngle operator()(const SireMol::AngleID &angle) const;

    QList<SireMol::MolViewPtr> toList() const;

    QList<SireMol::AngleID> IDs() const;

    SelectorAngle add(const Angle &Angle) const;

    int count() const;
    int size() const;
    int nViews() const;

    bool isEmpty() const;
    bool selectedAll() const;

    SireMol::MolViewPtr toSelector() const;

    SireMol::AtomSelection selection() const;

    SelectorAngle add(const SelectorAngle &other) const;

    SelectorAngle intersection(const SelectorAngle &other) const;

    SelectorAngle invert(const SireBase::PropertyMap &map) const;
    SelectorAngle invert() const;

    bool hasProperty(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key,
                     const SireBase::PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const SireBase::PropertyName &key) const;

    QList<SireBase::Properties> properties() const;

    SireMol::Mover<SelectorAngle> move() const;
    SireMol::Evaluator evaluate() const;

    QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key) const;
    QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key,
                                          const Property &default_value) const;

    QList<SireUnits::Dimension::Angle> sizes() const;
    QList<SireUnits::Dimension::Angle> sizes(const SireBase::PropertyMap &map) const;

    QList<SireUnits::Dimension::Angle> measures() const;
    QList<SireUnits::Dimension::Angle> measures(const SireBase::PropertyMap &map) const;

    QList<SireCAS::Expression> potentials() const;
    QList<SireCAS::Expression> potentials(const SireBase::PropertyMap &map) const;

    QList<SireUnits::Dimension::MolarEnergy> energies() const;
    QList<SireUnits::Dimension::MolarEnergy> energies(
                            const SireBase::PropertyMap &map) const;

    SireUnits::Dimension::MolarEnergy energy() const;
    SireUnits::Dimension::MolarEnergy energy(
                            const SireBase::PropertyMap &map) const;

    SelectorAngle::const_iterator constBegin() const;
    SelectorAngle::const_iterator begin() const;

    SelectorAngle::const_iterator constEnd() const;
    SelectorAngle::const_iterator end() const;

protected:
    /** The IDs of the angles (holding AtomIdx IDs) */
    QList<SireMol::AngleID> angs;
};

} // end of namespace SireMM

Q_DECLARE_METATYPE( SireMM::SelectorAngle )
Q_DECLARE_METATYPE( SireMol::Mover<SireMM::SelectorAngle> )

SIRE_EXPOSE_CLASS( SireMM::SelectorAngle )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMM::SelectorAngle>, SireMol::Mover_SelectorAngle_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "SireMol/mover.hpp"

namespace SireMol
{
    template class SireMol::Mover<SireMM::SelectorAngle>;
}

#endif

SIRE_END_HEADER

#endif
