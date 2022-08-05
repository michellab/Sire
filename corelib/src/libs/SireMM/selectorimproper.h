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

#ifndef SIREMM_SELECTORIMPROPER_H
#define SIREMM_SELECTORIMPROPER_H

#include "improper.h"

#include "SireMol/selector.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{
class SelectorImproper;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SelectorImproper&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SelectorImproper&);

namespace SireMM
{

/** This provides a Selector<T>-style interface for multiple impropers */
class SIREMM_EXPORT SelectorImproper :
    public SireBase::ConcreteProperty<SelectorImproper, SireMol::MoleculeView>
{

friend QDataStream& ::operator<<(QDataStream&, const SelectorImproper&);
friend QDataStream& ::operator>>(QDataStream&, SelectorImproper&);

public:
    SelectorImproper();
    SelectorImproper(const SireMol::MoleculeData &molecule,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorImproper(const MoleculeView &molecule,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorImproper(const SireMol::MoleculeView &molecule,
                     const QList<SireMol::ImproperID> &impropers);

    SelectorImproper(const SireMol::MoleculeData &molecule,
                     const SireMol::AtomID &atom,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SireMol::MoleculeData &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SireMol::MoleculeData &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireMol::AtomID &atom2,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SireMol::MoleculeData &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireMol::AtomID &atom2, const SireMol::AtomID &atom3,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SireMol::MoleculeView &molecule,
                     const SireMol::AtomID &atom,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SireMol::MoleculeView &molecule,
                     const SireMol::ImproperID &improper,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SireMol::MoleculeView &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SireMol::MoleculeView &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireMol::AtomID &atom2,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SireMol::MoleculeView &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireMol::AtomID &atom2, const SireMol::AtomID &atom3,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SireMol::Selector<SireMol::Atom> &atoms,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorImproper(const SireMol::Selector<SireMol::Atom> &atoms0,
                     const SireMol::Selector<SireMol::Atom> &atoms1,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorImproper(const SireMol::Selector<SireMol::Atom> &atoms0,
                     const SireMol::Selector<SireMol::Atom> &atoms1,
                     const SireMol::Selector<SireMol::Atom> &atoms2,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorImproper(const SireMol::Selector<SireMol::Atom> &atoms0,
                     const SireMol::Selector<SireMol::Atom> &atoms1,
                     const SireMol::Selector<SireMol::Atom> &atoms2,
                     const SireMol::Selector<SireMol::Atom> &atoms3,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorImproper(const SelectorImproper &other);

    virtual ~SelectorImproper();

    static const char* typeName();

    virtual const char* what() const
    {
        return SelectorImproper::typeName();
    }

    virtual SelectorImproper* clone() const
    {
        return new SelectorImproper(*this);
    }

    SelectorImproper& operator=(const SelectorImproper &Improper);

    bool operator==(const SelectorImproper &other) const;
    bool operator!=(const SelectorImproper &other) const;

    QString toString() const;

    SireMol::MolViewPtr operator[](int i) const;
    SireMol::MolViewPtr operator[](const SireBase::Slice &slice) const;
    SireMol::MolViewPtr operator[](const QList<qint64> &idxs) const;
    SireMol::MolViewPtr operator[](const SireMol::ImproperID &Improper) const;

    Improper operator()(int i) const;
    SelectorImproper operator()(int i, int j) const;
    SelectorImproper operator()(const SireBase::Slice &slice) const;
    SelectorImproper operator()(const QList<qint64> &idxs) const;
    SelectorImproper operator()(const SireMol::ImproperID &Improper) const;

    QList<SireMol::MolViewPtr> toList() const;

    QList<SireMol::ImproperID> IDs() const;

    SelectorImproper add(const Improper &Improper) const;

    int count() const;
    int size() const;
    int nViews() const;

    bool isEmpty() const;
    bool selectedAll() const;

    SireMol::MolViewPtr toSelector() const;

    SireMol::AtomSelection selection() const;

    SelectorImproper add(const SelectorImproper &other) const;

    SelectorImproper intersection(const SelectorImproper &other) const;

    SelectorImproper invert(const SireBase::PropertyMap &map) const;
    SelectorImproper invert() const;

    bool hasProperty(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key,
                     const SireBase::PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const SireBase::PropertyName &key) const;

    QList<SireBase::Properties> properties() const;

    SireMol::Mover<SelectorImproper> move() const;
    SireMol::Evaluator evaluate() const;

    QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key) const;
    QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key,
                                          const Property &default_value) const;

    QList<SireUnits::Dimension::Angle> sizes() const;
    QList<SireUnits::Dimension::Angle> sizes(const SireBase::PropertyMap &map) const;

    QList<SireCAS::Expression> potentials() const;
    QList<SireCAS::Expression> potentials(const SireBase::PropertyMap &map) const;

    QList<SireUnits::Dimension::MolarEnergy> energies() const;
    QList<SireUnits::Dimension::MolarEnergy> energies(
                            const SireBase::PropertyMap &map) const;

    SireUnits::Dimension::MolarEnergy energy() const;
    SireUnits::Dimension::MolarEnergy energy(
                            const SireBase::PropertyMap &map) const;

protected:
    /** The IDs of the Impropers (holding AtomIdx IDs) */
    QList<SireMol::ImproperID> imps;
};

} // end of namespace SireMM

Q_DECLARE_METATYPE( SireMM::SelectorImproper )
Q_DECLARE_METATYPE( SireMol::Mover<SireMM::SelectorImproper> )

SIRE_EXPOSE_CLASS( SireMM::SelectorImproper )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMM::SelectorImproper>, SireMol::Mover_SelectorImproper_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "SireMol/mover.hpp"

namespace SireMol
{
    template class SireMol::Mover<SireMM::SelectorImproper>;
}

#endif

SIRE_END_HEADER

#endif
