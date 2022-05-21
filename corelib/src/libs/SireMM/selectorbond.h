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

#ifndef SIREMM_SELECTORBOND_H
#define SIREMM_SELECTORBOND_H

#include "bond.h"

#include "SireMol/selector.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{
class SelectorBond;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SelectorBond&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SelectorBond&);

namespace SireMM
{

/** This provides a Selector<T>-style interface for multiple bonds */
class SIREMM_EXPORT SelectorBond :
    public SireBase::ConcreteProperty<SelectorBond, SireMol::MoleculeView>
{

friend QDataStream& ::operator<<(QDataStream&, const SelectorBond&);
friend QDataStream& ::operator>>(QDataStream&, SelectorBond&);

public:
    SelectorBond();
    SelectorBond(const SireMol::MoleculeData &molecule,
                 const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorBond(const MoleculeView &molecule,
                 const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorBond(const SireMol::MoleculeView &molecule,
                 const QList<SireMol::BondID> &bonds);

    SelectorBond(const SireMol::MoleculeData &molecule,
                 const SireMol::AtomID &atom,
                 const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorBond(const SireMol::MoleculeData &molecule,
                 const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                 const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorBond(const SireMol::MoleculeView &molecule,
                 const SireMol::AtomID &atom,
                 const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorBond(const SireMol::MoleculeView &molecule,
                 const SireMol::BondID &atom,
                 const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorBond(const SireMol::MoleculeView &molecule,
                 const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                 const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorBond(const SireMol::Selector<SireMol::Atom> &atoms,
                 const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorBond(const SireMol::Selector<SireMol::Atom> &atoms0,
                 const SireMol::Selector<SireMol::Atom> &atoms1,
                 const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorBond(const SelectorBond &other);

    virtual ~SelectorBond();

    static const char* typeName();

    virtual const char* what() const
    {
        return SelectorBond::typeName();
    }

    virtual SelectorBond* clone() const
    {
        return new SelectorBond(*this);
    }

    SelectorBond& operator=(const SelectorBond &bond);

    bool operator==(const SelectorBond &other) const;
    bool operator!=(const SelectorBond &other) const;

    QString toString() const;

    SireMol::MolViewPtr operator[](int i) const;
    SireMol::MolViewPtr operator[](const SireBase::Slice &slice) const;
    SireMol::MolViewPtr operator[](const QList<qint64> &idxs) const;
    SireMol::MolViewPtr operator[](const SireMol::BondID &bond) const;

    Bond operator()(int i) const;
    SelectorBond operator()(int i, int j) const;
    SelectorBond operator()(const SireBase::Slice &slice) const;
    SelectorBond operator()(const QList<qint64> &idxs) const;
    SelectorBond operator()(const SireMol::BondID &bond) const;

    QList<SireMol::MolViewPtr> toList() const;

    QList<SireMol::BondID> IDs() const;

    SelectorBond add(const Bond &bond) const;

    int count() const;
    int size() const;
    int nViews() const;

    bool isEmpty() const;
    bool selectedAll() const;

    SireMol::AtomSelection selection() const;

    bool hasProperty(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key,
                     const SireBase::PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const SireBase::PropertyName &key) const;

    QList<SireBase::Properties> properties() const;

    SireMol::Mover<SelectorBond> move() const;
    SireMol::Evaluator evaluate() const;

    QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key) const;
    QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key,
                                          const Property &default_value) const;

    QList<SireUnits::Dimension::Length> lengths() const;
    QList<SireUnits::Dimension::Length> lengths(const SireBase::PropertyMap &map) const;

    QList<SireCAS::Expression> potentials() const;
    QList<SireCAS::Expression> potentials(const SireBase::PropertyMap &map) const;

    QList<SireUnits::Dimension::MolarEnergy> energies() const;
    QList<SireUnits::Dimension::MolarEnergy> energies(
                            const SireBase::PropertyMap &map) const;

    SireUnits::Dimension::MolarEnergy energy() const;
    SireUnits::Dimension::MolarEnergy energy(
                            const SireBase::PropertyMap &map) const;

protected:
    /** The IDs of the bond (holding AtomIdx IDs) */
    QList<SireMol::BondID> bnds;
};

} // end of namespace SireMM

Q_DECLARE_METATYPE( SireMM::SelectorBond )
Q_DECLARE_METATYPE( SireMol::Mover<SireMM::SelectorBond> )

SIRE_EXPOSE_CLASS( SireMM::SelectorBond )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMM::SelectorBond>, SireMol::Mover_SelectorBond_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "SireMol/mover.hpp"

namespace SireMol
{
    template class SireMol::Mover<SireMM::SelectorBond>;
}

#endif

SIRE_END_HEADER

#endif
