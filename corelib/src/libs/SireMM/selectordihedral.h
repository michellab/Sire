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

#ifndef SIREMM_SELECTORDIHEDRAL_H
#define SIREMM_SELECTORDIHEDRAL_H

#include "dihedral.h"

#include "SireMol/selector.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{
class SelectorDihedral;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SelectorDihedral&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SelectorDihedral&);

namespace SireMM
{

/** This provides a Selector<T>-style interface for multiple dihedrals */
class SIREMM_EXPORT SelectorDihedral :
    public SireBase::ConcreteProperty<SelectorDihedral, SireMol::MoleculeView>
{

friend QDataStream& ::operator<<(QDataStream&, const SelectorDihedral&);
friend QDataStream& ::operator>>(QDataStream&, SelectorDihedral&);

public:
    SelectorDihedral();
    SelectorDihedral(const SireMol::MoleculeData &molecule,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorDihedral(const MoleculeView &molecule,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorDihedral(const SireMol::MoleculeView &molecule,
                     const QList<SireMol::DihedralID> &dihedrals);

    SelectorDihedral(const SireMol::MoleculeData &molecule,
                     const SireMol::AtomID &atom,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SireMol::MoleculeData &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SireMol::MoleculeData &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireMol::AtomID &atom2,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SireMol::MoleculeData &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireMol::AtomID &atom2, const SireMol::AtomID &atom3,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SireMol::MoleculeView &molecule,
                     const SireMol::AtomID &atom,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SireMol::MoleculeView &molecule,
                     const SireMol::DihedralID &dihedral,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SireMol::MoleculeView &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SireMol::MoleculeView &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireMol::AtomID &atom2,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SireMol::MoleculeView &molecule,
                     const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                     const SireMol::AtomID &atom2, const SireMol::AtomID &atom3,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SireMol::Selector<SireMol::Atom> &atoms,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorDihedral(const SireMol::Selector<SireMol::Atom> &atoms0,
                     const SireMol::Selector<SireMol::Atom> &atoms1,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorDihedral(const SireMol::Selector<SireMol::Atom> &atoms0,
                     const SireMol::Selector<SireMol::Atom> &atoms1,
                     const SireMol::Selector<SireMol::Atom> &atoms2,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorDihedral(const SireMol::Selector<SireMol::Atom> &atoms0,
                     const SireMol::Selector<SireMol::Atom> &atoms1,
                     const SireMol::Selector<SireMol::Atom> &atoms2,
                     const SireMol::Selector<SireMol::Atom> &atoms3,
                     const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorDihedral(const SelectorDihedral &other);

    virtual ~SelectorDihedral();

    static const char* typeName();

    virtual const char* what() const
    {
        return SelectorDihedral::typeName();
    }

    virtual SelectorDihedral* clone() const
    {
        return new SelectorDihedral(*this);
    }

    SelectorDihedral& operator=(const SelectorDihedral &dihedral);

    bool operator==(const SelectorDihedral &other) const;
    bool operator!=(const SelectorDihedral &other) const;

    QString toString() const;

    SireMol::MolViewPtr operator[](int i) const;
    SireMol::MolViewPtr operator[](const SireBase::Slice &slice) const;
    SireMol::MolViewPtr operator[](const QList<qint64> &idxs) const;
    SireMol::MolViewPtr operator[](const SireMol::DihedralID &dihedral) const;

    Dihedral operator()(int i) const;
    SelectorDihedral operator()(int i, int j) const;
    SelectorDihedral operator()(const SireBase::Slice &slice) const;
    SelectorDihedral operator()(const QList<qint64> &idxs) const;
    SelectorDihedral operator()(const SireMol::DihedralID &dihedral) const;

    QList<SireMol::MolViewPtr> toList() const;

    QList<SireMol::DihedralID> IDs() const;

    SelectorDihedral add(const Dihedral &dihedral) const;

    int count() const;
    int size() const;
    int nViews() const;

    bool isEmpty() const;
    bool selectedAll() const;

    SireMol::MolViewPtr toSelector() const;

    SireMol::AtomSelection selection() const;

    SelectorDihedral add(const SelectorDihedral &other) const;

    SelectorDihedral intersection(const SelectorDihedral &other) const;

    SelectorDihedral invert(const SireBase::PropertyMap &map) const;
    SelectorDihedral invert() const;

    bool hasProperty(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key,
                     const SireBase::PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const SireBase::PropertyName &key) const;

    QList<SireBase::Properties> properties() const;

    SireMol::Mover<SelectorDihedral> move() const;
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
    /** The IDs of the Dihedrals (holding AtomIdx IDs) */
    QList<SireMol::DihedralID> dihs;
};

} // end of namespace SireMM

Q_DECLARE_METATYPE( SireMM::SelectorDihedral )
Q_DECLARE_METATYPE( SireMol::Mover<SireMM::SelectorDihedral> )

SIRE_EXPOSE_CLASS( SireMM::SelectorDihedral )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMM::SelectorDihedral>, SireMol::Mover_SelectorDihedral_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "SireMol/mover.hpp"

namespace SireMol
{
    template class SireMol::Mover<SireMM::SelectorDihedral>;
}

#endif

SIRE_END_HEADER

#endif
