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

#ifndef SIREMOL_PARTIALMOLECULE_H
#define SIREMOL_PARTIALMOLECULE_H

#include "moleculeview.h"
#include "atomselection.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class PartialMolecule;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::PartialMolecule&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::PartialMolecule&);

namespace SireMol
{

template<class T>
class Mover;

template<class T>
class Selector;

class Evaluator;

class Molecule;
class Atom;

/** This class provides a view to an arbitrary part of a molecule
    (which can range from just a single atom all the way through to
    the entire molecule). As such, this class can be used to
    represent Molecule, Residue and Atom, as well as everything
    in between!

    @author Christopher Woods
*/
class SIREMOL_EXPORT PartialMolecule
        : public SireBase::ConcreteProperty<PartialMolecule,MoleculeView>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const PartialMolecule&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, PartialMolecule&);

public:
    PartialMolecule();
    PartialMolecule(const MoleculeView &molecule);

    PartialMolecule(const MoleculeData &moldata,
                    const AtomSelection &atoms);

    PartialMolecule(const PartialMolecule &other);

    ~PartialMolecule();

    static const char* typeName();

    PartialMolecule* clone() const;

    PartialMolecule& operator=(const MoleculeView &other);
    PartialMolecule& operator=(const PartialMolecule &other);

    bool operator==(const PartialMolecule &other) const;
    bool operator!=(const PartialMolecule &other) const;

    QString toString() const;

    bool isEmpty() const;
    bool selectedAll() const;

    MolViewPtr toSelector() const;

    const MolName& name() const;
    MolNum number() const;

    quint64 version() const;
    quint64 version(const PropertyName &key) const;

    int nAtoms() const;
    int nCutGroups() const;
    int nResidues() const;
    int nChains() const;
    int nSegments() const;

    Mover<PartialMolecule> move() const;

    Evaluator evaluate() const;

    AtomSelection selection() const;

    PartialMolecule extract() const;

    bool hasProperty(const PropertyName &key) const;

    bool hasMetadata(const PropertyName &metakey) const;

    bool hasMetadata(const PropertyName &key,
                     const PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const PropertyName &key) const;

    const Property& property(const PropertyName &key) const;

    const Property& metadata(const PropertyName &metakey) const;
    const Property& metadata(const PropertyName &key,
                             const PropertyName &metakey) const;

    MolViewPtr toUnit() const;

private:
    /** The atoms that are selected in this view */
    AtomSelection selected_atoms;
};

}

Q_DECLARE_METATYPE(SireMol::PartialMolecule);
Q_DECLARE_METATYPE(SireMol::Mover<SireMol::PartialMolecule>);

SIRE_EXPOSE_CLASS( SireMol::PartialMolecule )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::PartialMolecule>,
                   SireMol::Mover_PartialMolecule_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "mover.hpp"
#include "selector.hpp"

template class SireMol::Mover<SireMol::PartialMolecule>;

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
