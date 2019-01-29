/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREMOL_BEADS_H
#define SIREMOL_BEADS_H

#include "bead.h"
#include "beadproperty.hpp"
#include "atomselection.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Beads;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream &ds, const SireMol::Beads&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream &ds, SireMol::Beads&);

namespace SireMol
{

/** This class is a view of all of the beads (for a specific "Beading")
    in a molecule
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT Beads : public SireBase::ConcreteProperty<Beads,MoleculeView>
{

friend QDataStream& ::operator<<(QDataStream&, const Beads&);
friend QDataStream& ::operator>>(QDataStream&, Beads&);

public:
    Beads();
    Beads(const MoleculeData &moldata,
          const PropertyMap &map = PropertyMap());
          
    Beads(const Beads &other);
    
    ~Beads();
    
    Beads& operator=(const Beads &other);
    
    bool operator==(const Beads &other) const;
    bool operator!=(const Beads &other) const;
    
    static const char* typeName();
    
    Beads* clone() const;
    
    Bead operator[](BeadIdx beadidx) const;
    
    Bead at(BeadIdx beadidx) const;
    
    Bead bead(BeadIdx beadidx) const;
    
    int count() const;
    int size() const;
    
    QString toString() const;
    
    bool isEmpty() const;
    bool selectedAll() const;
    
    AtomSelection selection() const;
    
    int nBeads() const;
    int nAtoms() const;
    
    void update(const MoleculeData &moldata);

    Mover<Beads> move() const;
    Evaluator evaluate() const;

    QList<AtomIdx> atomIdxs() const;
    
    const Beading& beading() const;
    
    bool contains(AtomIdx atomidx) const;
    bool contains(const AtomID &atomid) const;
    bool intersects(const AtomID &atomid) const;

    SireBase::PropertyPtr atomProperty(const PropertyName &key) const;

    bool hasProperty(const PropertyName &key) const;
    bool hasMetadata(const PropertyName &metakey) const;
    bool hasMetadata(const PropertyName &key,
                     const PropertyName &metakey) const;
                     
    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const PropertyName &key) const;

private:
    /** The beading used to divide the molecule into beads */
    BeadingPtr bdng;
    
    /** The location of the beading property */
    PropertyName beading_property;
    
    /** The atoms that are part of the beads */
    AtomSelection selected_atoms;
};

}

Q_DECLARE_METATYPE(SireMol::Beads);
Q_DECLARE_METATYPE(SireMol::Mover<SireMol::Beads>);

SIRE_EXPOSE_CLASS( SireMol::Beads )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Beads>, SireMol::Mover_Beads_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "mover.hpp"

template class SireMol::Mover<SireMol::Beads>;

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
