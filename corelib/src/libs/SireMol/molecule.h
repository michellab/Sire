/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMOL_MOLECULE_H
#define SIREMOL_MOLECULE_H

#include "moleculeview.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Molecule;
}

QDataStream& operator<<(QDataStream&, const SireMol::Molecule&);
QDataStream& operator>>(QDataStream&, SireMol::Molecule&);

namespace SireMol
{

class MolID;
class MolName;
class MolNum;

class Evaluator;
class MolEditor;

template<class T>
class Mover;

template<class T>
class Selector;

class Atom;
class CutGroup;
class Residue;
class Chain;
class Segment;

/** A Molecule represents a complete molecule. 

    Most of the manipulation of a molecule is handled by the 'or/er' classes,
    e.g. Mover, Selector, Editer, Evaluator.

    These classes provide additional member functions, thereby allowing me
    to keep the API of Molecule small.

    Examples of use include;

    mol = mol.move().translate( Vector(1,2,3) )
    point = mol.evaluate().center()
    mass = mol.evaluate().mass()

    mol = mol.edit().rename( ResNum(43)[0], "ALA" ).commit()

    Equally, we can quickly select well-defined subgroups within the
    molecule, e.g. atom(s), residue(e), chain(s), CutGroup(s) and
    segment(s), via the 'select' functions, e.g.

    ala49 = mol.select( ResName("ala") + ResNum(49) );

    or if there is more than one residue designate ALA:49

    ala49_0 = mol.select( (ResName("ala")+ResNum(49))[0] );

    or to get all of these residues, do

    all_ala49 = mol.selectAll( ResName("ala") + ResNum(49) );

    @author Christopher Woods
*/
class SIREMOL_EXPORT Molecule : public SireBase::ConcreteProperty<Molecule,MoleculeView>
{

friend QDataStream& ::operator<<(QDataStream&, const Molecule&);
friend QDataStream& ::operator>>(QDataStream&, Molecule&);

public:
    Molecule();
    Molecule(const QString &molname);
    
    Molecule(const MoleculeData &moldata);

    Molecule(const Molecule &other);

    ~Molecule();

    Molecule& operator=(const Molecule &other);

    bool operator==(const Molecule &other) const;
    bool operator!=(const Molecule &other) const;

    static const char* typeName();
    
    Molecule* clone() const;

    QString toString() const;
    
    bool isEmpty() const;
    bool selectedAll() const;
    
    AtomSelection selection() const;
    
    const MolName& name() const;
    MolNum number() const;
    
    quint64 version() const;
    quint64 version(const PropertyName &key) const;
    
    int nAtoms() const;
    int nCutGroups() const;
    int nResidues() const;
    int nChains() const;
    int nSegments() const;
    
    MoleculeInfo info() const;
    
    Mover<Molecule> move() const;
    Evaluator evaluate() const;
    MolEditor edit() const;
    
    void update(const MoleculeData &moldata);
    
    bool hasProperty(const PropertyName &key) const;
    bool hasMetadata(const PropertyName &metakey) const;
    bool hasMetadata(const PropertyName &key,
                     const PropertyName &metakey) const;
                     
    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const PropertyName &key) const;
    
    const Properties& properties() const;
    
    const Property& property(const PropertyName &key) const;
                             
    const Property& metadata(const PropertyName &metakey) const;
    
    const Property& metadata(const PropertyName &key,
                             const PropertyName &metakey) const;

    void assertContainsProperty(const PropertyName &key) const;
    
    void assertContainsMetadata(const PropertyName &metakey) const;
    void assertContainsMetadata(const PropertyName &key,
                                const PropertyName &metakey) const;

protected:
    void setProperty(const QString &key, const Property &value);
    
    void setMetadata(const QString &metakey,
                     const Property &value);
                     
    void setMetadata(const QString &key, const QString &metakey,
                     const Property &value);
};

}

Q_DECLARE_METATYPE(SireMol::Molecule);
Q_DECLARE_METATYPE(SireMol::Mover<SireMol::Molecule>);

SIRE_EXPOSE_CLASS( SireMol::Molecule )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Molecule>, SireMol::Mover_Molecule_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "mover.hpp"

template class SireMol::Mover<SireMol::Molecule>;

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
