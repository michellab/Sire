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

#ifndef SIREMOL_MOLECULES_H
#define SIREMOL_MOLECULES_H

#include <boost/tuple/tuple.hpp>

#include "viewsofmol.h"
#include "molnum.h"

#include "SireBase/chunkedhash.hpp"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Molecules;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Molecules&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Molecules&);

namespace SireMol
{

class MolNumViewIdx;
class SelectResult;

/** This class provides a container for lots of molecules. This
    forms a general purpose molecule container, which is used as the argument
    to functions which expect to be passed lots of molecules or parts
    of molecules. This class holds the Molecules using the
    ViewsOfMol class, thereby allowing multiple arbitrary views of each 
    molecule to be held.

    @author Christopher Woods
*/
class SIREMOL_EXPORT Molecules 
            : public SireBase::ConcreteProperty<Molecules,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const Molecules&);
friend QDataStream& ::operator>>(QDataStream&, Molecules&);

public:

    typedef SireBase::ChunkedHash<MolNum,ViewsOfMol>::const_iterator const_iterator;
    typedef SireBase::ChunkedHash<MolNum,ViewsOfMol>::iterator iterator;

    Molecules();

    Molecules(const MoleculeView &molecule);
    Molecules(const ViewsOfMol &molviews);

    Molecules(const SelectResult &result);

    template<class T>
    explicit Molecules(const QList<T> &molecules);

    template<class T>
    explicit Molecules(const QVector<T> &molecules);

    Molecules(const Molecules &other);

    ~Molecules();

    static const char* typeName();

    const char* what() const
    {
        return Molecules::typeName();
    }

    Molecules* clone() const;

    Molecules& operator=(const Molecules &other);
    
    bool operator==(const Molecules &other) const;
    bool operator!=(const Molecules &other) const;
    
    QString toString() const;
    
    const ViewsOfMol& operator[](MolNum molnum) const;
    PartialMolecule operator[](const boost::tuple<MolNum,SireID::Index> &viewidx) const;
    
    const ViewsOfMol& at(MolNum molnum) const;
    PartialMolecule at(const boost::tuple<MolNum,SireID::Index> &viewidx) const;

    PartialMolecule at(MolNum molnum, int idx) const;

    Molecules operator+(const Molecules &other) const;
    Molecules operator-(const Molecules &other) const;

    Molecules& operator+=(const Molecules &other);
    Molecules& operator-=(const Molecules &other);

    const ViewsOfMol& molecule(MolNum molnum) const;

    bool isEmpty() const;

    bool contains(MolNum molnum) const;
    bool contains(const MoleculeView &molview) const;
    bool contains(const ViewsOfMol &molviews) const;
    bool contains(const Molecules &molecules) const;

    bool intersects(const MoleculeView &molview) const;
    bool intersects(const Molecules &other) const;

    SelectResult search(const QString &search_string) const;

    int count() const;
    int nMolecules() const;
    
    int nViews() const;
    int nViews(MolNum molnum) const;

    const ViewsOfMol& first() const;
    const ViewsOfMol& last() const;
    
    const ViewsOfMol& front() const;
    const ViewsOfMol& back() const;

    const_iterator begin() const;
    const_iterator end() const;

    const_iterator constBegin() const;
    const_iterator constEnd() const;

    const_iterator find(MolNum molnum) const;
    const_iterator constFind(MolNum molnum) const;

    QSet<MolNum> molNums() const;

    void assertContains(MolNum molnum) const;

    void add(const MoleculeView &molview);
    void add(const ViewsOfMol &molviews);
    void add(const Molecules &molecules);

    bool addIfUnique(const MoleculeView &molview);
    ViewsOfMol addIfUnique(const ViewsOfMol &molviews);
    QList<ViewsOfMol> addIfUnique(const Molecules &molecules);

    bool unite(const MoleculeView &molview);
    ViewsOfMol unite(const ViewsOfMol &molviews);
    QList<ViewsOfMol> unite(const Molecules &other);
    
    bool remove(const MoleculeView &molview);
    ViewsOfMol remove(const ViewsOfMol &molviews);
    QList<ViewsOfMol> remove(const Molecules &molecules);
    
    bool removeAll(const MoleculeView &molview);
    ViewsOfMol removeAll(const ViewsOfMol &molviews);
    QList<ViewsOfMol> removeAll(const Molecules &molecules);
    
    ViewsOfMol remove(MolNum molnum);
    bool removeAll();

    void clear();

    bool uniteViews();
    bool removeDuplicates();

    void reserve(int nmolecules);

    bool update(const MoleculeData &moldata);
    bool update(const MoleculeView &molview);
    QList<Molecule> update(const Molecules &molecules);

private:
    template<class T>
    static Molecules from(const T &molecules);

    /** Hash that contains all of the views of
        all of the molecules, indexed by 
        their molecule number */
    SireBase::ChunkedHash<MolNum,ViewsOfMol> mols;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Molecules Molecules::from(const T &molecules)
{
    Molecules mols;

    if (molecules.count() == 0)
        return mols;

    mols.mols.reserve(molecules.count());

    for (typename T::const_iterator it = molecules.begin();
         it != molecules.end();
         ++it)
    {
        SireBase::ChunkedHash<MolNum,ViewsOfMol>::iterator mol 
                                                    = mols.mols.find(it->number());
        
        if (mol != mols.mols.end())
        {
            mol->add(*it);
        }
        else
        {
            mols.mols.insert(it->number(), *it);
        }
    }

    return mols;
}

/** Converting constructor used to convert from general
    containers of molecules to this container */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Molecules::Molecules(const QList<T> &molecules)
          : SireBase::ConcreteProperty<Molecules,SireBase::Property>()
{
    *this = Molecules::from(molecules);
}

/** Converting constructor used to convert from general
    containers of molecules to this container */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Molecules::Molecules(const QVector<T> &molecules)
          : SireBase::ConcreteProperty<Molecules,SireBase::Property>()
{
    *this = Molecules::from(molecules);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMol::Molecules);

SIRE_EXPOSE_CLASS( SireMol::Molecules )

SIRE_END_HEADER

#endif
