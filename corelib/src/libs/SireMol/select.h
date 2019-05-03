/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Christopher Woods
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

#ifndef SIREMOL_SELECT_H
#define SIREMOL_SELECT_H

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireMol/viewsofmol.h"
#include "SireMol/mover.hpp"
#include "SireMol/atom.h"
#include "SireMol/cutgroup.h"
#include "SireMol/residue.h"
#include "SireMol/chain.h"
#include "SireMol/segment.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/moleculegroups.h"

#include <boost/shared_ptr.hpp>

#include <QList>

SIRE_BEGIN_HEADER

namespace SireMol
{
class Select;
class SelectResult;
class SelectResultMover;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Select&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Select&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::SelectResult&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::SelectResult&);

namespace SireMol
{

class MolGroupsBase;
class MoleculeGroup;
class Molecules;
class MoleculeView;
class ViewsOfMol;

namespace parser
{

class SelectEngine;

using SelectEnginePtr = boost::shared_ptr<SelectEngine>;
using SelectEngineWeakPtr = boost::weak_ptr<SelectEngine>;

/** This is the base class of all of the select objects. It is a private
    object that should only be used by Select

    @author Christopher Woods
*/
class SelectEngine
{

public:
    enum ObjType { COMPLEX = 0,
                   ATOM = 1,
                   CUTGROUP = 2,
                   RESIDUE = 3,
                   CHAIN = 4,
                   SEGMENT = 5,
                   MOLECULE = 6 };

    virtual ~SelectEngine();
    
    SelectResult operator()(const SelectResult &result,
                            const PropertyMap &map = PropertyMap()) const;
    
    SelectResult operator()(const MolGroupsBase &molgroups,
                            const PropertyMap &map = PropertyMap()) const;
    
    SelectResult operator()(const MoleculeGroup &molgroup,
                            const PropertyMap &map = PropertyMap()) const;
    
    SelectResult operator()(const Molecules &molecules,
                            const PropertyMap &map = PropertyMap()) const;
    
    SelectResult operator()(const MoleculeView &molecule,
                            const PropertyMap &map = PropertyMap()) const;
    
    virtual SelectEnginePtr simplify();
    
    virtual bool usesCoordinates() const;
    
    bool hasParent() const;
    
    SelectEnginePtr self();
    
    void setParent(SelectEnginePtr parent) const;
    
    virtual ObjType objectType() const=0;
    
protected:
    SelectEngine();

    virtual SelectResult select(const SelectResult &result, const PropertyMap &map) const=0;

    static SelectEnginePtr makePtr(SelectEngine *ptr);
    
    virtual ViewsOfMol expandMol(const ViewsOfMol &mol) const;
    virtual SelectResult expand(const SelectResult &result) const;
    
    /** The parent engine */
    SelectEngineWeakPtr parent;
    
    /** Weak pointer to self */
    SelectEngineWeakPtr selfptr;
};

} //end of namespace parser

/** This is the only publicly visible selector class. This provides a 
    front-end interface to selecting atoms and molecules
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT Select : public SireBase::ConcreteProperty<Select,SireBase::Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Select&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, Select&);

public:
    Select();
    Select(const QString &str);
    
    Select(const Select &other);
    
    ~Select();
    
    Select& operator=(const Select &other);
    
    bool operator==(const Select &other) const;
    bool operator!=(const Select &other) const;

    Select* clone() const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;

    QString objectType() const;

    SelectResult operator()(const MolGroupsBase &molgroups,
                            const PropertyMap &map = PropertyMap()) const;
    
    SelectResult operator()(const MoleculeGroup &molgroup,
                            const PropertyMap &map = PropertyMap()) const;
    
    SelectResult operator()(const Molecules &molecules,
                            const PropertyMap &map = PropertyMap()) const;
    
    SelectResult operator()(const MoleculeView &molecule,
                            const PropertyMap &map = PropertyMap()) const;

    SelectResult operator()(const SelectResult &result,
                            const PropertyMap &map = PropertyMap()) const;

    static void setToken(const QString &token, const QString &selection);
    static void resetTokens();

private:
    /** The actual search string */
    QString search_string;

    /** The underlying engine used to perform the selection */
    parser::SelectEnginePtr e;
};

/** This class holds the result of a Select 

    @author Christopher Woods
*/
class SIREMOL_EXPORT SelectResult
    : public SireBase::ConcreteProperty<SelectResult,SireBase::Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream &ds, const SelectResult&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream &ds, SelectResult&);

public:
    typedef QList<ViewsOfMol> Container;
    typedef Container::const_iterator const_iterator;
    typedef const_iterator iterator;

    SelectResult();
    
    SelectResult(const MolGroupsBase &molgroups);
    SelectResult(const MoleculeGroup &molgroup);
    SelectResult(const Molecules &molecules);
    SelectResult(const MoleculeView &molview);
    
    SelectResult(const QList<ViewsOfMol> molviews);
    
    SelectResult(const SelectResult &other);

    ~SelectResult();

    SelectResult& operator=(const SelectResult &other);
    
    bool operator==(const SelectResult &other) const;
    bool operator!=(const SelectResult &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    const_iterator begin() const;
    const_iterator end() const;
    
    const_iterator constBegin() const;
    const_iterator constEnd() const;
    
    SelectResult* clone() const;

    QString toString() const;

    MolViewPtr operator[](int i) const;
    MolViewPtr operator[](MolNum molnum) const;

    SelectResultMover move() const;
    
    bool isEmpty() const;
    
    int count() const;
    int size() const;

    bool contains(MolNum molnum) const;

    QList<ViewsOfMol> views() const;
    ViewsOfMol views(MolNum molnum) const;
    
    QList<MolNum> molNums() const;

    MoleculeGroup toGroup() const;
    MoleculeGroup toGroup(const QString &name) const;
    
    Molecules toMolecules() const;

    SelectResult search(const QString &search_term) const;

    SelectResult join() const;
    
    SelectResult atoms() const;
    SelectResult cutGroups() const;
    SelectResult residues() const;
    SelectResult chains() const;
    SelectResult segments() const;
    SelectResult molecules() const;

private:
    /** The list of all ViewsOfMol views */
    QList<ViewsOfMol> molviews;
};

/** This class provides a simple "move" interface to move all 
    views in a SelectResult 
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT SelectResultMover
    : public SireBase::ConcreteProperty<SelectResultMover,SireBase::Property>
{
public:
    SelectResultMover();
    SelectResultMover(const SelectResult &other);
    
    SelectResultMover(const SelectResultMover &other);
    
    ~SelectResultMover();

    SelectResultMover& operator=(const SelectResultMover &other);
    
    bool operator==(const SelectResultMover &other) const;
    bool operator!=(const SelectResultMover &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    SelectResultMover* clone() const;

    QString toString() const;
    
    SelectResultMover& translate(const Vector &delta);
    
    SelectResult commit() const;

private:
    /** All of the views that will be moved */
    QList< Mover<ViewsOfMol> > molviews;
};

} //end of namespace SireMol

Q_DECLARE_METATYPE( SireMol::Select )
Q_DECLARE_METATYPE( SireMol::SelectResult )
Q_DECLARE_METATYPE( SireMol::SelectResultMover )

SIRE_EXPOSE_CLASS( SireMol::Select )
SIRE_EXPOSE_CLASS( SireMol::SelectResult )
SIRE_EXPOSE_CLASS( SireMol::SelectResultMover )

SIRE_END_HEADER

#endif
