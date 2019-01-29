/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIREMOL_ATOMMATCHER_H
#define SIREMOL_ATOMMATCHER_H

#include <QHash>
#include <QList>
#include <QPair>

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireUnits/dimensions.h"

#include <boost/tuple/tuple.hpp>

namespace SireMol
{
class AtomMatcher;
class AtomIdxMatcher;
class AtomNameMatcher;
class AtomResultMatcher;
class AtomMatchInverter;
class AtomIDMatcher;
class AtomMultiMatcher;
class AtomMCSMatcher;
class ResIdxAtomNameMatcher;
class ResIdxAtomMCSMatcher;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomMatcher&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomResultMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomResultMatcher&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomMatchInverter&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomMatchInverter&);

namespace SireMol
{

class AtomIdx;
class AtomName;
class AtomIdentifier;
class MoleculeView;
class MoleculeInfoData;

using SireBase::PropertyMap;

/** Virtual base class of all of the functions used to match
    atoms in one molecule layout with atoms in another layout
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomMatcher : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const AtomMatcher&);
friend QDataStream& ::operator>>(QDataStream&, AtomMatcher&);

public:
    AtomMatcher();
    AtomMatcher(const AtomMatcher &other);
    
    virtual ~AtomMatcher();
    
    static const char* typeName()
    {
        return "SireMol::AtomMatcher";
    }
    
    virtual bool isNull() const;
    
    AtomMultiMatcher operator+(const AtomMatcher &other) const;
    
    AtomMultiMatcher add(const AtomMatcher &other) const;
    
    virtual AtomMatcher* clone() const=0;

    bool changesOrder(const MoleculeInfoData &molinfo0,
                      const MoleculeInfoData &molinfo1) const;

    bool changesOrder(const MoleculeView &molview0,
                      const MoleculeView &molview1) const;
    
    bool changesOrder(const MoleculeView &molview0,
                      const PropertyMap &map0,
                      const MoleculeView &molview1,
                      const PropertyMap &map1) const;
    
    bool changesOrder(const MoleculeView &molview0,
                      const MoleculeView &molview1,
                      const PropertyMap &map) const;
    
    QHash<AtomIdx,AtomIdx> match(const MoleculeInfoData &molinfo0,
                                 const MoleculeInfoData &molinfo1) const;
    
    QHash<AtomIdx,AtomIdx> match(const MoleculeView &molview0,
                                 const PropertyMap &map0,
                                 const MoleculeView &molview1,
                                 const PropertyMap &map1) const;

    QHash<AtomIdx,AtomIdx> match(const MoleculeView &molview0,
                                 const MoleculeView &molview1) const;
    
    QHash<AtomIdx,AtomIdx> match(const MoleculeView &molview0,
                                 const MoleculeView &molview1,
                                 const PropertyMap &map) const;

    static const AtomMultiMatcher& null();

protected:
    virtual bool pvt_changesOrder(const MoleculeInfoData &molinfo0,
                                  const MoleculeInfoData &molinfo1) const;

    virtual bool pvt_changesOrder(const MoleculeView &molview0,
                                  const MoleculeView &molview1) const;
    
    virtual bool pvt_changesOrder(const MoleculeView &molview0,
                                  const PropertyMap &map0,
                                  const MoleculeView &molview1,
                                  const PropertyMap &map1) const;
    
    virtual bool pvt_changesOrder(const MoleculeView &molview0,
                                  const MoleculeView &molview1,
                                  const PropertyMap &map) const;
    
    virtual QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeInfoData &molinfo0,
                                             const MoleculeInfoData &molinfo1) const;
    
    virtual QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                             const PropertyMap &map0,
                                             const MoleculeView &molview1,
                                             const PropertyMap &map1) const=0;

    virtual QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                             const MoleculeView &molview1) const;
    
    virtual QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                             const MoleculeView &molview1,
                                             const PropertyMap &map) const;
};

typedef SireBase::PropPtr<AtomMatcher> AtomMatcherPtr;

/** This is a simple atom matcher that can be used to repeat a match
    based on a prior result
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomResultMatcher
         : public SireBase::ConcreteProperty<AtomResultMatcher,AtomMatcher>
{

friend QDataStream& ::operator<<(QDataStream&, const AtomResultMatcher&);
friend QDataStream& ::operator>>(QDataStream&, AtomResultMatcher&);

public:
    AtomResultMatcher();
    AtomResultMatcher(const QHash<AtomIdx,AtomIdx> &results, bool invert=false);
    
    AtomResultMatcher(const AtomResultMatcher &other);
    
    ~AtomResultMatcher();
    
    static const char* typeName();
    
    const char* what() const
    {
        return AtomResultMatcher::typeName();
    }

    bool isNull() const;

    QString toString() const;

    AtomResultMatcher& operator=(const AtomResultMatcher &other);
    
    bool operator==(const AtomResultMatcher &other) const;
    bool operator!=(const AtomResultMatcher &other) const;

protected:
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeInfoData &molinfo0,
                                     const MoleculeInfoData &molinfo1) const;
    
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;

private:
    /** The result of matching using another AtomMatcher */
    QHash<AtomIdx,AtomIdx> m;
};

/** This is a atom matcher that inverts the match of the sub-matcher.
    This is useful when you want to match from molecule 1 to molecule 0
    as opposed to molecule 0 to molecule 1, and but don't want to change
    the order of the match at the calling site
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomMatchInverter
         : public SireBase::ConcreteProperty<AtomMatchInverter,AtomMatcher>
{

friend QDataStream& ::operator<<(QDataStream&, const AtomMatchInverter&);
friend QDataStream& ::operator>>(QDataStream&, AtomMatchInverter&);

public:
    AtomMatchInverter();
    AtomMatchInverter(const AtomMatcher &matcher);
    
    AtomMatchInverter(const AtomMatchInverter &other);
    
    ~AtomMatchInverter();
    
    static const char* typeName();
    
    const char* what() const
    {
        return AtomMatchInverter::typeName();
    }

    bool isNull() const;

    QString toString() const;

    AtomMatchInverter& operator=(const AtomMatchInverter &other);
    
    bool operator==(const AtomMatchInverter &other) const;
    bool operator!=(const AtomMatchInverter &other) const;

protected:
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeInfoData &molinfo0,
                                     const MoleculeInfoData &molinfo1) const;
    
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;

private:
    /** The matcher used for the forwards match */
    AtomMatcherPtr m;
};

}

Q_DECLARE_METATYPE( SireMol::AtomResultMatcher )
Q_DECLARE_METATYPE( SireMol::AtomMatchInverter )

SIRE_EXPOSE_CLASS( SireMol::AtomMatcher )
SIRE_EXPOSE_CLASS( SireMol::AtomResultMatcher )
SIRE_EXPOSE_CLASS( SireMol::AtomMatchInverter )

SIRE_EXPOSE_PROPERTY( SireMol::AtomMatcherPtr, SireMol::AtomMatcher )

#endif
