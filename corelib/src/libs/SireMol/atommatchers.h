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

#ifndef SIREMOL_ATOMMATCHERS_H
#define SIREMOL_ATOMMATCHERS_H

#include "atommatcher.h"

SIRE_BEGIN_HEADER


namespace SireMol
{
class AtomIdxMatcher;
class AtomNameMatcher;
class AtomIDMatcher;
class AtomMultiMatcher;
class AtomMCSMatcher;
class ResIdxAtomNameMatcher;
class ResIdxAtomMCSMatcher;
class ResIdxAtomCoordMatcher;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomIdxMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomIdxMatcher&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomNameMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomNameMatcher&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomIDMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomIDMatcher&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomMultiMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomMultiMatcher&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomMCSMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomMCSMatcher&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ResIdxAtomNameMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ResIdxAtomNameMatcher&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ResIdxAtomMCSMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ResIdxAtomMCSMatcher&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ResIdxAtomCoordMatcher&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ResIdxAtomCoordMatcher&);

namespace SireMol
{

class AtomIdx;
class AtomName;
class AtomIdentifier;
class MoleculeView;
class MoleculeInfoData;

using SireBase::PropertyMap;

/** This is a simple atom matcher that matches the atoms based
    on their index in the molecule - e.g. it matches the first
    atom in molinfo0 to the first atom in molinfo1, the second
    atom in molinfo0 to the second atom in molinfo1, and the
    nth atom in molinfo0 to the nth atom in molinfo1

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomIdxMatcher
         : public SireBase::ConcreteProperty<AtomIdxMatcher,AtomMatcher>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomIdxMatcher&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, AtomIdxMatcher&);

public:
    AtomIdxMatcher();
    AtomIdxMatcher(const AtomIdxMatcher&);

    ~AtomIdxMatcher();

    static const char* typeName();

    const char* what() const
    {
        return AtomIdxMatcher::typeName();
    }

    QString toString() const;

    AtomIdxMatcher& operator=(const AtomIdxMatcher &other);

    bool operator==(const AtomIdxMatcher &other) const;
    bool operator!=(const AtomIdxMatcher &other) const;

protected:
    bool pvt_changesOrder(const MoleculeInfoData &molinfo0,
                          const MoleculeInfoData &molinfo1) const;

    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeInfoData &molinfo0,
                                     const MoleculeInfoData &molinfo1) const;

    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;
};

/** This is a simple atom matcher that matches the atoms based
    on their names, so the atom called "CA1" in molinfo0 will
    be matched to the atom called "CA1" in molinfo1

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomNameMatcher
         : public SireBase::ConcreteProperty<AtomNameMatcher,AtomMatcher>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomNameMatcher&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, AtomNameMatcher&);

public:
    AtomNameMatcher();
    AtomNameMatcher(const AtomNameMatcher&);

    ~AtomNameMatcher();

    static const char* typeName();

    const char* what() const
    {
        return AtomNameMatcher::typeName();
    }

    AtomNameMatcher& operator=(const AtomNameMatcher &other);

    bool operator==(const AtomNameMatcher &other) const;
    bool operator!=(const AtomNameMatcher &other) const;

    QString toString() const;

protected:
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeInfoData &molinfo0,
                                     const MoleculeInfoData &molinfo1) const;

    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;
};

/** This is an atom matcher that allows the user to specify
    exactly how one atom matches another in the molecule
    by mapping one AtomIdentifier to another

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomIDMatcher
        : public SireBase::ConcreteProperty<AtomIDMatcher,AtomMatcher>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomIDMatcher&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, AtomIDMatcher&);

public:
    AtomIDMatcher();

    AtomIDMatcher(const QList< QPair<QString,QString> > &match_names);
    AtomIDMatcher(const QList< QPair<int,int> > &match_idxs);
    AtomIDMatcher(const QList< QPair<AtomIdentifier,AtomIdentifier> > &match_ids);

    AtomIDMatcher(const QList< boost::tuple<QString,QString> > &match_names);
    AtomIDMatcher(const QList< boost::tuple<int,int> > &match_idxs );
    AtomIDMatcher(const QList< boost::tuple<AtomIdentifier,AtomIdentifier> > &match_ids);

    AtomIDMatcher(const QHash<QString,QString> &match_names);
    AtomIDMatcher(const QHash<int,int> &match_idxs);
    AtomIDMatcher(const QHash<AtomIdentifier,AtomIdentifier> &match_ids);

    AtomIDMatcher(const QString &match_names);

    AtomIDMatcher(const AtomIDMatcher &other);

    ~AtomIDMatcher();

    static const char* typeName();

    bool isNull() const;

    const char* what() const
    {
        return AtomIDMatcher::typeName();
    }

    AtomIDMatcher& operator=(const AtomIDMatcher &other);

    bool operator==(const AtomIDMatcher &other) const;
    bool operator!=(const AtomIDMatcher &other) const;

    QString toString() const;

protected:
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeInfoData &molinfo0,
                                     const MoleculeInfoData &molinfo1) const;

    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;

private:
    /** The mapping from atom ID in molecule 0 to atom ID in molecule 1 */
    QList< QPair<AtomIdentifier,AtomIdentifier > > m;
};

/** This is an atom matcher that matches using the maximum common substructure
    of the two molecules

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomMCSMatcher
        : public SireBase::ConcreteProperty<AtomMCSMatcher,AtomMatcher>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomMCSMatcher&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, AtomMCSMatcher&);

public:
    AtomMCSMatcher();
    AtomMCSMatcher(bool verbose);
    AtomMCSMatcher(const SireUnits::Dimension::Time &timeout,
                   bool verbose);
    AtomMCSMatcher(const AtomMatcher &prematcher,
                   bool verbose);
    AtomMCSMatcher(const AtomMatcher &prematcher,
                   const SireUnits::Dimension::Time &timeout,
                   bool verbose);
    AtomMCSMatcher(bool match_light_atoms,
                   bool verbose);
    AtomMCSMatcher(const SireUnits::Dimension::Time &timeout,
                   bool match_light_atoms,
                   bool verbose);
    AtomMCSMatcher(const AtomMatcher &prematcher,
                   bool match_light_atoms,
                   bool verbose);
    AtomMCSMatcher(const AtomMatcher &prematcher,
                   const SireUnits::Dimension::Time &timeout,
                   bool match_light_atoms,
                   bool verbose);

    AtomMCSMatcher(const AtomMCSMatcher &other);

    ~AtomMCSMatcher();

    static const char* typeName();

    const char* what() const
    {
        return AtomMCSMatcher::typeName();
    }

    AtomMCSMatcher& operator=(const AtomMCSMatcher &other);

    bool operator==(const AtomMCSMatcher &other) const;
    bool operator!=(const AtomMCSMatcher &other) const;

    QString toString() const;

    const AtomMatcher& preMatcher() const;

    SireUnits::Dimension::Time timeout() const;

    bool matchingLightAtoms() const;

    bool isVerbose() const;

protected:
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;

private:
    /** The pre-matcher */
    AtomMatcherPtr prematcher;

    /** Timeout for the MCS match */
    SireUnits::Dimension::Time t;

    /** Whether or not to match light atoms */
    bool match_light;

    /** Whether or not to report progess to stdout */
    bool verbose;
};

/** This is an atom matcher combines several sub-AtomMatchers together

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomMultiMatcher
        : public SireBase::ConcreteProperty<AtomMultiMatcher,AtomMatcher>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomMultiMatcher&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, AtomMultiMatcher&);

public:
    AtomMultiMatcher();
    AtomMultiMatcher(const AtomMatcher &matcher);
    AtomMultiMatcher(const AtomMatcher &m0, const AtomMatcher &m1);

    AtomMultiMatcher(const AtomMultiMatcher &other);

    ~AtomMultiMatcher();

    static const char* typeName();

    bool isNull() const;

    const char* what() const
    {
        return AtomMultiMatcher::typeName();
    }

    AtomMultiMatcher& operator=(const AtomMultiMatcher &other);

    bool operator==(const AtomMultiMatcher &other) const;
    bool operator!=(const AtomMultiMatcher &other) const;

    QString toString() const;

protected:
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeInfoData &molinfo0,
                                     const MoleculeInfoData &molinfo1) const;

    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;

private:
    /** The set of matches, which are processed in order */
    QList<AtomMatcherPtr> m;
};

/** Match atoms by name within each residue.

    @author Lester Hedges
*/
class SIREMOL_EXPORT ResIdxAtomNameMatcher
         : public SireBase::ConcreteProperty<ResIdxAtomNameMatcher,AtomMatcher>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ResIdxAtomNameMatcher&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ResIdxAtomNameMatcher&);

public:
    ResIdxAtomNameMatcher();
    ResIdxAtomNameMatcher(const ResIdxAtomNameMatcher&);

    ~ResIdxAtomNameMatcher();

    static const char* typeName();

    const char* what() const
    {
        return ResIdxAtomNameMatcher::typeName();
    }

    ResIdxAtomNameMatcher& operator=(const ResIdxAtomNameMatcher &other);

    bool operator==(const ResIdxAtomNameMatcher &other) const;
    bool operator!=(const ResIdxAtomNameMatcher &other) const;

    QString toString() const;

protected:
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeInfoData &molinfo0,
                                     const MoleculeInfoData &molinfo1) const;

    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;
};

/** Match atoms by name MCS within each residue.

    @author Lester Hedges
*/
class SIREMOL_EXPORT ResIdxAtomMCSMatcher
         : public SireBase::ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ResIdxAtomMCSMatcher&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ResIdxAtomMCSMatcher&);

public:
    ResIdxAtomMCSMatcher();
    ResIdxAtomMCSMatcher(bool verbose);
    ResIdxAtomMCSMatcher(const SireUnits::Dimension::Time &timeout,
                         bool verbose);
    ResIdxAtomMCSMatcher(const AtomMatcher &prematcher,
                         bool verbose);
    ResIdxAtomMCSMatcher(const AtomMatcher &prematcher,
                         const SireUnits::Dimension::Time &timeout,
                         bool verbose);
    ResIdxAtomMCSMatcher(bool match_light_atoms,
                         bool verbose);
    ResIdxAtomMCSMatcher(const SireUnits::Dimension::Time &timeout,
                         bool match_light_atoms,
                         bool verbose);
    ResIdxAtomMCSMatcher(const AtomMatcher &prematcher,
                         bool match_light_atoms,
                         bool verbose);
    ResIdxAtomMCSMatcher(const AtomMatcher &prematcher,
                         const SireUnits::Dimension::Time &timeout,
                         bool match_light_atoms,
                         bool verbose);

    ResIdxAtomMCSMatcher(const ResIdxAtomMCSMatcher&);

    ~ResIdxAtomMCSMatcher();

    static const char* typeName();

    const char* what() const
    {
        return ResIdxAtomMCSMatcher::typeName();
    }

    ResIdxAtomMCSMatcher& operator=(const ResIdxAtomMCSMatcher &other);

    bool operator==(const ResIdxAtomMCSMatcher &other) const;
    bool operator!=(const ResIdxAtomMCSMatcher &other) const;

    QString toString() const;

    const AtomMatcher& preMatcher() const;

    SireUnits::Dimension::Time timeout() const;

    bool matchingLightAtoms() const;

    bool isVerbose() const;

protected:
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;

private:
    /** The pre-matcher */
    AtomMatcherPtr prematcher;

    /** Timeout for the MCS match */
    SireUnits::Dimension::Time t;

    /** Whether or not to match light atoms */
    bool match_light;

    /** Whether or not to report progess to stdout */
    bool verbose;
};

/** Match atoms by coordinates within each residue.

    @author Lester Hedges
*/
class SIREMOL_EXPORT ResIdxAtomCoordMatcher
         : public SireBase::ConcreteProperty<ResIdxAtomCoordMatcher,AtomMatcher>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ResIdxAtomCoordMatcher&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ResIdxAtomCoordMatcher&);

public:
    ResIdxAtomCoordMatcher();
    ResIdxAtomCoordMatcher(const ResIdxAtomCoordMatcher&);

    ~ResIdxAtomCoordMatcher();

    static const char* typeName();

    const char* what() const
    {
        return ResIdxAtomCoordMatcher::typeName();
    }

    ResIdxAtomCoordMatcher& operator=(const ResIdxAtomCoordMatcher &other);

    bool operator==(const ResIdxAtomCoordMatcher &other) const;
    bool operator!=(const ResIdxAtomCoordMatcher &other) const;

    QString toString() const;

protected:
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &molview0,
                                     const PropertyMap &map0,
                                     const MoleculeView &molview1,
                                     const PropertyMap &map1) const;
};

}

Q_DECLARE_METATYPE( SireMol::AtomIdxMatcher )
Q_DECLARE_METATYPE( SireMol::AtomNameMatcher )
Q_DECLARE_METATYPE( SireMol::AtomIDMatcher )
Q_DECLARE_METATYPE( SireMol::AtomMultiMatcher )
Q_DECLARE_METATYPE( SireMol::AtomMCSMatcher )
Q_DECLARE_METATYPE( SireMol::ResIdxAtomNameMatcher )
Q_DECLARE_METATYPE( SireMol::ResIdxAtomMCSMatcher )
Q_DECLARE_METATYPE( SireMol::ResIdxAtomCoordMatcher )

SIRE_EXPOSE_CLASS( SireMol::AtomIdxMatcher )
SIRE_EXPOSE_CLASS( SireMol::AtomNameMatcher )
SIRE_EXPOSE_CLASS( SireMol::AtomIDMatcher )
SIRE_EXPOSE_CLASS( SireMol::AtomMultiMatcher )
SIRE_EXPOSE_CLASS( SireMol::AtomMCSMatcher )
SIRE_EXPOSE_CLASS( SireMol::ResIdxAtomNameMatcher )
SIRE_EXPOSE_CLASS( SireMol::ResIdxAtomMCSMatcher )
SIRE_EXPOSE_CLASS( SireMol::ResIdxAtomCoordMatcher )

SIRE_END_HEADER

#endif
