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

#ifndef SIREMM_CLJDELTA_H
#define SIREMM_CLJDELTA_H

#include "cljatoms.h"
#include "cljboxes.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJDelta;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJDelta&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJDelta&);

namespace SireMM
{

class CLJBoxes;

using boost::tuple;

/** This class is used to hold the change in coordinates etc. of a set of atoms caused
    by e.g. a Monte Carlo move
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJDelta
{

friend QDataStream& ::operator<<(QDataStream&, const CLJDelta&);
friend QDataStream& ::operator>>(QDataStream&, CLJDelta&);

public:
    CLJDelta();
    CLJDelta(qint32 idnum, const CLJAtoms &oldatoms, const CLJAtoms &newatoms);
    
    CLJDelta(const CLJDelta &other);
    
    ~CLJDelta();
    
    CLJDelta& operator=(const CLJDelta &other);
    
    bool operator==(const CLJDelta &other) const;
    bool operator!=(const CLJDelta &other) const;
    
    static const char* typeName();
    
    const char* what() const;

    QString toString() const;

    bool isEmpty() const;
    bool isNull() const;
    
    qint32 ID() const;
    
    CLJAtoms newAtoms() const;
    CLJAtoms oldAtoms() const;

    CLJAtoms changedAtoms() const;
    
    void assertIdenticalTo(const CLJDelta &other) const;
    
    static CLJAtoms mergeChanged(const CLJDelta *deltas, int count);
    static CLJAtoms mergeChanged(const QVector<CLJDelta> &deltas);

    static CLJAtoms mergeNew(const CLJDelta *deltas, int count);
    static CLJAtoms mergeNew(const QVector<CLJDelta> &deltas);

    static CLJAtoms mergeOld(const CLJDelta *deltas, int count);
    static CLJAtoms mergeOld(const QVector<CLJDelta> &deltas);
    
    static tuple<CLJAtoms,CLJAtoms,CLJAtoms> merge(const CLJDelta *deltas, int count);
    static tuple<CLJAtoms,CLJAtoms,CLJAtoms> merge(const QVector<CLJDelta> &deltas);
    
private:
    /** The old atoms */
    CLJAtoms old_atoms;

    /** The new atoms */
    CLJAtoms new_atoms;
    
    /** The ID number of this delta in the CLJWorkspace that created
        and manages it */
    qint32 idnum;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the ID number of this delta. This can be used for book-keeping
    by the object that created this delta */
inline qint32 CLJDelta::ID() const
{
    return idnum;
}

/** Return the old version of the atoms in a format that should make
    them easy to rebox */
inline CLJAtoms CLJDelta::oldAtoms() const
{
    return old_atoms;
}

/** Return the new version of the atoms in a format that should make
    them easy to rebox */
inline CLJAtoms CLJDelta::newAtoms() const
{
    return new_atoms;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::CLJDelta )

SIRE_EXPOSE_CLASS( SireMM::CLJDelta )

SIRE_END_HEADER

#endif
