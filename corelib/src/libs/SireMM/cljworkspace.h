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

#ifndef SIREMM_CLJWORKSPACE_H
#define SIREMM_CLJWORKSPACE_H

#include "cljdelta.h"

#include "SireMol/moleculeview.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJWorkspace;
namespace detail{ class CLJWorkspaceData; }
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJWorkspace&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJWorkspace&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::detail::CLJWorkspaceData&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::detail::CLJWorkspaceData&);

namespace SireMM
{

using boost::tuple;

/** This class provides a workspace in which to hold the details of the changes
    that occur in a CLJ forcefield during a Monte Carlo move. The class is optimised
    to avoid copying or duplicating data during Sire's copy-on-write copying
    (e.g. the memory allocated in a workspace will always be available for the 
     new copy of a forcefield rather than the old, which, if things work correctly,
     will mean that there should be no memory allocation during simple MC moves...)
     
     @author Christopher Woods
*/
class SIREMM_EXPORT CLJWorkspace
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const CLJWorkspace&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, CLJWorkspace&);

public:
    CLJWorkspace();
    CLJWorkspace(const CLJWorkspace &other);
    
    ~CLJWorkspace();
    
    CLJWorkspace& operator=(const CLJWorkspace &other);
    
    bool operator==(const CLJWorkspace &other) const;
    bool operator!=(const CLJWorkspace &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    QString toString() const;
    
    const CLJDelta& operator[](int i) const;

    const CLJDelta& at(int i) const;
    
    CLJDelta getitem(int i) const;
    
    int count() const;
    int size() const;
    
    const CLJDelta* data() const;
    const CLJDelta* constData() const;
    
    CLJDelta push(CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
                  const CLJAtoms &new_atoms, const CLJDelta &old_delta);
    
    void removeSameIDAtoms(CLJBoxes &boxes);
    
    QVector<CLJBoxIndex> commit(CLJBoxes &boxes, const CLJDelta &delta);
    QVector<CLJBoxIndex> revert(CLJBoxes &boxes, const CLJDelta &delta);
    
    int nDeltas() const;
    
    bool isSingleID() const;
    
    bool isEmpty() const;
    
    tuple<CLJAtoms,CLJAtoms,CLJAtoms> merge() const;
    
    CLJAtoms changedAtoms() const;
    CLJAtoms oldAtoms() const;
    CLJAtoms newAtoms() const;

    void clear();

    bool needsAccepting() const;

    void accept(CLJBoxes &boxes);
    void mustRecalculateFromScratch(CLJBoxes &boxes);

    bool recalculatingFromScratch() const;

private:
    void returnToMemoryPool();
    void createFromMemoryPool();
    
    void detach();

    /** Implicitly shared pointer to the data */
    boost::shared_ptr<detail::CLJWorkspaceData> d;
    
    /** Whether or not we are recalculating everything from scratch */
    bool recalc_from_scratch;
};

}

Q_DECLARE_METATYPE( SireMM::CLJWorkspace )

SIRE_EXPOSE_CLASS( SireMM::CLJWorkspace )

SIRE_END_HEADER

#endif
