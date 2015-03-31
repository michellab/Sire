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

#include "cljdelta.h"

#include "SireMol/moleculeview.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<CLJDelta> r_delta(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJDelta &delta)
{
    writeHeader(ds, r_delta, 1);
    
    SharedDataStream sds(ds);
    
    sds << delta.old_atoms << delta.new_atoms << delta.idnum;
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJDelta &delta)
{
    VersionID v = readHeader(ds, r_delta);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> delta.old_atoms >> delta.new_atoms >> delta.idnum;
    }
    else
        throw version_error(v, "1", r_delta, CODELOC);
    
    return ds;
}

/** Null constructor */
CLJDelta::CLJDelta() : idnum(-1)
{}

/** Construct the delta that changes from 'oldatoms' to 'newatoms' */
CLJDelta::CLJDelta(qint32 num, const CLJAtoms &oldatoms, const CLJAtoms &newatoms)
         : old_atoms(oldatoms), new_atoms(newatoms), idnum(num)
{
    if (idnum < 0)
    {
        idnum = -1;
        old_atoms = CLJAtoms();
        new_atoms = CLJAtoms();
    }
}

/** Copy constructor */
CLJDelta::CLJDelta(const CLJDelta &other)
         : old_atoms(other.old_atoms), new_atoms(other.new_atoms), idnum(other.idnum)
{}

/** Destructor */
CLJDelta::~CLJDelta()
{}

/** Copy assignment operator */
CLJDelta& CLJDelta::operator=(const CLJDelta &other)
{
    if (this != &other)
    {
        old_atoms = other.old_atoms;
        new_atoms = other.new_atoms;
        idnum = other.idnum;
    }
    
    return *this;
}

/** Comparison operator */
bool CLJDelta::operator==(const CLJDelta &other) const
{
    return this == &other or
           (new_atoms == other.new_atoms and
            old_atoms == other.old_atoms and
            idnum == other.idnum);
}

/** Comparison operator */
bool CLJDelta::operator!=(const CLJDelta &other) const
{
    return not operator==(other);
}

const char* CLJDelta::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJDelta>() );
}

const char* CLJDelta::what() const
{
    return CLJDelta::typeName();
}

bool CLJDelta::isNull() const
{
    return idnum < 0;
}

/** Assert that this CLJDelta is equal to 'other' */
void CLJDelta::assertIdenticalTo(const CLJDelta &other) const
{
    if (this->operator!=(other))
    {
        QStringList differences;
        
        for (int i=0; i<qMin(old_atoms.count(),other.old_atoms.count()); ++i)
        {
            if (old_atoms.at(i) != other.old_atoms.at(i))
            {
                differences.append( QObject::tr(
                        "OLD %1: %2 vs. %3")
                            .arg(i).arg(old_atoms.at(i).toString())
                            .arg(other.old_atoms.at(i).toString()) );
            }
        }

        if (old_atoms.count() != other.old_atoms.count())
            differences.append( QObject::tr(
                    "OLD: number of atoms - %1 vs. %2")
                        .arg(old_atoms.count()).arg(other.old_atoms.count()) );
        
        for (int i=0; i<qMin(new_atoms.count(),other.new_atoms.count()); ++i)
        {
            if (new_atoms.at(i) != other.new_atoms.at(i))
            {
                differences.append( QObject::tr(
                        "NEW %1: %2 vs. %3")
                            .arg(i).arg(new_atoms.at(i).toString())
                            .arg(other.new_atoms.at(i).toString()) );
            }
        }

        if (new_atoms.count() != other.new_atoms.count())
            differences.append( QObject::tr(
                    "NEW: number of atoms - %1 vs. %2")
                        .arg(new_atoms.count()).arg(other.new_atoms.count()) );
    
        if (not differences.isEmpty())
            throw SireError::assertation_failed( QObject::tr(
                    "This CLJDelta is not equal to the other. Difference are:\n%1\n")
                        .arg(differences.join("\n")), CODELOC );
    }
}

/** Return whether or not this change is empty (has no atoms or no change) */
bool CLJDelta::isEmpty() const
{
    return old_atoms == new_atoms;
}

QString CLJDelta::toString() const
{
    if (isNull())
        return QObject::tr("CLJDelta::null");
    else
        return QObject::tr("CLJDelta( nChanged() = %1 )")
                    .arg(changedAtoms().count());
}

/** Return difference between the old and new atoms. This returns the change
    as only the atoms that have changed, with the parameters of the old atoms
    negated so that a delta energy can be calculated easily */
CLJAtoms CLJDelta::changedAtoms() const
{
    //work out which atoms have changed and which ones haven't...
    QVarLengthArray<CLJAtom> changed_atoms;
    
    for (int i=0; i<qMin(old_atoms.count(),new_atoms.count()); ++i)
    {
        CLJAtom old_atom = old_atoms.at(i);
        CLJAtom new_atom = new_atoms.at(i);
        
        if (old_atom != new_atom)
        {
            if (not old_atom.isDummy())
                changed_atoms.append( old_atom.negate() );
            
            if (not new_atom.isDummy())
                changed_atoms.append( new_atom );
        }
    }
    
    if (old_atoms.count() > new_atoms.count())
    {
        for (int i=new_atoms.count(); i<old_atoms.count(); ++i)
        {
            CLJAtom old_atom = old_atoms.at(i);
            
            if (not old_atom.isDummy())
                changed_atoms.append( old_atom.negate() );
        }
    }
    else if (new_atoms.count() > old_atoms.count())
    {
        for (int i=old_atoms.count(); i<new_atoms.count(); ++i)
        {
            CLJAtom new_atom = new_atoms.at(i);
            
            if (not new_atom.isDummy())
                changed_atoms.append( new_atom );
        }
    }

    return CLJAtoms(changed_atoms.constData(), changed_atoms.count());
}

/** Merge together the changed atoms from the 'n' deltas from the passed array 
    into a single changed atoms object. The resulting set of changed atoms will
    thus be able to be used to calculate energy changes from a lot of changed
    atoms */
CLJAtoms CLJDelta::mergeChanged(const CLJDelta *deltas, int n)
{
    if (n == 0)
        return CLJAtoms();
    
    else if (n == 1)
        return deltas[0].changedAtoms();

    else
    {
        //work out which atoms have changed and which ones haven't...
        QVarLengthArray<CLJAtom> changed_atoms;

        for (int l=0; l<n; ++l)
        {
            CLJAtoms old_atoms = deltas[l].oldAtoms();
            CLJAtoms new_atoms = deltas[l].newAtoms();
        
            for (int i=0; i<qMin(old_atoms.count(),new_atoms.count()); ++i)
            {
                CLJAtom old_atom = old_atoms.at(i);
                CLJAtom new_atom = new_atoms.at(i);
                
                if (old_atom != new_atom)
                {
                    if (not old_atom.isDummy())
                        changed_atoms.append( old_atom.negate() );
                    
                    if (not new_atom.isDummy())
                        changed_atoms.append( new_atom );
                }
            }
            
            if (old_atoms.count() > new_atoms.count())
            {
                for (int i=new_atoms.count(); i<old_atoms.count(); ++i)
                {
                    CLJAtom old_atom = old_atoms.at(i);
                    
                    if (not old_atom.isDummy())
                        changed_atoms.append( old_atom.negate() );
                }
            }
            else if (new_atoms.count() > old_atoms.count())
            {
                for (int i=old_atoms.count(); i<new_atoms.count(); ++i)
                {
                    CLJAtom new_atom = new_atoms.at(i);
                    
                    if (not new_atom.isDummy())
                        changed_atoms.append( new_atom );
                }
            }
        }
        
        return CLJAtoms(changed_atoms.constData(), changed_atoms.count());
    }
}

/** Merge together the changed atoms from the passed deltas
    into a single changed atoms object. The resulting set of changed atoms will
    thus be able to be used to calculate energy changes from a lot of changed
    atoms */
CLJAtoms CLJDelta::mergeChanged(const QVector<CLJDelta> &deltas)
{
    return mergeChanged(deltas.constData(), deltas.count());
}

/** Merge together the new atoms from the 'n' deltas from the passed array
    into a single new atoms object. The resulting set of new atoms will
    thus be able to be used to calculate energy changes from a lot of changed
    atoms */
CLJAtoms CLJDelta::mergeNew(const CLJDelta *deltas, int n)
{
    if (n == 0)
        return CLJAtoms();
    else
    {
        //work out which atoms have changed and which ones haven't...
        QVarLengthArray<CLJAtom> changed_atoms;

        for (int l=0; l<n; ++l)
        {
            CLJAtoms old_atoms = deltas[l].oldAtoms();
            CLJAtoms new_atoms = deltas[l].newAtoms();
        
            for (int i=0; i<qMin(old_atoms.count(),new_atoms.count()); ++i)
            {
                CLJAtom old_atom = old_atoms.at(i);
                CLJAtom new_atom = new_atoms.at(i);
                
                if (old_atom != new_atom)
                {
                    if (not new_atom.isDummy())
                        changed_atoms.append( new_atom );
                }
            }
            
            if (new_atoms.count() > old_atoms.count())
            {
                for (int i=old_atoms.count(); i<new_atoms.count(); ++i)
                {
                    CLJAtom new_atom = new_atoms.at(i);
                    
                    if (not new_atom.isDummy())
                        changed_atoms.append( new_atom );
                }
            }
        }
        
        return CLJAtoms(changed_atoms.constData(), changed_atoms.count());
    }
}

/** Merge together the new atoms from deltas
    into a single new atoms object. The resulting set of new atoms will
    thus be able to be used to calculate energy changes from a lot of changed
    atoms */
CLJAtoms CLJDelta::mergeNew(const QVector<CLJDelta> &deltas)
{
    return mergeNew(deltas.constData(), deltas.count());
}

/** Merge together the old atoms from the 'n' deltas from the passed array
    into a single old atoms object. The resulting set of old atoms will
    thus be able to be used to calculate energy changes from a lot of changed
    atoms */
CLJAtoms CLJDelta::mergeOld(const CLJDelta *deltas, int n)
{
    if (n == 0)
        return CLJAtoms();
    
    else
    {
        //work out which atoms have changed and which ones haven't...
        QVarLengthArray<CLJAtom> changed_atoms;

        for (int l=0; l<n; ++l)
        {
            CLJAtoms old_atoms = deltas[l].oldAtoms();
            CLJAtoms new_atoms = deltas[l].newAtoms();
        
            for (int i=0; i<qMin(old_atoms.count(),new_atoms.count()); ++i)
            {
                CLJAtom old_atom = old_atoms.at(i);
                CLJAtom new_atom = new_atoms.at(i);
                
                if (old_atom != new_atom)
                {
                    if (not old_atom.isDummy())
                        changed_atoms.append( old_atom );
                }
            }
            
            if (old_atoms.count() > new_atoms.count())
            {
                for (int i=new_atoms.count(); i<old_atoms.count(); ++i)
                {
                    CLJAtom old_atom = old_atoms.at(i);
                    
                    if (not old_atom.isDummy())
                        changed_atoms.append( old_atom );
                }
            }
        }
        
        return CLJAtoms(changed_atoms.constData(), changed_atoms.count());
    }
}

/** Merge together the old atoms from deltas
    into a single old atoms object. The resulting set of old atoms will
    thus be able to be used to calculate energy changes from a lot of changed
    atoms */
CLJAtoms CLJDelta::mergeOld(const QVector<CLJDelta> &deltas)
{
    return mergeOld(deltas.constData(), deltas.count());
}

/** Merge together the changed atoms from the 'n' deltas from the passed array 
    into a tuple of the changed, old and new atoms. The resulting set of changed atoms will
    thus be able to be used to calculate energy changes from a lot of changed
    atoms */
tuple<CLJAtoms,CLJAtoms,CLJAtoms> CLJDelta::merge(const CLJDelta *deltas, int n)
{
    if (n == 0)
        return CLJAtoms();
    else
    {
        //work out which atoms have changed and which ones haven't...
        QVarLengthArray<CLJAtom> changed_atoms;
        QVarLengthArray<CLJAtom> all_old_atoms;
        QVarLengthArray<CLJAtom> all_new_atoms;

        for (int l=0; l<n; ++l)
        {
            CLJAtoms old_atoms = deltas[l].oldAtoms();
            CLJAtoms new_atoms = deltas[l].newAtoms();
        
            for (int i=0; i<qMin(old_atoms.count(),new_atoms.count()); ++i)
            {
                CLJAtom old_atom = old_atoms.at(i);
                CLJAtom new_atom = new_atoms.at(i);
                
                if (old_atom != new_atom)
                {
                    if (not old_atom.isDummy())
                    {
                        changed_atoms.append( old_atom.negate() );
                        all_old_atoms.append( old_atom );
                    }
                    
                    if (not new_atom.isDummy())
                    {
                        changed_atoms.append( new_atom );
                        all_new_atoms.append( new_atom );
                    }
                }
            }
            
            if (old_atoms.count() > new_atoms.count())
            {
                for (int i=new_atoms.count(); i<old_atoms.count(); ++i)
                {
                    CLJAtom old_atom = old_atoms.at(i);
                    
                    if (not old_atom.isDummy())
                    {
                        changed_atoms.append( old_atom.negate() );
                        all_old_atoms.append( old_atom );
                    }
                }
            }
            else if (new_atoms.count() > old_atoms.count())
            {
                for (int i=old_atoms.count(); i<new_atoms.count(); ++i)
                {
                    CLJAtom new_atom = new_atoms.at(i);
                    
                    if (not new_atom.isDummy())
                    {
                        changed_atoms.append( new_atom );
                        all_new_atoms.append( new_atom );
                    }
                }
            }
        }
        
        return tuple<CLJAtoms,CLJAtoms,CLJAtoms>(
                        CLJAtoms(changed_atoms.constData(), changed_atoms.count()),
                        CLJAtoms(all_old_atoms.constData(), all_old_atoms.count()),
                        CLJAtoms(all_new_atoms.constData(), all_new_atoms.count()) );
    }
}

/** Merge together the changed atoms from the passed deltas
    into a tuple of changed, old and new atoms. The resulting set of changed atoms will
    thus be able to be used to calculate energy changes from a lot of changed
    atoms */
tuple<CLJAtoms,CLJAtoms,CLJAtoms> CLJDelta::merge(const QVector<CLJDelta> &deltas)
{
    return merge(deltas.constData(), deltas.count());
}
