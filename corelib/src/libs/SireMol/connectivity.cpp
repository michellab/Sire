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

#include <QDataStream>
#include <QElapsedTimer>

#include <boost/assert.hpp>

#include "atomselection.h"
#include "connectivity.h"
#include "moleculedata.h"
#include "moleculeinfodata.h"
#include "moleculeview.h"
#include "atommatcher.h"

#include "angleid.h"
#include "bondid.h"
#include "dihedralid.h"
#include "improperid.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <QDebug>

using namespace SireStream;
using namespace SireMol;
using namespace SireBase;

/////////
///////// Implementation of ConnectivityBase
/////////

static const RegisterMetaType<ConnectivityBase> r_conbase(MAGIC_ONLY,
                                                          "SireMol::ConnectivityBase");

/** Serialise ConnectivityBase */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, 
                                       const ConnectivityBase &conbase)
{
    writeHeader(ds, r_conbase, 2);

    SharedDataStream sds(ds);

    sds << conbase.connected_atoms << conbase.connected_res
        << conbase.d << static_cast<const MolViewProperty&>(conbase);

    return ds;
}

/** Deserialise a MoleculeBonds */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, 
                                       ConnectivityBase &conbase)
{
    VersionID v = readHeader(ds, r_conbase);

    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> conbase.connected_atoms >> conbase.connected_res
            >> conbase.d
            >> static_cast<MolViewProperty&>(conbase);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> conbase.connected_atoms >> conbase.connected_res
            >> conbase.d
            >> static_cast<Property&>(conbase);
    }
    else
        throw version_error(v, "1,2", r_conbase, CODELOC);

    return ds;
}

/** Null constructor */
ConnectivityBase::ConnectivityBase()
                 : MolViewProperty()
{
    MoleculeInfo null_info;
    d = null_info.data();
}

const MoleculeInfoData& ConnectivityBase::info() const
{
    return *d;
}

/** Construct the connectivity for molecule described by 
    the passed info object */
ConnectivityBase::ConnectivityBase(const MoleculeData &moldata)
                 : MolViewProperty(), d(moldata.info())
{
    if (info().nAtoms() > 0)
    {
        connected_atoms.resize(info().nAtoms());
        connected_atoms.squeeze();
    }
    
    if (info().nResidues() > 0)
    {
        connected_res.resize(info().nResidues());
        connected_res.squeeze();
    }
}
    
/** Copy constructor */
ConnectivityBase::ConnectivityBase(const ConnectivityBase &other)
                 : MolViewProperty(other),
                   connected_atoms(other.connected_atoms),
                   connected_res(other.connected_res),
                   d(other.d)
{}

/** Destructor */
ConnectivityBase::~ConnectivityBase()
{}

/** Copy assignment operator */
ConnectivityBase& ConnectivityBase::operator=(const ConnectivityBase &other)
{
    if (this != &other)
    {
        connected_atoms = other.connected_atoms;
        connected_res = other.connected_res;
        d = other.d;
    }
    
    return *this;
}

/** Comparison operator */
bool ConnectivityBase::operator==(const ConnectivityBase &other) const
{
    return (d == other.d or *d == *(other.d)) and 
           connected_atoms == other.connected_atoms;
}

/** Comparison operator */
bool ConnectivityBase::operator!=(const ConnectivityBase &other) const
{
    return (d != other.d and *d != *(other.d)) or
           connected_atoms != other.connected_atoms;
}

bool ConnectivityBase::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return molinfo == this->info();
}

PropertyPtr ConnectivityBase::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                      const AtomMatcher &atommatcher) const
{
    try
    {
        if (not atommatcher.changesOrder(this->info(), molinfo))
        {
            //the order of the atoms remains the same - this means that the 
            //AtomIdx indicies are still valid
            Connectivity ret;
            ret.connected_atoms = connected_atoms;
            ret.connected_res = connected_res;
            ret.d = molinfo;
            return ret;
        }

        QHash<AtomIdx,AtomIdx> matched_atoms = atommatcher.match(this->info(), molinfo);

        ConnectivityEditor editor;
        editor.d = molinfo;
        editor.connected_atoms = QVector< QSet<AtomIdx> >( molinfo.nAtoms() );
        editor.connected_res = QVector< QSet<ResIdx> >( molinfo.nResidues() );

        for (int i=0; i<connected_atoms.count(); ++i)
        {
            AtomIdx old_idx(i);
            
            AtomIdx new_idx = matched_atoms.value(old_idx, AtomIdx(-1));
            
            if (new_idx != -1)
            {
                foreach (AtomIdx old_bond, this->connectionsTo(old_idx))
                {
                    AtomIdx new_bond = matched_atoms.value(old_bond, AtomIdx(-1));
                    
                    if (new_bond != -1)
                    {
                        if (new_bond > new_idx)
                            editor.connect(new_idx, new_bond);
                    }
                }
            }
        }

        return editor.commit();
    }
    catch(const SireError::exception &e)
    {
        qDebug() << e.toString();
        throw;
        return Connectivity();
    }
}

static QString atomString(const MoleculeInfoData &molinfo, AtomIdx atom)
{
    if (molinfo.isWithinResidue(atom))
        return QString("%1:%2:%3")
                    .arg( molinfo.name(atom) )
                    .arg( molinfo.name( molinfo.parentResidue(atom) ) )
                    .arg( molinfo.number( molinfo.parentResidue(atom) ) );
    else
        return QString("%1:%2").arg(molinfo.name(atom)).arg(atom);
}

QString ConnectivityBase::toString() const
{
    QStringList lines;
    
    lines.append( QObject::tr("Connectivity: nConnections() == %1.")
                        .arg(this->nConnections()) );

    if (nConnections() == 0)
        return lines.at(0);
                        
    if (not connected_res.isEmpty())
    {
        lines.append( QObject::tr("Connected residues:") );
        
        for (int i=0; i<connected_res.count(); ++i)
        {
            QStringList resnums;
            
            foreach (ResIdx j, connected_res.at(i))
            {
                resnums.append( QString("%1:%2")
                                 .arg(d->name(j))
                                 .arg(d->number(j)) );
            }
            
            if (not connected_res.at(i).isEmpty())
                lines.append( QObject::tr("  * Residue %1:%2 bonded to %3.")
                        .arg(d->name(ResIdx(i)))
                        .arg(d->number(ResIdx(i)))
                        .arg(resnums.join(" ")) );
        }
    }
    
    if (not connected_atoms.isEmpty())
    {
        lines.append( QObject::tr("Connected atoms:") );
        
        for (int i=0; i<connected_atoms.count(); ++i)
        {
            QStringList atoms;
            
            foreach (AtomIdx j, connected_atoms.at(i))
            {
                atoms.append( ::atomString(info(),j) );
            }
            
            if (not connected_atoms.at(i).isEmpty())
                lines.append( QObject::tr("  * Atom %1 bonded to %2.")
                        .arg( ::atomString(info(),AtomIdx(i)) )
                        .arg( atoms.join(" ") ) );
        }
    }
    
    return lines.join("\n");
}

/** Return the total number of connections between atoms
    in this connectivity object */
int ConnectivityBase::nConnections() const
{
    int n = 0;
    
    const QSet<AtomIdx> *connected_atoms_array = connected_atoms.constData();
    int nats = connected_atoms.count();
    
    for (int i=0; i<nats; ++i)
    {
        foreach (AtomIdx j, connected_atoms_array[i])
        {
            if (i < j)
            {
                //only count connections when the i atom is less
                //than j - this avoids double counting the bond
                // i-j and j-i
                ++n;
            }
        }
    }
    
    return n;
}

/** Return the indicies of atoms connected to the atom at index 'atomidx'.
    This returns an empty set if there are no atoms connected to 
    this atom
    
    \throw SireError::invalid_index
*/
const QSet<AtomIdx>& ConnectivityBase::connectionsTo(AtomIdx atomidx) const
{
    return connected_atoms.constData()[ atomidx.map(connected_atoms.count()) ];
}

/** Return the indicies of atoms connected to the atom identified
    by 'resid' - this returns an empty set if there are no connections
    to this atom
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
const QSet<AtomIdx>& ConnectivityBase::connectionsTo(const AtomID &atomid) const
{
    return connected_atoms.constData()[ info().atomIdx(atomid) ];
}

/** Return the indicies of the residues connected to the residue at 
    index 'residx'. This returns an empty set if there are no residues
    connected to this residue
    
    \throw SireError::invalid_index
*/
const QSet<ResIdx>& ConnectivityBase::connectionsTo(ResIdx residx) const
{
    return connected_res.constData()[ residx.map(connected_res.count()) ];
}

/** Return the indicies of the residues connectd to the residue
    identified by 'resid'. This returns an empty set if there are
    no residues connected to this residue.
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
const QSet<ResIdx>& ConnectivityBase::connectionsTo(const ResID &resid) const
{
    return connected_res.constData()[ info().resIdx(resid) ];
}

/** Return the number of connections to the atom at index 'atomidx' 

    \throw SireError::index_error
*/
int ConnectivityBase::nConnections(AtomIdx atomidx) const
{
    return this->connectionsTo(atomidx).count();
}

/** Return the number of connections to the atom with ID 'atomid' 

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(const AtomID &atomid) const
{
    return this->connectionsTo(atomid).count();
}

/** Return the number of connections to the residue at index 'residx'

    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(ResIdx residx) const
{
    return this->connectionsTo(residx).count();
}

/** Return the number of connections to the residue identified
    by 'resid'
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(const ResID &resid) const
{
    return this->connectionsTo(resid).count();
}

/** Return the number of atom connections between the residues at
    indicies 'res0' and 'res1' 
    
    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(ResIdx res0, ResIdx res1) const
{
    if (not this->areConnected(res0, res1))
        return 0;
        
    int n = 0;
    
    QList<AtomIdx> atoms0 = info().getAtomsIn(res0);
    QList<AtomIdx> atoms1 = info().getAtomsIn(res1);
    
    foreach (AtomIdx atom0, atoms0)
    {
        const QSet<AtomIdx> &connected = this->connectionsTo(atom0);

        foreach (AtomIdx atom1, atoms1)
        {
            if (connected.contains(atom1))
                ++n;
        }
    }
    
    return n;
}

/** Return the number of atom connections between the residues
    identified by 'res0' and 'res1'
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(const ResID &res0, const ResID &res1) const
{
    return this->nConnections( info().resIdx(res0), info().resIdx(res1) );
}

/** Return whether or not the atoms at indicies 'atom0' and 'atom1'
    are connected
    
    \throw SireError::invalid_index
*/
bool ConnectivityBase::areConnected(AtomIdx atom0, AtomIdx atom1) const
{
    return this->connectionsTo(atom0).contains( 
                                AtomIdx(atom1.map(connected_atoms.count())) );
}

/** Return whether or not the atoms identified by 'atom0' and 'atom1'
    are connected
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ConnectivityBase::areConnected(const AtomID &atom0, const AtomID &atom1) const
{
    return this->connectionsTo(atom0).contains( info().atomIdx(atom1) );
}

/** Return whether or not the two atoms are bonded together */
bool ConnectivityBase::areBonded(AtomIdx atom0, AtomIdx atom1) const
{
    return areConnected(atom0, atom1);
}

/** Return whether or not the two atoms are angled together */
bool ConnectivityBase::areAngled(AtomIdx atom0, AtomIdx atom2) const
{
    atom0 = d.read().atomIdx(atom0);
    atom2 = d.read().atomIdx(atom2);

    if (atom0 == atom2)
        return false;

    foreach (const AtomIdx &atom1, connected_atoms[atom0.value()])
    {
        if (connected_atoms[atom2.value()].contains(atom1))
            return true;
    }
    
    return false;
}

/** Return whether or not the two atoms are dihedraled together */
bool ConnectivityBase::areDihedraled(AtomIdx atom0, AtomIdx atom3) const
{
    atom0 = d.read().atomIdx(atom0);
    atom3 = d.read().atomIdx(atom3);
    
    if (atom0 == atom3)
        return false;
    
    foreach (const AtomIdx &atom1, connected_atoms[atom0.value()])
    {
        if (atom1.value() != atom3.value())
        {
            foreach (const AtomIdx &atom2, connected_atoms[atom3.value()])
            {
                if (atom2.value() != atom0.value())
                {
                    if (connected_atoms[atom1.value()].contains(atom2))
                        return true;
                }
            }
        }
    }
    
    return false;
}

/** Return whether or not the two atoms are bonded together */
bool ConnectivityBase::areBonded(const AtomID &atom0, const AtomID &atom1) const
{
    return areBonded(d.read().atomIdx(atom0), d.read().atomIdx(atom1));
}

/** Return whether or not the two atoms are angled together */
bool ConnectivityBase::areAngled(const AtomID &atom0, const AtomID &atom2) const
{
    return areAngled(d.read().atomIdx(atom0), d.read().atomIdx(atom2));
}

/** Return whether or not the two atoms are bonded together */
bool ConnectivityBase::areDihedraled(const AtomID &atom0, const AtomID &atom3) const
{
    return areDihedraled(d.read().atomIdx(atom0), d.read().atomIdx(atom3));
}

/** Return whether or not the residues at indicies 'res0' and 'res1' 
    are connected
    
    \throw SireError::invalid_index
*/
bool ConnectivityBase::areConnected(ResIdx res0, ResIdx res1) const
{
    return this->connectionsTo(res0).contains( 
                                        ResIdx(res1.map(connected_res.count())) );
}

/** Return whether the residues identified by 'res0' and 'res1' are connected */
bool ConnectivityBase::areConnected(const ResID &res0, const ResID &res1) const
{
    return this->connectionsTo(res0).contains( info().resIdx(res1) );
}

/** Non-checking version of Connectivity::connectedTo(AtomIdx) */
const QSet<AtomIdx>& ConnectivityBase::_pvt_connectedTo(AtomIdx atom) const
{
    return connected_atoms.constData()[atom];
}

/** Internal recursive function used to find all paths between two atoms */
QList< QList<AtomIdx> > ConnectivityBase::_pvt_findPaths(AtomIdx cursor, const AtomIdx end_atom,
                                                         QSet<AtomIdx> &done) const
{
    //create the list containing all paths from the cursor atom to the end atom
    QList< QList<AtomIdx> > all_paths;

    if (not done.contains(cursor))
    {
        //we have not traced through this atom before...
        done.insert(cursor);

        //loop through all atoms bonded to the cursor
        foreach (const AtomIdx &bonded_to_cursor, this->_pvt_connectedTo(cursor))
        {
            if (bonded_to_cursor == end_atom)
            {
                //we have found a path to the end atom. Return a single list containing
                //cursor and end_atom, so that the functions that call this can then add their
                //atoms to create all of the paths
                QList<AtomIdx> path;
                path.append(cursor);
                path.append(end_atom);
                all_paths.append(path);
            }
            else
            {
                QSet<AtomIdx> new_done = done;
            
                QList< QList<AtomIdx> > paths = this->_pvt_findPaths(bonded_to_cursor,
                                                                     end_atom, new_done);
        
                if (not paths.isEmpty())
                {
                    for (QList< QList<AtomIdx> >::iterator it = paths.begin();
                         it != paths.end();
                         ++it)
                    {
                        (*it).prepend(cursor);
                        all_paths.append( *it );
                    }
                }
            }
        }
    }
    
    return all_paths;
}

/** Return all possible bonded paths between two atoms. This returns an empty
    list if there are no bonded paths between the two atoms */
QList< QList<AtomIdx> > ConnectivityBase::findPaths(AtomIdx atom0, AtomIdx atom1) const
{
    atom0 = atom0.map( d->nAtoms() );
    atom1 = atom1.map( d->nAtoms() );
    
    if (atom0 == atom1)
        return QList< QList<AtomIdx> >();
    
    QSet<AtomIdx> done;
    done.reserve(d->nAtoms());
    
    return this->_pvt_findPaths(atom0, atom1, done);
}

/** Find the shortest bonded path between two atoms. This returns an empty
    list if there is no bonded path between these two atoms */
QList<AtomIdx> ConnectivityBase::findPath(AtomIdx atom0, AtomIdx atom1) const
{
    QList< QList<AtomIdx> > paths = findPaths(atom0, atom1);
    
    QList<AtomIdx> shortest;
    
    foreach (const QList<AtomIdx> &path, paths)
    {
        if (shortest.isEmpty())
            shortest = path;
        
        else if (shortest.count() > path.count())
            shortest = path;
    }
    
    return shortest;
}

/** Return all possible bonded paths between two atoms. This returns an empty
    list if there are no bonded paths between the two atoms */
QList<AtomIdx> ConnectivityBase::findPath(const AtomID &atom0, const AtomID &atom1) const
{
    return this->findPath( d->atomIdx(atom0), d->atomIdx(atom1) );
}

/** Find the shortest bonded path between two atoms. This returns an empty
    list if there is no bonded path between these two atoms */
QList< QList<AtomIdx> > ConnectivityBase::findPaths(const AtomID &atom0, const AtomID &atom1) const
{
    return this->findPaths( d->atomIdx(atom0), d->atomIdx(atom1) );
}

/** This function returns whether or not the two passed atoms are part of
    the same ring */
bool ConnectivityBase::inRing(AtomIdx atom0, AtomIdx atom1) const
{
    QList< QList<AtomIdx> > paths = findPaths(atom0, atom1);
    
    //if there is more than one path between the atoms then they must
    //be part of a ring
    return (paths.count() > 1);
}

/** This function returns whether or not the three passed atoms are all part of
    the same ring */
bool ConnectivityBase::inRing(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2) const
{
    QList< QList<AtomIdx> > paths = findPaths(atom0, atom2);
    
    atom1 = atom1.map(d->nAtoms());
    
    if (paths.count() > 1)
    {
        //the three are part of the same ring if any of the paths contains
        //atom1
        foreach (const QList<AtomIdx> &path, paths)
        {
            if (path.contains(atom1))
                return true;
        }
    }
    
    return false;
}

/** This function returns whether or not the four passed atoms are part of
    the same ring */
bool ConnectivityBase::inRing(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3) const
{
    QList< QList<AtomIdx> > paths = findPaths(atom0, atom3);
    
    atom1 = atom1.map(d->nAtoms());
    atom2 = atom2.map(d->nAtoms());
    
    if (paths.count() > 1)
    {
        //the four are part of the same ring if atom1 and atom2 are contained
        //in any of the paths
        bool have_atom1 = false;
        bool have_atom2 = false;
        
        foreach (const QList<AtomIdx> &path, paths)
        {
            if (path.contains(atom1))
            {
                have_atom1 = true;
                if (have_atom2)
                    return this->inRing(atom1,atom2);
            }
            
            if (path.contains(atom2))
            {
                have_atom2 = true;
                if (have_atom1)
                    return this->inRing(atom1,atom2);
            }
        }
    }
    
    return false;
}

/** This function returns whether or not the two passed atoms are part of
    the same ring */
bool ConnectivityBase::inRing(const AtomID &atom0, const AtomID &atom1) const
{
    return this->inRing( info().atomIdx(atom0), info().atomIdx(atom1) );
}

/** This function returns whether or not the three passed atoms are all part of
    the same ring */
bool ConnectivityBase::inRing(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2) const
{
    return this->inRing( info().atomIdx(atom0), info().atomIdx(atom1),
                         info().atomIdx(atom2) );
}

/** This function returns whether or not the two passed atoms are part of
    the same ring */
bool ConnectivityBase::inRing(const AtomID &atom0, const AtomID &atom1,
                              const AtomID &atom2, const AtomID &atom3) const
{
    return this->inRing( info().atomIdx(atom0), info().atomIdx(atom1),
                         info().atomIdx(atom2), info().atomIdx(atom3) );
}

/** This function returns whether or not the two atoms in the passed bond
    are both part of the same ring */
bool ConnectivityBase::inRing(const BondID &bond) const
{
    return this->inRing(bond.atom0(), bond.atom1());
}

/** This function returns whether or not the three atoms in the passed angle
    are all part of the same ring */
bool ConnectivityBase::inRing(const AngleID &angle) const
{
    return this->inRing(angle.atom0(), angle.atom1(), angle.atom2());
}

/** This function returns whether or not the four atoms in the passed dihedral
    are all part of the same ring */
bool ConnectivityBase::inRing(const DihedralID &dihedral) const
{
    return this->inRing(dihedral.atom0(), dihedral.atom1(),
                        dihedral.atom2(), dihedral.atom3());
}

/** This is a recursive function that traces all atoms that can trace their bonding
    to 'strt' but that don't go through any of the atoms in 'root',
    and to add those atoms to 'group'.

    \throw SireMol::ring_error
*/
void ConnectivityBase::traceRoute(AtomIdx start, QSet<AtomIdx> &root,
                                  QSet<AtomIdx> &group) const
{
    //add this atom to the group
    group.insert(start);
    root.insert(start);

    //now see if any of its bonded atoms need to be added
    const QSet<AtomIdx> &bonded_atoms = this->_pvt_connectedTo(start);

    //loop over every bond that involves the 'start' atom
    for (QSet<AtomIdx>::const_iterator it = bonded_atoms.constBegin();
         it != bonded_atoms.constEnd();
         ++it)
    {
        //if this is a root atom then ignore it, as we don't
        //want to build a path through this atom
        if (root.contains(*it))
        {
            continue;
        }
        //has this atom already been selected?
        else if (group.contains(*it))
        {
            //yes, this atom is already included!
            continue;
        }
        else
        {
            //now we can trace the atoms from the 'other' atom...
            this->traceRoute(*it, root, group);
        }
    }

    //ok - we have added all of the atoms that are connected to this atom. We
    //have finished with this atom, so we can return.
    return;
}
    
/** This is a recursive function that traces all atoms that can trace their bonding
    to 'strt' but that don't go through any of the atoms in 'root',
    and to add those atoms to 'group'. This traces only atoms that are contained
    in 'selected_atoms'

    \throw SireMol::ring_error
*/
void ConnectivityBase::traceRoute(const AtomSelection &selected_atoms,
                                  AtomIdx start, QSet<AtomIdx> &root,
                                  QSet<AtomIdx> &group) const
{
    if (not selected_atoms.selected(start))
        //the atom is not one of the selected atoms
        return;

    //add this atom to the group
    group.insert(start);
    root.insert(start);

    //now see if any of its bonded atoms need to be added
    const QSet<AtomIdx> &bonded_atoms = this->_pvt_connectedTo(start);

    //loop over every bond that involves the 'start' atom
    for (QSet<AtomIdx>::const_iterator it = bonded_atoms.constBegin();
         it != bonded_atoms.constEnd();
         ++it)
    {
        //if this is a root atom then ignore it, as we don't
        //want to move backwards!
        if (root.contains(*it))
            continue;
    
        //has this atom or residue already been selected?
        else if (group.contains(*it))
            //yes, this atom or residue is already included!
            continue;

        else
            //now we can trace the atoms from the 'other' atom...
            this->traceRoute(selected_atoms, *it, root, group);
    }

    //ok - we have added all of the atoms that are connected to this atom. We
    //have finished with this atom, so we can return.
    return;
}

/** Return the two AtomSelection objects corresponding to the atoms
    selected in 'group0' and 'group1' */
tuple<AtomSelection,AtomSelection>
ConnectivityBase::selectGroups(const QSet<AtomIdx> &group0,
                               const QSet<AtomIdx> &group1) const
{
    AtomSelection grp0(*d);
    AtomSelection grp1(*d);

    tuple<AtomSelection,AtomSelection> groups( grp0.selectOnly(group0),
                                               grp1.selectOnly(group1) );
    
    return groups;
}

/** Split this molecule into two parts about the atoms
    atom0 and atom1. For example;
    
    For example;

          C1--C2--C3--C4--C5--C6--C7--C8

    Splitting C3 and C4 would result in two groups, {C1,C2,C3} and {C4,C5,C6,C7,C8}

          C1\
          C2-C4--C5
          C3/

    Splitting C4 and C5 would result in two groups, {C1,C2,C3,C4} and {C5}

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/   
tuple<AtomSelection,AtomSelection> 
ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1) const
{
    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    //map the atoms
    int nats = d->nAtoms();
    atom0 = atom0.map(nats);
    atom1 = atom1.map(nats);

    if (atom0 == atom1)
        throw SireMol::ring_error( QObject::tr(
            "You cannot split a molecule into two parts using the same atom! (%1).")
                .arg( ::atomString(info(),atom0) ), CODELOC );

    //make sure that there is sufficient space for the
    //selections - this prevents mallocs while tracing
    //the bonds
    group0.reserve(nats);
    group1.reserve(nats);
    root0.reserve(nats);
    root1.reserve(nats);
    
    //add the two atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom1);
    root0.insert(atom0);
    root0.insert(atom1);
    root1.insert(atom1);
    root1.insert(atom0);
    
    //add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if (bonded_atom != atom1)
        {
            this->traceRoute(bonded_atom, root0, group0);
        }
    }
    
    //remove atom1 from group0, in case it was found as part of a ring
    group0.remove(atom1);
    
    //now add the atoms bonded to atom1 to group1
    bool has_rings = false;
    
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom1))
    {
        if (bonded_atom != atom0)
        {
            if (group0.contains(bonded_atom))
            {
                has_rings = true;
            }

            this->traceRoute(bonded_atom, root1, group1);
        }
    }
    
    group1.remove(atom0);
    
    //if there is any overlap in the two sets then that means that
    //the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (bonded_atom != atom1 and group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                //qDebug() << "DISCONNECTING(0)" << d->name(bonded_atom) << d->name(atom0);
            }
        }
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom1))
        {
            if (bonded_atom != atom0 and group0.contains(bonded_atom))
            {
                editor.disconnect(atom1, bonded_atom);
                //qDebug() << "DISCONNECTING(1)" << d->name(bonded_atom) << d->name(atom1);
            }
        }
        
        //release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();
        
        //split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;
        
        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( d->name(atom) );
        }
        
        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( d->name(atom) );
        }
    
        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/
        
        return this->selectGroups(group0, group1);
    }
}

/** Split the molecule into two parts about the bond between atom0 and atom1.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/   
tuple<AtomSelection,AtomSelection> 
ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1) const
{
    return this->split( d->atomIdx(atom0), d->atomIdx(atom1) );
}

/** Split the molecule into two parts about the bond 'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/   
tuple<AtomSelection,AtomSelection> 
ConnectivityBase::split(const BondID &bond) const
{
    return this->split( bond.atom0(), bond.atom1() );
}

/** Split the selected atoms of this molecule into two parts about the atoms
    atom0 and atom1. For example;
    
    For example;

          C1--C2--C3--C4--C5--C6--C7--C8

    Splitting C3 and C4 would result in two groups, {C1,C2,C3} and {C4,C5,C6,C7,C8}

          C1\
          C2-C4--C5
          C3/

    Splitting C4 and C5 would result in two groups, {C1,C2,C3,C4} and {C5}

    However splitting C1 and C5 would add a bond between C1 and C5. This would mean
    than C1-C4-C5 would form a ring, so an exception would be thrown.

    Note that both atom0 and atom1 *must* be selected as part of 'selected_atoms'
    or else a missing_atom exception will be thrown.

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/   
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1, 
                        const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*d);
    
    if (selected_atoms.selectedAll())
        return this->split(atom0, atom1);
 
    selected_atoms.assertSelected(atom0);
    selected_atoms.assertSelected(atom1);

    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    //make sure that there is sufficient space for the
    //selections - this prevents mallocs while tracing
    //the bonds
    group0.reserve(selected_atoms.nSelected());
    group1.reserve(selected_atoms.nSelected());
    root0.reserve(selected_atoms.nSelected());
    root1.reserve(selected_atoms.nSelected());
    
    //map the atoms
    atom0 = atom0.map(d->nAtoms());
    atom1 = atom1.map(d->nAtoms());

    if (atom0 == atom1)
        throw SireMol::ring_error( QObject::tr(
            "You cannot split a molecule into two parts using the same atom! (%1).")
                .arg( ::atomString(info(),atom0) ), CODELOC );

    //add the two atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom1);
    root0.insert(atom0);
    root0.insert(atom1);
    root1.insert(atom1);
    root1.insert(atom0);
    
    //add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if ( (bonded_atom != atom1) and 
             selected_atoms.selected(bonded_atom) )
        {
            this->traceRoute(selected_atoms, bonded_atom, root0, group0);
        }
    }
    
    //now add the atoms bonded to atom1 to group1
    bool has_rings = false;
    
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom1))
    {
        if ( (bonded_atom != atom0) and
             selected_atoms.selected(bonded_atom) )
        {
            if (group0.contains(bonded_atom))
                has_rings = true;
                    
            this->traceRoute(selected_atoms, bonded_atom, root1, group1);
        }
    }
    
    //if there is any overlap in the two sets then that means that
    //the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (bonded_atom != atom1 and group1.contains(bonded_atom))
                editor.disconnect(atom0, bonded_atom);
        }
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom1))
        {
            if (bonded_atom != atom0 and group0.contains(bonded_atom))
                editor.disconnect(atom1, bonded_atom);
        }
        
        //release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();
        
        //split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1);
    }
    else
        return this->selectGroups(group0, group1);
}

/** Split the selected atoms of this molecule about the atoms 
    'atom0' and 'atom1'

    \throw SireMol::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/   

tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1,
                        const AtomSelection &selected_atoms) const
{
    return this->split( d->atomIdx(atom0), d->atomIdx(atom1),
                        selected_atoms );
}

/** Split the selected atoms of this molecule into two parts
    about the bond 'bond'

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/   
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const BondID &bond, const AtomSelection &selected_atoms) const
{
    return this->split( bond.atom0(), bond.atom1(), selected_atoms );
}

/** Split this molecule into three parts about the atoms
    'atom0', 'atom1' and 'atom2'.

    For example;

      C1   C3--C8
        \ /
        C2
       /  \
     C4   C5--C6-C7

    Splitting C5,C2,C3 would return {C5,C6,C7} in one group, and {C3,C8} 
    in the other group.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/   
tuple<AtomSelection,AtomSelection> 
ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2) const
{
    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    //map the atoms
    int nats = d->nAtoms();
    atom0 = atom0.map(nats);
    atom1 = atom1.map(nats);
    atom2 = atom2.map(nats);

    if (atom0 == atom1 or atom0 == atom2 or atom1 == atom2)
        throw SireMol::ring_error( QObject::tr(
            "You cannot split a molecule into two parts using the same atoms! "
            "(%1, %2, %3).")
                .arg(::atomString(info(),atom0),
                     ::atomString(info(),atom1),
                     ::atomString(info(),atom2)), CODELOC );

    //make sure that there is sufficient space for the
    //selections - this prevents mallocs while tracing
    //the bonds
    group0.reserve(nats);
    group1.reserve(nats);
    root0.reserve(nats);
    root1.reserve(nats);
    
    //add the end atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom2);
    root0.insert(atom0);
    root0.insert(atom1);
    root0.insert(atom2);
    root1.insert(atom2);
    root1.insert(atom1);
    root1.insert(atom0);
    
    //add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if (bonded_atom != atom1)
        {
            this->traceRoute(bonded_atom, root0, group0);
        }
    }
    
    //now add the atoms bonded to atom1 to group1
    bool has_rings = false;
    
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom2))
    {
        if (bonded_atom != atom1)
        {
            if (group0.contains(bonded_atom))
                has_rings = true;

            this->traceRoute(bonded_atom, root1, group1);
        }
    }
    
    //if there is any overlap in the two sets then that means that
    //the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                //qDebug() << "DISCONNECTING(0)" << d->name(bonded_atom) << d->name(atom0);
            }
        }
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom2))
        {
            if (group0.contains(bonded_atom))
            {
                editor.disconnect(atom2, bonded_atom);
                //qDebug() << "DISCONNECTING(1)" << d->name(bonded_atom) << d->name(atom2);
            }
        }
        
        //release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();
        
        //split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1, atom2);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;
        
        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( d->name(atom) );
        }
        
        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( d->name(atom) );
        }
    
        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/
        
        return this->selectGroups(group0, group1);
    }
}

/** Split the molecule into two parts based on the three supplied atoms

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/   
tuple<AtomSelection,AtomSelection> 
ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1,
                        const AtomID &atom2) const
{
    return this->split( d->atomIdx(atom0), d->atomIdx(atom1),
                        d->atomIdx(atom2) );
}

/** Split the molecule into two parts based on the supplied angle

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/   
tuple<AtomSelection,AtomSelection> 
ConnectivityBase::split(const AngleID &angle) const
{
    return this->split( angle.atom0(), angle.atom1(), angle.atom2() );
}

/** Split the selected atoms of this molecule into three parts about the atoms
    'atom0', 'atom1' and 'atom2'.
    
    Note that all three atoms must be contained in the selection or else
    a missing_atom exception will be thrown
    
    An exception will be thrown if it is not possible to split the molecule
    unambiguously in two, as the angle is part of a ring.

    For example;

      C1   C3--C8
        \ /
        C2
       /  \
     C4   C5--C6-C7

    Splitting C5,C2,C3 would return {C5,C6,C7} in one group, and {C3,C8} in the other group.

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/   
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2,
                        const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*d);
    
    if (selected_atoms.selectedAll())
        return this->split(atom0, atom1, atom2);
 
    selected_atoms.assertSelected(atom0);
    selected_atoms.assertSelected(atom1);
    selected_atoms.assertSelected(atom2);

    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    //make sure that there is sufficient space for the
    //selections - this prevents mallocs while tracing
    //the bonds
    group0.reserve(selected_atoms.nSelected());
    group1.reserve(selected_atoms.nSelected());
    root0.reserve(selected_atoms.nSelected());
    root1.reserve(selected_atoms.nSelected());
    
    //map the atoms
    atom0 = atom0.map(d->nAtoms());
    atom1 = atom1.map(d->nAtoms());
    atom2 = atom2.map(d->nAtoms());

    if (atom0 == atom1 or atom0 == atom2 or atom1 == atom2)
        throw SireMol::ring_error( QObject::tr(
            "You cannot split a molecule into two parts using the same atoms! "
            "(%1, %2, %3).")
                .arg(::atomString(info(),atom0),
                     ::atomString(info(),atom1),
                     ::atomString(info(),atom2)), CODELOC );

    //add the two end atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom2);
    
    root0.insert(atom0);
    root0.insert(atom1);
    root0.insert(atom2);
    root1.insert(atom2);
    root1.insert(atom1);
    root1.insert(atom0);
    
    //add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if ( (bonded_atom != atom1) and 
             selected_atoms.selected(bonded_atom) )
        {
            this->traceRoute(selected_atoms, bonded_atom, root0, group0);
        }
    }
    
    //now add the atoms bonded to atom1 to group1
    bool has_rings = false;
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom2))
    {
        if ( (bonded_atom != atom1) and
             selected_atoms.selected(bonded_atom) )
        {
            if (group0.contains(bonded_atom))
                has_rings = true;
                
            this->traceRoute(selected_atoms, bonded_atom, root1, group1);
        }
    }
    
    //if there is any overlap in the two sets then that means that
    //the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                //qDebug() << "DISCONNECTING(0)" << d->name(bonded_atom) << d->name(atom0);
            }
        }
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom2))
        {
            if (group0.contains(bonded_atom))
            {
                editor.disconnect(atom2, bonded_atom);
                //qDebug() << "DISCONNECTING(1)" << d->name(bonded_atom) << d->name(atom2);
            }
        }
        
        //release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();
        
        //split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1, atom2);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;
        
        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( d->name(atom) );
        }
        
        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( d->name(atom) );
        }
    
        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/
        
        return this->selectGroups(group0, group1);
    }
}

/** Split the selected atoms of the molecule into two groups around the 
    three supplied atoms

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2,
                        const AtomSelection &selected_atoms) const
{
    return this->split( d->atomIdx(atom0), d->atomIdx(atom1),
                        d->atomIdx(atom2), selected_atoms );
}
      
/** Split the selected atoms 'selected_atoms' of this molecule  
    into two parts based on the angle identified in 
    'angle'. This splits the molecule about atom0() and atom2()
    of the angle, ignoring atom atom1().
    
    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/      
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const AngleID &angle, 
                        const AtomSelection &selected_atoms) const
{
    return this->split( angle.atom0(), angle.atom1(),
                        angle.atom2(), selected_atoms );
}

/** Split this molecule into two parts based on the passed atoms. 
    This splits the molecule between atom0 and atom3, ignoring 
    atom1 and atom2.

    C1   C4--C5--C6
      \ /
      C2    C8--C9
     /  \  /
    C3   C7
           \
            C10--C11

    Splitting C4,C2,C7,C10 will return {C4,C5,C6} and {C10,C11}. 
    If this molecule had been split by just Bond(C2,C7) using the above 
    function, then the first returned group would
    be {C1,C2,C3,C4,C5,C6}, while the second group would be {C7,C8,C9,C10,C11}.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/      
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1, 
                        AtomIdx atom2, AtomIdx atom3) const
{
    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    //map the atoms
    int nats = d->nAtoms();
    atom0 = atom0.map(nats);
    atom1 = atom1.map(nats);
    atom2 = atom2.map(nats);
    atom3 = atom3.map(nats);

    if (atom0 == atom1 or atom0 == atom2 or atom0 == atom3 or
        atom1 == atom2 or atom1 == atom3 or
        atom2 == atom3)
        throw SireMol::ring_error( QObject::tr(
            "You cannot split a molecule into two parts using the same atoms! "
            "(%1, %2, %3, %4).")
                .arg(::atomString(info(),atom0),
                     ::atomString(info(),atom1),
                     ::atomString(info(),atom2),
                     ::atomString(info(),atom3)), CODELOC );

    //make sure that there is sufficient space for the
    //selections - this prevents mallocs while tracing
    //the bonds
    group0.reserve(nats);
    group1.reserve(nats);
    root0.reserve(nats);
    root1.reserve(nats);
    
    //add the end atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom3);
    
    root0.insert(atom0);
    root0.insert(atom1);
    root0.insert(atom2);
    root0.insert(atom3);
    
    root1.insert(atom3);
    root1.insert(atom2);
    root1.insert(atom1);
    root1.insert(atom0);
    
    //add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if (bonded_atom != atom1)
        {
            this->traceRoute(bonded_atom, root0, group0);
        }
    }
    
    //now add the atoms bonded to atom1 to group1
    bool has_rings = false;
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom3))
    {
        if (bonded_atom != atom2)
        {
            if (group0.contains(bonded_atom))
                has_rings = true;

            this->traceRoute(bonded_atom, root1, group1);
        }
    }

    //if there is any overlap in the two sets then that means that
    //the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                //qDebug() << "DISCONNECTING(0)" << d->name(bonded_atom) << d->name(atom0);
            }
        }
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom3))
        {
            if (group0.contains(bonded_atom))
            {
                editor.disconnect(atom3, bonded_atom);
                //qDebug() << "DISCONNECTING(1)" << d->name(bonded_atom) << d->name(atom3);
            }
        }
        
        //release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();
        
        //split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1, atom2, atom3);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;
        
        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( d->name(atom) );
        }
        
        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( d->name(atom) );
        }
    
        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/
        
        return this->selectGroups(group0, group1);
    }
}

/** Split this molecule into two parts based on the passed atoms. 
    This splits the molecule between atom0 and atom3, ignoring 
    atom1 and atom2.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/      
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1, 
                    const AtomID &atom2, const AtomID &atom3) const
{
    return this->split( d->atomIdx(atom0), d->atomIdx(atom1),
                        d->atomIdx(atom2), d->atomIdx(atom3) );
}
      
/** Split this molecule into two parts based on the dihedral identified in 
    'dihedral'. This splits the molecule about atom0() and atom3()
    of the dihedral, ignoring atoms atom1() and atom2().
    
    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/      
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const DihedralID &dihedral) const
{
    return this->split( dihedral.atom0(), dihedral.atom1(),
                        dihedral.atom2(), dihedral.atom3() );
}

/** Split the selected atoms of this molecule into two parts 
    based on the passed atoms. 
    
    This splits the molecule between atom0 and atom3, ignoring 
    atom1 and atom2.

    All four atoms must be selected in 'selected_atoms' or else
    a missing_atom exception will be thrown

    C1   C4--C5--C6
      \ /
      C2    C8--C9
     /  \  /
    C3   C7
           \
            C10--C11

    Splitting C4,C2,C7,C10 will return {C4,C5,C6} and {C10,C11}. 
    If this molecule had been split by just Bond(C2,C7) using the above 
    function, then the first returned group would
    be {C1,C2,C3,C4,C5,C6}, while the second group would be {C7,C8,C9,C10,C11}.

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/      
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1, 
                        AtomIdx atom2, AtomIdx atom3,
                        const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*d);
    
    if (selected_atoms.selectedAll())
        return this->split(atom0, atom1, atom2, atom3);
 
    selected_atoms.assertSelected(atom0);
    selected_atoms.assertSelected(atom1);
    selected_atoms.assertSelected(atom2);
    selected_atoms.assertSelected(atom3);

    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    //make sure that there is sufficient space for the
    //selections - this prevents mallocs while tracing
    //the bonds
    group0.reserve(selected_atoms.nSelected());
    group1.reserve(selected_atoms.nSelected());
    root0.reserve(selected_atoms.nSelected());
    root1.reserve(selected_atoms.nSelected());
    
    //map the atoms
    atom0 = atom0.map(d->nAtoms());
    atom1 = atom1.map(d->nAtoms());
    atom2 = atom2.map(d->nAtoms());
    atom3 = atom3.map(d->nAtoms());

    if (atom0 == atom1 or atom0 == atom2 or atom0 == atom3 or
        atom1 == atom2 or atom1 == atom3 or
        atom2 == atom3)
        throw SireMol::ring_error( QObject::tr(
            "You cannot split a molecule into two parts using the same atoms! "
            "(%1, %2, %3, %4).")
                .arg(::atomString(info(),atom0),
                     ::atomString(info(),atom1),
                     ::atomString(info(),atom2),
                     ::atomString(info(),atom3)), CODELOC );

    //add the two end atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom3);
    
    root0.insert(atom0);
    root0.insert(atom1);
    root0.insert(atom2);
    root0.insert(atom3);
    
    root1.insert(atom3);
    root1.insert(atom2);
    root1.insert(atom1);
    root1.insert(atom0);
    
    //add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if ( (bonded_atom != atom1) and 
             selected_atoms.selected(bonded_atom) )
        {
            this->traceRoute(selected_atoms, bonded_atom, root0, group0);
        }
    }
    
    //now add the atoms bonded to atom1 to group1
    bool has_rings = true;
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom3))
    {
        if ( (bonded_atom != atom2) and
             selected_atoms.selected(bonded_atom) )
        {
            if (group0.contains(bonded_atom))
                has_rings = true;
                
            this->traceRoute(selected_atoms, bonded_atom, root1, group1);
        }
    }
    
    //if there is any overlap in the two sets then that means that
    //the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                //qDebug() << "DISCONNECTING(0)" << d->name(bonded_atom) << d->name(atom0);
            }
        }
        
        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom3))
        {
            if (group0.contains(bonded_atom))
            {
                editor.disconnect(atom3, bonded_atom);
                //qDebug() << "DISCONNECTING(1)" << d->name(bonded_atom) << d->name(atom3);
            }
        }
        
        //release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();
        
        //split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1, atom2, atom3);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;
        
        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( d->name(atom) );
        }
        
        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( d->name(atom) );
        }
    
        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/
        
        return this->selectGroups(group0, group1);
    }
}
      
/** Split the selected atoms 'selected_atoms' of this molecule
    into two parts based on the passed atoms. This splits
    the molecule between atom0 and atom3, ignoring atom1 and 
    atom2.
    
    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/      
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1,
                    const AtomID &atom2, const AtomID &atom3,
                    const AtomSelection &selected_atoms) const
{
    return this->split( d->atomIdx(atom0), d->atomIdx(atom1),
                        d->atomIdx(atom2), d->atomIdx(atom3),
                        selected_atoms );
}
      
/** Split the selected atoms 'selected_atoms' of this molecule  
    into two parts based on the dihedral identified in 
    'dihedral'. This splits the molecule about atom0() and atom3()
    of the dihedral, ignoring atoms atom1() and atom2().
    
    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/      
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const DihedralID &dihedral,
                    const AtomSelection &selected_atoms) const
{
    return this->split( dihedral.atom0(), dihedral.atom1(),
                        dihedral.atom2(), dihedral.atom3(),
                        selected_atoms );
}

/** Split this molecule into two parts based on the improper angle
    identified by 'improper'. This splits the molecule about
    bond between atom0() and atom1() of the improper
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/      
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const ImproperID &improper) const
{
    return this->split( improper.atom0(), improper.atom1() );
}

/** Split the selected atoms in 'selected_atoms' in this molecule
    into two parts based on the improper angle
    identified by 'improper'. This splits the molecule about
    bond between atom0() and atom1() of the improper
    
    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/      
tuple<AtomSelection,AtomSelection>
ConnectivityBase::split(const ImproperID &improper,
                    const AtomSelection &selected_atoms) const
{
    return this->split( improper.atom0(), improper.atom1(),
                        selected_atoms );
}

/** Return the list of bonds present in this connectivity*/
QList<BondID> ConnectivityBase::getBonds() const
{
  QList<BondID> bonds;
  int nats = connected_atoms.count();
  for (int i=0; i < nats; ++i)
    {
      AtomIdx atomidx = AtomIdx(i);
      QSet<AtomIdx> neighbors = this->connectionsTo(atomidx);
      foreach (AtomIdx neighbor, neighbors)
	{
	  // Do not add the bond to the list if already found. Note that we 
	  // have to check if the bond has been defined in reverse order too
	  // and we do that using the mirror() method
	  BondID bond = BondID(atomidx, neighbor);
	  if ( not ( bonds.contains(bond) or bonds.contains(bond.mirror()) ) )
	    bonds.append(bond);
	}
    }
  return bonds;
}
/** Return the list of bonds in the connectivity containing atom */
QList<BondID> ConnectivityBase::getBonds(const AtomID &atom) const
{
   QList<BondID> bonds;

   // the "connectionsTo()" function will throw exceptions
   // if no atom matches the ID

   foreach (AtomIdx bonded_atom, this->connectionsTo(atom))
   {
       bonds.append( BondID(atom, bonded_atom) );
   }
   // It is ok to return an empty list of bonds
   //if (bonds.isEmpty())
   //    throw SireMol::missing_bond( QObject::tr(
   //       "There are no bonds to the atom with ID %1.")
   //           .arg(atom.toString()), CODELOC );

   return bonds;
}
/** Return a list of angles defined by the connectivity*/
QList<AngleID> ConnectivityBase::getAngles() const
{
  QList<AngleID> angles;
  int nats = connected_atoms.count();
  for (int i=0; i < nats; ++i)
    {
      AtomIdx atom0idx = AtomIdx(i);
      foreach (AtomIdx atom1idx, this->connectionsTo(atom0idx))
	{
	  foreach (AtomIdx atom2idx, this->connectionsTo(atom1idx))
	    {
	      if (atom2idx != atom0idx)
		{
		  AngleID angle = AngleID( atom0idx, atom1idx, atom2idx );
		  if ( not ( angles.contains(angle) or angles.contains(angle.mirror()) ) ) 
		    angles.append(angle);
		}
	    }
	}
    }
  return angles;
}
/** Return a list of angles defined by the connectivity that involve atom0 and atom1*/
QList<AngleID> ConnectivityBase::getAngles(const AtomID &atom0, const AtomID &atom1) const
{
  QList<AngleID> angles;
  // Is this the best way to do this? What if map returns multiple atoms?
  AtomIdx atom0idx = atom0.map(this->info())[0];
  AtomIdx atom1idx = atom1.map(this->info())[0];
  // check that atom0 and atom1 are bonded
  QSet<AtomIdx> bonded_atoms0 = this->connectionsTo(atom0);
  if (not bonded_atoms0.contains(atom1idx))
    throw SireMol::missing_bond( QObject::tr(
					     "There is no bond between atoms with ID %1 and %2.")
				 .arg(atom0.toString(),atom1.toString()), CODELOC );  
  // Get all the neighbors of atom1 that are not atom0 (to get the set of atom0-atom1-atom2)
  foreach (AtomIdx atom2idx, this->connectionsTo(atom1idx))
    {
      if ( atom2idx != atom0idx )
	{
	  AngleID angle = AngleID( atom0idx, atom1idx, atom2idx );
	  if ( not ( angles.contains(angle) or angles.contains(angle.mirror()) ) )
	    angles.append(angle);
	}
    }
  // And now the neighbors of atom0 so we get the set (atom1-atom0-atom2)
  foreach (AtomIdx atom2idx, this->connectionsTo(atom0idx))
    {
      if (atom2idx != atom1idx)
	{
	  AngleID angle = AngleID( atom1idx, atom0idx, atom2idx );
	  if ( not ( angles.contains(angle) or angles.contains(angle.mirror()) ) )
	    angles.append(angle);	  
	}
    }
  // It is ok not to return an empty list of angles
  //if (angles.isEmpty())
  //  throw SireMol::missing_angle( QObject::tr(
  //					      "There are no angles to the atoms with ID %1 and %2.")
  //				  .arg( atom0.toString(),atom1.toString() ), CODELOC );
  
  return angles;
}

/** Return a list of angles defined by the connectivity that involve atom0*/
QList<AngleID> ConnectivityBase::getAngles(const AtomID &atom0) const
{
  QList<AngleID> angles;
  AtomIdx atom0idx = atom0.map(this->info())[0];
  foreach (AtomIdx atom1idx, this->connectionsTo(atom0idx))
    {
      QList<AngleID> angles01 = this->getAngles(atom0idx, atom1idx);
      foreach (AngleID angle, angles01)	
	{
	  if ( not ( angles.contains(angle) or angles.contains(angle.mirror()) ) )
	    angles.append(angle);
	}
    }
  return angles;
}

/** Return a list of dihedrals defined by the connectivity*/
QList<DihedralID> ConnectivityBase::getDihedrals() const
{
  QList<DihedralID> dihedrals;
  int nats = connected_atoms.count();
  for (int i=0; i < nats ; ++i)
    {
      AtomIdx atom0idx = AtomIdx(i);
      foreach (AtomIdx atom1idx, this->connectionsTo(atom0idx))
	{
	  foreach (AtomIdx atom2idx, this->connectionsTo(atom1idx))
	    {
	      if (atom2idx != atom0idx)
		{
		  foreach (AtomIdx atom3idx, this->connectionsTo(atom2idx))
		    {
		      if (atom3idx != atom1idx)
			{
			  DihedralID dihedral = DihedralID( atom0idx, atom1idx, atom2idx, atom3idx);
			  if ( not ( dihedrals.contains(dihedral) or dihedrals.contains(dihedral.mirror()) ) )
			    dihedrals.append(dihedral);
			}
		    } 
		}
	    }
	}
    }
  return dihedrals;
}
/** Return a list of dihedrals defined by the connectivity that involve atom0, atom1 and atom2*/
QList<DihedralID> ConnectivityBase::getDihedrals(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2) const
{
  QList<DihedralID> dihedrals;
  AtomIdx atom0idx = atom0.map(this->info())[0];
  AtomIdx atom1idx = atom1.map(this->info())[0];
  AtomIdx atom2idx = atom2.map(this->info())[0];
  //Check that atom0, atom1 & atom2 form an angle
  QSet<AtomIdx> bonded_atoms0 = this->connectionsTo(atom0);
  if (not bonded_atoms0.contains(atom1idx))
        throw SireMol::missing_bond( QObject::tr(
					     "There is no bond between atoms with ID %1 and %2.")
				     .arg(atom0.toString(),atom1.toString()), CODELOC );  
  QSet<AtomIdx> bonded_atoms1 = this->connectionsTo(atom1);
  if (not bonded_atoms1.contains(atom2idx))
        throw SireMol::missing_bond( QObject::tr(
					     "There is no bond between atoms with ID %1 and %2.")
				     .arg(atom0.toString(),atom2.toString()), CODELOC );  
  // Get all the neighbors of atom2, that are not atom1 or atom0 so we build atom0-atom1-atom2-atom3
  foreach (AtomIdx atom3idx, this->connectionsTo(atom2idx))
    {
      if ( (atom3idx != atom1idx) and (atom3idx != atom0idx) )
	{
	  DihedralID dihedral = DihedralID( atom0idx, atom1idx, atom2idx, atom3idx );
	  if ( not ( dihedrals.contains(dihedral) or dihedrals.contains(dihedral.mirror()) ) )
	    dihedrals.append(dihedral);
	}
    }
  // And now building atom2-atom1-atom0-atom3
  foreach (AtomIdx atom3idx, this->connectionsTo(atom0idx))
    {
      if ( (atom3idx != atom1idx) and (atom3idx != atom2idx) )
	{
	  DihedralID dihedral = DihedralID( atom2idx, atom1idx, atom0idx, atom3idx );
	  if ( not ( dihedrals.contains(dihedral) or dihedrals.contains(dihedral.mirror()) ) )
	    dihedrals.append(dihedral);
	}
    }  
  return dihedrals;
}
/** Return a list of dihedrals defined by the connectivity that involve atom0 and atom1*/
QList<DihedralID> ConnectivityBase::getDihedrals(const AtomID &atom0, const AtomID &atom1) const
{
  QList<DihedralID> dihedrals;
  AtomIdx atom0idx = atom0.map(this->info())[0];
  AtomIdx atom1idx = atom1.map(this->info())[0];
  //Check that atom0 and atom1 are bonded
  QSet<AtomIdx> bonded_atoms0 = this->connectionsTo(atom0);
  if (not bonded_atoms0.contains(atom1idx))
    throw SireMol::missing_bond( QObject::tr(
					     "There is no bond between atoms with ID %1 and %2.")
				 .arg(atom0.toString(),atom1.toString()), CODELOC );  
  foreach (AtomIdx atom2idx, this->connectionsTo(atom1idx))
    {
      if (atom2idx != atom0idx) 
	{
	  // the dihedrals that are atom0-atom1-atom2-XXX
	  QList<DihedralID> dihedrals012 = this->getDihedrals(atom0idx, atom1idx, atom2idx);
	  foreach (DihedralID dihedral,dihedrals012)
	    {
	      if ( not ( dihedrals.contains(dihedral) or dihedrals.contains(dihedral.mirror()) ) )
		dihedrals.append(dihedral);
	    }
	  // the dihedrals that are atom2-atom1-atom0-XXX
	  QList<DihedralID> dihedrals210 = this->getDihedrals(atom2idx, atom1idx, atom0idx);
	  foreach (DihedralID dihedral,dihedrals210)
	    {
	      if ( not ( dihedrals.contains(dihedral) or dihedrals.contains(dihedral.mirror()) ) )
		dihedrals.append(dihedral);
	    }	  
	}
    }
  return dihedrals;
}
/** Return a list of dihedrals defined by the connectivity that involve atom0*/
QList<DihedralID> ConnectivityBase::getDihedrals(const AtomID &atom0) const
{
  QList<DihedralID> dihedrals;
  AtomIdx atom0idx = atom0.map(this->info())[0];
  foreach (AtomIdx atom1idx, this->connectionsTo(atom0idx))
    {
      // Dihedrals going along atom0-atom1-XXX
      QList<DihedralID> dihedrals01 = this->getDihedrals(atom0idx, atom1idx);
      foreach (DihedralID dihedral,dihedrals01)
	{
	  if ( not ( dihedrals.contains(dihedral) or dihedrals.contains(dihedral.mirror()) ) )
	    dihedrals.append(dihedral);
	}
      // Dihedrals going along atom1-atom0-XXX
      QList<DihedralID> dihedrals10 = this->getDihedrals(atom1idx, atom0idx);
      foreach (DihedralID dihedral,dihedrals10)
	{
	  if ( not ( dihedrals.contains(dihedral) or dihedrals.contains(dihedral.mirror()) ) )
	    dihedrals.append(dihedral);
	}     
    }
  return dihedrals;
}

/** Return a matrix (organised by AtomIdx) that says which atoms are bonded between
    order 'start' and order 'end' (e.g. if order is two, it returns true for each atom pair that
    are bonded together, if order is three, then true for each atom pair that are
    bonded or angled together, if order is four, then true for each atom pair
    that are bonded, angled or dihedraled) */
QVector< QVector<bool> > ConnectivityBase::getBondMatrix(int start, int end) const
{
    if (start < 0)
        start = 0;
    
    if (end < 0)
        end = 0;
    
    if (start > end)
        qSwap(start, end);
    
    QVector< QVector<bool> > ret;
    
    const int nats = d.read().nAtoms();
    
    if (nats == 0)
        return ret;
    
    else
    {
        ret = QVector< QVector<bool> >(nats);
        ret.squeeze();
     
        QVector<bool> row;
        
        if (start == 0)
            row = QVector<bool>(nats, true);
        else
            row = QVector<bool>(nats, false);
        
        row.squeeze();
        
        for (int i=0; i<nats; ++i)
        {
            ret.data()[i] = row;
        }
    }
    
    if (start <= 0)
        return ret;
    
    for (int order = start; order <= end; ++order)
    {
        if (order == 1)
        {
            for (int i=0; i<nats; ++i)
            {
                ret[i][i] = true;
            }
        }
    
        if (order == 2)
        {
            for (int i=0; i<nats; ++i)
            {
                QVector<bool> &row = ret.data()[i];
            
                for (QSet<AtomIdx>::const_iterator it = connected_atoms[i].constBegin();
                     it != connected_atoms[i].constEnd();
                     ++it)
                {
                    row[it->value()] = true;
                }
            }
        }
        
        if (order == 3)
        {
            for (int atm0=0; atm0<nats; ++atm0)
            {
                QVector<bool> &row = ret.data()[atm0];
                
                for (QSet<AtomIdx>::const_iterator it = connected_atoms[atm0].constBegin();
                     it != connected_atoms[atm0].constEnd();
                     ++it)
                {
                    const int atm1 = it->value();
                
                    for (QSet<AtomIdx>::const_iterator it2 = connected_atoms[atm1].constBegin();
                         it2 != connected_atoms[atm1].constEnd();
                         ++it2)
                    {
                        const int atm2 = it2->value();
                        
                        if (atm2 != atm0)
                            row[atm2] = true;
                    }
                }
            }
        }
        
        if (order == 4)
        {
            for (int atm0=0; atm0<nats; ++atm0)
            {
                QVector<bool> &row = ret.data()[atm0];
                
                for (QSet<AtomIdx>::const_iterator it = connected_atoms[atm0].constBegin();
                     it != connected_atoms[atm0].constEnd();
                     ++it)
                {
                    const int atm1 = it->value();
                
                    for (QSet<AtomIdx>::const_iterator it2 = connected_atoms[atm1].constBegin();
                         it2 != connected_atoms[atm1].constEnd();
                         ++it2)
                    {
                        const int atm2 = it2->value();
                        
                        if (atm2 != atm0)
                        {
                            for (QSet<AtomIdx>::const_iterator
                                                it3 = connected_atoms[atm2].constBegin();
                                 it3 != connected_atoms[atm2].constEnd();
                                 ++it3)
                            {
                                const int atm3 = it3->value();
                                
                                if (atm3 != atm0 and atm3 != atm1 and atm3 != atm2)
                                {
                                    row[atm3] = true;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if (order > 4)
        {
            qDebug() << "Cannot build a bond matrix for values greater than 4";
            break;
        }
    }
    
    return ret;
}

/** Return a matrix (organised by AtomIdx) that says which atoms are bonded up to 
    order 'order' (e.g. if order is two, it returns true for each atom pair that
    are bonded together, if order is three, then true for each atom pair that are
    bonded or angled together, if order is four, then true for each atom pair
    that are bonded, angled or dihedraled) */
QVector< QVector<bool> > ConnectivityBase::getBondMatrix(int order) const
{
    if (order <= 0)
        return getBondMatrix(0,0);
    else
        return getBondMatrix(1,order);
}

/////////
///////// Implementation of Connectivity
/////////

static const RegisterMetaType<Connectivity> r_connectivity;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const Connectivity &conn)
{
    writeHeader(ds, r_connectivity, 1);
    
    ds << static_cast<const ConnectivityBase&>(conn);
    
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, Connectivity &conn)
{
    VersionID v = readHeader(ds, r_connectivity);
    
    if (v == 1)
    {
        ds >> static_cast<ConnectivityBase&>(conn);
    }
    else
        throw version_error(v, "1", r_connectivity, CODELOC);
        
    return ds;
}

/** Reduce the memory usage of this object to a minimum */
void Connectivity::squeeze()
{
    connected_atoms.squeeze();
    
    QSet<AtomIdx> *connected_atoms_array = connected_atoms.data();
    int nats = connected_atoms.count();
    
    for (int i=0; i<nats; ++i)
    {
        connected_atoms_array[i].squeeze();
    }
    
    connected_res.squeeze();
    
    QSet<ResIdx> *connected_res_array = connected_res.data();
    int nres = connected_res.count();
    
    for (int i=0; i<nres; ++i)
    {
        connected_res_array[i].squeeze();
    }
}

/** Null constructor */
Connectivity::Connectivity() 
             : ConcreteProperty<Connectivity,ConnectivityBase>()
{}

/** Construct the connectivity for the molecule whose data
    is in 'moldata' */
Connectivity::Connectivity(const MoleculeData &moldata)
             : ConcreteProperty<Connectivity,ConnectivityBase>(moldata)
{}

    
/** Construct the connectivity for the molecule viewed in the 
    passed view. This automatically uses the bond hunting 
    function to add all of the bonds for the atoms in this view */
Connectivity::Connectivity(const MoleculeView &molview,
                           const BondHunter &bondhunter,
                           const PropertyMap &map)
             : ConcreteProperty<Connectivity,ConnectivityBase>()
{
    this->operator=( bondhunter(molview, map) );
}

/** Construct the connectivity from the passed editor */
Connectivity::Connectivity(const ConnectivityEditor &editor)
             : ConcreteProperty<Connectivity,ConnectivityBase>(editor)
{
    this->squeeze();
}

/** Private constructor allowing a ConnectivityBase to become a Connectivity */
Connectivity::Connectivity(const ConnectivityBase &base)
             : ConcreteProperty<Connectivity,ConnectivityBase>(base)
{}

/** Copy constructor */
Connectivity::Connectivity(const Connectivity &other)
             : ConcreteProperty<Connectivity,ConnectivityBase>(other)
{}

/** Destructor */
Connectivity::~Connectivity()
{}

/** Copy assignment from another Connectivity object */
Connectivity& Connectivity::operator=(const Connectivity &other)
{
    ConnectivityBase::operator=(other);
    return *this;
}

/** Copy assignment from a ConnectivityEditor */
Connectivity& Connectivity::operator=(const ConnectivityEditor &editor)
{
    ConnectivityBase::operator=(editor);
    this->squeeze();
    return *this;
}

/** Comparison operator */
bool Connectivity::operator==(const Connectivity &other) const
{
    return ConnectivityBase::operator==(other);
}

/** Comparison operator */
bool Connectivity::operator!=(const Connectivity &other) const
{
    return ConnectivityBase::operator!=(other);
}

/** Return an editor that can edit a copy of this connectivity */
ConnectivityEditor Connectivity::edit() const
{
    return ConnectivityEditor(*this);
}

const char* Connectivity::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Connectivity>() );
}

/////////
///////// Implementation of ConnectivityEditor
/////////

static const RegisterMetaType<ConnectivityEditor> r_conneditor;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const ConnectivityEditor &conn)
{
    writeHeader(ds, r_conneditor, 1);
    
    ds << static_cast<const ConnectivityBase&>(conn);
    
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, ConnectivityEditor &conn)
{
    VersionID v = readHeader(ds, r_conneditor);
    
    if (v == 1)
    {
        ds >> static_cast<ConnectivityBase&>(conn);
    }
    else
        throw version_error(v, "1", r_conneditor, CODELOC);
        
    return ds;
}

/** Null constructor */
ConnectivityEditor::ConnectivityEditor()
                   : ConcreteProperty<ConnectivityEditor,ConnectivityBase>()
{}

/** Construct an editor to edit a copy of the passed 
    Connectivity object */
ConnectivityEditor::ConnectivityEditor(const Connectivity &connectivity)
                   : ConcreteProperty<ConnectivityEditor,ConnectivityBase>(connectivity)
{}

/** Copy constructor */
ConnectivityEditor::ConnectivityEditor(const ConnectivityEditor &other)
                   : ConcreteProperty<ConnectivityEditor,ConnectivityBase>(other)
{}

/** Destructor */
ConnectivityEditor::~ConnectivityEditor()
{}

/** Copy assignment operator */
ConnectivityEditor& ConnectivityEditor::operator=(const ConnectivityBase &other)
{
    ConnectivityBase::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool ConnectivityEditor::operator==(const ConnectivityEditor &other) const
{
    return ConnectivityBase::operator==(other);
}

/** Comparison operator */
bool ConnectivityEditor::operator!=(const ConnectivityEditor &other) const
{
    return ConnectivityBase::operator!=(other);
}

/** Record the connection between the atoms at indicies 'atom0' 
    and 'atom1'
    
    \throw SireError::invalid_index
*/
ConnectivityEditor& ConnectivityEditor::connect(AtomIdx atom0, AtomIdx atom1)
{
    AtomIdx atomidx0 = AtomIdx( atom0.map(connected_atoms.count()) );
    AtomIdx atomidx1 = AtomIdx( atom1.map(connected_atoms.count()) );
    
    if (atomidx0 == atomidx1)
        return *this;
    
    QSet<AtomIdx> *connected_atoms_array = connected_atoms.data();
    
    connected_atoms_array[atomidx0].insert(atomidx1);
    connected_atoms_array[atomidx1].insert(atomidx0);
    
    if (info().isWithinResidue(atomidx0) and
        info().isWithinResidue(atomidx1))
    {
        QSet<ResIdx> *connected_res_array = connected_res.data();

        ResIdx residx0 = info().parentResidue(atomidx0);
        ResIdx residx1 = info().parentResidue(atomidx1);
        
        if (residx0 != residx1)
        {
            connected_res_array[residx0].insert(residx1);
            connected_res_array[residx1].insert(residx0);
        }
    }
    
    return *this;
}

/** Record a connection between the atom identified by 'atom0' and 
    the atom identified by 'atom1'
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
ConnectivityEditor& ConnectivityEditor::connect(const AtomID &atom0,
                                                const AtomID &atom1)
{
    return this->connect( info().atomIdx(atom0), info().atomIdx(atom1) );
}

/** Remove the connection between the atoms at indicies 'atom0'
    and 'atom1' - this does nothing if there isn't already a connection
    
    \throw SireError::invalid_index
*/
ConnectivityEditor& ConnectivityEditor::disconnect(AtomIdx atom0, AtomIdx atom1)
{
    if (this->areConnected(atom0, atom1))
    {
        AtomIdx atomidx0 = AtomIdx( atom0.map( connected_atoms.count() ) );
        AtomIdx atomidx1 = AtomIdx( atom1.map( connected_atoms.count() ) );
        
        connected_atoms[atomidx0].remove(atomidx1);
        connected_atoms[atomidx1].remove(atomidx0);
        
        //now check to see if the residues are still connected
        ResIdx residx0 = info().parentResidue(atomidx0);
        ResIdx residx1 = info().parentResidue(atomidx1);

        if (this->nConnections(residx0, residx1) == 0)
        {
            connected_res[residx0].remove(residx1);
            connected_res[residx1].remove(residx0);
        }
    }

    return *this;
}

/** Disconnect the atoms that are identified by 'atom0' and 'atom1' - 
    this does nothing if there isn't a connection between these atoms
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
ConnectivityEditor& ConnectivityEditor::disconnect(const AtomID &atom0, 
                                                   const AtomID &atom1)
{
    return this->disconnect( info().atomIdx(atom0), info().atomIdx(atom1) );
}

/** Remove all of the connections to the atom at index 'atomidx' 

    \throw SireError::invalid_index
*/
ConnectivityEditor& ConnectivityEditor::disconnectAll(AtomIdx atomidx)
{
    QSet<AtomIdx> connected = this->connectionsTo(atomidx);
    
    foreach (AtomIdx atom1, connected)
    {
        this->disconnect(atomidx, atom1);
    }
    
    return *this;
}

/** Remove all bonds from this molecule */
ConnectivityEditor& ConnectivityEditor::disconnectAll()
{
    connected_atoms.clear();
    connected_res.clear();

    if (info().nAtoms() > 0)
    {
        connected_atoms.resize(info().nAtoms());
        connected_atoms.squeeze();
    }
    
    if (info().nResidues() > 0)
    {
        connected_res.resize(info().nResidues());
        connected_res.squeeze();
    }
    
    return *this;
}

/** Remove all of the connections to the atom identified by 'atomid' 

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
ConnectivityEditor& ConnectivityEditor::disconnectAll(const AtomID &atomid)
{
    return this->disconnectAll(info().atomIdx(atomid));
}

/** Remove all of the connections that involve any of the atoms
    in the residue at index 'residx'
    
    \throw SireError::invalid_index
*/
ConnectivityEditor& ConnectivityEditor::disconnectAll(ResIdx residx)
{
    foreach (AtomIdx atomidx, info().getAtomsIn(residx))
    {
        this->disconnectAll(atomidx);
    }
    
    return *this;
}

/** Remove all of the connections that involve any of the atoms
    in the residue identified by 'resid'
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ConnectivityEditor& ConnectivityEditor::disconnectAll(const ResID &resid)
{
    return this->disconnectAll( info().resIdx(resid) );
}

/** Return the editied connectivity */
Connectivity ConnectivityEditor::commit() const
{
    return Connectivity(*this);
}

const char* ConnectivityEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ConnectivityEditor>() );
}
