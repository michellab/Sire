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

#include "patches.h"

#include "SireError/errors.h"

#include "tostring.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireFF;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

////////////
//////////// Implementation of FFBead
////////////

static const RegisterMetaType<FFBead> r_ffbead(NO_ROOT);

QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const FFBead &ffbead)
{
    writeHeader(ds, r_ffbead, 1);
    
    SharedDataStream sds(ds);
    
    sds << ffbead.coords << ffbead.params;
    
    return ds;
}

QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, FFBead &ffbead)
{
    VersionID v = readHeader(ds, r_ffbead);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> ffbead.coords >> ffbead.params;
    }
    else
        throw version_error(v, "1", r_ffbead, CODELOC);
        
    return ds;
}

/** Null Constructor */
FFBead::FFBead()
{}

/** Construct a bead from the passed coordinates and parameters */
FFBead::FFBead(const CoordGroup &coordinates, const FFParameters &parameters)
       : coords(coordinates), params(parameters)
{}

/** Copy constructor */
FFBead::FFBead(const FFBead &other)
       : coords(other.coords), params(other.params)
{}

/** Destructor */
FFBead::~FFBead()
{}

/** Copy assignment operator */
FFBead& FFBead::operator=(const FFBead &other)
{
    coords = other.coords;
    params = other.params;
    return *this;
}

/** Comparison operator */
bool FFBead::operator==(const FFBead &other) const
{
    return coords.constData() == other.coords.constData() and
           (params.constData() == other.params.constData() or
            params.read().equals(other.params.read()));
}

/** Comparison operator */
bool FFBead::operator!=(const FFBead &other) const
{
    return not FFBead::operator==(other);
}

const char* FFBead::typeName()
{
    return "SireFF::FFBead";
}

////////////
//////////// Implementation of FFBeadChange
////////////

static const RegisterMetaType<FFBeadChange> r_ffbeadchange(NO_ROOT);

QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const FFBeadChange &bead)
{
    writeHeader(ds, r_ffbeadchange, 1);
    
    SharedDataStream sds(ds);
    
    sds << bead.old_bead << bead.new_bead;
    
    return ds;
}

QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, FFBeadChange &bead)
{
    VersionID v = readHeader(ds, r_ffbeadchange);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> bead.old_bead >> bead.new_bead;
    }
    else
        throw version_error(v, "1", r_ffbeadchange, CODELOC);
        
    return ds;
}

/** Null constructor */
FFBeadChange::FFBeadChange()
{}

/** Construct for the change from 'old_bead' to 'new_bead' */
FFBeadChange::FFBeadChange(const FFBead &oldbead, const FFBead &newbead)
             : old_bead(oldbead), new_bead(newbead)
{
    if (old_bead == new_bead)
    {
        old_bead = FFBead();
        new_bead = FFBead();
    }
}

/** Copy constructor */
FFBeadChange::FFBeadChange(const FFBeadChange &other)
             : old_bead(other.old_bead), new_bead(other.new_bead)
{}

/** Destructor */
FFBeadChange::~FFBeadChange()
{}

/** Copy assignment operator */
FFBeadChange& FFBeadChange::operator=(const FFBeadChange &other)
{
    old_bead = other.old_bead;
    new_bead = other.new_bead;
    return *this;
}

/** Comparison operator */
bool FFBeadChange::operator==(const FFBeadChange &other) const
{
    return old_bead == other.old_bead and
           new_bead == other.new_bead;
}

/** Comparison operator */
bool FFBeadChange::operator!=(const FFBeadChange &other) const
{
    return not FFBeadChange::operator==(other);
}

/** Return the change from this->oldBead() to 'newer_bead' */
FFBeadChange FFBeadChange::update(const FFBead &newer_bead) const
{
    return FFBeadChange(old_bead, newer_bead);
}

const char* FFBeadChange::typeName()
{
    return "SireFF::FFBeadChange";
}

////////////
//////////// Implementation of Patches
////////////

static const RegisterMetaType<Patches> r_patches;

QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const Patches &patches)
{
    writeHeader(ds, r_patches, 1);
    
    SharedDataStream sds(ds);
    
    sds << patches.ptchng << patches.ptchs
        << patches.beadid_to_patch << patches.last_beadid
        << static_cast<const Property&>(patches);
        
    return ds;
}

QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, Patches &patches)
{
    VersionID v = readHeader(ds, r_patches);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> patches.ptchng >> patches.ptchs
            >> patches.beadid_to_patch >> patches.last_beadid
            >> static_cast<Property&>(patches);
    }
    else
        throw version_error(v, "1", r_patches, CODELOC);
        
    return ds;
}

/** Construct the patches for a single infinite space (just 1 patch for all beads) */
Patches::Patches() : ConcreteProperty<Patches,Property>(), last_beadid(0)
{
    ptchs = QVector<Patch>(1);
    ptchs.squeeze();
}

/** Construct the patches for the passed space, using the passed
    patching */
Patches::Patches(const Space &space, const Patching &patchng)
        : ConcreteProperty<Patches,Property>(),
          ptchng( patchng.repatch(space) ), last_beadid(0)
{
    ptchs = QVector<Patch>(ptchng.read().nPatches());
    ptchs.squeeze();
}

/** Construct for the passed patching */
Patches::Patches(const Patching &patchng)
        : ConcreteProperty<Patches,Property>(),
          ptchng(patchng), last_beadid(0)
{
    ptchs = QVector<Patch>(ptchng.read().nPatches());
    ptchs.squeeze();
}
      
/** Copy constructor */
Patches::Patches(const Patches &other)
        : ConcreteProperty<Patches,Property>(other),
          ptchng(other.ptchng), ptchs(other.ptchs),
          beadid_to_patch(other.beadid_to_patch),
          last_beadid(other.last_beadid)
{}

/** Destructor */
Patches::~Patches()
{}

/** Copy assignment operator */
Patches& Patches::operator=(const Patches &other)
{
    if (this != &other)
    {
        ptchng = other.ptchng;
        ptchs = other.ptchs;
        beadid_to_patch = other.beadid_to_patch;
        last_beadid = other.last_beadid;
    }
    
    return *this;
}

/** Comparison operator */
bool Patches::operator==(const Patches &other) const
{
    return this == &other or
           (ptchng == other.ptchng and
            ptchs == other.ptchs and
            beadid_to_patch == other.beadid_to_patch and
            last_beadid == other.last_beadid);
}

/** Comparison operator */
bool Patches::operator!=(const Patches &other) const
{
    return not Patches::operator==(other);
}

/** Return a string describing the patches */
QString Patches::toString() const
{
    if (this->isEmpty())
        return QObject::tr("Patches::null");
        
    QStringList lines;
    
    for (int i=0; i<this->nPatches(); ++i)
    {
        lines.append( this->constData()[i].toString() );
    }
    
    return QObject::tr("Patches( patching == %1\n"
                       "         nBeads() == %2\n    %3\n        )")
                .arg(ptchng.read().toString())
                .arg(this->nBeads())
                .arg(lines.join("\n    "));
                       
}

/** Return the patching function used to divide the beads between patches */
const Patching& Patches::patching() const
{
    return ptchng.read();
}

/** Return the space in which these patches exist */
const Space& Patches::space() const
{
    return ptchng.read().space();
}

/** Repatch the space and all of the beads */
void Patches::repatch(const Patching &patching)
{
    if (ptchng.read().equals(patching))
        return;

    QVector<Patch> new_patches( patching.nPatches() );
    new_patches.squeeze();
    
    QHash<quint32,quint32> new_beadid_to_patch;
    new_beadid_to_patch.reserve(beadid_to_patch.count());
    
    const Space &space = patching.space();
    
    for (int i=0; i<ptchs.count(); ++i)
    {
        const Patch &patch = ptchs.constData()[i];

        for (int j=0; j<patch.nBeads(); ++j)
        {
            quint32 beadid = patch.beadIDs().constData()[j];
            const CoordGroup &coords = patch.coordinates().constData()[j];
            FFParametersPtr params = patch.parameters()[j];
            
            QPair<int,Vector> idx = patching.patchIndexAndCenter(coords.aaBox().center());
            
            BOOST_ASSERT(idx.first >= 0 and idx.first < new_patches.count());
            
            new_patches.data()[idx.first].add(beadid, 
                                              space.getMinimumImage(coords, idx.second),
                                              params.read());
                                       
            new_beadid_to_patch.insert(beadid, idx.first);
        }
    }
    
    ptchs = new_patches;
    beadid_to_patch = new_beadid_to_patch;
    ptchng = patching;
}

/** Repatch for a new space */
void Patches::repatch(const Space &space)
{
    this->repatch( ptchng.read().repatch(space) );
}

/** Return the number of patches */
int Patches::nPatches() const
{
    return ptchs.count();
}

/** Return the number of patches */
int Patches::size() const
{
    return this->nPatches();
}

/** Return the number of patches */
int Patches::count() const
{
    return this->nPatches();
}

/** Return whether or not this this empty */
bool Patches::isEmpty() const
{
    return ptchs.isEmpty();
}

/** Return the total number of beads in the patches */
int Patches::nBeads() const
{
    return beadid_to_patch.count();
}

/** Return the ith patch 

    \throw SireError::invalid_index
*/
const Patch& Patches::operator[](int i) const
{
    if (i < 0 or i >= ptchs.count())
        throw SireError::invalid_index( QObject::tr(
                "There is no patch with index %1 - the number of patches is %2.")
                    .arg(i).arg(this->nPatches()), CODELOC );
                    
    return ptchs.constData()[i];
}

/** Return the ith patch 

    \throw SireError::invalid_index
*/
const Patch& Patches::at(int i) const
{
    return this->operator[](i);
}

/** Return a pointer to the array of patches */
const Patch* Patches::data() const
{
    return ptchs.constData();
}

/** Return a pointer to the array of patches */
const Patch* Patches::constData() const
{
    return ptchs.constData();
}

/** Return the patch index, and index in that patch, of the bead
    with the passed beadid
    
    \throw SireError::invalid_key
*/
QPair<int,int> Patches::getLocation(quint32 beadid) const
{
    int idx = beadid_to_patch.value(beadid, -1);
    
    if (idx == -1)
        throw SireError::invalid_key( QObject::tr(
                "There is no bead with ID == %1. Available IDs are %2.")
                    .arg(beadid).arg(Sire::toString(beadid_to_patch.keys())),
                        CODELOC );
                        
    return QPair<int,int>(idx, ptchs.constData()[idx].getLocation(beadid));
}

/** Add a new bead (with specified coordinates and parameters) to these patches,
    returning the ID of the bead in these patches */
QPair<quint32,FFBead> Patches::add(const CoordGroup &coords, const FFParameters &params)
{
    QPair<int,Vector> idx = ptchng.read().patchIndexAndCenter(
                                                coords.aaBox().center() );
                                                
    BOOST_ASSERT( idx.first >= 0 and idx.first < ptchs.count() );
    
    FFBead bead = ptchs.data()[idx.first].add( last_beadid + 1,
                                                ptchng.read().space()
                                                .getMinimumImage(coords, idx.second),
                                                params );

    ++last_beadid;
    
    beadid_to_patch.insert(last_beadid, idx.first);
    
    return QPair<quint32,FFBead>(last_beadid, bead);
}

/** Add the passed bead, returning the ID of the bead in these patches */
QPair<quint32,FFBead> Patches::add(const FFBead &bead)
{
    return this->add(bead.coordinates(), bead.parameters());
}

/** Add the passed beads, returning the IDs of the beads in these patches */
QHash<quint32,FFBead> Patches::add(const CoordGroupArray &coords, 
                                   const FFParametersArray &params)
{
    if (coords.count() != params.count())
        throw SireError::program_bug( QObject::tr(
                "The number of bead parameters (%1) must equal the number of "
                "coordinate groups (%2).")
                    .arg(params.count()).arg(coords.count()), CODELOC );

    Patches new_patches(*this);
    QHash<quint32,FFBead> beads;
    beads.reserve(coords.count());
    
    for (int i=0; i<coords.count(); ++i)
    {
        QPair<quint32,FFBead> bead = new_patches.add(coords.constData()[i], params[i]);
        beads.insert( bead.first, bead.second );
    }
    
    this->operator=(new_patches);
    return beads;
}

/** Add the passed beads, returning the IDs of the beads in these patches */
QHash<quint32,FFBead> Patches::add(const QVector<FFBead> &beads)
{
    Patches new_patches(*this);
        
    QHash<quint32,FFBead> new_beads;
    new_beads.reserve(beads.count());
    
    for (int i=0; i<beads.count(); ++i)
    {
        QPair<quint32,FFBead> bead = new_patches.add(beads[i].coordinates(), 
                                                     beads[i].parameters());
                                                     
        new_beads.insert(bead.first, bead.second);
    }
    
    this->operator=(new_patches);
    return new_beads;
}

/** Add the passed beads, returning the IDs of the beads in these patches */
QHash<quint32,FFBead> Patches::add(const QVector<CoordGroup> &coords,
                                   const QVector<FFParametersPtr> &params)
{
    if (coords.count() != params.count())
        throw SireError::program_bug( QObject::tr(
                "The number of bead parameters (%1) must equal the number of "
                "coordinate groups (%2).")
                    .arg(params.count()).arg(coords.count()), CODELOC );

    Patches new_patches(*this);

    QHash<quint32,FFBead> beads;
    beads.reserve(coords.count());
    
    for (int i=0; i<coords.count(); ++i)
    {
        QPair<quint32,FFBead> bead = new_patches.add(coords.constData()[i], 
                                                     params.constData()[i]);
                                                    
        beads.insert(bead.first, bead.second);
    }
    
    this->operator=(new_patches);
    return beads;
}

int Patches::getIdx(quint32 beadid) const
{
    int idx = beadid_to_patch.value(beadid,-1);
    
    if (idx == -1)
        throw SireError::invalid_key( QObject::tr(
                "There is no bead with ID == %1 in these patches. "
                "Available bead IDs are %2.")
                    .arg(beadid).arg(Sire::toString(beadid_to_patch.keys())), 
                        CODELOC );
                        
    return idx;
}

/** Update the patch with bead ID 'beadid' with the new coordinates 'coords'.
    Note that the number of coordinates must now change!
    
    \throw SireError::invalid_key
    \throw SireError::incompatible_error
*/
FFBeadChange Patches::update(quint32 beadid, const CoordGroup &coords)
{
    int old_idx = getIdx(beadid);

    //is this changing patch?
    QPair<int,Vector> new_idx = patching()
                                    .patchIndexAndCenter(coords.aaBox().center());

    if (new_idx.first == old_idx)
    {
        //just update the coordinates
        return ptchs.data()[old_idx].update(beadid, 
                                 space().getMinimumImage(coords,new_idx.second));
    }
    else
    {
        FFParametersPtr old_params = ptchs.constData()[old_idx]
                                      .parameters()[ptchs.constData()[old_idx]
                                                         .getLocation(beadid)];
    
        FFBead added = ptchs.data()[new_idx.first].add(beadid,
                                      space().getMinimumImage(coords,new_idx.second),
                                      old_params.read());
                             
        return FFBeadChange(ptchs.data()[old_idx].remove(beadid), added);
    }
}

/** Update the bead with ID 'beadid' with the parameters in 'params'.

    \throw SireError::invalid_key
    \throw SireError::incompatible_error
*/
FFBeadChange Patches::update(quint32 beadid, const FFParameters &params)
{
    int idx = getIdx(beadid);
    return ptchs.data()[idx].update(beadid, params);
}

/** Update the bead with ID 'beadid' with the passed coordinates and parameters

    \throw SireError:::invalid_key
    \throw SireError::incompatible_error
*/
FFBeadChange Patches::update(quint32 beadid, const CoordGroup &coords,
                             const FFParameters &params)
{
    int old_idx = getIdx(beadid);
    
    //is this changing patch?
    QPair<int,Vector> new_idx = patching().patchIndexAndCenter(coords.aaBox().center());
    
    if (new_idx.first == old_idx)
    {
        //just update the coordinates and parameters in place
        return ptchs.data()[old_idx].update(beadid, 
                                     space().getMinimumImage(coords,new_idx.second),
                                     params);
    }
    else
    {
        //add the bead to its new patch, and remove from the old patch
        FFBead added = ptchs.data()[new_idx.first].add(beadid,
                                       space().getMinimumImage(coords,new_idx.second),
                                       params);
                                 
        return FFBeadChange(ptchs.data()[old_idx].remove(beadid), added);
    }
}

/** Update the beads with the passed bead IDs with the coordinates in the
    passed array
    
    \throw SireError:::invalid_key
    \throw SireError::incompatible_error
*/
QHash<quint32,FFBeadChange> Patches::update(const QVector<quint32> &beadids, 
                                            const CoordGroupArray &coords)
{
    if (beadids.count() != coords.count())
        throw SireError::program_bug( QObject::tr(
                "You must pass in the same number of coordinate groups (%1) "
                "as you pass in bead IDs (%2).")
                    .arg(coords.count()).arg(beadids.count()), CODELOC );
                

    Patches new_patches(*this);
    
    QHash<quint32,FFBeadChange> deltas;
    deltas.reserve(beadids.count());
    
    for (int i=0; i<beadids.count(); ++i)
    {
        FFBeadChange delta = new_patches
                            .update( beadids.constData()[i], coords.constData()[i] );
                            
        if (not delta.isEmpty())
            deltas.insert(beadids.constData()[i], delta);
    }
    
    this->operator=(new_patches);
    
    return deltas;
}

/** Update the beads with the passed bead IDs with the parameters in the
    passed array
    
    \throw SireError:::invalid_key
    \throw SireError::incompatible_error
*/
QHash<quint32,FFBeadChange> Patches::update(const QVector<quint32> &beadids, 
                                            const FFParametersArray &params)
{
    if (beadids.count() != params.count())
        throw SireError::program_bug( QObject::tr(
                "You must pass in the same number of parameter groups (%1) "
                "as you pass in bead IDs (%2).")
                    .arg(params.count()).arg(beadids.count()), CODELOC );
                

    Patches new_patches(*this);
    
    QHash<quint32,FFBeadChange> deltas;
    deltas.reserve(beadids.count());
    
    for (int i=0; i<beadids.count(); ++i)
    {
        FFBeadChange delta = new_patches.update( beadids.constData()[i], 
                                                 params[i].read() );
        
        if (not delta.isEmpty())
            deltas.insert(beadids.constData()[i], delta);
    }
    
    this->operator=(new_patches);
    
    return deltas;
}

/** Update the beads with the passed bead IDs with the new
    coordinates and parameters in the passed arrays
    
    \throw SireError:::invalid_key
    \throw SireError::incompatible_error
*/
QHash<quint32,FFBeadChange> Patches::update(const QVector<quint32> &beadids,
                                            const CoordGroupArray &coords,
                                      const FFParametersArray &params)
{
    if (beadids.count() != params.count() or
        beadids.count() != coords.count())
        throw SireError::program_bug( QObject::tr(
                "You must pass in the same number of parameter groups (%1) "
                "and coordinates groups (%2) as you pass in bead IDs (%3).")
                    .arg(params.count()).arg(coords.count())
                    .arg(beadids.count()), CODELOC );
                

    Patches new_patches(*this);
    
    QHash<quint32,FFBeadChange> deltas;
    deltas.reserve(beadids.count());
    
    for (int i=0; i<beadids.count(); ++i)
    {
        FFBeadChange delta = new_patches.update( beadids.constData()[i], 
                                                 coords.constData()[i],
                                                 params[i].read() );
                                              
        if (not delta.isEmpty())
            deltas.insert(beadids.constData()[i], delta);
    }
    
    this->operator=(new_patches);
    
    return deltas;
}

/** Remove the bead with ID 'beadid' */
FFBeadChange Patches::remove(quint32 beadid)
{
    int idx = beadid_to_patch.value(beadid, -1);
    
    if (idx != -1)
    {
        FFBead old_bead = ptchs.data()[idx].remove(beadid);
        beadid_to_patch.remove(beadid);
        return FFBeadChange(old_bead, FFBead());
    }
    else
        return FFBeadChange();
}

/** Remove the beads whose indicies are in 'beadids' */
QHash<quint32,FFBeadChange> Patches::remove(const QVector<quint32> &beadids)
{
    Patches new_patches(*this);
    
    QHash<quint32,FFBeadChange> deltas;
    deltas.reserve(beadids.count());
    
    for (int i=0; i<beadids.count(); ++i)
    {
        FFBeadChange delta = new_patches.remove( beadids.constData()[i] );
        
        if (not delta.isEmpty())
            deltas.insert(beadids.constData()[i], delta);
    }
    
    return deltas;
}

/** Remove all of the beads from the patches */
void Patches::removeAll()
{
    for (int i=0; i<ptchs.count(); ++i)
    {
        ptchs.data()[i].removeAll();
    }
    
    beadid_to_patch.clear();
    last_beadid = 0;
}
