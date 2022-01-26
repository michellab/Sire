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

#include "patch.h"
#include "patches.h"

#include "SireError/errors.h"

#include "tostring.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireFF;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Patch> r_patch;

QDataStream &operator<<(QDataStream &ds, const Patch &patch)
{
    writeHeader(ds, r_patch, 1);

    SharedDataStream sds(ds);

    sds << patch.coords << patch.params << patch.idx_to_beadid
        << patch.beadid_to_idx << static_cast<const Property&>(patch);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Patch &patch)
{
    VersionID v = readHeader(ds, r_patch);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> patch.coords >> patch.params >> patch.idx_to_beadid
            >> patch.beadid_to_idx >> static_cast<Property&>(patch);

        patch.aabox = patch.coords.aaBox();
    }
    else
        throw version_error(v, "1", r_patch, CODELOC);

    return ds;
}

/** Constructor */
Patch::Patch() : ConcreteProperty<Patch,Property>()
{}

/** Copy constructor */
Patch::Patch(const Patch &other)
      : ConcreteProperty<Patch,Property>(other),
        coords(other.coords), params(other.params),
        idx_to_beadid(other.idx_to_beadid), beadid_to_idx(other.beadid_to_idx),
        aabox(other.aabox)
{}

/** Destructor */
Patch::~Patch()
{}

const char* Patch::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Patch>() );
}

/** Copy assignment operator */
Patch& Patch::operator=(const Patch &other)
{
    if (this != &other)
    {
        coords = other.coords;
        params = other.params;
        idx_to_beadid = other.idx_to_beadid;
        beadid_to_idx = other.beadid_to_idx;
        Property::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool Patch::operator==(const Patch &other) const
{
    return this == &other or
           (coords == other.coords and params == other.params and
            idx_to_beadid == other.idx_to_beadid and
            beadid_to_idx == other.beadid_to_idx and
            Property::operator==(other));
}

/** Comparison operator */
bool Patch::operator!=(const Patch &other) const
{
    return not Patch::operator==(other);
}

/** Return whether or not this patch is empty */
bool Patch::isEmpty() const
{
    return coords.count() == 0;
}

/** Return the AABox that completely surrounds this patch */
const AABox& Patch::aaBox() const
{
    return aabox;
}

/** Return the index (location) of the patch with ID 'beadid'

    \throw SireError::invalid_key
*/
int Patch::getLocation(quint32 beadid) const
{
    int idx = beadid_to_idx.value(beadid, -1);

    if (idx == -1)
        throw SireError::invalid_key( QObject::tr(
                "There is no bead with ID %1 in this patch. Available beads "
                "are %1.")
                    .arg(beadid).arg(Sire::toString(beadid_to_idx.keys())),
                        CODELOC );

    return idx;
}

/** Return the internal IDs of each of the beads in this patch,
    in the order that they appear in the patch */
const QVector<quint32> Patch::beadIDs() const
{
    return idx_to_beadid;
}

/** Return the number of beads in this patch */
int Patch::nBeads() const
{
    return coords.nCoordGroups();
}

/** Return a string representation of this patch */
QString Patch::toString() const
{
    if (this->isEmpty())
        return QObject::tr("Patch::null");

    return QObject::tr("Patch( aaBox() == %1, nBeads() == %2 )")
                .arg(this->aaBox().toString()).arg(this->nBeads());
}

/** Return the array of coordinates of the beads in this patch */
const CoordGroupArray& Patch::coordinates() const
{
    return coords;
}

/** Return the array of parameters for each bead in the patch */
const FFParametersArray& Patch::parameters() const
{
    return params.read();
}

/** Internal function used to add the bead with internal ID 'beadid'
    and specified coordinates and parameters to this patch */
FFBead Patch::add(quint32 beadid, const CoordGroup &groupcoords,
                  const FFParameters &groupparams)
{
    if (beadid_to_idx.contains(beadid))
        //this bead already exists in this patch - remove it as it is
        //an old version
        this->remove(beadid);

    beadid_to_idx.insert(beadid, idx_to_beadid.count());
    idx_to_beadid.append(beadid);

    if (params.read().isEmpty())
    {
        //make this parameter array into an array of the right type
        params = groupparams.toArray();
    }
    else
    {
        params.edit().append(groupparams);
    }

    coords.append(groupcoords);

    aabox = coords.aaBox();

    return FFBead(groupcoords, groupparams);
}

/** Internal function used to add the array of beads with the passed
    array of beadids */
QHash<quint32,FFBead> Patch::add(const QVarLengthArray<quint32> &beadids,
                                 const CoordGroupArray &groupcoords,
                                 const FFParametersArray &groupparams)
{
    if (beadids.count() != groupcoords.count() or
        beadids.count() != groupparams.count())
    {
        throw SireError::incompatible_error( QObject::tr(
                "You must ensure that the number of CoordGroups (%1) "
                "is equal to the number of parameter groups (%2) and "
                "that is equal to the number of specified beads (%3).")
                    .arg(groupcoords.count())
                    .arg(groupparams.count())
                    .arg(beadids.count()), CODELOC );
    }

    for (int i=0; i<beadids.count(); ++i)
    {
        if (beadid_to_idx.contains(beadids.constData()[i]))
            this->remove(beadids.constData()[i]);
    }

    const int old_n_beads = idx_to_beadid.count();
    idx_to_beadid.resize( idx_to_beadid.count() + beadids.count() );
    quint32 *idx_to_beadid_data = idx_to_beadid.data() + old_n_beads;

    for (int i=0; i<beadids.count(); ++i)
    {
        beadid_to_idx.insert(beadids.constData()[i], beadids.count() + i);
        idx_to_beadid_data[i] = beadids.constData()[i];
    }

    if (params.read().isEmpty())
    {
        //make an array of the right type
        params = groupparams;
    }
    else
    {
        params.edit().append(groupparams);
    }

    coords.append(groupcoords);
    aabox = coords.aaBox();

    QHash<quint32,FFBead> beads;
    beads.reserve(beadids.count());

    for (int i=0; i<beadids.count(); ++i)
    {
        beads.insert(beadids[i], FFBead(groupcoords.at(i), groupparams.at(i)));
    }

    return beads;
}

/** Return the index of the bead with passed beadid

    \throw SireError::program_bug
*/
int Patch::getBeadIdx(quint32 beadid) const
{
    int idx = beadid_to_idx.value(beadid, -1);

    if (idx == -1)
        throw SireError::program_bug( QObject::tr(
                "There is a bug as there is no bead with ID %1 in this Patch!")
                    .arg(beadid), CODELOC );

    return idx;
}

/** Update the bead with ID 'beadid' with the new set of passed coordinates

    \throw SireError::incompatible_error
*/
FFBeadChange Patch::update(quint32 beadid, const CoordGroup &new_coords)
{
    int idx = getBeadIdx(beadid);

    const CoordGroup &old_coords = coords.at(idx);

    FFBeadChange delta;

    if (old_coords != new_coords)
    {
        FFParametersPtr beadparams = params.read().at(idx);

        delta = FFBeadChange( FFBead(old_coords,beadparams),
                              FFBead(new_coords,beadparams) );

        coords.update( getBeadIdx(beadid), new_coords );
        aabox = coords.aaBox();
    }

    return delta;
}

/** Update the bead with ID 'beadid' with the new set of passed parameters

    \throw SireError::incompatible_error
*/
FFBeadChange Patch::update(quint32 beadid, const FFParameters &new_params)
{
    int idx = getBeadIdx(beadid);

    FFBead old_bead( coords.at(idx), params.read().at(idx) );

    params.edit().update( getBeadIdx(beadid), new_params );

    return FFBeadChange(old_bead, FFBead(coords.at(idx),new_params));
}

/** Update the bead with ID 'beadid' with the new set of passed coordinates
    and parameters

    \throw SireError::incompatible_error
*/
FFBeadChange Patch::update(quint32 beadid, const CoordGroup &new_coords,
                                           const FFParameters &new_params)
{
    int idx = getBeadIdx(beadid);

    FFBead old_bead( coords.at(idx), params.read().at(idx) );

    //update a copy, so that the original is safe
    //if the parameters update causes an exception
    CoordGroupArray ncoords(coords);

    ncoords.update(idx, new_coords);
    params.edit().update(idx, new_params);
    coords = ncoords;
    aabox = coords.aaBox();

    return FFBeadChange(old_bead, FFBead(new_coords,new_params));
}

/** Update the specified beads with the specified coords

    \throw SireError::incompatible_error
*/
QHash<quint32,FFBeadChange> Patch::update(const QVarLengthArray<quint32> &beadids,
                                          const CoordGroupArray &new_coords)
{
    if (beadids.count() != new_coords.count())
        throw SireError::invalid_arg( QObject::tr(
                "You must pass in the same number of CoordGroups (%1) as there are "
                "beads specified in the array (%2).")
                    .arg(new_coords.count()).arg(beadids.count()), CODELOC );

    QHash<quint32,FFBeadChange> deltas;
    deltas.reserve(beadids.count());

    if (beadids.isEmpty())
        return deltas;

    else if (beadids.count() == 1)
    {
        FFBeadChange delta = this->update( beadids.constData()[0],
                                           new_coords.constData()[0] );

        if (not delta.isEmpty())
            deltas.insert( beadids.constData()[0], delta );
    }
    else
    {
        CoordGroupArray ncoords(coords);

        for (int i=0; i<beadids.count(); ++i)
        {
            int idx = getBeadIdx(beadids.constData()[i]);

            if (coords.at(idx) != new_coords.constData()[i])
            {
                ncoords.update( idx, new_coords.constData()[i] );

                FFParametersPtr beadparams = params.read().at(idx);

                deltas.insert( beadids.constData()[i], FFBeadChange(
                                  FFBead(coords.at(idx),beadparams),
                                  FFBead(new_coords.at(i),beadparams) ) );
            }
        }

        if (not deltas.isEmpty())
        {
            coords = ncoords;
            aabox = coords.aaBox();
        }
    }

    return deltas;
}

/** Update the specified beads with the specified forcefield parameters

    \throw SireError::incompatible_error
*/
QHash<quint32,FFBeadChange> Patch::update(const QVarLengthArray<quint32> &beadids,
                                          const FFParametersArray &new_params)
{
    if (beadids.count() != new_params.count())
        throw SireError::invalid_arg( QObject::tr(
                "You must pass in the same number of parameters (%1) as there are "
                "beads specified in the array (%2).")
                    .arg(new_params.count()).arg(beadids.count()), CODELOC );

    QHash<quint32,FFBeadChange> deltas;
    deltas.reserve(beadids.count());

    if (beadids.isEmpty())
        return deltas;

    else
    {
        QVarLengthArray<int> idxs(beadids.count());

        int *idxs_array = idxs.data();

        for (int i=0; i<beadids.count(); ++i)
        {
            int idx = getBeadIdx(beadids.constData()[i]);
            idxs_array[i] = idx;

            deltas.insert( beadids.constData()[i], FFBeadChange(
                              FFBead(coords.at(idx), params.read().at(idx)),
                              FFBead(coords.at(idx), new_params.at(i)) ) );
        }

        params.edit().update(idxs, new_params);

        return deltas;
    }
}

/** Update the coordinates and parameters of the specified beads

    \throw SireError::incompatible_error
*/
QHash<quint32,FFBeadChange> Patch::update(const QVarLengthArray<quint32> &beadids,
                                          const CoordGroupArray &new_coords,
                                          const FFParametersArray &new_params)
{
    if (beadids.count() != new_coords.count() or
        beadids.count() != new_params.count())
    {
        throw SireError::incompatible_error( QObject::tr(
                "You must ensure that the number of CoordGroups (%1) "
                "is equal to the number of parameter groups (%2) and "
                "that is equal to the number of specified beads (%3).")
                    .arg(new_coords.count())
                    .arg(new_params.count())
                    .arg(beadids.count()), CODELOC );
    }

    QHash<quint32,FFBeadChange> deltas;
    deltas.reserve(beadids.count());

    if (beadids.isEmpty())
        return deltas;

    else
    {
        QVarLengthArray<int> idxs(beadids.count());

        int *idxs_array = idxs.data();

        for (int i=0; i<beadids.count(); ++i)
        {
            int idx = getBeadIdx(beadids.constData()[i]);
            idxs_array[i] = idx;

            FFParametersPtr old_params = params.read().at(idx);

            deltas.insert( beadids.constData()[i],
                              FFBeadChange( FFBead(coords.at(idx),old_params),
                                            FFBead(new_coords.at(i),new_params.at(i)) ) );
        }

        CoordGroupArray ncoords(coords);

        for (int i=0; i<idxs.count(); ++i)
        {
            ncoords.update(idxs_array[i], new_coords.constData()[i]);
        }

        params.edit().update(idxs, new_params);
        coords = ncoords;
        aabox = coords.aaBox();

        return deltas;
    }
}

/** Remove the bead with ID 'beadid' */
FFBead Patch::remove(quint32 beadid)
{
    FFBead old_bead;

    int idx = getBeadIdx(beadid);
    {
        old_bead = FFBead(coords.at(idx), params.read().at(idx));

        CoordGroupArray ncoords(coords);
        ncoords.remove(idx);
        params.edit().remove(idx);
        coords = ncoords;
    }

    //update the index
    if (coords.isEmpty())
    {
        idx_to_beadid = QVector<quint32>();
        beadid_to_idx = QHash<quint32,int>();
        aabox = AABox();
    }
    else
    {
        beadid_to_idx.remove(beadid);
        idx_to_beadid.remove(idx);

        for (int i=0; i<idx_to_beadid.count(); ++i)
        {
            beadid_to_idx.insert( idx_to_beadid.constData()[i], i );
        }

        aabox = coords.aaBox();
    }

    return old_bead;
}

/** Remove all of the beads whose IDs are in 'beadids' */
QHash<quint32,FFBead> Patch::remove(const QVarLengthArray<quint32> &beadids)
{
    if (beadids.isEmpty())
        return QHash<quint32,FFBead>();

    QHash<quint32,FFBead> old_beads;
    old_beads.reserve(beadids.count());

    QVarLengthArray<int> idxs(beadids.count());

    int *idxs_array = idxs.data();

    for (int i=0; i<beadids.count(); ++i)
    {
        int idx = getBeadIdx(beadids.constData()[i]);
        old_beads.insert(beadids.constData()[i],
                         FFBead(coords.at(idx), params.read().at(idx)));

        idxs_array[i] = idx;
    }

    std::sort( idxs_array, idxs_array + idxs.count() );
    {
        CoordGroupArray ncoords(coords);

        QVarLengthArray<quint32> uidxs32(idxs.count());
        quint32 *uidxs32_array = uidxs32.data();

        for (int i=0; i<idxs.count(); ++i)
        {
            uidxs32_array[i] = idxs_array[i];
        }

        ncoords.remove(uidxs32);
        params.edit().remove(idxs);
        coords = ncoords;
    }

    //update the index
    if (coords.isEmpty())
    {
        idx_to_beadid = QVector<quint32>();
        beadid_to_idx = QHash<quint32,int>();
        aabox = AABox();
    }
    else
    {
        for (int i=0; i<beadids.count(); ++i)
        {
            beadid_to_idx.remove( beadids.constData()[i] );
        }

        int last_idx = -1;

        //have to remove from the end so that the indicies
        //to remove remain valid
        for (int i=idxs.count()-1; i>=0; --i)
        {
            if (idxs_array[i] != last_idx)
            {
                idx_to_beadid.remove(idxs_array[i]);
                last_idx = idxs_array[i];
            }
        }

        for (int i=0; i<idx_to_beadid.count(); ++i)
        {
            beadid_to_idx.insert( idx_to_beadid.constData()[i], i );
        }

        aabox = coords.aaBox();
    }

    return old_beads;
}

/** Remove all groups from this patch */
void Patch::removeAll()
{
    coords = CoordGroupArray();
    params.edit().removeAll();

    beadid_to_idx = QHash<quint32,int>();
    idx_to_beadid = QVector<quint32>();

    aabox = AABox();
}
