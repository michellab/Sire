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

#include "SireMM/cljboxes.h"
#include "SireMM/cljdelta.h"

#include "SireVol/aabox.h"
#include "SireVol/space.h"
#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include <QElapsedTimer>
#include <QDebug>

#include "tostring.h"

using namespace SireMM;
using namespace SireMaths;
using namespace SireBase;
using namespace SireVol;
using namespace SireStream;

///////////
/////////// Implementation of CLJBox
///////////

static const RegisterMetaType<CLJBox> r_cljbox(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const CLJBox &box)
{
    writeHeader(ds, r_cljbox, 1);

    SharedDataStream sds(ds);
    sds << box.atms << box.box_index << box.box_length;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJBox &box)
{
    VersionID v = readHeader(ds, r_cljbox);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> box.atms >> box.box_index >> box.box_length;
        box.findGaps();
    }
    else
        throw version_error(v, "1", r_cljbox, CODELOC);

    return ds;
}

/** Null constructor */
CLJBox::CLJBox() : RefCountData(), box_length(0)
{}

/** Construct an empty box at a specific location */
CLJBox::CLJBox(const CLJBoxIndex &index, Length length)
       : RefCountData(), box_index(index), box_length(length.value())
{}

/** Construct a box that holds the passed atoms */
CLJBox::CLJBox(const CLJBoxIndex &index, Length length, const CLJAtoms &atoms)
       : RefCountData(), atms(atoms), box_index(index), box_length(length.value())
{
    findGaps();
}

/** Copy constructor */
CLJBox::CLJBox(const CLJBox &other)
       : RefCountData(), atms(other.atms), gaps(other.gaps),
         box_index(other.box_index), box_length(other.box_length)
{}

/** Destructor */
CLJBox::~CLJBox()
{}

/** Search for all of the gaps (dummy atoms) in the box */
void CLJBox::findGaps()
{
    gaps.clear();

    const QVector<MultiInt> &ids = atms.ID();

    const quint32 dummy_id = CLJAtoms::idOfDummy()[0];

    for (int i=ids.count()-1; i>=0; --i)
    {
        const MultiInt &id = ids[i];

        for (int j=MultiInt::count()-1; j>=0; --j)
        {
            if (id[j] == dummy_id)
            {
                gaps.push( i*MultiInt::count() + j );
            }
        }
    }
}

QString CLJBox::toString() const
{
    return QObject::tr("CLJBox( nAtoms() == %1 )").arg(nAtoms());
}

/** Remove the atom at index 'atom' from the box */
void CLJBox::remove(int atom)
{
    if (atom < 0 or atom >= atms.count())
    {
        //this is an invalid atom
        return;
    }

    if (not atms.isDummy(atom))
    {
        atms.makeDummy(atom);
        gaps.push(atom);
    }
}

/** Remove the atoms whose indicies are in 'atoms' */
void CLJBox::remove(const QList<int> &atoms)
{
    foreach (int atom, atoms)
    {
        this->remove(atom);
    }
}

/** Remove the atoms at the specified indicies */
void CLJBox::remove(const QVector<CLJBoxIndex> &atoms)
{
    for (int i=0; i<atoms.count(); ++i)
    {
        const CLJBoxIndex &atom = atoms.constData()[i];

        if (atom.sameBox(box_index))
            this->remove(atom.index());
    }
}

/** Remove the atom at index 'atom', returning the atom removed */
CLJAtom CLJBox::take(int atom)
{
    if (atom < 0 or atom >= atms.count())
    {
        return CLJAtom();
    }

    if (not atms.isDummy(atom))
    {
        CLJAtom atm = atms[atom];
        atms.makeDummy(atom);
        gaps.push(atom);
        return atm;
    }
    else
        return CLJAtom();
}

/** Add a single passed atom into this box. This returns the index
    of the added atom */
CLJBoxIndex CLJBox::add(const CLJAtom &atom)
{
    if (atom.isDummy())
        return CLJBoxIndex();

    if (not gaps.isEmpty())
    {
        //sort the list of gaps so that we add from the bottom first
        std::sort(gaps.begin(), gaps.end(), std::greater<int>());

        int gap = gaps.pop();

        atms.set(gap, atom);
        return CLJBoxIndex( box_index.i(), box_index.j(), box_index.k(), gap );
    }
    else
    {
        //we need to append a new vector onto the box
        int idx = atms.count();
        atms.append(atom);

        for (int i=1; i<MultiFloat::count(); ++i)
        {
            gaps.push(idx + i);
        }

        return CLJBoxIndex( box_index.i(), box_index.j(), box_index.k(), idx );
    }
}

/** Add the passed atoms into this box. This returns the indicies
    of each added atom (in the same order as they were added) */
QVector<CLJBoxIndex> CLJBox::add(const CLJAtoms &atoms)
{
    if (atoms.isEmpty())
        return QVector<CLJBoxIndex>();

    //sort the list of gaps so that we add from the bottom first
    std::sort(gaps.begin(), gaps.end(), std::greater<int>());

    int n = atoms.nAtoms();

    if (n == 0)
        //nothing to add
        return QVector<CLJBoxIndex>();

    QVector<CLJBoxIndex> indicies(atoms.count(), CLJBoxIndex());

    //now we use 'n' as an index
    n -= 1;

    while (not gaps.isEmpty())
    {
        //take atoms from the top to fill in the gaps
        while (atoms.isDummy(n))
        {
            n -= 1;
            if (n < 0)
            {
                //all remaining atoms are dummies - nothing left to add
                return indicies;
            }
        }

        //this atom is not a dummy. Copy it and put it into the gap
        int gap = gaps.pop();
        atms.set(gap, atoms.at(n));
        indicies[n] = CLJBoxIndex(box_index.i(), box_index.j(), box_index.k(), gap);

        n -= 1;

        if (n < 0)
        {
            //all of the atoms have been added
            return indicies;
        }
    }

    //n was used to refer to index, so need to add 1 as we now
    //use n to refer to the number of atoms that need to still be added
    n += 1;

    //there are still atoms to add and there are no gaps to add them
    //Just add them onto the end of the vector
    int start = atms.count();
    atms.append( atoms, n );

    if (atms.isPadded())
    {
        //there is some padding
        for (int i=atms.nAtoms(); i<atms.count(); ++i)
        {
            gaps.push(i);
        }
    }

    for (int i = 0; i<n; ++i)
    {
        indicies[i] = CLJBoxIndex(box_index.i(), box_index.j(), box_index.k(), start + i);
    }

    return indicies;
}

/** Combine the atoms of the two passed boxes */
CLJBox CLJBox::operator+(const CLJBox &other) const
{
    if (box_index != other.box_index or box_length != other.box_length)
        throw SireError::incompatible_error( QObject::tr(
                "You cannot add together two boxes that are at different points in space! "
                "Box0 = %1, dimensions %2, while Box1 = %3, dimensions %4.")
                    .arg(index().toString()).arg(dimensions().toString())
                    .arg(other.index().toString()).arg(other.dimensions().toString()),
                        CODELOC );

    if (this->isEmpty())
        return other;
    else if (other.isEmpty())
        return *this;
    else
    {
        CLJBox ret(*this);
        ret.atms += other.atms;
        ret.findGaps();
        return ret;
    }
}

/** Copy assignment operator */
CLJBox& CLJBox::operator=(const CLJBox &other)
{
    atms = other.atms;
    box_index = other.box_index;
    box_length = other.box_length;
    return *this;
}

/** Comparison operator */
bool CLJBox::operator==(const CLJBox &other) const
{
    return box_length == other.box_length and box_index == other.box_index and atms == other.atms;
}

/** Comparison operator */
bool CLJBox::operator!=(const CLJBox &other) const
{
    return not operator==(other);
}

/** Return the number of atoms in the box. This is equal to the
    number of actual atoms (i.e. not including padding or dummy atoms) */
int CLJBox::nAtoms() const
{
    return atms.nAtoms();
}

/** Return the count of the box - this includes dummy atoms and padding */
int CLJBox::count() const
{
    return atms.count();
}

/** Return the count of the box - this includes dummy atoms and padding */
int CLJBox::size() const
{
    return atms.size();
}

/** Return the ith CLJAtom - this uses the CLJAtoms index, i.e. includes
    dummy atoms and padding */
CLJAtom CLJBox::operator[](int i) const
{
    return atms[i];
}

/** Return the ith CLJAtom - this uses the CLJAtoms index, i.e. includes
    dummy atoms and padding */
CLJAtom CLJBox::at(int i) const
{
    return atms.at(i);
}

/** Return the ith CLJAtom - this uses the CLJAtoms index, i.e. includes
    dummy atoms and padding */
CLJAtom CLJBox::getitem(int i) const
{
    return atms.getitem(i);
}

const char* CLJBox::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJBox>() );
}

const char* CLJBox::what() const
{
    return CLJBox::typeName();
}

/** Return the dimensions of this box */
AABox CLJBox::dimensions() const
{
    return box_index.box( Length(box_length) );
}

/** Return a copy of this box where the CLJAtoms are squeezed */
CLJBox CLJBox::squeeze() const
{
    if (gaps.isEmpty())
        return *this;

    //see if all of the gaps are at the end of the vector
    //(because of padding the MultiFloat values)
    const int natoms = atms.count() - atms.nPadded();

    bool all_padding = true;

    foreach (const int gap, gaps)
    {
        if (gap < natoms)
        {
            all_padding = false;
            break;
        }
    }

    if (all_padding)
        //nothing needs to be done
        return *this;

    //we need to squeeze out the gaps
    CLJBox ret(*this);
    ret.atms = atms.squeeze();
    ret.findGaps();
    return ret;
}

///////////
/////////// Implementation of CLJBoxPtr
///////////

QDataStream &operator<<(QDataStream &ds, const CLJBoxPtr &ptr)
{
    SharedDataStream sds(ds);
    sds << ptr.d;
    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJBoxPtr &ptr)
{
    SharedDataStream sds(ds);
    sds >> ptr.d;
    return ds;
}

static const CLJBoxPtr *null_box = 0;

/** Constructor */
CLJBoxPtr::CLJBoxPtr()
{
    if (null_box == 0)
    {
        null_box = new CLJBoxPtr( new CLJBox() );
    }

    this->operator=(*null_box);
}

/** Construct from the passed box. This takes over ownership of the pointer,
    and will delete it once there are no more CLJBoxPtrs pointing to this box */
CLJBoxPtr::CLJBoxPtr(CLJBox *box)
{
    if (box == 0)
    {
        if (null_box == 0)
        {
            null_box = new CLJBoxPtr( new CLJBox() );
        }

        this->operator=(*null_box);
    }
    else
    {
        d = box;
    }
}

/** Construct from the passed box */
CLJBoxPtr::CLJBoxPtr(const CLJBox &box) : d(new CLJBox(box))
{}

/** Copy constructor */
CLJBoxPtr::CLJBoxPtr(const CLJBoxPtr &other) : d(other.d)
{}

/** Destructor */
CLJBoxPtr::~CLJBoxPtr()
{
    //this will delete the box if it is the last
    //reference to the box
}

/** Copy assignment operator */
CLJBoxPtr& CLJBoxPtr::operator=(const CLJBoxPtr &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
CLJBoxPtr& CLJBoxPtr::operator=(const CLJBox &other)
{
    if (d.constData() == &other)
    {
        return *this;
    }

    if (d.constData() == 0)
    {
        d = new CLJBox(other);
    }
    else
    {
        d.data()->operator=(other);
    }

    return *this;
}


/** Comparison operator */
bool CLJBoxPtr::operator==(const CLJBoxPtr &other) const
{
    return d.constData() == other.d.constData() or
           *(d.constData()) == *(other.d.constData());
}

/** Comparison operator */
bool CLJBoxPtr::operator!=(const CLJBoxPtr &other) const
{
    return not operator==(other);
}

///////////
/////////// Implementation of CLJBoxIndex
///////////

static const RegisterMetaType<CLJBoxIndex> r_cljboxindex(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const CLJBoxIndex &index)
{
    writeHeader(ds, r_cljboxindex, 1);
    ds << index.v.index.ii << index.v.index.jj << index.v.index.kk << index.v.index.idx;
    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJBoxIndex &index)
{
    VersionID v = readHeader(ds, r_cljboxindex);

    if (v == 1)
    {
        ds >> index.v.index.ii >> index.v.index.jj
           >> index.v.index.kk >> index.v.index.idx;
    }
    else
        throw version_error(v, "1", r_cljboxindex, CODELOC);

    return ds;
}

/** Null constructor */
CLJBoxIndex::CLJBoxIndex()
{
    v.index.ii = std::numeric_limits<qint16>::min();
    v.index.jj = v.index.ii;
    v.index.kk = v.index.ii;
    v.index.idx = v.index.ii;
}

/** Construct the index of the box at index i,j,k with (optionally supplied)
    index of a particular atom in the box */
CLJBoxIndex::CLJBoxIndex(qint16 ii, qint16 jj, qint16 kk, qint16 atom_idx)
{
    v.index.ii = ii;
    v.index.jj = jj;
    v.index.kk = kk;
    v.index.idx = atom_idx;
}

/** Copy constructor */
CLJBoxIndex::CLJBoxIndex(const CLJBoxIndex &other)
{
    v.val = other.v.val;
}

/** Destructor */
CLJBoxIndex::~CLJBoxIndex()
{}

/** Copy assignment operator */
CLJBoxIndex& CLJBoxIndex::operator=(const CLJBoxIndex &other)
{
    v.val = other.v.val;
    return *this;
}

/** Return a null CLJBoxIndex */
CLJBoxIndex CLJBoxIndex::null()
{
    return CLJBoxIndex();
}

/** Return whether or not this is null */
bool CLJBoxIndex::isNull() const
{
    return this->operator==( CLJBoxIndex() );
}

const char* CLJBoxIndex::typeName()
{
    return "SireMM::CLJBoxIndex";
}

const char* CLJBoxIndex::what() const
{
    return CLJBoxIndex::typeName();
}

/** Return the minimum box indicies of the two passed boxes */
CLJBoxIndex CLJBoxIndex::min(const CLJBoxIndex &other) const
{
    return CLJBoxIndex( qMin(i(), other.i()), qMin(j(), other.j()),
                        qMin(k(), other.k()), qMin(index(), other.index()) );
}

/** Return the maximum box indicies of the two passed boxes */
CLJBoxIndex CLJBoxIndex::max(const CLJBoxIndex &other) const
{
    return CLJBoxIndex( qMax(i(), other.i()), qMax(j(), other.j()),
                        qMax(k(), other.k()), qMax(index(), other.index()) );
}

/** Return the AABox that describes this box (for a given box length of 'box_length').
    The boxes are arranged so that the box at (0,0,0) has its center at (0,0,0) and
    extends to (-0.5*length,-0.5*length,-0.5*length) to (0.5*length,0.5*length,0.5*length) */
AABox CLJBoxIndex::box(Length box_length) const
{
    Vector length(box_length.value());

    Vector origin( i() * length.x(),
                   j() * length.y(),
                   k() * length.z() );

    return AABox(origin, 0.5*length);
}

/** Return a copy of this index that contains only the box index (not the atom index) */
CLJBoxIndex CLJBoxIndex::boxOnly() const
{
    return CLJBoxIndex( i(), j(), k() );
}

QString CLJBoxIndex::toString() const
{
    if (isNull())
        return QObject::tr("CLJBoxIndex::null()");
    else
        return QObject::tr("CLJBoxIndex( %1, %2, %3 : %4 )")
                    .arg(v.index.ii)
                    .arg(v.index.jj)
                    .arg(v.index.kk)
                    .arg(v.index.idx);
}

/** Create the index for the box that contains the point 'x,y,z' in a set of boxes
    of length 1 / inv_box_length */
CLJBoxIndex CLJBoxIndex::createWithInverseBoxLength(float x, float y, float z, float inv_length)
{
    int i = int( std::floor(x * inv_length + 0.5) );  //std::floor will round the floating point
    int j = int( std::floor(y * inv_length + 0.5) );  //down to the nearest integer. Adding 0.5
    int k = int( std::floor(z * inv_length + 0.5) );  //ensures we instead round to nearest int

    const int min16 = std::numeric_limits<qint16>::min();
    const int max16 = std::numeric_limits<qint16>::max();

    if (i < min16 or i > max16 or
        j < min16 or j > max16 or
        k < min16 or k > max16)
    {
        throw SireError::program_bug( QObject::tr(
                "It is not possible to get the index of the box containing the point "
                "(%1,%2,%3) when the boxes have length %4 A, as the index lies outside "
                "the bounds of a 16bit integer (%5 to %6), (%7,%8,%9).")
                    .arg(x).arg(y).arg(z)
                    .arg(1 / inv_length)
                    .arg(min16).arg(max16)
                    .arg(i).arg(j).arg(k), CODELOC );
    }

    return CLJBoxIndex(i, j, k);
}

/** Create the index for the box that contains the point 'x,y,z' in a set of boxes
    of length 1 / inv_box_length */
CLJBoxIndex CLJBoxIndex::createWithInverseBoxLength(const Vector &coords, float inv_length)
{
    return createWithInverseBoxLength(coords.x(), coords.y(), coords.z(), inv_length);
}

/** Create the index for the box that contains the point 'x,y,z' in a set of boxes
    of length 'box_length' */
CLJBoxIndex CLJBoxIndex::createWithBoxLength(float x, float y, float z, Length box_length)
{
    return CLJBoxIndex::createWithInverseBoxLength(x, y, z, 1.0/box_length.value());
}

/** Create the index for the box that contains the point 'x,y,z' in a set of boxes
    of length 'box_length' */
CLJBoxIndex CLJBoxIndex::createWithBoxLength(const Vector &coords, Length box_length)
{
    return CLJBoxIndex::createWithInverseBoxLength(coords.x(), coords.y(),
                                                   coords.z(), 1.0/box_length.value());
}

/** Return the number of non-dummy indicies in the passed array */
int CLJBoxIndex::countNonDummies(const QVector<CLJBoxIndex> &indicies)
{
    int n = 0;

    for (int i=0; i<indicies.count(); ++i)
    {
        if (not indicies.constData()[i].isNull())
            n += 1;
    }

    return n;
}

///////////
/////////// Implementation of CLJBoxDistance
///////////

static const RegisterMetaType<CLJBoxDistance> r_boxdist(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const CLJBoxDistance &boxdist)
{
    writeHeader(ds, r_boxdist, 1);

    ds << boxdist.b0 << boxdist.b1 << boxdist.dist;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJBoxDistance &boxdist)
{
    VersionID v = readHeader(ds, r_boxdist);

    if (v == 1)
    {
        ds >> boxdist.b0 >> boxdist.b1 >> boxdist.dist;
    }
    else
        throw version_error(v, "1", r_boxdist, CODELOC);

    return ds;
}

/** Constructor */
CLJBoxDistance::CLJBoxDistance() : dist(0)
{}

/** Construct saying that the minimum distance between the box with index 'box0'
    and the box with index 'box1' is 'distance' */
CLJBoxDistance::CLJBoxDistance(quint32 box0, quint32 box1, float distance)
               : dist(distance), b0(box0), b1(box1)
{}

/** Copy constructor */
CLJBoxDistance::CLJBoxDistance(const CLJBoxDistance &other)
               : dist(other.dist), b0(other.b0), b1(other.b1)
{}

/** Destructor */
CLJBoxDistance::~CLJBoxDistance()
{}

/** Copy assignment operator */
CLJBoxDistance& CLJBoxDistance::operator=(const CLJBoxDistance &other)
{
    if (this != &other)
    {
        b0 = other.b0;
        b1 = other.b1;
        dist = other.dist;
    }

    return *this;
}

/** Comparison operator. Note that this will compare equal only if
    the distances are the same and if either, box0() == other.box() and
    box1() == other.box1() */
bool CLJBoxDistance::operator==(const CLJBoxDistance &other) const
{
    return dist == other.dist and
           b0 == other.b0 and b1 == other.b1;
}

/** Comparison operator */
bool CLJBoxDistance::operator!=(const CLJBoxDistance &other) const
{
    return not operator==(other);
}

/** Distance comparison operator - compares less if the distance
    is less than 'other'. This is used to allow CLJBoxDistance objects
    to be sorted according to their minimum distances */
bool CLJBoxDistance::operator<(const CLJBoxDistance &other) const
{
    return dist < other.dist;
}

/** Distance comparison operator - compares greater if the distance
    is greater than 'other'. This is used to allow CLJBoxDistance objects
    to be sorted according to their minimum distances */
bool CLJBoxDistance::operator>(const CLJBoxDistance &other) const
{
    return dist > other.dist;
}

const char* CLJBoxDistance::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJBoxDistance>() );
}

QString CLJBoxDistance::toString() const
{
    return QObject::tr("CLJBoxDistance( %1->%2, %3 A )")
                .arg(b0).arg(b1).arg(dist);
}

const char* CLJBoxDistance::what() const
{
    return CLJBoxDistance::typeName();
}


///////////
/////////// Implementation of CLJBoxes
///////////

static const RegisterMetaType<CLJBoxes> r_cljboxes(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const CLJBoxes &boxes)
{
    writeHeader(ds, r_cljboxes, 1);

    SharedDataStream sds(ds);

    sds << boxes.bxs << boxes.box_length;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJBoxes &boxes)
{
    VersionID v = readHeader(ds, r_cljboxes);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> boxes.bxs >> boxes.box_length;

        //reconstruct the index of box_index to integer
        boxes.box_to_idx.clear();

        for (int i=0; i<boxes.bxs.count(); ++i)
        {
            boxes.box_to_idx.insert( boxes.bxs.constData()[i].read().index(), i );
        }
    }
    else
        throw version_error(v, "1", r_cljboxes, CODELOC);

    return ds;
}

const double default_box_length = 10;

/** Null constructor */
CLJBoxes::CLJBoxes() : box_length(default_box_length)
{}

/** Construct, specifying the box length */
CLJBoxes::CLJBoxes(Length size) : box_length(size)
{
    if (size.value() <= 2)
        //don't be silly!
        box_length = 2.0;
}

/** Construct from the passed set of atoms */
void CLJBoxes::constructFrom(const CLJAtoms &atoms0, const CLJAtoms &atoms1)
{
    /*QElapsedTimer t;
    t.start();*/

    const float inv_length = 1.0 / box_length;

    //first check if all of the atoms are in the same box...
    if (atoms0.count() < 50 and atoms1.count() < 50)
    {
        CLJBoxIndex same_box;
        bool in_same_box = true;

        for (int i=0; i<atoms0.count(); ++i)
        {
            const CLJAtom atom = atoms0[i];

            if (atom.ID() != 0)
            {
                CLJBoxIndex index = CLJBoxIndex::createWithInverseBoxLength(
                                                            atom.coordinates(), inv_length);

                if (same_box.isNull())
                    same_box = index;
                else if (same_box != index)
                {
                    in_same_box = false;
                    break;
                }
            }
        }

        if (in_same_box)
        {
            for (int i=0; i<atoms1.count(); ++i)
            {
                const CLJAtom atom = atoms1[i];

                if (atom.ID() != 0)
                {
                    CLJBoxIndex index = CLJBoxIndex::createWithInverseBoxLength(
                                                                atom.coordinates(), inv_length);

                    if (same_box.isNull())
                        same_box = index;
                    else if (same_box != index)
                    {
                        in_same_box = false;
                        break;
                    }
                }
            }
        }

        if (in_same_box)
        {
            bxs.resize(1);
            box_to_idx.clear();
            box_to_idx.insert(same_box.boxOnly(), 0);

            const Length l(box_length);

            //just copy both sets of atoms into the same box
            if (atoms1.isEmpty())
                bxs[0] = CLJBox(same_box, l, atoms0);
            else if (atoms0.isEmpty())
                bxs[0] = CLJBox(same_box, l, atoms1);
            else
                bxs[0] = CLJBox(same_box, l, atoms0+atoms1);

            return;
        }
    }

    QHash< CLJBoxIndex,QList<CLJAtom> > boxed_atoms;

    for (int i=0; i<atoms0.count(); ++i)
    {
        const CLJAtom atom = atoms0[i];

        if (atom.ID() != 0)
        {
            CLJBoxIndex cljindex = CLJBoxIndex::createWithInverseBoxLength(
                                                        atom.coordinates(), inv_length);

            boxed_atoms[cljindex].append(atom);
        }
    }

    for (int i=0; i<atoms1.count(); ++i)
    {
        const CLJAtom atom = atoms1[i];

        if (atom.ID() != 0)
        {
            CLJBoxIndex cljindex = CLJBoxIndex::createWithInverseBoxLength(
                                                        atom.coordinates(), inv_length);

            boxed_atoms[cljindex].append(atom);
        }
    }

    //now build the CLJAtoms for each box
    box_to_idx.clear();
    bxs.resize(boxed_atoms.count());
    int i = 0;

    const Length l(box_length);

    for (QHash< CLJBoxIndex,QList<CLJAtom> >::const_iterator it = boxed_atoms.constBegin();
         it != boxed_atoms.constEnd();
         ++it)
    {
        //qDebug() << "Box" << it.key().box(Length(box_length)).toString() << it.value().count();
        bxs[i] = CLJBox( it.key(), l, CLJAtoms(it.value()) );
        box_to_idx.insert( it.key(), i );
        i += 1;
    }

    /*quint64 ns = t.nsecsElapsed();
    qDebug() << "Boxing up" << this->nAtoms() << "atoms took"
             << (0.000001*ns) << "ms";

    qDebug() << "number of boxes ==" << this->nOccupiedBoxes();*/
}

/** Add together two boxes - both boxes must have the same box length */
CLJBoxes CLJBoxes::operator+(const CLJBoxes &other) const
{
    if (box_length != other.box_length)
    {
        throw SireError::incompatible_error( QObject::tr(
                "You cannot add together two CLJBoxes objects that have different "
                "box lengths (%1 A vs. %2 A)")
                    .arg(box_length).arg(other.box_length), CODELOC );
    }

    if (this->isEmpty())
        return other;
    else if (other.isEmpty())
        return *this;

    CLJBoxes ret(*this);

    for (ContainerMap::const_iterator it = other.box_to_idx.constBegin();
         it != other.box_to_idx.constEnd();
         ++it)
    {
        int idx = box_to_idx.value(it.key(), -1);

        if (idx >= 0)
        {
            ret.bxs[idx] = bxs[idx].read() + other.bxs[it.value()].read();
        }
        else
        {
            ret.bxs.append( other.bxs[it.value()] );
            ret.box_to_idx.insert( it.key(), ret.bxs.count()-1 );
        }
    }

    return ret;
}

/** Box up the passed set of atoms into boxes of the default box size
    (5 angstroms) */
CLJBoxes::CLJBoxes(const CLJAtoms &atoms) : box_length(default_box_length)
{
    constructFrom(atoms, CLJAtoms());
}

/** Box up the passed set of atoms into boxes of specified box size */
CLJBoxes::CLJBoxes(const CLJAtoms &atoms, Length box_size)
         : box_length(box_size.value())
{
    if (box_length < 2)
        //don't be silly
        box_length = 2.0;

    constructFrom(atoms, CLJAtoms());
}

/** Box up the passed set of atoms into boxes of the default box size
    (5 angstroms) */
CLJBoxes::CLJBoxes(const CLJAtoms &atoms0, const CLJAtoms &atoms1) : box_length(default_box_length)
{
    constructFrom(atoms0, atoms1);
}

/** Box up the passed set of atoms into boxes of specified box size */
CLJBoxes::CLJBoxes(const CLJAtoms &atoms0, const CLJAtoms &atoms1, Length box_size)
         : box_length(box_size.value())
{
    if (box_length < 2)
        //don't be silly
        box_length = 2.0;

    constructFrom(atoms0, atoms1);
}

/** Copy constructor */
CLJBoxes::CLJBoxes(const CLJBoxes &other)
         : box_to_idx(other.box_to_idx), bxs(other.bxs), box_length(other.box_length)
{}

/** Destructor */
CLJBoxes::~CLJBoxes()
{}

/** Copy assignment operator */
CLJBoxes& CLJBoxes::operator=(const CLJBoxes &other)
{
    if (this != &other)
    {
        box_to_idx = other.box_to_idx;
        bxs = other.bxs;
        box_length = other.box_length;
    }

    return *this;
}

/** Comparison operator */
bool CLJBoxes::operator==(const CLJBoxes &other) const
{
    return box_length == other.box_length and bxs == other.bxs;
}

/** Comparison operator */
bool CLJBoxes::operator!=(const CLJBoxes &other) const
{
    return not operator==(other);
}

/** Return the atom at the specified index, or a null atom if
    none such atom exists */
CLJAtom CLJBoxes::operator[](const CLJBoxIndex &idx) const
{
    if (idx.isNull())
        return CLJAtom();

    ContainerMap::const_iterator it = box_to_idx.constFind(idx.boxOnly());

    if (it == box_to_idx.constEnd())
        return CLJAtom();

    const CLJBox &box = bxs.constData()[it.value()].read();

    if (idx.index() < 0 or idx.index() >= box.count())
        return CLJAtom();
    else
        return box.at( idx.index() );
}

/** Return the atom at the specified index, or a null atom if
    none such atom exists */
CLJAtom CLJBoxes::at(const CLJBoxIndex &idx) const
{
    return this->operator[](idx);
}

/** Return the atom at the specified index, or a null atom if
    none such atom exists */
CLJAtom CLJBoxes::getitem(const CLJBoxIndex &idx) const
{
    return this->operator[](idx);
}

QString CLJBoxes::toString() const
{
    return QObject::tr("CLJBoxes( nAtoms() == %1, nOccupiedBoxes() == %2 )")
                .arg(nAtoms()).arg(nOccupiedBoxes());
}

/** Return the indicies of all occupied boxes */
QVector<CLJBoxIndex> CLJBoxes::occupiedBoxIndicies() const
{
    return box_to_idx.keys().toVector();
}

/** Return the ith box */
CLJBox CLJBoxes::boxAt(int i) const
{
    if (i >= 0 and i < bxs.count())
        return bxs.constData()[i].read();
    else
        return CLJBox();
}

/** Return the box at index 'index' */
CLJBox CLJBoxes::boxAt(const CLJBoxIndex &index) const
{
    int idx = box_to_idx.value(index.boxOnly());

    if (idx >= 0)
        return bxs.constData()[idx].read();
    else
        return CLJBox();
}

/** Return the dimensions of the ith box */
AABox CLJBoxes::boxDimensionsAt(int i) const
{
    return this->boxAt(i).dimensions();
}

/** Return the AABox that describes the boundary of the box at index 'index' */
AABox CLJBoxes::boxDimensionsAt(const CLJBoxIndex &index) const
{
    return index.box( Length(box_length) );
}

/** Return the box that contains the point with coordinates 'coords' */
CLJBox CLJBoxes::boxAt(const Vector &coords) const
{
    return this->boxAt( CLJBoxIndex::createWithInverseBoxLength(coords, 1.0/box_length) );
}

/** Return the dimensions of the box that contains the point with coordinates 'coords' */
AABox CLJBoxes::boxDimensionsAt(const Vector &coords) const
{
    return this->boxDimensionsAt( CLJBoxIndex::createWithInverseBoxLength(coords, 1.0/box_length));
}

/** Return an array containing all occupied boxes. This is in the same order
    as the box dimensions returned by 'boxDimensions' */
QVector<CLJBox> CLJBoxes::boxes() const
{
    QVector<CLJBox> b;
    b.reserve(bxs.count());

    for (int i=0; i<bxs.count(); ++i)
    {
        b.append( bxs.constData()[i].read() );
    }

    return b;
}

/** Return an array containing all occupied box dimensions. This is in the same order
    as the boxes returned by 'boxes' */
QVector<AABox> CLJBoxes::boxDimensions() const
{
    if (bxs.isEmpty())
        return QVector<AABox>();

    QVector<AABox> b( bxs.count() );
    AABox *ba = b.data();

    for (int i=0; i<bxs.count(); ++i)
    {
        ba[i] = bxs.constData()[i].read().dimensions();
    }

    return b;
}

/** Return the number of occupied boxes */
int CLJBoxes::nOccupiedBoxes() const
{
    return bxs.count();
}

/** Return the number of atoms in the boxes */
int CLJBoxes::nAtoms() const
{
    int n = 0;

    for (const_iterator it = bxs.constBegin();
         it != bxs.constEnd();
         ++it)
    {
        n += it->read().nAtoms();
    }

    return n;
}

/** Return all of the atoms in all of the boxes (these may
    be returned with a lot of padding) */
CLJAtoms CLJBoxes::atoms() const
{
    CLJAtoms atms;

    for (const_iterator it = bxs.constBegin();
         it != bxs.constEnd();
         ++it)
    {
        atms += it->read().atoms();
    }

    return atms;
}

/** Return all of the atoms whose indicies are in 'idxs'. The atoms are returned
    in the same order as they appear in 'idxs' */
CLJAtoms CLJBoxes::atoms(const QVector<CLJBoxIndex> &idxs) const
{
    CLJAtoms ret;
    ret.resize(idxs.count());

    for (int i=0; i<idxs.count(); ++i)
    {
        ret.set(i, this->at(idxs.constData()[i]));
    }

    return ret;
}

/** Return all of the atoms whose indicies are in 'idxs'. The atoms are returned
    in the same order as they appear in 'idxs' */
CLJAtoms CLJBoxes::get(const QVector<CLJBoxIndex> &idxs) const
{
    return this->atoms(idxs);
}

/** Add a set of CLJAtoms to the box, returning the indicies of each added atom */
QVector<CLJBoxIndex> CLJBoxes::add(const CLJAtoms &atoms)
{
    if (atoms.isEmpty())
        return QVector<CLJBoxIndex>();

    QVector<CLJBoxIndex> indicies(atoms.count(), CLJBoxIndex::null());

    const float inv_length = 1.0 / box_length;

    QHash< CLJBoxIndex, tuple< QList<CLJAtom>, QList<int> > > boxed_atoms;

    for (int i=0; i<atoms.count(); ++i)
    {
        const CLJAtom atom = atoms[i];

        if (atom.ID() != 0)
        {
            CLJBoxIndex cljindex = CLJBoxIndex::createWithInverseBoxLength(
                                                        atom.coordinates(), inv_length);

            boxed_atoms[cljindex].get<0>().append(atom);
            boxed_atoms[cljindex].get<1>().append(i);
        }
    }

    //now build the CLJAtoms for each box
    for (QHash< CLJBoxIndex, tuple< QList<CLJAtom>,QList<int> > >::iterator
                                                it = boxed_atoms.begin();
         it != boxed_atoms.end();
         ++it)
    {
        int idx = box_to_idx.value(it.key(), -1);

        if (idx >= 0)
        {
            QVector<CLJBoxIndex> box_idxs = bxs[idx].write().add( CLJAtoms(it.value().get<0>()) );
            const QList<int> &atom_idxs = it.value().get<1>();

            for (int i=0; i<atom_idxs.count(); ++i)
            {
                indicies[ atom_idxs.at(i) ] = box_idxs.at(i);
            }
        }
        else
        {
            const QList<int> &atom_idxs = it.value().get<1>();
            bxs.append( new CLJBox(it.key(), Length(box_length), CLJAtoms(it.value().get<0>())) );
            box_to_idx.insert( it.key(), bxs.count() - 1 );

            for (int i=0; i<atom_idxs.count(); ++i)
            {
                indicies[ atom_idxs.at(i) ] = CLJBoxIndex( it.key().i(), it.key().j(),
                                                           it.key().k(), i );
            }
        }
    }

    return indicies;
}

/** Remove the atoms at the specified indicies. This does a rapid remove, i.e.
    it just turns the specified atoms into dummies (which may be overwritten by
    subsequent "add" operations). If you want to completely remove the atoms then
    use "remove" followed by "squeeze". This will turn the atoms into dummies and will
    then remove all dummy atoms from the boxes */
void CLJBoxes::remove(const QVector<CLJBoxIndex> &atoms)
{
    foreach (const CLJBoxIndex &atom, atoms)
    {
        if (atom.hasAtomIndex())
        {
            int idx = box_to_idx.value( atom.boxOnly(), -1 );

            if (idx >= 0)
            {
                bxs[idx].write().remove(atom.index());
            }
        }
    }
}

/** Remove the atoms at the specified indicies, returning the atoms that are
    removed. This does a rapid remove, i.e. it just turns the specified atoms into
    dummies (which may be overwritten by subsequent "add" operations). If you want
    to completely remove the atoms, use "take" followed by "squeeze" */
CLJAtoms CLJBoxes::take(const QVector<CLJBoxIndex> &atoms)
{
    QList<CLJAtom> atms;

    foreach (const CLJBoxIndex &atom, atoms)
    {
        if (atom.hasAtomIndex())
        {
            int idx = box_to_idx.value(atom.boxOnly(), -1);

            if (idx >= 0)
            {
                CLJAtom cljatom = bxs[idx].write().take(atom.index());

                if (not cljatom.isDummy())
                    atms.append(cljatom);
            }
        }
    }

    return CLJAtoms(atms);
}

/** Return a copy of the boxes where all of the CLJAtoms objects have been squeezed,
    and all empty boxes have been removed */
CLJBoxes CLJBoxes::squeeze() const
{
    CLJBoxes ret(*this);

    ret.bxs.clear();
    ret.box_to_idx.clear();

    for (const_iterator it = bxs.constBegin();
         it != bxs.constEnd();
         ++it)
    {
        if (it->read().nAtoms() != 0)
        {
            ret.bxs.append( it->read().squeeze() );
            ret.box_to_idx.insert( ret.bxs.last().read().index(), ret.bxs.count() - 1 );
        }
    }

    return ret;
}

/** Return the length of each side of each box */
Length CLJBoxes::length() const
{
    return Length(box_length);
}

const char* CLJBoxes::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJBoxes>() );
}

const char* CLJBoxes::what() const
{
    return CLJBoxes::typeName();
}

// this is the collection of square roots of integers between 0 and 31
// This is used as a lookup table for the square_root function, which
// we know will always have integer distances
static const int nsqrts = 32;
static const float sqrts[] = { 0, 1, 1.41421356237, 1.73205080757,
                               2.0, 2.2360679775, 2.44948974278,
                               2.64575131106, 2.82842712475,
                               3.0, 3.16227766017, 3.31662479036,
                               3.46410161514, 3.60555127546,
                               3.74165738677, 3.87298334621,
                               4.0, 4.12310562562, 4.24264068712,
                               4.35889894354, 4.472135955, 4.58257569496,
                               4.69041575982, 4.79583152331, 4.89897948557,
                               5.0, 5.09901951359, 5.19615242271, 5.29150262213,
                               5.38516480713, 5.47722557505, 5.56776436283 };

/** Calculate the square root of the passed integer */
SIRE_ALWAYS_INLINE float square_root(int delta2)
{
    if (delta2 < nsqrts)
    {
        return sqrts[delta2];
    }
    else
    {
        return std::sqrt( float(delta2) );
    }
}

/** Get the number of box lengths separating boxes 0 and 1 at indicies
    i0 and i1 */
SIRE_ALWAYS_INLINE int getDelta( const qint16 i0, const qint16 i1 )
{
    // the below single-line expression is doing
    // if (i0 == i1)
    // {
    //     return 0;
    // }
    // else if (i0 < i1)
    // {
    //     return i1 - i0 - 1;
    // }
    // else
    // {
    //     return i0 - i1 - 1;
    // }
    //
    // I found on my mac that the single line compiled to faster code
    // than the above if statements.

    return (i0 == i1) ? 0 : (i0 < i1) ? (i1-i0-1) : (i0-i1-1);
}

/** Get the number of box lengths separating box 1 from the set of boxes
    from minbox0 to maxbox0 */
SIRE_ALWAYS_INLINE int getDelta( const qint16 min0, const qint16 max0, const qint16 i1 )
{
    // the below single-line expression is doing
    // if (i1 >= min0 and i1 <= max0)
    // {
    //     //i1 is contained within the box
    //     return 0;
    // }
    // else
    //     return qMin( getDelta(min0,i1), getDelta(max0,i1) );

    return (i1 >= min0 and i1 <= max0) ? 0 : qMin( getDelta(min0,i1), getDelta(max0,i1) );
}

SIRE_ALWAYS_INLINE int getDelta2(const CLJBoxIndex &box0, const CLJBoxIndex &box1)
{
    const int dx = getDelta(box0.i(), box1.i());
    const int dy = getDelta(box0.j(), box1.j());
    const int dz = getDelta(box0.k(), box1.k());
    return dx*dx + dy*dy + dz*dz;
}

SIRE_ALWAYS_INLINE int getDelta2(const CLJBoxIndex &minbox0, const CLJBoxIndex &maxbox0,
                     const CLJBoxIndex &box1)
{
    const int dx = getDelta(minbox0.i(), maxbox0.i(), box1.i());
    const int dy = getDelta(minbox0.j(), maxbox0.j(), box1.j());
    const int dz = getDelta(minbox0.k(), maxbox0.k(), box1.k());
    return dx*dx + dy*dy + dz*dz;
}

/** Get the number of box lengths separating boxes 0 and 1 at indicies
    i0 and i1 */
SIRE_ALWAYS_INLINE int getBoxDelta( const qint16 i0, const qint16 i1 )
{
    // the below single-line expression is doing
    // if (i0 == i1)
    // {
    //     return 0;
    // }
    // else if (i0 < i1)
    // {
    //     return i1 - i0;
    // }
    // else
    // {
    //     return i0 - i1;
    // }
    //
    // I found on my mac that the single line compiled to faster code
    // than the above if statements.

    return (i0 == i1) ? 0 : (i0 < i1) ? (i1-i0) : (i0-i1);
}

/** Get the number of box lengths separating box 1 from the set of boxes
    from minbox0 to maxbox0 */
SIRE_ALWAYS_INLINE int getBoxDelta( const qint16 min0, const qint16 max0, const qint16 i1 )
{
    // the below single-line expression is doing
    // if (i1 >= min0 and i1 <= max0)
    // {
    //     //i1 is contained within the box
    //     return 0;
    // }
    // else
    //     return qMin( getBoxDelta(min0,i1), getBoxDelta(max0,i1) );

    return (i1 >= min0 and i1 <= max0) ? 0 : qMin( getBoxDelta(min0,i1), getBoxDelta(max0,i1) );
}

/** Return the distance between the two boxes, assuming they are in an infinite cartesian space */
float CLJBoxes::getDistance(const CLJBoxIndex &box0, const CLJBoxIndex &box1) const
{
    return box_length * square_root( getDelta2(box0,box1) );
}

/** Return the distance between the two boxes based on the space 'space' */
float CLJBoxes::getDistance(const Space &space, const CLJBoxIndex &box0,
                            const CLJBoxIndex &box1) const
{
    if (space.isCartesian() and not space.isPeriodic())
    {
        return box_length * square_root( getDelta2(box0,box1) );
    }
    else
    {
        const Length l(box_length);
        return space.minimumDistance(box0.box(l), box1.box(l));
    }
}

/** Return the distance between the two boxes based on the space 'space', assuming that
    one of the boxes covers a volume of nx,ny,nz box lengths above its minimum dimensions */
float CLJBoxes::getDistance(const Space &space, const CLJBoxIndex &box0,
                            const CLJBoxIndex &box1, quint32 nx, quint32 ny, quint32 nz) const
{
    const Length l(box_length);

    AABox min_box = box0.box(l);
    AABox max_box = CLJBoxIndex(box0.i()+nx-1, box0.j()+ny-1, box0.k()+ny-1).box(l);

    return space.minimumDistance( AABox::from(min_box.minCoords(),max_box.maxCoords()),
                                  box1.box(l));
}

/** Return the distances between all of the occupied boxes in 'boxes'
    based on the space 'space' */
QVector<CLJBoxDistance> CLJBoxes::getDistances(const Space &space, const CLJBoxes &boxes)
{
    QVector<CLJBoxDistance> dists;

    const quint32 nboxes = boxes.bxs.count();
    const CLJBoxPtr *b = boxes.bxs.constData();

    dists.reserve((nboxes*nboxes) / 2);

    if (space.isCartesian() and not space.isPeriodic())
    {

        for (quint32 i=0; i<nboxes; ++i)
        {
            const CLJBoxIndex &box0 = b[i].read().index();

            for (quint32 j=i; j<nboxes; ++j)
            {
                const CLJBoxIndex &box1 = b[j].read().index();

                const float dist = square_root( getDelta2(box0,box1) );
                dists.append( CLJBoxDistance(i, j, dist*boxes.box_length) );
            }
        }
    }
    else
    {
        for (quint32 i=0; i<nboxes; ++i)
        {
            const AABox box0 = b[i].read().dimensions();

            for (quint32 j=i; j<nboxes; ++j)
            {
                const AABox box1 = b[j].read().dimensions();

                dists.append( CLJBoxDistance(i, j, space.minimumDistance(box0,box1)) );
            }
        }
    }

    return dists;
}

/** Return the distances between all of the occupied boxes in 'boxes'
    based on the space 'space', only returning boxes that are separated
    by distances of less than 'cutoff' */
QVector<CLJBoxDistance> CLJBoxes::getDistances(const Space &space, const CLJBoxes &boxes,
                                               Length cutoff)
{
    //QElapsedTimer t;
    //t.start();

    QVector<CLJBoxDistance> dists;

    const quint32 nboxes = boxes.bxs.count();
    const CLJBoxPtr *b = boxes.bxs.constData();

    dists.reserve((nboxes*nboxes) / 2);

    if (space.isCartesian())
    {
        if (space.isPeriodic())
        {
            float box_x, box_y, box_z;
            try
            {
                // Attempt to cast as a PeriodicBox.
                Vector dimensions = space.asA<PeriodicBox>().dimensions();

                box_x = dimensions.x() / boxes.box_length;
                box_y = dimensions.y() / boxes.box_length;
                box_z = dimensions.z() / boxes.box_length;
            }
            catch(...)
            {
                // A cubic TriclinicBox is Cartesian.
                Vector v0 = space.asA<TriclinicBox>().vector0();
                Vector v1 = space.asA<TriclinicBox>().vector1();
                Vector v2 = space.asA<TriclinicBox>().vector2();
                box_x = v0.x() / boxes.box_length;
                box_y = v1.y() / boxes.box_length;
                box_z = v2.z() / boxes.box_length;
            }

            const float half_box_x = 0.5 * box_x;
            const float half_box_y = 0.5 * box_y;
            const float half_box_z = 0.5 * box_z;

            const float box_cutoff = cutoff.value() / boxes.box_length;

            //std::ceil rounds up to the neares integer, so we round up the cutoff to
            //the nearest integer to allow us to use integer distance math
            const int int_box_cutoff2 = int( std::ceil( box_cutoff*box_cutoff ) );

            for (quint32 i=0; i<nboxes; ++i)
            {
                const CLJBoxIndex &box0 = b[i].read().index();

                for (quint32 j=i; j<nboxes; ++j)
                {
                    const CLJBoxIndex &box1 = b[j].read().index();

                    int idx = getBoxDelta(box0.i(), box1.i());
                    int idy = getBoxDelta(box0.j(), box1.j());
                    int idz = getBoxDelta(box0.k(), box1.k());

                    if (idx > half_box_x or idy > half_box_y or idz > half_box_z)
                    {
                        float dx = idx;
                        float dy = idy;
                        float dz = idz;

                        while (dx > half_box_x)
                            dx -= box_x;

                        while (dy > half_box_y)
                            dy -= box_y;

                        while (dz > half_box_z)
                            dz -= box_z;

                        dx -= 1;
                        dy -= 1;
                        dz -= 1;

                        dx = qMax(0.0f, dx);
                        dy = qMax(0.0f, dy);
                        dz = qMax(0.0f, dz);

                        const float dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                        if (dist < box_cutoff)
                            dists.append( CLJBoxDistance(i, j, dist*boxes.box_length) );
                    }
                    else
                    {
                        idx = qMax(0, idx-1);
                        idy = qMax(0, idy-1);
                        idz = qMax(0, idz-1);

                        const int d2 = idx*idx + idy*idy + idz*idz;

                        if (d2 <= int_box_cutoff2)
                        {
                            //the box-pair are within cutoff
                            const float dist = square_root(d2);

                            if (dist < box_cutoff)
                                dists.append( CLJBoxDistance(i, j, dist*boxes.box_length) );
                        }
                    }
                }
            }
        }
        else
        {
            const float box_cutoff = cutoff.value() / boxes.box_length;

            //std::ceil rounds up to the neares integer, so we round up the cutoff to
            //the nearest integer to allow us to use integer distance math
            const int int_box_cutoff2 = int( std::ceil( box_cutoff*box_cutoff ) );

            for (quint32 i=0; i<nboxes; ++i)
            {
                const CLJBoxIndex &box0 = b[i].read().index();

                for (quint32 j=i; j<nboxes; ++j)
                {
                    const CLJBoxIndex &box1 = b[j].read().index();

                    const int d2 = getDelta2(box0, box1);

                    if (d2 <= int_box_cutoff2)
                    {
                        //the box-pair are within cutoff
                        const float dist = square_root(d2);

                        if (dist < box_cutoff)
                            dists.append( CLJBoxDistance(i, j, dist*boxes.box_length) );
                    }
                }
            }
        }
    }
    else
    {
        for (quint32 i=0; i<nboxes; ++i)
        {
            const AABox box0 = b[i].read().dimensions();

            for (quint32 j=i; j<nboxes; ++j)
            {
                const AABox box1 = b[j].read().dimensions();

                const float dist = space.minimumDistance(box0,box1);

                if (dist < cutoff.value())
                    dists.append( CLJBoxDistance(i, j, dist) );
            }
        }
    }

    //quint64 ns = t.nsecsElapsed();
    //qDebug() << "Getting box distances took" << (0.000001*ns) << "ms";

    return dists;
}

/** Return the distances between all pairs of occupied boxes between the boxes in
    'boxes0' and the boxes in 'boxes1' */
QVector<CLJBoxDistance> CLJBoxes::getDistances(const Space &space, const CLJBoxes &boxes0,
                                               const CLJBoxes &boxes1)
{
    QVector<CLJBoxDistance> dists;

    const quint32 n0 = boxes0.bxs.count();
    const quint32 n1 = boxes1.bxs.count();
    const CLJBoxPtr *b0 = boxes0.bxs.constData();
    const CLJBoxPtr *b1 = boxes1.bxs.constData();

    dists.reserve((n0*n1)/2);

    if (space.isCartesian() and (boxes0.box_length == boxes1.box_length)
        and not space.isPeriodic())
    {
        for (quint32 i=0; i<n0; ++i)
        {
            const CLJBoxIndex &box0 = b0[i].read().index();

            for (quint32 j=0; j<n1; ++j)
            {
                const CLJBoxIndex &box1 = b1[j].read().index();

                const float dist = square_root( getDelta2(box0,box1) );
                dists.append( CLJBoxDistance(i, j, dist*boxes0.box_length) );
            }
        }
    }
    else
    {
        for (quint32 i=0; i<n0; ++i)
        {
            const AABox box0 = b0[i].read().dimensions();

            for (quint32 j=0; j<n1; ++j)
            {
                const AABox box1 = b1[j].read().dimensions();

                dists.append( CLJBoxDistance(i, j, space.minimumDistance(box0,box1)) );
            }
        }
    }

    return dists;
}

/** Return the distances between all pairs of occupied boxes between the boxes in
    'boxes0' and the boxes in 'boxes1' */
QVector<CLJBoxDistance> CLJBoxes::getDistances(const Space &space, const CLJBoxes &boxes0,
                                               const CLJBoxes &boxes1, Length cutoff)
{
    //QElapsedTimer t;
    //t.start();

    QVector<CLJBoxDistance> dists;

    const quint32 n0 = boxes0.bxs.count();
    const quint32 n1 = boxes1.bxs.count();
    const CLJBoxPtr *b0 = boxes0.bxs.constData();
    const CLJBoxPtr *b1 = boxes1.bxs.constData();

    dists.reserve((n0*n1)/2);

    if (space.isCartesian() and (boxes0.box_length == boxes1.box_length))
    {
        if (space.isPeriodic())
        {
            const float box_cutoff = cutoff.value() / boxes0.box_length;

            //std::ceil rounds up to the neares integer, so we round up the cutoff to
            //the nearest integer to allow us to use integer distance math
            const int int_box_cutoff2 = std::ceil( box_cutoff*box_cutoff );

            Vector dimensions = space.asA<PeriodicBox>().dimensions();

            const float box_x = dimensions.x() / boxes0.box_length;
            const float box_y = dimensions.y() / boxes0.box_length;
            const float box_z = dimensions.z() / boxes0.box_length;

            const float half_box_x = 0.5 * box_x;
            const float half_box_y = 0.5 * box_y;
            const float half_box_z = 0.5 * box_z;

            for (quint32 i=0; i<n0; ++i)
            {
                const CLJBoxIndex &box0 = b0[i].read().index();

                for (quint32 j=0; j<n1; ++j)
                {
                    const CLJBoxIndex &box1 = b1[j].read().index();

                    int idx = getBoxDelta(box0.i(), box1.i());
                    int idy = getBoxDelta(box0.j(), box1.j());
                    int idz = getBoxDelta(box0.k(), box1.k());

                    if (idx > half_box_x or idy > half_box_y or idz > half_box_z)
                    {
                        float dx = idx;
                        float dy = idy;
                        float dz = idz;

                        while (dx > half_box_x)
                            dx -= box_x;

                        while (dy > half_box_y)
                            dy -= box_y;

                        while (dz > half_box_z)
                            dz -= box_z;

                        dx -= 1;
                        dy -= 1;
                        dz -= 1;

                        dx = qMax(0.0f, dx);
                        dy = qMax(0.0f, dy);
                        dz = qMax(0.0f, dz);

                        const float dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                        if (dist < box_cutoff)
                            dists.append( CLJBoxDistance(i, j, dist*boxes0.box_length) );
                    }
                    else
                    {
                        idx = qMax( 0, idx-1 );
                        idy = qMax( 0, idy-1 );
                        idz = qMax( 0, idz-1 );

                        const int d2 = idx*idx + idy*idy + idz*idz;

                        if (d2 <= int_box_cutoff2)
                        {
                            //the box-pair are within cutoff
                            const float dist = square_root(d2);

                            if (dist < box_cutoff)
                                dists.append( CLJBoxDistance(i, j, dist*boxes0.box_length) );
                        }
                    }
                }
            }
        }
        else
        {
            const float box_cutoff = cutoff.value() / boxes0.box_length;

            //std::ceil rounds up to the neares integer, so we round up the cutoff to
            //the nearest integer to allow us to use integer distance math
            const int int_box_cutoff2 = std::ceil( box_cutoff*box_cutoff );

            for (quint32 i=0; i<n0; ++i)
            {
                const CLJBoxIndex &box0 = b0[i].read().index();

                for (quint32 j=0; j<n1; ++j)
                {
                    const CLJBoxIndex &box1 = b1[j].read().index();

                    int d2 = getDelta2(box0,box1);

                    if (d2 <= int_box_cutoff2)
                    {
                        //the box-pair are within cutoff
                        const float dist = square_root(d2);

                        if (dist < box_cutoff)
                            dists.append( CLJBoxDistance(i, j, dist*boxes0.box_length) );
                    }
                }
            }
        }
    }
    else
    {
        for (quint32 i=0; i<n0; ++i)
        {
            const AABox box0 = b0[i].read().dimensions();

            for (quint32 j=0; j<n1; ++j)
            {
                const AABox box1 = b1[j].read().dimensions();

                const float dist = space.minimumDistance(box0, box1);

                if (dist < cutoff.value())
                    dists.append( CLJBoxDistance(i, j, dist) );
            }
        }
    }

    //quint64 ns = t.nsecsElapsed();
    //qDebug() << "Getting box distances took" << (0.000001*ns) << "ms";

    return dists;
}

/** Return the distances between all the passed atoms and all of the occupied
    boxes in 'boxes1', returning a set of CLJBoxDistance objects where box0() is 0 */
QVector<CLJBoxDistance> CLJBoxes::getDistances(const Space &space, const CLJAtoms &atoms0,
                                               const CLJBoxes &boxes1, Length cutoff)
{
    //QElapsedTimer t;
    //t.start();

    //get the minimum and maximum coordinates of the atoms
    Vector mincoords0 = atoms0.minCoords();
    Vector maxcoords0 = atoms0.maxCoords();

    QVector<CLJBoxDistance> dists;

    const quint32 n1 = boxes1.bxs.count();
    const CLJBoxPtr *b1 = boxes1.bxs.constData();

    dists.reserve(n1);

    //having timed this, it is just as quick to just convert to AABox and use
    //that to calculate the minimum distance
    const AABox box0 = AABox::from(mincoords0, maxcoords0);

    for (quint32 j=0; j<n1; ++j)
    {
        const AABox box1 = b1[j].read().dimensions();

        const float dist = space.minimumDistance(box0, box1);

        if (dist < cutoff.value())
            dists.append( CLJBoxDistance(0, j, dist) );
    }

    //qint64 ns = t.nsecsElapsed();
    //qDebug() << "Getting box distances took" << (0.000001*ns) << "ms";

    return dists;
}

/** Return the distances between all the passed atoms and all of the occupied
    boxes in 'boxes1', returning a set of CLJBoxDistance objects where box0() is 0 */
QVector<CLJBoxDistance> CLJBoxes::getDistances(const Space &space, const CLJAtoms &atoms0,
                                               const CLJBoxes &boxes1)
{
    //QElapsedTimer t;
    //t.start();

    //get the minimum and maximum coordinates of the atoms
    Vector mincoords0 = atoms0.minCoords();
    Vector maxcoords0 = atoms0.maxCoords();

    QVector<CLJBoxDistance> dists;

    const quint32 n1 = boxes1.bxs.count();
    const CLJBoxPtr *b1 = boxes1.bxs.constData();

    dists.reserve(n1);

    //having timed this, it is just as quick to just convert to AABox and use
    //that to calculate the minimum distance
    const AABox box0 = AABox::from(mincoords0, maxcoords0);

    for (quint32 j=0; j<n1; ++j)
    {
        const AABox box1 = b1[j].read().dimensions();

        const float dist = space.minimumDistance(box0, box1);
        dists.append( CLJBoxDistance(0, j, dist) );
    }

    //qint64 ns = t.nsecsElapsed();
    //qDebug() << "Getting box distances took" << (0.000001*ns) << "ms";

    return dists;
}
