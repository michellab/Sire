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

#ifndef SIREMM_CLJATOMS_H
#define SIREMM_CLJATOMS_H

#include <QVector>

#include "ljparameter.h"

#include "SireUnits/dimensions.h"

#include "SireBase/propertymap.h"

#include "SireMaths/vector.h"
#include "SireMaths/multifloat.h"
#include "SireMaths/multiint.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJAtom;
class CLJAtoms;
}

QDataStream& operator<<(QDataStream&, const SireMM::CLJAtom&);
QDataStream& operator>>(QDataStream&, SireMM::CLJAtom&);

QDataStream& operator<<(QDataStream&, const SireMM::CLJAtoms&);
QDataStream& operator>>(QDataStream&, SireMM::CLJAtoms&);

namespace SireMol
{
class Molecule;
class Molecules;
class MoleculeGroup;
class MoleculeView;
}

namespace SireMM
{

using SireBase::PropertyMap;

using SireUnits::Dimension::Charge;

using SireMaths::Vector;
using SireMaths::MultiFloat;
using SireMaths::MultiInt;

using SireMol::MoleculeGroup;
using SireMol::Molecules;
using SireMol::MoleculeView;

/** This class holds everything about a single CLJAtom */
class SIREMM_EXPORT CLJAtom
{

friend class CLJAtoms;

friend QDataStream& ::operator<<(QDataStream&, const CLJAtom&);
friend QDataStream& ::operator>>(QDataStream&, CLJAtom&);

public:
    CLJAtom();
    CLJAtom(Vector coords, Charge charge, LJParameter ljparam, qint32 idnum=1);
    
    CLJAtom(const CLJAtom &other);
    
    ~CLJAtom();
    
    CLJAtom& operator=(const CLJAtom &other);
    
    bool operator==(const CLJAtom &other) const;
    bool operator!=(const CLJAtom &other) const;

    static const char* typeName();
    
    const char* what() const;
    
    static QVector<CLJAtom> buildFrom(const MoleculeView &molecule,
                                      const PropertyMap &map = PropertyMap());

    QString toString() const;
    
    Vector coordinates() const;
    Charge charge() const;
    LJParameter ljParameter() const;
    qint32 ID() const;
    
    bool isDummy() const;
    bool isNull() const;
    
    CLJAtom negate() const;
    
private:
    /** The coordinates of the atom */
    float x;
    float y;
    float z;
    
    /** The reduced charge of the atom */
    float chg;
    
    /** The reduced sigma parameter for the atom */
    float sig;
    
    /** The reduced epsilon parameter for the atom */
    float eps;
    
    /** The ID number for the atom - atoms with ID 0 are dummies,
        while atoms with the same ID are in the same molecule */
    qint32 idnum;
};

/** This class holds vectorised arrays of the coordinates,
    reduced charges and reduced Lennard Jones parameters of 
    a set of atoms. This class is intended only to be used
    in the fast functions used to calculated coulomb and LJ
    energies, and is not intended for general use.
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJAtoms
{

friend QDataStream& ::operator<<(QDataStream&, const CLJAtoms&);
friend QDataStream& ::operator>>(QDataStream&, CLJAtoms&);

public:
    enum ID_SOURCE
    {
        USE_MOLNUM = 0,
        USE_ATOMIDX = 1
    };

    CLJAtoms();
    CLJAtoms(const QVector<CLJAtom> &atoms);
    CLJAtoms(const QList<CLJAtom> &atoms);
    CLJAtoms(const CLJAtom *atoms, int natoms);
    
    CLJAtoms(const QVector<Vector> &coordinates,
             const QVector<Charge> &charges,
             const QVector<LJParameter> &ljparams,
             qint32 atomid=1);

    CLJAtoms(const QVector<Vector> &coordinates,
             const QVector<Charge> &charges,
             const QVector<LJParameter> &ljparams,
             const QVector<qint32> &ids);

    CLJAtoms(const MoleculeView &molecule,
             const PropertyMap &map = PropertyMap());

    CLJAtoms(const MoleculeView &molecule,
             ID_SOURCE id_source,
             const PropertyMap &map = PropertyMap());

    CLJAtoms(const Molecules &molecules,
             const PropertyMap &map = PropertyMap());

    CLJAtoms(const MoleculeGroup &molecules,
             const PropertyMap &map = PropertyMap());

    CLJAtoms(const Molecules &molecules,
             ID_SOURCE id_source,
             const PropertyMap &map = PropertyMap());

    CLJAtoms(const MoleculeGroup &molecules,
             ID_SOURCE id_source,
             const PropertyMap &map = PropertyMap());

    CLJAtoms(const CLJAtoms &other);

    ~CLJAtoms();
    
    void reconstruct(const MoleculeView &molecule,
                     const PropertyMap &map = PropertyMap());
    
    void reconstruct(const MoleculeView &molecule,
                     ID_SOURCE id_source,
                     const PropertyMap &map = PropertyMap());
    
    CLJAtoms& operator=(const CLJAtoms &other);
    
    bool operator==(const CLJAtoms &other) const;
    bool operator!=(const CLJAtoms &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    QString toString() const;

    int count() const;
    int size() const;
    
    int nAtoms() const;
    
    bool isEmpty() const;
    
    CLJAtom operator[](int i) const;
    
    CLJAtom at(int i) const;
    CLJAtom getitem(int i) const;
    
    void resize(int new_size);
    void copyIn(const CLJAtoms &other);
    
    CLJAtoms operator+(const CLJAtoms &other) const;
    
    CLJAtoms operator+(const CLJAtom &atom) const;
    CLJAtoms operator+(const QVector<CLJAtom> &atoms) const;
    
    CLJAtoms& operator+=(const CLJAtoms &other);
    
    CLJAtoms& operator+=(const CLJAtom &atom);
    CLJAtoms& operator+=(const QVector<CLJAtom> &atoms);
    
    void append(const CLJAtom &atom);
    void append(const CLJAtoms &atoms, int n=-1);
    
    void set(int i, const CLJAtom &atom);
    
    void setCoordinates(int i, Vector coords);
    void setCharge(int i, Charge charge);
    void setLJParameter(int i, LJParameter ljparam);
    void setID(int i, qint32 idnum);

    void setAllID(qint32 idnum);

    void makeDummy(int i);
    bool isDummy(int i) const;
    
    Vector minCoords() const;
    Vector maxCoords() const;
    
    CLJAtoms negate() const;
    
    QVector<CLJAtom> atoms() const;
    
    QVector<Vector> coordinates() const;
    QVector<Charge> charges() const;
    QVector<LJParameter> ljParameters() const;
    QVector<qint32> IDs() const;
    
    const QVector<MultiFloat>& x() const;
    const QVector<MultiFloat>& y() const;
    const QVector<MultiFloat>& z() const;
    
    const QVector<MultiFloat>& q() const;
    const QVector<MultiFloat>& sigma() const;
    const QVector<MultiFloat>& epsilon() const;
    
    const QVector<MultiInt>& ID() const;
    
    static MultiInt idOfDummy();
    
    CLJAtoms squeeze() const;
    
    bool isPadded() const;
    int nPadded() const;

    bool hasDummies() const;
    int nDummies() const;

private:
    void constructFrom(const MoleculeGroup &molecules,
                       const ID_SOURCE id_source, const PropertyMap &map);
    void constructFrom(const Molecules &molecules,
                       const ID_SOURCE id_source, const PropertyMap &map);
    void constructFrom(const MoleculeView &molecule,
                       const ID_SOURCE id_source, const PropertyMap &map);

    /** Vector of the x-coordinates of the atoms */
    QVector<MultiFloat> _x;
    
    /** Vector of the y-coordinates of the atoms */
    QVector<MultiFloat> _y;
    
    /** Vector of the z-coordinates of the atoms */
    QVector<MultiFloat> _z;
    
    /** Vector of the reduced partial charges of the 
        atoms (charge divided by sqrt(4 pi epsilon0) ) */
    QVector<MultiFloat> _q;
    
    /** Vector of the square root of the sigma LJ parameters */
    QVector<MultiFloat> _sig;
    
    /** Vector of the square root of the epsilon LJ parameters */
    QVector<MultiFloat> _eps;

    /** The molecule number for each atom - atoms with the same
        number are part of the same molecule. Also, if this number is
        zero, then this is a dummy atom */
    QVector<MultiInt> _id;
    
    /** ID number given to all dummy atoms */
    static qint32 id_of_dummy;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the vector of the vectorised x coordinates */
inline const QVector<MultiFloat>& CLJAtoms::x() const
{
    return _x;
}

/** Return the vector of the vectorised y coordinates */
inline const QVector<MultiFloat>& CLJAtoms::y() const
{
    return _y;
}

/** Return the vector of the vectorised z coordinates */
inline const QVector<MultiFloat>& CLJAtoms::z() const
{
    return _z;
}

/** Return the vector of the vectorised reduced charges
    (square root of charge divided by 4 pi epsilon 0) */
inline const QVector<MultiFloat>& CLJAtoms::q() const
{
    return _q;
}

/** Return the vector of the vectorised reduced (square root) sigma parameters */
inline const QVector<MultiFloat>& CLJAtoms::sigma() const
{
    return _sig;
}

/** Return the vector of the vectorised reduced (square root) epsilon parameters */
inline const QVector<MultiFloat>& CLJAtoms::epsilon() const
{
    return _eps;
}

/** Return the vector of vectorised atom IDs */
inline const QVector<MultiInt>& CLJAtoms::ID() const
{
    return _id;
}

/** Return whether or not the ith atom is a dummy atom (has an ID of 0) */
inline bool CLJAtoms::isDummy(int i) const
{
    int idx = i / MultiFloat::count();
    int sub_idx = i % MultiFloat::count();

    return (_id[idx][sub_idx] == 0);
}

/** Return whether or not this set of atoms is padded
    (has dummy atoms at the end of the array to pad out
     the MultiFloat/MultiInt vectors) */
inline bool CLJAtoms::isPadded() const
{
    if (_id.isEmpty())
        return false;
    else
        return _id.constData()[ _id.count()-1 ][ MultiInt::count() -1 ] == id_of_dummy;
}

/** Return the number of padded elements on the end of this set of atoms */
inline int CLJAtoms::nPadded() const
{
    if (_id.isEmpty())
        return 0;
    else
    {
        const MultiInt &last = _id.constData()[ _id.count() - 1 ];
        
        int n = 0;
        
        for (int i=MultiInt::count()-1; i>=0; --i)
        {
            if (last[i] != id_of_dummy)
                return n;
            
            n += 1;
        }
        
        return n;
    }
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMM::CLJAtom);
Q_DECLARE_METATYPE(SireMM::CLJAtoms);

Q_DECLARE_TYPEINFO(SireMM::CLJAtoms, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMM::CLJAtom )
SIRE_EXPOSE_CLASS( SireMM::CLJAtoms )

SIRE_END_HEADER

#endif

