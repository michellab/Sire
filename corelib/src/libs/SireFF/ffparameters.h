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

#ifndef SIREFF_FFPARAMETERS_H
#define SIREFF_FFPARAMETERS_H

#include "SireBase/property.h"

#include <QVarLengthArray>

SIRE_BEGIN_HEADER

namespace SireFF
{
class FFParameters;
class FFParametersArray;

class NullFFParameters;
class NullFFParametersArray;

}

QDataStream& operator<<(QDataStream&, const SireFF::FFParameters&);
QDataStream& operator>>(QDataStream&, SireFF::FFParameters&);

QDataStream& operator<<(QDataStream&, const SireFF::FFParametersArray&);
QDataStream& operator>>(QDataStream&, SireFF::FFParametersArray&);

QDataStream& operator<<(QDataStream&, const SireFF::NullFFParameters&);
QDataStream& operator>>(QDataStream&, SireFF::NullFFParameters&);

QDataStream& operator<<(QDataStream&, const SireFF::NullFFParametersArray&);
QDataStream& operator>>(QDataStream&, SireFF::NullFFParametersArray&);

namespace SireFF
{

typedef SireBase::PropPtr<FFParameters> FFParametersPtr;
typedef SireBase::PropPtr<FFParametersArray> FFParametersArrayPtr;

/** This is the virtual base class of all FFParameters. These classes
    are used to hold the extracted parameters for beads in a forcefield
    in a format that allows for rapid processing by the forcefield.

    The FFParameters holds all of the forcefield parameters for a single
    bead in the forcefield. All of the parameters for all of the beads
    are held in the corresponding FFParametersArray-derived class.
    
    These classes are designed to be used internally by the forcefields
    and are not designed to be part of the public API
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFParameters : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const FFParameters&);
friend QDataStream& ::operator>>(QDataStream&, FFParameters&);

public:
    FFParameters();
    FFParameters(const FFParameters &other);
    
    virtual ~FFParameters();
    
    static const char* typeName();
    
    virtual FFParameters* clone() const=0;

    virtual FFParametersArrayPtr toArray() const=0;

    static NullFFParameters null();

protected:
    FFParameters& operator=(const FFParameters &other);
    
    bool operator==(const FFParameters &other) const;
    bool operator!=(const FFParameters &other) const;
};

/** This is the virtual base class of all FFParametersArray objects. These
    classes are used to hold all of the parameters for all of the beads
    in the forcefield, in a format that allows for rapid computation by
    the forcefield.
    
    The parameters for an individual bead are held in the corresponding
    FFParameters-derived class.
    
    These classes are designed to be used internally by the forcefields
    and are not designed to be part of the public API
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFParametersArray : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const FFParametersArray&);
friend QDataStream& ::operator>>(QDataStream&, FFParametersArray&);

public:
    FFParametersArray();
    FFParametersArray(const FFParametersArray &other);
    
    virtual ~FFParametersArray();
    
    static const char* typeName();
    
    virtual FFParametersArray* clone() const=0;
    
    static NullFFParametersArray null();

    virtual FFParametersPtr operator[](int i) const=0;
    
    FFParametersPtr at(int i) const;

    virtual int count() const=0;
    virtual bool isEmpty() const=0;

    virtual void append(const FFParameters &params)=0;
    virtual void append(const FFParametersArray &params)=0;
    
    virtual void update(int idx, const FFParameters &params)=0;
    virtual void update(const QVarLengthArray<int> &idxs, 
                        const FFParametersArray &params)=0;

    virtual void remove(int idx)=0;
    virtual void remove(const QVarLengthArray<int> &idxs)=0;
    
    virtual void removeAll()=0;

protected:
    FFParametersArray& operator=(const FFParametersArray &other);
    
    bool operator==(const FFParametersArray &other) const;
    bool operator!=(const FFParametersArray &other) const;
};

/** Null FFParameters */
class SIREFF_EXPORT NullFFParameters 
            : public SireBase::ConcreteProperty<NullFFParameters,FFParameters>
{

friend QDataStream& ::operator<<(QDataStream&, const NullFFParameters&);
friend QDataStream& ::operator>>(QDataStream&, NullFFParameters&);

public:
    NullFFParameters();
    NullFFParameters(const NullFFParameters &other);
    
    ~NullFFParameters();
    
    static const char* typeName();
    
    NullFFParameters& operator=(const NullFFParameters &other);
    
    bool operator==(const NullFFParameters &other) const;
    bool operator!=(const NullFFParameters &other) const;
    
    FFParametersArrayPtr toArray() const;
};

/** Null FFParameters */
class SIREFF_EXPORT NullFFParametersArray
            : public SireBase::ConcreteProperty<NullFFParametersArray,FFParametersArray>
{

friend QDataStream& ::operator<<(QDataStream&, const NullFFParametersArray&);
friend QDataStream& ::operator>>(QDataStream&, NullFFParametersArray&);

public:
    NullFFParametersArray();
    NullFFParametersArray(const NullFFParametersArray &other);
    
    ~NullFFParametersArray();
    
    static const char* typeName();
    
    NullFFParametersArray& operator=(const NullFFParametersArray &other);
    
    bool operator==(const NullFFParametersArray &other) const;
    bool operator!=(const NullFFParametersArray &other) const;
    
    FFParametersPtr operator[](int i) const;

    int count() const;
    bool isEmpty() const;

    void append(const FFParameters &params);
    void append(const FFParametersArray &params);
    
    void update(int idx, const FFParameters &params);
    void update(const QVarLengthArray<int> &idxs, const FFParametersArray &params);

    void remove(int idx);
    void remove(const QVarLengthArray<int> &idxs);
    
    void removeAll();
};

} // end of namespace SireFF

Q_DECLARE_METATYPE( SireFF::NullFFParameters )
Q_DECLARE_METATYPE( SireFF::NullFFParametersArray )

SIRE_END_HEADER

#endif
