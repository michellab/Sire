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

#ifndef SIREFF_ATOMICFFPARAMETERS_HPP
#define SIREFF_ATOMICFFPARAMETERS_HPP

#include "ffparameters.h"

#include "SireBase/packedarray2d.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"
#include "tostring.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
template<class T>
class AtomicFFParameters;

template<class T>
class AtomicFFParametersArray;
}

template<class T>
QDataStream& operator<<(QDataStream&, const SireFF::AtomicFFParameters<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireFF::AtomicFFParameters<T>&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireFF::AtomicFFParametersArray<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireFF::AtomicFFParametersArray<T>&);

namespace SireFF
{

/** This template class holds the FFParameters of type T, where
    there is a parameter for each atom in the bead
    
    @author Christopher Woods
*/
template<class T>
class AtomicFFParameters
        : public SireBase::ConcreteProperty<AtomicFFParameters<T>,FFParameters>
{

friend SIREFF_EXPORT QDataStream& ::operator<<<>(QDataStream&, const AtomicFFParameters<T>&);
friend SIREFF_EXPORT QDataStream& ::operator>><>(QDataStream&, AtomicFFParameters<T>&);

public:
    typedef typename SireBase::PackedArray2D<T>::Array Array;

    AtomicFFParameters();
    AtomicFFParameters(const QVector<T> &params);
    AtomicFFParameters(const typename SireBase::PackedArray2D<T>::Array &params);
    
    AtomicFFParameters(const AtomicFFParameters<T> &other);
    
    ~AtomicFFParameters();
    
    static const char* typeName();
    
    AtomicFFParameters<T>& operator=(const AtomicFFParameters<T> &other);
    
    bool operator==(const AtomicFFParameters<T> &other) const;
    bool operator!=(const AtomicFFParameters<T> &other) const;
    
    FFParametersArrayPtr toArray() const;
    
    const T* constData() const;
    
    int count() const;
    
    const typename SireBase::PackedArray2D<T>::Array& array() const;
    
    operator typename SireBase::PackedArray2D<T>::Array() const;
    
private:
    /** The actual parameters */
    typename SireBase::PackedArray2D<T>::Array params;
};

/** This class holds the array of AtomicFFParameters<T> 

    @author Christopher Woods
*/
template<class T>
class AtomicFFParametersArray
    : public SireBase::ConcreteProperty<AtomicFFParametersArray<T>,FFParametersArray>
{

friend SIREFF_EXPORT QDataStream& ::operator<<<>(QDataStream&, const AtomicFFParametersArray<T>&);
friend SIREFF_EXPORT QDataStream& ::operator>><>(QDataStream&, AtomicFFParametersArray<T>&);

public:
    typedef typename SireBase::PackedArray2D<T>::Array Array;

    AtomicFFParametersArray();
    AtomicFFParametersArray(const SireBase::PackedArray2D<T> &params);
    AtomicFFParametersArray(const QVector<T> &params);
    AtomicFFParametersArray(const QVector< QVector<T> > &params);
    AtomicFFParametersArray(const AtomicFFParameters<T> &params);
    
    AtomicFFParametersArray(const AtomicFFParametersArray<T> &other);
    
    ~AtomicFFParametersArray();
    
    static const char* typeName();
    
    AtomicFFParametersArray<T>& operator=(const AtomicFFParametersArray<T> &other);
    
    bool operator==(const AtomicFFParametersArray<T> &other) const;
    bool operator!=(const AtomicFFParametersArray<T> &other) const;
    
    FFParametersPtr operator[](int i) const;
    
    int count() const;
    
    bool isEmpty() const;
    
    const SireBase::PackedArray2D<T>& array() const;
    
    operator SireBase::PackedArray2D<T>() const;
    
    const Array* constData() const;
    
    void append(const QVector<T> &params);
    void append(const QVector< QVector<T> > &params);
    
    void append(const AtomicFFParameters<T> &params);
    void append(const AtomicFFParametersArray<T> &params);
    
    void append(const FFParameters &params);
    void append(const FFParametersArray &params);
    
    void update(int idx, const QVector<T> &params);
    void update(const QVarLengthArray<int> &idxs, const QVector< QVector<T> > &params);
    
    void update(int idx, const AtomicFFParameters<T> &params);
    void update(const QVarLengthArray<int> &idxs, 
                const AtomicFFParametersArray<T> &params);
    
    void update(int idx, const FFParameters &params);
    void update(const QVarLengthArray<int> &idxs, const FFParametersArray &params);
    
    void remove(int idx);
    void remove(const QVarLengthArray<int> &idxs);
    
    void removeAll();
    
private:
    /** The packed array of all of the bead's parameters */
    SireBase::PackedArray2D<T> params;
};

/////////
///////// Implementation of AtomicFFParameters
/////////

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParameters<T>::AtomicFFParameters()
                      : SireBase::ConcreteProperty<AtomicFFParameters<T>,FFParameters>()
{}

/** Construct to hold the parameters in 'params' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParameters<T>::AtomicFFParameters(const QVector<T> &parameters)
                      : SireBase::ConcreteProperty<AtomicFFParameters<T>,FFParameters>(),
                        params(parameters)
{}

/** Construct to hold the parameters in 'params' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParameters<T>::AtomicFFParameters(
                        const typename SireBase::PackedArray2D<T>::Array &parameters)
                      : SireBase::ConcreteProperty<AtomicFFParameters<T>,FFParameters>(),
                        params(parameters)
{}
                        
/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParameters<T>::AtomicFFParameters(const AtomicFFParameters<T> &other)
                 : SireBase::ConcreteProperty<AtomicFFParameters<T>,FFParameters>(other),
                   params(other.params)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParameters<T>::~AtomicFFParameters()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* AtomicFFParameters<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< AtomicFFParameters<T> >() );
}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParameters<T>& AtomicFFParameters<T>::operator=(
                                                const AtomicFFParameters<T> &other)
{
    params = other.params;
    FFParameters::operator=(other);
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomicFFParameters<T>::operator==(const AtomicFFParameters<T> &other) const
{
    return params == other.params and FFParameters::operator==(other);
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomicFFParameters<T>::operator!=(const AtomicFFParameters<T> &other) const
{
    return not AtomicFFParameters<T>::operator==(other);
}

/** Return these parameters in an FFParameterArray */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
FFParametersArrayPtr AtomicFFParameters<T>::toArray() const
{
    return AtomicFFParametersArray<T>(*this);
}

/** Return the number of parameters in the array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomicFFParameters<T>::count() const
{
    return params.count();
}

/** Return a raw pointer to the parameters */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* AtomicFFParameters<T>::constData() const
{
    return params.constData();
}

/** Return the raw array of parameters */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename SireBase::PackedArray2D<T>::Array& AtomicFFParameters<T>::array() const
{
    return params;
}

/** Allow automatic casting to a SireBase::PackedArray2D<T>::Array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParameters<T>::operator typename SireBase::PackedArray2D<T>::Array() const
{
    return params;
}

/////////
///////// Implementation of AtomicFFParametersArray
/////////

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParametersArray<T>::AtomicFFParametersArray()
    : SireBase::ConcreteProperty<AtomicFFParametersArray<T>,FFParametersArray>()
{}

/** Construct to hold the passed parameters */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParametersArray<T>::AtomicFFParametersArray(
                                    const SireBase::PackedArray2D<T> &parameters)
    : SireBase::ConcreteProperty<AtomicFFParametersArray<T>,FFParametersArray>(),
      params(parameters)
{}

/** Construct to hold the passed parameters */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParametersArray<T>::AtomicFFParametersArray(const QVector<T> &parameters)
    : SireBase::ConcreteProperty<AtomicFFParametersArray<T>,FFParametersArray>(),
      params(parameters)
{}

/** Construct to hold the passed parameters */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParametersArray<T>::AtomicFFParametersArray(
                                        const QVector< QVector<T> > &parameters)
    : SireBase::ConcreteProperty<AtomicFFParametersArray<T>,FFParametersArray>(),
      params(parameters)
{}

/** Construct from the passed AtomicFFParameters */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParametersArray<T>::AtomicFFParametersArray(const AtomicFFParameters<T> &prms)
    : SireBase::ConcreteProperty<AtomicFFParametersArray<T>,FFParametersArray>(),
      params(prms)
{}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParametersArray<T>::AtomicFFParametersArray(
                                    const AtomicFFParametersArray<T> &other)
    : SireBase::ConcreteProperty<AtomicFFParametersArray<T>,FFParametersArray>(other),
      params(other.params)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParametersArray<T>::~AtomicFFParametersArray()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* AtomicFFParametersArray<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< AtomicFFParametersArray<T> >() );
}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParametersArray<T>& AtomicFFParametersArray<T>::operator=(
                                            const AtomicFFParametersArray<T> &other)
{
    params = other.params;
    FFParametersArray::operator=(other);
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomicFFParametersArray<T>::operator==(const AtomicFFParametersArray<T> &other) const
{
    return params == other.params and FFParametersArray::operator==(other);
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomicFFParametersArray<T>::operator!=(const AtomicFFParametersArray<T> &other) const
{
    return not FFParametersArray::operator==(other);
}

/** Return the ith group of atomic parameters

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
FFParametersPtr AtomicFFParametersArray<T>::operator[](int i) const
{
    return FFParametersPtr( AtomicFFParameters<T>(params[i]) );
}

/** Return the raw parameter array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const SireBase::PackedArray2D<T>& AtomicFFParametersArray<T>::array() const
{
    return params;
}

/** Return the raw pointer to the array of parameters for a bead */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename SireBase::PackedArray2D<T>::Array*
AtomicFFParametersArray<T>::constData() const
{
    return params.constData();
}

/** Return the number of arrays of parameters */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomicFFParametersArray<T>::count() const
{
    return params.count();
}

/** Return whether or not this array is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomicFFParametersArray<T>::isEmpty() const
{
    return params.isEmpty();
}

/** Allow for automatic casting to a raw parameter array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomicFFParametersArray<T>::operator SireBase::PackedArray2D<T>() const
{
    return params;
}

/** Append the passed parameters onto the end of the array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::append(const QVector<T> &parameters)
{
    params.append(parameters);
}

/** Append the passed arrays of parameters onto the end of the array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::append(const QVector< QVector<T> > &parameters)
{
    params.append(parameters);
}

/** Append the passed parameters onto the end of the array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::append(const AtomicFFParameters<T> &parameters)
{
    params.append( parameters.array() );
}

/** Append the passed array of parameters onto the end of the array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::append(const AtomicFFParametersArray<T> &parameters)
{
    params.append( parameters.array() );
}

/** Append the passed parameters onto the end of the array

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::append(const FFParameters &params)
{
    AtomicFFParametersArray<T>::append( params.asA< AtomicFFParameters<T> >() );
}

/** Append the passed array of parameters onto the end of this array

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::append(const FFParametersArray &params)
{
    AtomicFFParametersArray<T>::append( params.asA< AtomicFFParametersArray<T> >() );
}

/** Update the parameters at index 'idx' with the passed values

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::update(int idx, const QVector<T> &parameters)
{
    params.update(idx, parameters);
}

/** Update the parameters at the passed indicies with the passed values

    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::update(const QVarLengthArray<int> &idxs, 
                                        const QVector< QVector<T> > &parameters)
{
    params.updateAll(idxs, AtomicFFParametersArray<T>(parameters));
}

/** Update the parameters at index 'idx' with the passed values

    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::update(int idx, const AtomicFFParameters<T> &parameters)
{
    params.update(idx, parameters.array());
}

/** Update the parameters at the passed indicies with the passed values

    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::update(const QVarLengthArray<int> &idxs, 
                                        const AtomicFFParametersArray<T> &parameters)
{
    params.updateAll(idxs, parameters.array());
}

/** Update the parameters at index 'idx' with the passed values

    \throw SireError::invalid_cast
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::update(int idx, const FFParameters &parameters)
{
    AtomicFFParametersArray<T>::update(idx, parameters.asA< AtomicFFParameters<T> >());
}

/** Update the parameters at the passed indicies with the passed values

    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::update(const QVarLengthArray<int> &idxs, 
                                        const FFParametersArray &parameters)
{
    AtomicFFParametersArray<T>::update(idxs, 
                                       parameters.asA< AtomicFFParametersArray<T> >());
}

/** Remove all of the parameters at index 'idx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::remove(int idx)
{
    params.remove(idx);
}

/** Remove all of the parameters at the passed indicies

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::remove(const QVarLengthArray<int> &idxs)
{
    params.removeAll(idxs);
}

/** Remove all of the parameters */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomicFFParametersArray<T>::removeAll()
{
    params = SireBase::PackedArray2D<T>();
}

/////////
///////// Implementation of streaming operators
/////////

namespace detail
{

template<class T>
struct AtomicFFParametersMT
{
    static const RegisterMetaType< SireFF::AtomicFFParameters<T> > r_atomffparams;
};

template<class T>
struct AtomicFFParametersArrayMT
{
    static const RegisterMetaType< SireFF::AtomicFFParametersArray<T> > r_atomffparams;
};

template<class T>
const RegisterMetaType< SireFF::AtomicFFParameters<T> > 
                                AtomicFFParametersMT<T>::r_atomffparams;

template<class T>
const RegisterMetaType< SireFF::AtomicFFParametersArray<T> >
                                AtomicFFParametersArrayMT<T>::r_atomffparams;

} // end of namespace detail

} // end of namespace SireFF

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireFF::AtomicFFParameters<T> &ffparams)
{
    SireStream::writeHeader(ds, 
                            SireFF::detail::AtomicFFParametersMT<T>::r_atomffparams,
                            1);
                            
    SireStream::SharedDataStream sds(ds);
    
    sds << ffparams.params;
    
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireFF::AtomicFFParameters<T> &ffparams)
{
    SireStream::VersionID v = SireStream::readHeader(ds,
                                SireFF::detail::AtomicFFParametersMT<T>::r_atomffparams);
                                
    if (v == 1)
    {
        SireStream::SharedDataStream sds(ds);
        sds >> ffparams.params;
    }
    else
        throw SireStream::version_error(v, "1", 
                               SireFF::detail::AtomicFFParametersMT<T>::r_atomffparams,
                               CODELOC);
                               
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, 
                        const SireFF::AtomicFFParametersArray<T> &ffparams)
{
    SireStream::writeHeader(ds, 
                            SireFF::detail::AtomicFFParametersArrayMT<T>::r_atomffparams,
                            1);
                            
    SireStream::SharedDataStream sds(ds);
    
    sds << ffparams.params;
    
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireFF::AtomicFFParametersArray<T> &ffparams)
{
    SireStream::VersionID v = SireStream::readHeader(ds,
                    SireFF::detail::AtomicFFParametersArrayMT<T>::r_atomffparams);
                                
    if (v == 1)
    {
        SireStream::SharedDataStream sds(ds);
        sds >> ffparams.params;
    }
    else
        throw SireStream::version_error(v, "1", 
                            SireFF::detail::AtomicFFParametersArrayMT<T>::r_atomffparams,
                            CODELOC);
                               
    return ds;
}

SIRE_END_HEADER

#endif
