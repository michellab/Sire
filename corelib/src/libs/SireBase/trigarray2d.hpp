/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREBASE_TRIGARRAY2D_HPP
#define SIREBASE_TRIGARRAY2D_HPP

#include <QStringList>

#include "trigarray2d.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
template<class T>
class TrigArray2D;
}

template<class T>
QDataStream& operator<<(QDataStream&, const SireBase::TrigArray2D<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireBase::TrigArray2D<T>&);

namespace SireBase
{

/** This provides a 2D symmetric matrix of objects of type T. All objects
    are packed together in memory, in packed format
    
        j0   j1   j2
    i0   0    1    2
    i1   1    3    4
    i2   2    4    5
    
    Arranged in memory as;
    
    [ (i0,j0), (i0,j1), (i0,j2), (i1,j1), (i1,j2), (i2,j2) ]
    [    0   ,    1   ,    2   ,   3    ,    4,       5    ]
    
    @author Christopher Woods
*/
template<class T>
class TrigArray2D : public TrigArray2DBase
{

friend SIREBASE_EXPORT QDataStream& ::operator<<<>(QDataStream&, const TrigArray2D<T>&);
friend SIREBASE_EXPORT QDataStream& ::operator>><>(QDataStream&, TrigArray2D<T>&);

public:
    TrigArray2D();

    TrigArray2D(int dimension);
    TrigArray2D(int dimension, const T &default_value);

    TrigArray2D(const TrigArray2D<T> &other);

    ~TrigArray2D();

    TrigArray2D<T>& operator=(const TrigArray2D<T> &other);

    bool operator==(const TrigArray2D<T> &other) const;
    bool operator!=(const TrigArray2D<T> &other) const;

    const T& operator()(int i, int j) const;
    T& operator()(int i, int j);
    const T& at(int i, int j) const;

    QString toString() const;

    void set(int i, int j, const T &value);
    void setAll(const T &value);

    const T& get(int i, int j) const;

    void redimension(int dimension);

    const T* data() const;
    T* data();
    const T* constData() const;

    TrigArray2D<T> transpose() const;

private:
    /** The 1D array of entries in this TrigArray2D */
    QVector<T> array;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Construct a null TrigArray2D */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
TrigArray2D<T>::TrigArray2D() : TrigArray2DBase()
{}

/** Construct a TrigArray2D that holds a [dimension,dimension] square
    symmetric matrix of values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
TrigArray2D<T>::TrigArray2D(int dimension) : TrigArray2DBase(dimension)
{
    if (TrigArray2DBase::nRows() > 0)
    {
        const int dim = TrigArray2D::nRows();

        array = QVector<T>( (dim*dim + dim)/2 );
        array.squeeze();
    }
}

/** Construct a TrigArray2D that holds a square symmetric matrix  
    of dimension [dimension,dimension] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
TrigArray2D<T>::TrigArray2D(int dimension, const T &default_value)
               : TrigArray2DBase(dimension)
{
    if (TrigArray2DBase::nRows() > 0)
    {
        const int dim = TrigArray2D::nRows();
    
        array = QVector<T>( (dim*dim + dim)/2, default_value );
        array.squeeze();
    }
}

/** Copy constructor. Fast as this class is implicitly shared */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
TrigArray2D<T>::TrigArray2D(const TrigArray2D<T> &other)
               : TrigArray2DBase(other), array(other.array)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
TrigArray2D<T>::~TrigArray2D()
{}

/** Copy assignment operator - fast as this class is implicitly 
    shared */
template<class T>
SIRE_INLINE_TEMPLATE
TrigArray2D<T>& TrigArray2D<T>::operator=(const TrigArray2D<T> &other)
{
    TrigArray2DBase::operator=(other);
    array = other.array;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool TrigArray2D<T>::operator==(const TrigArray2D<T> &other) const
{
    return TrigArray2DBase::operator==(other) and array == other.array;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool TrigArray2D<T>::operator!=(const TrigArray2D<T> &other) const
{
    return TrigArray2DBase::operator!=(other) or array != other.array;
}

/** Return a reference to the object in the ith row and
    the jth column 
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_INLINE_TEMPLATE
const T& TrigArray2D<T>::operator()(int i, int j) const
{
    return array.constData()[ TrigArray2DBase::checkedOffset(i,j) ];
}

/** Return a reference to the object in the ith row and
    the jth column 
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_INLINE_TEMPLATE
T& TrigArray2D<T>::operator()(int i, int j)
{
    return array.data()[ TrigArray2DBase::checkedOffset(i,j) ];
}

/** Return a reference to the object in the ith row and
    the jth column 
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_INLINE_TEMPLATE
const T& TrigArray2D<T>::at(int i, int j) const
{
    return TrigArray2D<T>::operator()(i,j);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void TrigArray2D<T>::set(int i, int j, const T &value)
{
    this->operator()(i,j) = value;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& TrigArray2D<T>::get(int i, int j) const
{
    return this->operator()(i,j);
}

/** Redimension this TrigArray2D to [dimension,dimension] objects.
    This will keep any existing data in the top left of the
    matrix if the matrix gets bigger, or it will crop any
    extra data to the bottom right if the matrix gets smaller 
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void TrigArray2D<T>::redimension(int dimension)
{
    if (TrigArray2DBase::nRows() == dimension)
        return;

    TrigArray2D<T> new_array(dimension);
    
    //copy the data...
    dimension = qMin(dimension, TrigArray2DBase::nRows());
    
    T *data = new_array.data();
    const T *old_data = array.constData();
    
    for (int i=0; i<TrigArray2DBase::nRows(); ++i)
    {
        for (int j=0; j<TrigArray2DBase::nRows(); ++j)
        {
            data[this->map(i,j)] = old_data[this->map(i,j)];
        }
    }
    
    this->operator=(new_array);
}

/** Return a raw pointer into the data of the array. Remember
    that the data is in row-major order */
template<class T>
SIRE_INLINE_TEMPLATE
const T* TrigArray2D<T>::data() const
{
    return array.data();
}

/** Return a raw pointer into the data of the array. Remember
    that the data is in row-major order */
template<class T>
SIRE_INLINE_TEMPLATE
T* TrigArray2D<T>::data()
{
    return array.data();
}

/** Return a raw pointer into the data of the array. Remember
    that the data is in row-major order */
template<class T>
SIRE_INLINE_TEMPLATE
const T* TrigArray2D<T>::constData() const
{
    return array.constData();
}

/** Return the transpose of this matrix */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
TrigArray2D<T> TrigArray2D<T>::transpose() const
{
    return *this;
}

/** Set all values in this array equal to 'value' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void TrigArray2D<T>::setAll(const T &value)
{
    T *array_data = this->array.data();
    int count = this->array.count();
    
    for (int i=0; i<count; ++i)
        array_data[i] = value;
}

/** Return a string representation of this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString TrigArray2D<T>::toString() const
{
    if (this->nRows() == 0)
        return "( )";
        
    else if (this->nRows() == 1)
    {
        QStringList row;
        for (int j=0; j<this->nColumns(); ++j)
        {
            row.append( Sire::toString( this->operator()(0,j) ) );
        }
        
        return QString("( %1 )").arg( row.join(", ") );
    }

    QStringList rows;
    
    for (int i=0; i<this->nRows(); ++i)
    {
        QStringList row;
        
        for (int j=0; j<this->nColumns(); ++j)
        {
            row.append( Sire::toString( this->operator()(i,j) ) );
        }
        
        if (i == 0)
            rows.append( QString("/ %1 \\").arg( row.join(", ") ) );
        else if (i == this->nRows() - 1)
            rows.append( QString("\\ %1 /").arg( row.join(", ") ) );
        else
            rows.append( QString("| %1 |").arg( row.join(", ") ) );
    }
    
    return rows.join("\n");
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Serialise to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireBase::TrigArray2D<T> &array)
{
    SireStream::SharedDataStream sds(ds);
    
    sds << static_cast<const SireBase::TrigArray2DBase&>(array)
        << array.array;
        
    return ds;
}

/** Extract from a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireBase::TrigArray2D<T> &array)
{
    SireStream::SharedDataStream sds(ds);
    
    sds >> static_cast<SireBase::TrigArray2DBase&>(array)
        >> array.array;
        
    return ds;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

SIRE_EXPOSE_ALIAS( SireBase::TrigArray2D<double>, SireBase::TrigArray2D_double_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireBase::TrigArray2D<double>;
#endif

SIRE_END_HEADER

#endif
