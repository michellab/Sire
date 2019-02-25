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

#ifndef SIREBASE_ARRAY2D_HPP
#define SIREBASE_ARRAY2D_HPP

#include <QStringList>

#include "array2d.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
template<class T>
class Array2D;
}

template<class T>
QDataStream& operator<<(QDataStream&, const SireBase::Array2D<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireBase::Array2D<T>&);

namespace SireBase
{

/** This provides a 2D matrix of objects of type T. All objects
    are packed together in memory in row-major order, e.g.
    
        j0   j1   j2
    i0   0    1    2
    i1   3    4    5
    i2   6    7    8
    
    Arranged in memory as;
    
    [ (i0,j0), (i0,j1), (i0,j2), (i1,j0),  ... etc. ]
    [    0   ,    1   ,    2   ,   3    ,  ... etc. ]
    
    @author Christopher Woods
*/
template<class T>
class Array2D : public Array2DBase
{

friend SIREBASE_EXPORT QDataStream& ::operator<<<>(QDataStream&, const Array2D<T>&);
friend SIREBASE_EXPORT QDataStream& ::operator>><>(QDataStream&, Array2D<T>&);

public:
    Array2D();

    Array2D(int nrows, int ncolumns);
    
    Array2D(int nrows, int ncolumns, const T &default_value);

    Array2D(const Array2D<T> &other);

    ~Array2D();

    Array2D<T>& operator=(const Array2D<T> &other);

    bool operator==(const Array2D<T> &other) const;
    bool operator!=(const Array2D<T> &other) const;

    const T& operator()(int i, int j) const;
    T& operator()(int i, int j);
    const T& at(int i, int j) const;

    QString toString() const;

    void set(int i, int j, const T &value);
    void setAll(const T &value);

    const T& get(int i, int j) const;

    void redimension(int nrows, int ncolumns);

    const T* data() const;
    T* data();
    const T* constData() const;

    const T* row(int i) const;
    T* row(int i);
    const T* constRow(int i) const;

    Array2D<T> transpose() const;

private:
    /** The 1D array of entries in this Array2D */
    QVector<T> array;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Construct a null Array2D */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Array2D<T>::Array2D() : Array2DBase()
{}

/** Construct a Array2D that holds nrow rows of ncolumn columns. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Array2D<T>::Array2D(int nrows, int ncolumns)
           : Array2DBase(nrows,ncolumns)
{
    if (this->nRows() * this->nColumns() > 0)
    {
        array = QVector<T>(this->nRows() * this->nColumns());
        array.squeeze();
    }
}

/** Construct a Array2D that holds nrow rows of ncolumn columns. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Array2D<T>::Array2D(int nrows, int ncolumns, const T &default_value)
           : Array2DBase(nrows,ncolumns)
{
    if (this->nRows() * this->nColumns() > 0)
    {
        array = QVector<T>(this->nRows() * this->nColumns(), default_value);
        array.squeeze();
    }
}

/** Copy constructor. Fast as this class is implicitly shared */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Array2D<T>::Array2D(const Array2D<T> &other)
           : Array2DBase(other), array(other.array)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Array2D<T>::~Array2D()
{}

/** Copy assignment operator - fast as this class is implicitly 
    shared */
template<class T>
SIRE_INLINE_TEMPLATE
Array2D<T>& Array2D<T>::operator=(const Array2D<T> &other)
{
    Array2DBase::operator=(other);
    array = other.array;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Array2D<T>::operator==(const Array2D<T> &other) const
{
    return Array2DBase::operator==(other) and array == other.array;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Array2D<T>::operator!=(const Array2D<T> &other) const
{
    return Array2DBase::operator!=(other) or array != other.array;
}

/** Return a reference to the object in the ith row and
    the jth column 
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_INLINE_TEMPLATE
const T& Array2D<T>::operator()(int i, int j) const
{
    Array2DBase::assertValidIndex(i,j);
    return array.constData()[ Array2DBase::map(i,j) ];
}

/** Return a reference to the object in the ith row and
    the jth column 
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_INLINE_TEMPLATE
T& Array2D<T>::operator()(int i, int j)
{
    Array2DBase::assertValidIndex(i,j);
    return array.data()[ Array2DBase::map(i,j) ];
}

/** Return a reference to the object in the ith row and
    the jth column 
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_INLINE_TEMPLATE
const T& Array2D<T>::at(int i, int j) const
{
    return Array2D<T>::operator()(i,j);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Array2D<T>::set(int i, int j, const T &value)
{
    this->operator()(i,j) = value;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Array2D<T>::get(int i, int j) const
{
    return this->operator()(i,j);
}

/** Redimension this Array2D to nrows by ncolumns objects.
    This will keep any existing data in the top left of the
    matrix if the matrix gets bigger, or it will crop any
    extra data to the bottom right if the matrix gets smaller 
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Array2D<T>::redimension(int nrows, int ncolumns)
{
    if (nrows < 0)
        nrows = 0;
        
    if (ncolumns < 0)
        ncolumns = 0;

    if (nrows == this->nRows() and ncolumns == this->nColumns())
        return;

    Array2D<T> new_array(nrows, ncolumns);
    
    //copy the data...
    nrows = qMin(nrows, this->nRows());
    ncolumns = qMin(ncolumns, this->nColumns());
    
    T *data = new_array.data();
    const T *old_data = array.constData();
    
    for (int i=0; i<nrows; ++i)
    {
        for (int j=0; j<ncolumns; ++j)
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
const T* Array2D<T>::data() const
{
    return array.data();
}

/** Return a raw pointer into the data of the array. Remember
    that the data is in row-major order */
template<class T>
SIRE_INLINE_TEMPLATE
T* Array2D<T>::data()
{
    return array.data();
}

/** Return a raw pointer into the data of the array. Remember
    that the data is in row-major order */
template<class T>
SIRE_INLINE_TEMPLATE
const T* Array2D<T>::constData() const
{
    return array.constData();
}

/** Return a raw pointer to the first item in row 'i' */
template<class T>
SIRE_INLINE_TEMPLATE
const T* Array2D<T>::row(int i) const
{
    return array.data() + Array2DBase::map(i,0);
}

/** Return a raw pointer to the first item in row 'i' */
template<class T>
SIRE_INLINE_TEMPLATE
T* Array2D<T>::row(int i)
{
    return array.data() + Array2DBase::map(i,0);
}

/** Return a raw pointer to the first item in row 'i' */
template<class T>
SIRE_INLINE_TEMPLATE
const T* Array2D<T>::constRow(int i) const
{
    return array.constData() + Array2DBase::map(i,0);
}

/** Return the transpose of this matrix */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Array2D<T> Array2D<T>::transpose() const
{
    if (array.isEmpty())
        return Array2D<T>();
    
    Array2D trans( this->nColumns(), this->nRows() );
    
    T *new_array = trans.data();
    const T *old_array = array.constData();
    
    for (int i=0; i < this->nRows(); ++i)
    {
        for (int j=0; j < this->nColumns(); ++j)
        {
            new_array[trans.map(j,i)] = old_array[this->map(i,j)];
        }
    }
    
    return trans;
}

/** Set all values in this array equal to 'value' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Array2D<T>::setAll(const T &value)
{
    T *array_data = this->array.data();
    int count = this->array.count();
    
    for (int i=0; i<count; ++i)
        array_data[i] = value;
}

/** Return a string representation of this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString Array2D<T>::toString() const
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
QDataStream& operator<<(QDataStream &ds, const SireBase::Array2D<T> &array)
{
    SireStream::SharedDataStream sds(ds);
    
    sds << static_cast<const SireBase::Array2DBase&>(array)
        << array.array;
        
    return ds;
}

/** Extract from a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireBase::Array2D<T> &array)
{
    SireStream::SharedDataStream sds(ds);
    
    sds >> static_cast<SireBase::Array2DBase&>(array)
        >> array.array;
        
    return ds;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

SIRE_EXPOSE_ALIAS( SireBase::Array2D<double>, SireBase::Array2D_double_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireBase::Array2D<double>;
#endif

SIRE_END_HEADER

#endif
