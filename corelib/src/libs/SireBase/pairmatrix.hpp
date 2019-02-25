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

#ifndef SIREBASE_PAIRMATRIX_HPP
#define SIREBASE_PAIRMATRIX_HPP

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireBase
{

/** This provides a 2D matrix that contains information about all pairs of 
    two groups. This matrix is designed to be used in a pair-pair loop, 
    where each element is accessed sequentially, and it is thus highly optimised
    for such use. It does provide random access via indicies, but this will be slightly
    slower (though still quite quick!).

    If you want a general purpose 2D, implicitly shared array, then 
    use Array2D

    @author Christopher Woods
*/
template<class T>
class PairMatrix
{
public:
    /** Construct a null PairMatrix */
    PairMatrix();

    /** Construct a PairMatrix that holds n_outer elements on the outer loop,
        and n_inner elements on the inner loop. If either or these are 0 then
        an invalid PairMatrix will be created. Also note that you can redimension
        this PairMatrix at any time using a potentially very quick function. */
    PairMatrix(unsigned int n_outer, unsigned int n_inner);

    /** Copy constructor. This is not very fast... */
    PairMatrix(const PairMatrix<T> &other);

    ~PairMatrix();

    /** Return a reference to the element at index 'i,j', where
        'i' is the outer loop, and j is the inner loop. As this is optimised
        for speed, there is no checking to ensure that this array position
        is valid! */
    const T& operator()(unsigned int i, unsigned int j) const;
    T& operator()(unsigned int i, unsigned int j);

    /** Return a reference to the element at index 'j', where 'j'
        refers to the inner loop, and the outer loop index has been
        set via 'setOuterIndex()' */
    const T& operator()(unsigned int j) const;
    T& operator()(unsigned int j);

    const T& operator[](unsigned int j) const;
    T& operator[](unsigned int j);

    /** Return a const reference to the element at index 'i,j'. */
    const T& at(unsigned int i, unsigned int j) const;

    /** Redimension this PairMatrix to n_outer by n_inner elements. This
        function is very quick if the total number of elements does not increase
        as the existing array is just reinterpreted. Reallocation will only
        occur if the number of elements becomes greater than n_elements. Note
        that you should repopulate the PairMatrix after it has been
        redimensioned. Also note that both i and j should be greater than 0,
        and that any index set using 'setOuterIndex' will now be invalid. */
    void redimension(unsigned int i, unsigned int j);

    /** Return the dimensions of the array */
    unsigned int nOuter() const;
    unsigned int nInner() const;

    /** Set the index of the outer loop. This is used in combination with
        operator[](unsigned int j) which then only varies the index of the
        inner loop */
    void setOuterIndex(unsigned int i);

private:
    /** Pointer to the array of entries in this PairMatrix */
    T *array;

    /** The current outer index (multiplied by n_outer) */
    unsigned int outer_index;

    /** The dimensions of the array */
    unsigned int n_outer, n_inner;

    /** The number of elements in the array */
    unsigned int n_elements;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PairMatrix<T>::PairMatrix() : array(0), n_outer(0), n_inner(0), n_elements(0)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PairMatrix<T>::PairMatrix(unsigned int i, unsigned int j)
              : n_outer(i), n_inner(j), n_elements(0)
{
    n_elements = n_outer * n_inner;
    array = new T[n_elements];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PairMatrix<T>::PairMatrix(const PairMatrix<T> &other)
              : n_outer(other.n_outer), n_inner(other.n_inner),
                n_elements(other.n_elements)
{
    if (n_elements > 0)
    {
        array = new T[n_elements];
        
        //have to manual copy as T may not be memcpy'able
        for (quint32 i=0; i<n_elements; ++i)
        {
            array[i] = other.array[i];
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PairMatrix<T>::~PairMatrix()
{
    delete[] array;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
unsigned int PairMatrix<T>::nOuter() const
{
    return n_outer;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
unsigned int PairMatrix<T>::nInner() const
{
    return n_inner;
}

template<class T>
SIRE_INLINE_TEMPLATE
const T& PairMatrix<T>::operator()(unsigned int i, unsigned int j) const
{
    return array[i*n_inner + j];
}

template<class T>
SIRE_INLINE_TEMPLATE
T& PairMatrix<T>::operator()(unsigned int i, unsigned int j)
{
    return array[i*n_inner + j];
}

template<class T>
SIRE_INLINE_TEMPLATE
const T& PairMatrix<T>::operator()(unsigned int j) const
{
    return array[outer_index + j];
}

template<class T>
SIRE_INLINE_TEMPLATE
T& PairMatrix<T>::operator()(unsigned int j)
{
    return array[outer_index + j];
}

template<class T>
SIRE_INLINE_TEMPLATE
const T& PairMatrix<T>::operator[](unsigned int j) const
{
    return array[outer_index + j];
}

template<class T>
SIRE_INLINE_TEMPLATE
T& PairMatrix<T>::operator[](unsigned int j)
{
    return array[outer_index + j];
}

template<class T>
SIRE_INLINE_TEMPLATE
void PairMatrix<T>::setOuterIndex(unsigned int i)
{
    outer_index = n_inner * i;
}

template<class T>
SIRE_INLINE_TEMPLATE
const T& PairMatrix<T>::at(unsigned int i, unsigned int j) const
{
    return PairMatrix::operator()(i,j);
}

template<class T>
SIRE_INLINE_TEMPLATE
void PairMatrix<T>::redimension(unsigned int i, unsigned int j)
{
    n_outer = i;
    n_inner = j;

    outer_index = 0;

    unsigned int new_n_elements = n_inner * n_outer;
    if (new_n_elements > n_elements)
    {
        n_elements = new_n_elements;

        delete[] array;
        array = new T[n_elements];
    }
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_END_HEADER

#endif
