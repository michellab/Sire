/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREBASE_QUICKCOPY_HPP
#define SIREBASE_QUICKCOPY_HPP

#include <cstring>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireBase
{

/** This function copies 'nvalues' from 'source' to 'destination'. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* quickCopy(T *destination, const T *source, int nvalues)
{
    //if (QTypeInfo<T>::isComplex)
    //{
        //we need to copy the values
        for (int i=0; i<nvalues; ++i)
        {
            destination[i] = source[i];
        }
        
        return destination;
    //}
    //else
    //{
    //    //we can use memcpy
    //    void *output = std::memcpy(destination, source, nvalues*sizeof(T));
    //    
    //    return static_cast<T*>(output);
    //}
}

}

SIRE_END_HEADER

#endif
