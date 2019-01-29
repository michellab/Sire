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

#include "ffparameters.h"
#include "atomicffparameters.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireFF;
using namespace SireBase;
using namespace SireStream;

Q_DECLARE_METATYPE( AtomicFFParameters<double> );
Q_DECLARE_METATYPE( AtomicFFParametersArray<double> );

namespace SireFF
{
    template class AtomicFFParameters<double>;
    template class AtomicFFParametersArray<double>;
}

////////////
//////////// Implementation of FFParameters
////////////

static const RegisterMetaType<FFParameters> r_ffparams( MAGIC_ONLY,
                                                        FFParameters::typeName() );
                                                        
QDataStream &operator<<(QDataStream &ds, const FFParameters &ffparams)
{
    writeHeader(ds, r_ffparams, 1);
    
    ds << static_cast<const Property&>(ffparams);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, FFParameters &ffparams)
{
    VersionID v = readHeader(ds, r_ffparams);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(ffparams);
    }
    else
        throw version_error(v, "1", r_ffparams, CODELOC);
        
    return ds;
}

/** Constructor */
FFParameters::FFParameters() : Property()
{}

/** Copy constructor */
FFParameters::FFParameters(const FFParameters &other) : Property(other)
{}

/** Destructor */
FFParameters::~FFParameters()
{}

/** Copy assignment operator */
FFParameters& FFParameters::operator=(const FFParameters &other)
{
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool FFParameters::operator==(const FFParameters &other) const
{
    return Property::operator==(other);
}

/** Comparison operator */
bool FFParameters::operator!=(const FFParameters &other) const
{
    return not FFParameters::operator==(other);
}

const char* FFParameters::typeName()
{
    return "SireFF::FFParameters";
}

NullFFParameters FFParameters::null()
{
    return NullFFParameters();
}

////////////
//////////// Implementation of FFParametersArray
////////////

static const RegisterMetaType<FFParametersArray> r_ffparamsarray( MAGIC_ONLY,
                                                    FFParametersArray::typeName() );

QDataStream &operator<<(QDataStream &ds,
                                      const FFParametersArray &ffparams)
{
    writeHeader(ds, r_ffparamsarray, 1);
    
    ds << static_cast<const Property&>(ffparams);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, FFParametersArray &ffparams)
{
    VersionID v = readHeader(ds, r_ffparamsarray);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(ffparams);
    }
    else
        throw version_error(v, "1", r_ffparamsarray, CODELOC);
        
    return ds;
}

/** Constructor */
FFParametersArray::FFParametersArray() : Property()
{}

/** Copy constructor */
FFParametersArray::FFParametersArray(const FFParametersArray &other)
                  : Property(other)
{}

/** Destructor */
FFParametersArray::~FFParametersArray()
{}

/** Copy assignment operator */
FFParametersArray& FFParametersArray::operator=(const FFParametersArray &other)
{
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool FFParametersArray::operator==(const FFParametersArray &other) const
{
    return Property::operator==(other);
}

/** Comparison operator */
bool FFParametersArray::operator!=(const FFParametersArray &other) const
{
    return not FFParametersArray::operator==(other);
}

const char* FFParametersArray::typeName()
{
    return "SireFF::FFParametersArray";
}

/** Return the ith set of parameters */
FFParametersPtr FFParametersArray::at(int i) const
{
    return this->operator[](i);
}

NullFFParametersArray FFParametersArray::null()
{
    return NullFFParametersArray();
}

////////////
//////////// Implementation of NullFFParameters
////////////

static const RegisterMetaType<NullFFParameters> r_nullffparams;

QDataStream &operator<<(QDataStream &ds, 
                                      const NullFFParameters &nullffparams)
{
    writeHeader(ds, r_nullffparams, 1);
    
    ds << static_cast<const FFParameters&>(nullffparams);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, NullFFParameters &nullffparams)
{
    VersionID v = readHeader(ds, r_nullffparams);
    
    if (v == 1)
    {
        ds >> static_cast<FFParameters&>(nullffparams);
    }
    else
        throw version_error(v, "1", r_nullffparams, CODELOC);
        
    return ds;
}

/** Constructor */
NullFFParameters::NullFFParameters() : ConcreteProperty<NullFFParameters,FFParameters>()
{}

/** Copy constructor */
NullFFParameters::NullFFParameters(const NullFFParameters &other)
                 : ConcreteProperty<NullFFParameters,FFParameters>(other)
{}

/** Destructor */
NullFFParameters::~NullFFParameters()
{}

const char* NullFFParameters::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullFFParameters>() );
}

/** Copy assignment operator */
NullFFParameters& NullFFParameters::operator=(const NullFFParameters &other)
{
    FFParameters::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullFFParameters::operator==(const NullFFParameters &other) const
{
    return FFParameters::operator==(other);
}

/** Comparison operator */
bool NullFFParameters::operator!=(const NullFFParameters &other) const
{
    return not NullFFParameters::operator==(other);
}

/** Return a null array */
FFParametersArrayPtr NullFFParameters::toArray() const
{
    return NullFFParametersArray();
}

////////////
//////////// Implementation of NullFFParametersArray
////////////

static const RegisterMetaType<NullFFParametersArray> r_nullffarray;


QDataStream &operator<<(QDataStream &ds, 
                                      const NullFFParametersArray &nullffparams)
{
    writeHeader(ds, r_nullffarray, 1);
    
    ds << static_cast<const FFParametersArray&>(nullffparams);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, 
                                      NullFFParametersArray &nullffparams)
{
    VersionID v = readHeader(ds, r_nullffarray);
    
    if (v == 1)
    {
        ds >> static_cast<FFParametersArray&>(nullffparams);
    }
    else
        throw version_error(v, "1", r_nullffarray, CODELOC);
        
    return ds;
}

/** Constructor */
NullFFParametersArray::NullFFParametersArray()
                      : ConcreteProperty<NullFFParametersArray,FFParametersArray>()
{}

/** Copy constructor */
NullFFParametersArray::NullFFParametersArray(const NullFFParametersArray &other)
                : ConcreteProperty<NullFFParametersArray,FFParametersArray>(other)
{}

/** Destructor */
NullFFParametersArray::~NullFFParametersArray()
{}

const char* NullFFParametersArray::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullFFParametersArray>() );
}

/** Copy assignment operator */
NullFFParametersArray& NullFFParametersArray::operator=(
                                                    const NullFFParametersArray &other)
{
    FFParametersArray::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullFFParametersArray::operator==(const NullFFParametersArray &other) const
{
    return FFParametersArray::operator==(other);
}

/** Comparison operator */
bool NullFFParametersArray::operator!=(const NullFFParametersArray &other) const
{
    return not NullFFParametersArray::operator==(other);
}

/** Return the number of groups of parameters in this array */
int NullFFParametersArray::count() const
{
    return 0;
}

/** Return whether or not this array is empty */
bool NullFFParametersArray::isEmpty() const
{
    return true;
}

/** Return the ith set of parameters 

    \throw SireError::invalid_index
*/
FFParametersPtr NullFFParametersArray::operator[](int i) const
{
    throw SireError::invalid_index( QObject::tr(
            "Cannot access element %1 from a NullFFParametersArray.")
                .arg(i), CODELOC );
                
    return FFParametersPtr();
}

/** Append the passed set of parameters onto the end of this array */
void NullFFParametersArray::append(const FFParameters &params)
{}

/** Append the passed array of parameters onto the end of this array */
void NullFFParametersArray::append(const FFParametersArray &params)
{}

/** Update the parameters at index 'idx' so that they equal 'params' */
void NullFFParametersArray::update(int idx, const FFParameters &params)
{}

/** Update the parameters with the passed indicies so that they equal 
    the values in 'params' */
void NullFFParametersArray::update(const QVarLengthArray<int> &idxs, 
                                   const FFParametersArray &params)
{}

/** Remove the parameters at index 'idx' */
void NullFFParametersArray::remove(int idx)
{}

/** Remove the parameters whose indicies are in 'idxs' */
void NullFFParametersArray::remove(const QVarLengthArray<int> &idxs)
{}

/** Remove all of the parameters from this array */
void NullFFParametersArray::removeAll()
{}

