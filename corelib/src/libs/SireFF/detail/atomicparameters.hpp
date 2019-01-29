/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREFF_DETAIL_ATOMICPARAMETERS_HPP
#define SIREFF_DETAIL_ATOMICPARAMETERS_HPP

#include <QSet>

#include "SireBase/packedarray2d.hpp"

SIRE_BEGIN_HEADER

namespace SireFF
{
namespace detail
{
template<class PARAM>
class AtomicParameters;
}
}

template<class PARAM>
QDataStream& operator<<(QDataStream&, const SireFF::detail::AtomicParameters<PARAM>&);
template<class PARAM>
QDataStream& operator>>(QDataStream&, SireFF::detail::AtomicParameters<PARAM>&);

namespace SireFF
{

namespace detail
{

using SireBase::PackedArray2D;

SIREFF_EXPORT bool selectedAll(const QSet<quint32> &idxs, quint32 n);

/** This class hold parameters of type PARAM, one for each atom

    @author Christopher Woods
*/
template<class PARAM>
class AtomicParameters
{

friend QDataStream& ::operator<<<>(QDataStream&, const AtomicParameters<PARAM>&);
friend QDataStream& ::operator>><>(QDataStream&, AtomicParameters<PARAM>&);

public:

    typedef PARAM Parameter;
    typedef typename SireBase::PackedArray2D<PARAM> Parameters;
    typedef typename Parameters::Array Array;

    AtomicParameters();
    
    AtomicParameters(const Parameters &parameters);
    
    AtomicParameters(const AtomicParameters<PARAM> &other);
    
    ~AtomicParameters();
    
    AtomicParameters<PARAM>& operator=(const AtomicParameters<PARAM> &other);
    
    bool operator==(const AtomicParameters<PARAM> &other) const;
    bool operator!=(const AtomicParameters<PARAM> &other) const;

    int nGroups() const;

    const Parameters& atomicParameters() const;

    void setAtomicParameters(const AtomicParameters<PARAM> &other);

    bool changedAllGroups(const AtomicParameters<PARAM> &params) const;
    
    QSet<quint32> getChangedGroups(const AtomicParameters<PARAM> &params) const;
    
    void addChangedGroups(const AtomicParameters<PARAM> &params,
                          QSet<quint32> &changed_groups) const;
    
    AtomicParameters<PARAM> applyMask(const QSet<quint32> &idxs) const;

protected:
    /** The atomic parameters, arranged by CutGroup */
    Parameters params;

};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
AtomicParameters<PARAM>::AtomicParameters()
{}

/** Construct from the passed parameters */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
AtomicParameters<PARAM>::AtomicParameters(
                  const typename AtomicParameters<PARAM>::Parameters &parameters)
              : params(parameters)
{}

/** Copy constructor */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
AtomicParameters<PARAM>::AtomicParameters(const AtomicParameters<PARAM> &other)
                        : params(other.params)
{}

/** Destructor */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
AtomicParameters<PARAM>::~AtomicParameters()
{}

/** Copy assignment operator */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
AtomicParameters<PARAM>& AtomicParameters<PARAM>::operator=(
                                        const AtomicParameters<PARAM> &other)
{
    params = other.params;
    
    return *this;
}

/** Comparison operator */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
bool AtomicParameters<PARAM>::operator==(const AtomicParameters<PARAM> &other) const
{
    return params == other.params;
}

/** Comparison operator */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
bool AtomicParameters<PARAM>::operator!=(const AtomicParameters<PARAM> &other) const
{
    return params != other.params;
}

/** Return the number of groups */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
int AtomicParameters<PARAM>::nGroups() const
{
    return params.count();
}

/** Return the actual atomic parameters */
template<class PARAM>
SIRE_INLINE_TEMPLATE
const typename AtomicParameters<PARAM>::Parameters& 
AtomicParameters<PARAM>::atomicParameters() const
{
    return params;
}

/** Set the atomic parameters */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
void AtomicParameters<PARAM>::setAtomicParameters(const AtomicParameters<PARAM> &other)
{
    params = other.params;
}

/** Return whether or not all of the groups have changed parameters
    compared to 'params' */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
bool AtomicParameters<PARAM>::changedAllGroups(
                                  const AtomicParameters<PARAM> &other) const
{
    int ngroups = qMin(params.count(), other.params.count());
        
    const typename Parameters::Array *this_array = params.constData();
    const typename Parameters::Array *other_array = other.params.constData();
                                                    
    for (int i=0; i<ngroups; ++i)
    {
        if (this_array[i] == other_array[i])
            return false;
    }
        
    return true;
}

/** Add the changed groups from 'params' to the set 'changed_groups' */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
void AtomicParameters<PARAM>::addChangedGroups(
                                   const AtomicParameters<PARAM> &other,
                                   QSet<quint32> &changed_groups) const
{
    //now add on the groups that have changed parameters
    quint32 ngroups = this->nGroups();
    
    if (SireFF::detail::selectedAll(changed_groups, ngroups))
        return;
        
    quint32 nsharedgroups = qMin(ngroups, quint32(other.params.count()));
   
    const typename Parameters::Array *this_array = params.constData();
    const typename Parameters::Array *other_array = other.params.constData();
    
    for (quint32 i=0; i<nsharedgroups; ++i)
    {
        if (not changed_groups.contains(i) 
            and (this_array[i] != other_array[i]))
        {
            changed_groups.insert(i);
            
            if (SireFF::detail::selectedAll(changed_groups, ngroups))
                return;
        }
    }
}

/** Return the set of groups that have changed parameters */
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
QSet<quint32> AtomicParameters<PARAM>::getChangedGroups(
                                    const AtomicParameters<PARAM> &other) const
{
    QSet<quint32> changed_groups;
    changed_groups.reserve(params.count());
    
    this->addChangedGroups(other, changed_groups);
    
    return changed_groups;
}

/** Mask this set of parameters so that only the CutGroups at indicies 
    'idxs' are present */ 
template<class PARAM>
SIRE_OUTOFLINE_TEMPLATE
AtomicParameters<PARAM> AtomicParameters<PARAM>::applyMask(
                                          const QSet<quint32> &idxs) const
{
    if (SireFF::detail::selectedAll(idxs, params.count()))
        //all of the groups have been selected
        return *this;
    
    else if (idxs.isEmpty())
        //definitely nothing has been selected
        return AtomicParameters<PARAM>();
    
    else if (idxs.count() == 1)
        //return the single array
        return AtomicParameters<PARAM>( params[ *(idxs.constBegin()) ] );
    
    //otherwise, some are marked - apply the mask
    //mask by the indicies
    quint32 ngroups = params.count();
    
    QVector< typename Parameters::Array > group_params;
    group_params.reserve(idxs.count());
    
    const typename Parameters::Array *this_array = params.constData();
    
    for (quint32 i=0; i<ngroups; ++i)
    {
        if (idxs.contains(i))
            group_params.append( this_array[i] );
    }
    
    Parameters new_params(group_params);
    
    return AtomicParameters<PARAM>(new_params);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} //end of namespace detail

} //end of namespace SireFF

/** Serialise to a binary datastream */
template<class PARAM>
QDataStream& operator<<(QDataStream &ds,
                        const SireFF::detail::AtomicParameters<PARAM> &params)
{
    ds << params.params;

    return ds;
}

/** Extract from a binary datastream */
template<class PARAM>
QDataStream& operator>>(QDataStream &ds,
                        SireFF::detail::AtomicParameters<PARAM> &params)
{
    ds >> params.params;

    return ds;
}

SIRE_END_HEADER

#endif
