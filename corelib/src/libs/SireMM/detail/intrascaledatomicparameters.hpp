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

#ifndef SIREMM_DETAIL_INTRASCALEDATOMICPARAMETERS_HPP
#define SIREMM_DETAIL_INTRASCALEDATOMICPARAMETERS_HPP

#include "SireBase/propertymap.h"

#include "SireMol/partialmolecule.h"

#include "SireFF/detail/atomicparameters.hpp"
#include "SireFF/detail/atomicparameters3d.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{
namespace detail
{
template<class SCALE_FACTORS>
class IntraScaledParameters;

template<class ATOMPARAM, class INTRASCALE>
class IntraScaledAtomicParameters;
}
}

template<class SCALE_FACTORS>
QDataStream& operator<<(QDataStream&, 
                        const SireMM::detail::IntraScaledParameters<SCALE_FACTORS>&);
template<class SCALE_FACTORS>
QDataStream& operator>>(QDataStream&, 
                        SireMM::detail::IntraScaledParameters<SCALE_FACTORS>&);

template<class ATOMPARAM, class INTRASCALE>
QDataStream& operator<<(QDataStream&,
              const SireMM::detail::IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>&);
template<class ATOMPARAM, class INTRASCALE>
QDataStream& operator>>(QDataStream&,
              SireMM::detail::IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>&);

namespace SireMM
{

namespace detail
{

using SireBase::PropertyName;

using SireMol::PartialMolecule;
using SireMol::CGIdx;

class SIREMM_EXPORT IntraScaleParameterName
{
public:
    IntraScaleParameterName()
    {}
    
    ~IntraScaleParameterName()
    {}
    
    const QString& intraScaleFactors() const
    {
        return nbscl_param;
    }

private:
    static QString nbscl_param;
};

/** This class represents parameters that are scaled using information
    stored in the object of type 'SCALE_FACTORS'
    
    @author Christopher Woods
*/ 
template<class SCALE_FACTORS>
class IntraScaledParameters
{

friend QDataStream& ::operator<<<>(QDataStream&, 
                                   const IntraScaledParameters<SCALE_FACTORS>&);
friend QDataStream& ::operator>><>(QDataStream&, 
                                   IntraScaledParameters<SCALE_FACTORS>&);

public:
    typedef SCALE_FACTORS ScaleFactors;

    IntraScaledParameters();
    
    IntraScaledParameters(const PartialMolecule &molecule,
                          const PropertyName &scale_property);
                          
    IntraScaledParameters(const SCALE_FACTORS &scale_factors);
                          
    IntraScaledParameters(const IntraScaledParameters<SCALE_FACTORS> &other);
    
    ~IntraScaledParameters();
    
    IntraScaledParameters<SCALE_FACTORS> operator=(
                        const IntraScaledParameters<SCALE_FACTORS> &other);
                                
    bool operator==(const IntraScaledParameters<SCALE_FACTORS> &other) const;
    bool operator!=(const IntraScaledParameters<SCALE_FACTORS> &other) const;
    
    int nGroups() const;
    
    const ScaleFactors& intraScaleFactors() const;
    
    void setIntraScaleFactors(const IntraScaledParameters<SCALE_FACTORS> &other);
    
    bool changedAllGroups(const IntraScaledParameters<SCALE_FACTORS> &other) const;
    
    QSet<quint32> getChangedGroups(
                    const IntraScaledParameters<SCALE_FACTORS> &params) const;
    
    void addChangedGroups(
            const IntraScaledParameters<SCALE_FACTORS> &params,
            QSet<quint32> &changed_groups) const;
    
    IntraScaledParameters<SCALE_FACTORS> applyMask(const QSet<quint32> &idxs) const;

protected:
    /** The intramolecular inter-atomic scale factors to apply to 
        the potential between intramolecular pairs of atoms */
    ScaleFactors sclfactors;
};

/** Atomic parameters that use a scale factor
    for intramolecular atom pairs
    
    @author Christopher Woods
*/
template<class ATOMPARAM, class INTRASCALE>
class IntraScaledAtomicParameters : public ATOMPARAM, public INTRASCALE
{

friend QDataStream& ::operator<<<>(QDataStream&, 
                         const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>&);
friend QDataStream& ::operator>><>(QDataStream&, 
                         IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>&);

public:
    typedef typename ATOMPARAM::Parameter Parameter;
    typedef typename ATOMPARAM::Parameters Parameters;
    typedef typename INTRASCALE::ScaleFactors ScaleFactors;

    IntraScaledAtomicParameters();
                          
    IntraScaledAtomicParameters(const ATOMPARAM &parameters,
                                const INTRASCALE &sclfactors);
                          
    IntraScaledAtomicParameters(
            const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other);
    
    ~IntraScaledAtomicParameters();
    
    IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> operator=(
                        const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other);
                                
    bool operator==(
            const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other) const;
    bool operator!=(
            const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other) const;
    
    int nGroups() const;
    
    bool changedAllGroups(
            const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other) const;
    
    QSet<quint32> getChangedGroups(
            const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &params) const;
    
    void addChangedGroups(
            const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &params,
            QSet<quint32> &changed_groups) const;
    
    IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> 
    applyMask(const QSet<quint32> &idxs) const;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/////////
///////// Implementation of IntraScaledParameters
/////////

/** Null constructor */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledParameters<SCALE_FACTORS>::IntraScaledParameters()
{}

/** Construct for the molecule 'molecule' using the specified properties
    to find the coordinates and scale factor properties
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledParameters<SCALE_FACTORS>::IntraScaledParameters(
                              const PartialMolecule &molecule,
                              const PropertyName &scale_property)
{
    const SireBase::Property &property = molecule.property(scale_property);
    sclfactors = property.asA<SCALE_FACTORS>();
}

/** Construct by combining some AtomicParameters3D with some scale factors */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledParameters<SCALE_FACTORS>::IntraScaledParameters(
                              const SCALE_FACTORS &scale_factors)
                : sclfactors(scale_factors)
{}
                              
/** Copy constructor */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledParameters<SCALE_FACTORS>::IntraScaledParameters(
                              const IntraScaledParameters<SCALE_FACTORS> &other)
                : sclfactors(other.sclfactors)
{}

/** Destructor */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledParameters<SCALE_FACTORS>::~IntraScaledParameters()
{}

/** Copy assignment operator */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledParameters<SCALE_FACTORS> IntraScaledParameters<SCALE_FACTORS>::operator=(
                            const IntraScaledParameters<SCALE_FACTORS> &other)
{
    sclfactors = other.sclfactors;
    return *this;
}
                 
/** Comparison operator */           
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
bool IntraScaledParameters<SCALE_FACTORS>::operator==(
                            const IntraScaledParameters<SCALE_FACTORS> &other) const
{
    return sclfactors == other.sclfactors;
}

/** Comparison operator */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
bool IntraScaledParameters<SCALE_FACTORS>::operator!=(
                            const IntraScaledParameters<SCALE_FACTORS> &other) const
{
    return sclfactors != other.sclfactors;
}

/** Return the number of groups in the molecule */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
int IntraScaledParameters<SCALE_FACTORS>::nGroups() const
{
    return sclfactors.nGroups();
}

/** Return the inter-atomic intramolecular scale factors for the 
    intramolecular atom-atom interactions */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
const SCALE_FACTORS& IntraScaledParameters<SCALE_FACTORS>::intraScaleFactors() const
{
    return sclfactors;
}

/** Set the scale factors */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
void IntraScaledParameters<SCALE_FACTORS>::setIntraScaleFactors(
                                    const IntraScaledParameters<SCALE_FACTORS> &other)
{
    sclfactors = other.sclfactors;
}

/** Return whether or not all CutGroups in this molecule have changed compared
    to 'other' */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
bool IntraScaledParameters<SCALE_FACTORS>::changedAllGroups(
                        const IntraScaledParameters<SCALE_FACTORS> &other) const
{
    return sclfactors != other.sclfactors;
}

/** Add to 'changed_groups' the indicies of groups that have changed
    compared to 'other' */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
void IntraScaledParameters<SCALE_FACTORS>::addChangedGroups(
                      const IntraScaledParameters<SCALE_FACTORS> &other,
                      QSet<quint32> &changed_groups) const
{
    if (sclfactors != other.sclfactors)
    {
        //all groups have changed
        for (CGIdx i(0); i<this->nGroups(); ++i)
        {
            changed_groups.insert(i);
        }
        
        return;
    }
}

/** Return the indicies of CutGroups that have changed compared to 'other' */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
QSet<quint32> IntraScaledParameters<SCALE_FACTORS>::getChangedGroups(
                        const IntraScaledParameters<SCALE_FACTORS> &other) const
{
    QSet<quint32> changed_groups;
 
    if (sclfactors != other.sclfactors)
    {
        changed_groups.reserve(this->nGroups());
        this->addChangedGroups(other, changed_groups);
    }
    
    return changed_groups;
}

/** Mask these parameters so that only the atomic parameters for the 
    CutGroups whose indicies are in 'cgidxs' are present. */
template<class SCALE_FACTORS>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledParameters<SCALE_FACTORS> 
IntraScaledParameters<SCALE_FACTORS>::applyMask(const QSet<quint32>&) const
{
    //there is no way to mask the intramolecular scale factors
    return *this;
}

/////////
///////// Implementation of IntraScaledAtomicParameters
/////////

/** Constructor */
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::IntraScaledAtomicParameters()
              : ATOMPARAM(), INTRASCALE()
{}

/** Construct using the passed parameters */                    
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::IntraScaledAtomicParameters(
                            const ATOMPARAM &parameters,
                            const INTRASCALE &sclfactors)
              : ATOMPARAM(parameters),
                INTRASCALE(sclfactors)
{}
 
/** Copy constructor */                     
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::IntraScaledAtomicParameters(
        const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other)
              : ATOMPARAM(other), INTRASCALE(other)
{}

/** Destructor */
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::~IntraScaledAtomicParameters()
{}

/** Copy assignment operator */
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> 
IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::operator=(
                    const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other)
{
    ATOMPARAM::operator=(other);
    INTRASCALE::operator=(other);
    
    return *this;
}

/** Comparison operator */           
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
bool IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::operator==(
        const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other) const
{
    return ATOMPARAM::operator==(other) and
           INTRASCALE::operator==(other);
}

/** Comparison operator */
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
bool IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::operator!=(
        const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other) const
{
    return ATOMPARAM::operator!=(other) or
           INTRASCALE::operator!=(other);
}

/** Return the number of groups in the molecule */
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
int IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::nGroups() const
{
    BOOST_ASSERT(ATOMPARAM::nGroups() == INTRASCALE::nGroups());
                 
    return ATOMPARAM::nGroups();
}

/** Return whether all groups have changed compared to 'other' */
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
bool IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::changedAllGroups(
        const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other) const
{
    return INTRASCALE::changedAllGroups(other) or
           ATOMPARAM::changedAllGroups(other);
}

/** Return the set of all groups that have changed */
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
QSet<quint32> IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::getChangedGroups(
        const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other) const
{
    QSet<quint32> changed_groups 
                = INTRASCALE::getChangedGroups(other);
                
    ATOMPARAM::addChangedGroups(other, changed_groups);
    
    return changed_groups;
}

/** Add the changed group to the set 'changed_groups' */
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
void IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::addChangedGroups(
        const IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &other,
        QSet<quint32> &changed_groups) const
{
    INTRASCALE::addChangedGroups(other, changed_groups);
    ATOMPARAM::addChangedGroups(other, changed_groups);
}

/** Mask these parameters so only those that are in the groups whose
    indicies are 'idxs' are present */
template<class ATOMPARAM, class INTRASCALE>
SIRE_OUTOFLINE_TEMPLATE
IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>
IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>::applyMask(
                                            const QSet<quint32> &idxs) const
{
    return IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE>(
                                ATOMPARAM::applyMask(idxs),
                                INTRASCALE::applyMask(idxs) );
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace detail

} // end of namespace SireMM

/** Serialise to a binary datastream */
template<class INTRASCALE>
QDataStream& operator<<(QDataStream &ds, 
                        const SireMM::detail::IntraScaledParameters<INTRASCALE> &params)
{
    ds << params.sclfactors;
    return ds;
}

/** Extract from a binary datastream */
template<class INTRASCALE>
QDataStream& operator>>(QDataStream &ds, 
                        SireMM::detail::IntraScaledParameters<INTRASCALE> &params)
{
    ds >> params.sclfactors;
    return ds;
}

/** Serialise to a binary datastream */
template<class ATOMPARAM, class INTRASCALE>
QDataStream& operator<<(QDataStream &ds, 
      const SireMM::detail::IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &params)
{
    ds << static_cast<const ATOMPARAM&>(params)
       << static_cast<const INTRASCALE&>(params);
       
    return ds;
}

/** Extract from a binary datastream */
template<class ATOMPARAM, class INTRASCALE>
QDataStream& operator>>(QDataStream &ds, 
      SireMM::detail::IntraScaledAtomicParameters<ATOMPARAM,INTRASCALE> &params)
{
    ds >> static_cast<ATOMPARAM&>(params)
       >> static_cast<INTRASCALE&>(params);
       
    return ds;
}

SIRE_END_HEADER

#endif
