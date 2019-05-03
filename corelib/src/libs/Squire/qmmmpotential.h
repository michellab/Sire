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

#ifndef SQUIRE_QMMMPOTENTIAL_H
#define SQUIRE_QMMMPOTENTIAL_H

#include "SireBase/propertymap.h"

SIRE_BEGIN_HEADER

namespace Squire
{
template<class QM, class MM>
class QMMMPotential;
}

template<class QM, class MM>
QDataStream& operator<<(QDataStream&, const Squire::QMMMPotential<QM,MM>&);

template<class QM, class MM>
QDataStream& operator>>(QDataStream&, Squire::QMMMPotential<QM,MM>&);

namespace SireMol
{
class PartialMolecule;
class MoleculeGroup;
}

namespace Squire
{

using SireBase::PropertyMap;

using SireMol::PartialMolecule;
using SireMol::MoleculeGroup;

/** This is a QM/MM forcefield that calculates the energy
    of QM molecules in a field of MM point charges
    
    @author Christopher Woods
*/
template<class QM, class MM>
class SQUIRE_EXPORT QMMMPotential : protected QM, protected MM
{

friend SQUIRE_EXPORT QDataStream& ::operator<<<>(QDataStream&, const QMMMPotential<QM,MM>&);

friend SQUIRE_EXPORT QDataStream& ::operator>><>(QDataStream&, QMMMPotential<QM,MM>&);

public:
    typedef QM QMPotential;
    typedef MM MMPotential;

    typedef typename QM::Parameter QMParameter;
    typedef typename QM::Parameters QMParameters;
    
    typedef typename MM::Parameter MMParameter;
    typedef typename MM::Parameters MMParameters;
    
    typedef typename QM::Molecule QMMolecule;
    typedef typename QM::Molecules QMMolecules;
    typedef typename QM::ChangedMolecule ChangedQMMolecule;
    
    typedef typename MM::Molecule MMMolecule;
    typedef typename MM::Molecules MMMolecules;
    typedef typename MM::ChangedMolecule ChangedMMMolecule;
    
    QMMMPotential();
    
    QMMMPotential(const QMMMPotential &other);
    
    ~QMMMPotential();

    QMMMPotential<QM,MM>& operator=(const QMMMPotential<QM,MM> &other);

    typename QMMMPotential<QM,MM>::QMParameters 
    getQMParameters(const PartialMolecule &mol,
                    const PropertyMap &map = PropertyMap());
                  
    typename QMMMPotential<QM,MM>::QMParameters
    updateQMParameters(const typename QMMMPotential<QM,MM>::QMParameters &old_params,
                       const PartialMolecule &old_mol,
                       const PartialMolecule &new_mol,
                       const PropertyMap &map = PropertyMap());
                     
    typename QMMMPotential<QM,MM>::QMParameters
    updateQMParameters(const typename QMMMPotential<QM,MM>::QMParameters &old_params,
                       const PartialMolecule &old_molecule,
                       const PartialMolecule &new_molecule,
                       const PropertyMap &old_map, const PropertyMap &new_map);
                     
    typename QMMMPotential<QM,MM>::QMMolecule
    parameteriseQM(const PartialMolecule &molecule,
                   const PropertyMap &map = PropertyMap());
    
    typename QMMMPotential<QM,MM>::QMMolecules 
    parameteriseQM(const MoleculeGroup &molecules,
                   const PropertyMap &map = PropertyMap());

    typename QMMMPotential<QM,MM>::MMParameters 
    getMMParameters(const PartialMolecule &mol,
                    const PropertyMap &map = PropertyMap());
                  
    typename QMMMPotential<QM,MM>::MMParameters
    updateMMParameters(const typename QMMMPotential<QM,MM>::MMParameters &old_params,
                       const PartialMolecule &old_mol,
                       const PartialMolecule &new_mol,
                       const PropertyMap &map = PropertyMap());
                     
    typename QMMMPotential<QM,MM>::MMParameters
    updateMMParameters(const typename QMMMPotential<QM,MM>::MMParameters &old_params,
                       const PartialMolecule &old_molecule,
                       const PartialMolecule &new_molecule,
                       const PropertyMap &old_map, const PropertyMap &new_map);
                     
    typename QMMMPotential<QM,MM>::MMMolecule
    parameteriseMM(const PartialMolecule &molecule,
                   const PropertyMap &map = PropertyMap());
    
    typename QMMMPotential<QM,MM>::MMMolecules 
    parameteriseMM(const MoleculeGroup &molecules,
                   const PropertyMap &map = PropertyMap());

};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor */
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
QMMMPotential<QM,MM>::QMMMPotential() : QM(), MM()
{}

/** Copy constructor */
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
QMMMPotential<QM,MM>::QMMMPotential(const QMMMPotential &other)
                     : QM(other), MM(other)
{}

/** Destructor */
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
QMMMPotential<QM,MM>::~QMMMPotential()
{}

/** Copy assignment operator */
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
QMMMPotential<QM,MM>& QMMMPotential<QM,MM>::operator=(const QMMMPotential<QM,MM> &other)
{
    if (this != &other)
    {
        QM::operator=(other);
        MM::operator=(other);
    }
    
    return *this;
}

/** Parameterise the passed molecule for the QM part of this potential,
    using the passed property map to require the properties that contain
    the required parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::QMParameters 
QMMMPotential<QM,MM>::getQMParameters(const PartialMolecule &mol,
                                      const PropertyMap &map)
{
    return QM::getParameters(mol, map);
}
              
/** Update the QM parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::QMParameters
QMMMPotential<QM,MM>::updateQMParameters(
                    const typename QMMMPotential<QM,MM>::QMParameters &old_params,
                    const PartialMolecule &old_molecule,
                    const PartialMolecule &new_molecule,
                    const PropertyMap &map)
{
    return QM::updateParameters(old_params, old_molecule,
                                new_molecule, map);
}
                 
/** Update the QM parameters for the molecule going from 'old_molecule' to 
    'new_molecule', also while the parameters of 'old_molecule'
    where found in 'old_map', now get the parameters using 'new_map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::QMParameters
QMMMPotential<QM,MM>::updateQMParameters(
                    const typename QMMMPotential<QM,MM>::QMParameters &old_params,
                    const PartialMolecule &old_molecule,
                    const PartialMolecule &new_molecule,
                    const PropertyMap &old_map, 
                    const PropertyMap &new_map)
{
    return QM::updateParameters(old_params, old_molecule,
                                new_molecule, old_map, new_map);
}
                 
/** Return the QMMMPotential::QMMolecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::QMMolecule
QMMMPotential<QM,MM>::parameteriseQM(const PartialMolecule &molecule,
                                     const PropertyMap &map)
{
    return QM::parameterise(molecule, map);
}

/** Convert the passed group of molecules into QMMMPotential::QMMolecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::QMMolecules 
QMMMPotential<QM,MM>::parameteriseQM(const MoleculeGroup &molecules,
                                     const PropertyMap &map)
{
    return QM::parameterise(molecules, map);
}

/** Parameterise the passed molecule for the MM part of this potential,
    using the passed property map to require the properties that contain
    the required parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::MMParameters 
QMMMPotential<QM,MM>::getMMParameters(const PartialMolecule &mol,
                                      const PropertyMap &map)
{
    return MM::getParameters(mol, map);
}
              
/** Update the MM parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::MMParameters
QMMMPotential<QM,MM>::updateMMParameters(
                    const typename QMMMPotential<QM,MM>::MMParameters &old_params,
                    const PartialMolecule &old_molecule,
                    const PartialMolecule &new_molecule,
                    const PropertyMap &map)
{
    return MM::updateParameters(old_params, old_molecule,
                                new_molecule, map);
}
                 
/** Update the MM parameters for the molecule going from 'old_molecule' to 
    'new_molecule', also while the parameters of 'old_molecule'
    where found in 'old_map', now get the parameters using 'new_map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::MMParameters
QMMMPotential<QM,MM>::updateMMParameters(
                    const typename QMMMPotential<QM,MM>::MMParameters &old_params,
                    const PartialMolecule &old_molecule,
                    const PartialMolecule &new_molecule,
                    const PropertyMap &old_map, 
                    const PropertyMap &new_map)
{
    return MM::updateParameters(old_params, old_molecule,
                                new_molecule, old_map, new_map);
}
                 
/** Return the QMMMPotential::MMMolecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::MMMolecule
QMMMPotential<QM,MM>::parameteriseMM(const PartialMolecule &molecule,
                                     const PropertyMap &map)
{
    return MM::parameterise(molecule, map);
}

/** Convert the passed group of molecules into QMMMPotential::MMMolecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
 
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class QM, class MM>
SIRE_OUTOFLINE_TEMPLATE
typename QMMMPotential<QM,MM>::MMMolecules 
QMMMPotential<QM,MM>::parameteriseMM(const MoleculeGroup &molecules,
                                     const PropertyMap &map)
{
    return MM::parameterise(molecules, map);
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace Squire

/** Serialise a QMMMPotential to a binary datastream */
template<class QM, class MM>
QDataStream& operator<<(QDataStream &ds, 
                        const Squire::QMMMPotential<QM,MM> &qmmmpot)
{
    ds << static_cast<const QM&>(qmmmpot)
       << static_cast<const MM&>(qmmmpot);
       
    return ds;
}

/** Extract a QMMMPotential from a binary datastream */
template<class QM, class MM>
QDataStream& operator>>(QDataStream &ds, 
                        Squire::QMMMPotential<QM,MM> &qmmmpot)
{
    ds >> static_cast<QM&>(qmmmpot) 
       >> static_cast<MM&>(qmmmpot);
       
    return ds;
}

SIRE_END_HEADER

#endif
