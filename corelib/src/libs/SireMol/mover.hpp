/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIREMOL_MOVER_HPP
#define SIREMOL_MOVER_HPP

#include "mover.h"
#include "evaluator.h"
#include "atommatcher.h"

#include "SireBase/propertymap.h"

#include "SireMaths/axisset.h"
#include "SireMaths/align.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class Space;
}

namespace SireMol
{

class AtomMatcher;

using SireVol::Space;

/** This manipulator class is used to move the atoms in 
    a molecule view object. 
    
    @author Christopher Woods
*/
template<class T>
class SIREMOL_EXPORT Mover : public SireBase::ConcreteProperty<Mover<T>,T>, public MoverBase
{
public:
    Mover();

    Mover(const T &view);
    Mover(const T &view, const AtomSelection &movable_atoms);
    
    Mover(const Mover<T> &other);
    
    ~Mover();
    
    static const char* typeName();

    Mover<T>* clone() const;
    
    Mover<T>& operator=(const Mover<T> &other);
    Mover<T>& operator=(const T &other);
    
    T commit() const;
    
    QString toString() const;
    
    Mover<T>& mapInto(const AxisSet &axes,
                      const PropertyMap &map = PropertyMap());
    
    Mover<T>& changeFrame(const AxisSet &from_frame,
                          const AxisSet &to_frame,
                          const PropertyMap &map = PropertyMap());
                      
    Mover<T>& translate(const Vector &delta,
                        const PropertyMap &map = PropertyMap());
    
    Mover<T>& rotate(const Quaternion &quat, const Vector &point,
                     const PropertyMap &map = PropertyMap());
                     
    Mover<T>& rotate(const Matrix &rotmat, const Vector &point,
                     const PropertyMap &map = PropertyMap());
    
    Mover<T>& transform(const Transform &transform,
                        const PropertyMap &map = PropertyMap());
    
    Mover<T>& change(const BondID &bond, SireUnits::Dimension::Length delta,
                     const PropertyMap &map = PropertyMap());
                     
    Mover<T>& change(const AngleID &angle, SireUnits::Dimension::Angle delta,
                     const PropertyMap &map = PropertyMap());
                     
    Mover<T>& change(const DihedralID &dihedral, SireUnits::Dimension::Angle delta,
                     const PropertyMap &map = PropertyMap());
    
    Mover<T>& change(const BondID &bond, SireUnits::Dimension::Angle delta,
                     const PropertyMap &map = PropertyMap());
                                                       
    Mover<T>& change(const ImproperID &improper, SireUnits::Dimension::Angle delta,
                     const PropertyMap &map = PropertyMap());
                     
    Mover<T>& set(const BondID &bond, SireUnits::Dimension::Length value,
                  const PropertyMap &map = PropertyMap());
                  
    Mover<T>& set(const AngleID &angle, SireUnits::Dimension::Angle value,
                  const PropertyMap &map = PropertyMap());
                  
    Mover<T>& set(const DihedralID &dihedral, SireUnits::Dimension::Angle value,
                  const PropertyMap &map = PropertyMap());
                  
    Mover<T>& setAll(const DihedralID &dihedral, SireUnits::Dimension::Angle value,
                     const PropertyMap &map = PropertyMap());
                     
    Mover<T>& set(const ImproperID &improper, SireUnits::Dimension::Angle value,
                  const PropertyMap &map = PropertyMap());

    Mover<T>& alignTo(const MoleculeView &other, 
                      const AtomMatcher &matcher,
                      const PropertyMap &map = PropertyMap());
                      
    Mover<T>& alignTo(const MoleculeView &other,
                      const AtomMatcher &matcher,
                      const PropertyMap &map0,
                      const PropertyMap &map1);
                    
    Mover<T>& alignTo(const MoleculeView &other,
                      const AtomSelection &aligning_atoms,
                      const AtomMatcher &matcher,
                      const PropertyMap &map = PropertyMap());
                      
    Mover<T>& alignTo(const MoleculeView &other,
                      const AtomSelection &aligning_atoms,
                      const AtomMatcher &matcher,
                      const PropertyMap &map0,
                      const PropertyMap &map1);

    Mover<T>& align(const MoleculeView &other,
                    const PropertyMap &map = PropertyMap());

    Mover<T>& align(const MoleculeView &other,
                    const PropertyMap &map0,
                    const PropertyMap &map1);

    Mover<T>& align(const MoleculeView &other,
                    const AtomMatcher &matcher,
                    const PropertyMap &map = PropertyMap());

    Mover<T>& align(const MoleculeView &other,
                    const AtomMatcher &matcher,
                    const PropertyMap &map0,
                    const PropertyMap &map1);
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>::Mover() : SireBase::ConcreteProperty<Mover<T>,T>(), MoverBase()
{}

/** Construct a mover that can move all of the atoms
    in the view 'view' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>::Mover(const T &view) 
         : SireBase::ConcreteProperty<Mover<T>,T>(view), MoverBase(view.selection())
{}

/** Construct a mover that can move the 'movable_atoms' of the 
    view 'view' (note that movable_atoms really should be a subset
    of 'view', but that this class does not check that this is the case) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>::Mover(const T &view, const AtomSelection &movable_atoms)
         : SireBase::ConcreteProperty<Mover<T>,T>(view), MoverBase(movable_atoms)
{
    movable_atoms.assertCompatibleWith(view);
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>::Mover(const Mover<T> &other)
         : SireBase::ConcreteProperty<Mover<T>,T>(other), MoverBase(other)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>::~Mover()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::operator=(const Mover<T> &other)
{
    T::operator=(other);
    MoverBase::operator=(other);
    return *this;
}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::operator=(const T &other)
{
    T::operator=(other);
    MoverBase::setMovableAtoms(other.selection());
    
    return *this;
}
    
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* Mover<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Mover<T> >() );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>* Mover<T>::clone() const
{
    return new Mover<T>(*this);
}

/** Commit this Mover - this returns the new copy of
    the Molecule after it has been moved */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T Mover<T>::commit() const
{
    return *this;
}

/** Return a string representation of this mover */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString Mover<T>::toString() const
{
    return QObject::tr( "Mover{ %1 }" ).arg( T::toString() );
}

/** Map the atoms in this view into the axis set described in 'axes',
    using the supplied property map to find the coordinates to be mapped */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::mapInto(const AxisSet &axes, const PropertyMap &map)
{
    MoverBase::mapInto(*(this->d), axes, map);
    return *this;
}

/** Change the coordinate frame of the atoms in this view from the 
    frame 'from_frame' to the frame 'to_frame' using the supplied
    property map to find the coordinates to be mapped */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::changeFrame(const AxisSet &from_frame,
                                const AxisSet &to_frame,
                                const PropertyMap &map)
{
    MoverBase::changeFrame(*(this->d), from_frame, to_frame, map);
    return *this;
}

/** Translate the movable atoms by 'delta', using the supplied
    property map to find the necessary properties */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::translate(const Vector &delta, const PropertyMap &map)
{
    MoverBase::translate(*(this->d), delta, map);
    return *this;
}

/** Rotate the movable atoms using the quaternion 'quat' about the 
    point 'point', using the supplied map to find the necessary 
    properties */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::rotate(const Quaternion &quat, const Vector &point,
                           const PropertyMap &map)
{
    MoverBase::rotate(*(this->d), quat, point, map);
    return *this;
}
                 
/** Rotate the movable atoms using the rotation matrix 'rotmat' about
    the point 'point', using the supplied map to find the necessary
    properties */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::rotate(const Matrix &rotmat, const Vector &point,
                           const PropertyMap &map)
{
    MoverBase::rotate(*(this->d), rotmat, point, map);
    return *this;
}

/** Transform the movable atoms using the transformation 't',
    using the supplied map to find the necessary
    properties */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::transform(const Transform &t, const PropertyMap &map)
{
    MoverBase::transform(*(this->d), t, map);
    return *this;
}

/** Change the bond identified by 'bond' by the amount 'delta',
    by only moving the movable atoms in this view, using the 
    supplied map to find the necessary properties.
    
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/                         
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::change(const BondID &bond, 
                           SireUnits::Dimension::Length delta,
                           const PropertyMap &map)
{
    MoverBase::change(*(this->d), bond, delta, map);
    return *this;
}
                 
/** Change the angle identified by 'angle' by the amount 'delta',
    by only moving the movable atoms in this view, using the 
    supplied map to find the necessary properties
    
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::change(const AngleID &angle, 
                           SireUnits::Dimension::Angle delta,
                           const PropertyMap &map)
{
    MoverBase::change(*(this->d), angle, delta, map);
    return *this;
}
                 
/** Change the dihedral identified by 'dihedral' by the amount 'delta',
    by only moving the movable atoms in this view, using the 
    supplied map to find the necessary properties
    
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::change(const DihedralID &dihedral, 
                           SireUnits::Dimension::Angle delta,
                           const PropertyMap &map)
{
    MoverBase::change(*(this->d), dihedral, delta, map);
    return *this;
}

/** Rotate all of the atoms around the bond identified by 'bond' 
    by the amount 'delta', by only moving the movable atoms in 
    this view, using the supplied map to find the necessary properties
    
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::change(const BondID &bond, 
                           SireUnits::Dimension::Angle delta,
                           const PropertyMap &map)
{
    MoverBase::change(*(this->d), bond, delta, map);
    return *this;
}

/** Change the improper angle by the amount 'delta'

    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::change(const ImproperID &improper, 
                           SireUnits::Dimension::Angle delta,
                           const PropertyMap &map)
{
    MoverBase::change(*(this->d), improper, delta, map);
    return *this;
}
                
/** Set the specified bond to the specified value 

    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::set(const BondID &bond, SireUnits::Dimension::Length value,
                        const PropertyMap &map)
{
    MoverBase::set(*(this->d), bond, value, map);
    return *this;
}
              
/** Set the specified angle to the specified value

    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::set(const AngleID &angle, SireUnits::Dimension::Angle value,
                        const PropertyMap &map)
{
    MoverBase::set(*(this->d), angle, value, map);
    return *this;
}
              
/** Set the specified dihedral to the specified value

    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::set(const DihedralID &dihedral, 
                        SireUnits::Dimension::Angle value,
                        const PropertyMap &map)
{
    MoverBase::set(*(this->d), dihedral, value, map);
    return *this;
}
              
/** Set the specified dihedral to the specified value

    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::setAll(const DihedralID &dihedral, 
                           SireUnits::Dimension::Angle value,
                           const PropertyMap &map)
{
    MoverBase::setAll(*(this->d), dihedral, value, map);
    return *this;
}
                 
/** Set the specified improper to the specified value

    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::set(const ImproperID &improper, SireUnits::Dimension::Angle value,
                        const PropertyMap &map)
{
    MoverBase::set(*(this->d), improper, value, map);
    return *this;
}

/** Align all of the atoms in this view against their equivalent atoms
    in 'other', using 'matcher' to match atoms in this molecule to
    atoms in 'other'. If an atom can't be found, then it is ignored. If
    no atoms can be found then an exception is raised.
    
    Property map 'map0' is used to find the coordinates property 
    in this molecule, and map1 is used to find the coordinates
    property in 'map1'
    
    \throw SireMol::missing_atom
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::alignTo(const MoleculeView &other, const AtomMatcher &matcher,
                            const PropertyMap &map0, const PropertyMap &map1)
{
    MoverBase::mapInto(*(this->d), this->evaluate().alignmentAxes(
                                        other, matcher, map0, map1), map0 );

    return *this;
}

/** Align all of the atoms in this view against their equivalent atoms
    in 'other', using 'matcher' to match atoms in this molecule to
    atoms in 'other'. If an atom can't be found, then it is ignored. If
    no atoms can be found then an exception is raised.
    
    Property map 'map' is used to find the coordinates property 
    in both molecules
    
    \throw SireMol::missing_atom
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::alignTo(const MoleculeView &other, const AtomMatcher &matcher,
                            const PropertyMap &map)
{
    return this->alignTo(other, matcher, map, map);
}
                
/** Align all of the atoms by matching 'aligning_atoms' in this molecule
    against their equivalent atoms in 'other', using 'matcher' to match 
    atoms in this molecule to atoms in 'other'. If an atom can't be found, 
    then it is ignored. If no atoms can be found then an exception is raised.
    
    Property map 'map0' is used to find the coordinates property 
    in this molecule, and map1 is used to find the coordinates
    property in 'map1'
    
    \throw SireMol::missing_atom
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::alignTo(const MoleculeView &other, 
                            const AtomSelection &aligning_atoms,
                            const AtomMatcher &matcher,
                            const PropertyMap &map0, const PropertyMap &map1)
{
    MoverBase::mapInto(*(this->d), Evaluator(*this, aligning_atoms)
                                        .alignmentAxes(other, matcher,
                                                       map0, map1), map0 );

    return *this;
}
                
/** Align all of the atoms by matching 'aligning_atoms' in this molecule
    against their equivalent atoms in 'other', using 'matcher' to match 
    atoms in this molecule to atoms in 'other'. If an atom can't be found, 
    then it is ignored. If no atoms can be found then an exception is raised.
    
    Property map 'map' is used to find the coordinates property 
    in both molecules
    
    \throw SireMol::missing_atom
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::alignTo(const MoleculeView &other, 
                            const AtomSelection &aligning_atoms,
                            const AtomMatcher &matcher,
                            const PropertyMap &map)
{
    return this->alignTo(other, aligning_atoms, matcher, map, map);
}

/** Align this molecule view against 'other' using the supplied AtomMatcher
    to match atoms between the two views, and the supplied property maps to
    find the required properties in each view (this, then other). */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::align(const MoleculeView &other,
                          const AtomMatcher &matcher,
                          const PropertyMap &map0,
                          const PropertyMap &map1)
{
    MoverBase::transform(*(this->d),
                         SireMol::getAlignment(other, map1,
                                               *this, map0, AtomMatchInverter(matcher)), map0);
    
    return *this;
}

/** Align this molecule view against 'other' using the optionally supplied property
    map to find the required properties in both molecules. This matches atoms using their
    AtomIdx, which may not be what you want. If you want more control on matching,
    then supply an AtomMatcher */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::align(const MoleculeView &other,
                          const PropertyMap &map)
{
    MoverBase::transform(*(this->d),
                         SireMol::getAlignment(other, *this, map), map);
    
    return *this;
}

/** Align this molecule view against 'other' using the supplied property
    maps to find the required properties in each molecule respectively. 
    This matches atoms using their
    AtomIdx, which may not be what you want. If you want more control on matching,
    then supply an AtomMatcher */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::align(const MoleculeView &other,
                          const PropertyMap &map0,
                          const PropertyMap &map1)
{
    MoverBase::transform(*(this->d),
                         SireMol::getAlignment(other, map1, *this, map0), map0);
    
    return *this;
}

/** Align this molecule view against 'other' using the supplied AtomMatcher
    to match atoms between the two views, and the supplied property map to
    find the required properties in both views. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover<T>& Mover<T>::align(const MoleculeView &other,
                          const AtomMatcher &matcher,
                          const PropertyMap &map)
{
    MoverBase::transform(*(this->d),
                         SireMol::getAlignment(other, *this, AtomMatchInverter(matcher), map), map);
    
    return *this;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_END_HEADER

#endif
