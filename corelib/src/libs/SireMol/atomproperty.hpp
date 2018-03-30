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

#ifndef SIREMOL_ATOMPROPERTY_HPP
#define SIREMOL_ATOMPROPERTY_HPP

#include <QVector>

#include "SireBase/qvariant_metatype.h"

#include "molviewproperty.h"
#include "moleculeinfo.h"
#include "moleculeinfodata.h"
#include "atomselection.h"

#include "SireBase/packedarray2d.hpp"
#include "SireBase/quickcopy.hpp"

#include "SireError/errors.h"

#include "tostring.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class AtomProp;

template<class T>
class AtomProperty;
}

QDataStream& operator<<(QDataStream&, const SireMol::AtomProp&);
QDataStream& operator>>(QDataStream&, SireMol::AtomProp&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireMol::AtomProperty<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireMol::AtomProperty<T>&);

namespace SireMol
{

//typedef the basic types
typedef AtomProperty<QString>  AtomStringProperty;
typedef AtomProperty<qint64>   AtomIntProperty;
typedef AtomProperty<double>   AtomFloatProperty;
typedef AtomProperty<QVariant> AtomVariantProperty;

using SireBase::Property;
using SireBase::PropertyPtr;
using SireBase::PackedArray2D;

/** Small class used to give a common base to all
    AtomProperty classes */
class SIREMOL_EXPORT AtomProp : public MolViewProperty
{
public:
    AtomProp();
    AtomProp(const AtomProp &other);
    
    virtual ~AtomProp();
    
    AtomProp& operator=(const AtomVariantProperty &other)
    {
        this->assignFrom(other);
        return *this;
    }
    
    virtual bool canConvert(const QVariant &value) const=0;
    
    virtual void assignFrom(const AtomVariantProperty &values)=0;
    
    virtual AtomVariantProperty toVariant() const=0;
    
    virtual void assertCanConvert(const QVariant &value) const=0;

    virtual PropertyPtr merge(const MoleculeInfoData &molinfo) const=0;
    virtual PropertyPtr divide(const QVector<AtomSelection> &beads) const=0;
    virtual PropertyPtr divideByResidue(const MoleculeInfoData &molinfo) const=0;

protected:
    void throwIncorrectNumberOfAtoms(int nats, int ntotal) const;
    void throwIncorrectNumberOfSelectedAtoms(int nats, int nselected) const;
};

/** This is a property that can hold one value for each
    atom in the molecule. The properties are held in 
    a packed array, with each array corresponding to a CutGroup,
    and the order of the properties stored according 
    to CutGroup
    
    mol.setProperty( "charge", AtomCharges( [....] ) )
    mol.setProperty( "lj", AtomLJs( [ [...], [...] ] ) )

    atom.setProperty( "charge", 0.0 * mod_electron )
    
    @author Christopher Woods
*/
template<class T>
class SIREMOL_EXPORT AtomProperty 
    : public SireBase::ConcreteProperty<AtomProperty<T>, AtomProp>
{

friend QDataStream& ::operator<<<>(QDataStream&, const AtomProperty<T>&);
friend QDataStream& ::operator>><>(QDataStream&, AtomProperty<T>&);

public:
    typedef typename PackedArray2D<T>::Array Array;

    AtomProperty();

    AtomProperty(const MoleculeInfo &molinfo);
    AtomProperty(const MoleculeInfo &molinfo, const T &default_value);
    
    AtomProperty(const MoleculeView &molview);
    AtomProperty(const MoleculeView &molview, const T &default_value);
    
    AtomProperty(const MoleculeInfoData &molinfo);
    AtomProperty(const MoleculeInfoData &molinfo, const T &default_value);
    
    AtomProperty(const T &value);
    AtomProperty(const PackedArray2D<T> &values);
    
    AtomProperty(const AtomProperty<T> &other);
    
    ~AtomProperty();
    
    AtomProperty<T>& operator=(const AtomProperty<T> &other);
    
    static const char* typeName();
    
    AtomProperty<T>* clone() const;
 
    bool operator==(const AtomProperty<T> &other) const;
    bool operator!=(const AtomProperty<T> &other) const;

    AtomProperty<QVariant> toVariant() const;
    static AtomProperty<T> fromVariant(const AtomProperty<QVariant> &variant);
    
    void assignFrom(const AtomProperty<QVariant> &values);

    const typename PackedArray2D<T>::Array& operator[](CGIdx cgidx) const;
    const typename PackedArray2D<T>::Array& at(CGIdx cgidx) const;
    const typename PackedArray2D<T>::Array& get(CGIdx cgidx) const;

    const T& operator[](const CGAtomIdx &cgatomidx) const;
    const T& at(const CGAtomIdx &cgatomidx) const;
    const T& get(const CGAtomIdx &cgatomidx) const;

    AtomProperty<T>& set(const CGAtomIdx &cgatomidx, const T &value);

    AtomProperty<T>& set(CGIdx cgidx, const QVector<T> &values);

    const PackedArray2D<T>& array() const;

    const typename PackedArray2D<T>::Array* data() const;
    const typename PackedArray2D<T>::Array* constData() const;

    const T* data(CGIdx cgidx) const;
    const T* constData(CGIdx cgidx) const;

    QString toString() const;

    bool isEmpty() const;

    int size() const;
    int count() const;
    
    int nCutGroups() const;
    
    int nAtoms() const;
    int nAtoms(CGIdx cgidx) const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;
    bool isCompatibleWith(const MoleculeInfo &molinfo) const;

    AtomProperty<T> matchToSelection(const AtomSelection &selection) const;

    QVector<T> toVector() const;
    QVector<T> toVector(const AtomSelection &selection) const;

    PropertyPtr merge(const MoleculeInfoData &molinfo) const;
    PropertyPtr divide(const QVector<AtomSelection> &beads) const;
    PropertyPtr divideByResidue(const MoleculeInfoData &molinfo) const;
    
    void copyFrom(const QVector<T> &values);
    void copyFrom(const QVector<T> &values, const AtomSelection &selection);

    bool canConvert(const QVariant &value) const;

    void assertCanConvert(const QVariant &value) const;

private:
    /** The actual atomic property values */
    SireBase::PackedArray2D<T> props;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return whether or not the variant 'value' can be converted to be
    held as an AtomProperty<T> */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::canConvert(const QVariant &value) const
{
    return value.isNull() or value.canConvert<T>();
}

/** Assert that the passed variant value can be converted to be
    held within this property
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomProperty<T>::assertCanConvert(const QVariant &value) const
{
    if (not (value.isNull() or value.canConvert<T>()))
    {
        throw SireError::invalid_cast( QObject::tr(
            "It is not possible to convert the value of type %1 to "
            "type %2, as is required for storing in the AtomProperty %3.")
                .arg(value.typeName()).arg( QMetaType::typeName(qMetaTypeId<T>()) )
                .arg(this->what()), CODELOC );
    }
}

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty()
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{}

/** Create an AtomProperty that holds one value for each 
    atom described in 'molinfo'. Each atom starts with
    a default-constructed value of the property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const MoleculeInfoData &molinfo)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{   
    int ncg = molinfo.nCutGroups();

    if (ncg > 0)
    {
        //create space for each CutGroup
        QVector< QVector<T> > tmp_props = QVector< QVector<T> >(ncg);
        QVector<T> *tmp_props_array = tmp_props.data();
        
        for (CGIdx i(0); i<ncg; ++i)
        {
            //now create space for all of the atoms
            tmp_props_array[i] = QVector<T>(molinfo.nAtoms(i));
        }

        //now copy this into the PackedArray
        props = PackedArray2D<T>(tmp_props);
    }
}

/** Create an AtomProperty that holds one value for each 
    atom described in 'molinfo'. Each atom starts with
    a default-constructed value of the property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const MoleculeInfo &molinfo)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{
    this->operator=( AtomProperty<T>(molinfo.data()) );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const MoleculeInfo &molinfo, const T &default_value)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{
    this->operator=( AtomProperty<T>(molinfo.data(), default_value) );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const MoleculeView &molview)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{
    this->operator=( AtomProperty<T>(MoleculeInfo(molview)) );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const MoleculeView &molview, const T &default_value)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{
    this->operator=( AtomProperty<T>(MoleculeInfo(molview), default_value) );
}

/** Create an AtomProperty that holds one value for each 
    atom described in 'molinfo'. Each atom starts with
    the value 'default_value' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const MoleculeInfoData &molinfo,
                              const T &default_value)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{   
    int ncg = molinfo.nCutGroups();

    if (ncg > 0)
    {
        //create space for each CutGroup
        QVector< QVector<T> > tmp_props = QVector< QVector<T> >(ncg);
        QVector<T> *tmp_props_array = tmp_props.data();
        
        for (CGIdx i(0); i<ncg; ++i)
        {
            //now create space for all of the atoms
            tmp_props_array[i] = QVector<T>(molinfo.nAtoms(i), default_value);
        }

        //now copy this into the PackedArray
        props = PackedArray2D<T>(tmp_props);
    }
}

/** Construct an Atom property that holds a single value (only
    suitable for a molecule that has just one atom in just one CutGroup) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const T &value)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{
    QVector<T> tmp_props(1, value);
    props = PackedArray2D<T>(tmp_props);
}

/** Construct the Atom property from the PackedArray2D of values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const PackedArray2D<T> &values)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>(),
                  props(values)
{}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const AtomProperty<T> &other)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>(), 
                  props(other.props)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::~AtomProperty()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>& AtomProperty<T>::operator=(const AtomProperty<T> &other)
{
    props = other.props;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::operator==(const AtomProperty<T> &other) const
{
    return props == other.props;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::operator!=(const AtomProperty<T> &other) const
{
    return props != other.props;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>* AtomProperty<T>::clone() const
{
    return new AtomProperty<T>(*this);
}
    
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* AtomProperty<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< AtomProperty<T> >() );
}

/** Return the array of properties for the atoms in the CutGroup
    identified by index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& AtomProperty<T>::operator[](CGIdx cgidx) const
{
    return props.constData()[ cgidx.map(props.count()) ];
}

/** Convert the contained properties into an array of arrays of QVariants.
    There is one array per CutGroup */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<QVariant> AtomProperty<T>::toVariant() const
{
    return AtomProperty<QVariant>( props.toVariant() );
}

/** Assign the values of the properties from the array of QVariants

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomProperty<T>::assignFrom(const AtomProperty<QVariant> &variant) 
{
    props = SireBase::PackedArray2D<T>::fromVariant(variant.array());
}

/** Return an AtomProperty constructed from an array of QVariants */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T> AtomProperty<T>::fromVariant(const AtomProperty<QVariant> &variant)
{
    return AtomProperty<T>( SireBase::PackedArray2D<T>::fromVariant(variant.array()) );
}

/** Return the array of properties for the atoms in the CutGroup
    identified by index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& AtomProperty<T>::at(CGIdx cgidx) const
{
    return this->operator[](cgidx);
}

/** Return the array of properties for the atoms in the CutGroup
    identified by index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& AtomProperty<T>::get(CGIdx cgidx) const
{
    return this->operator[](cgidx);
}

/** Return the property for the atom at index 'cgatomidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& AtomProperty<T>::operator[](const CGAtomIdx &cgatomidx) const
{
    const typename PackedArray2D<T>::Array &group_props 
                            = this->operator[](cgatomidx.cutGroup());

    return group_props.constData()[ cgatomidx.atom().map(group_props.count()) ];
}

/** Return the property for the atom at index 'cgatomidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& AtomProperty<T>::at(const CGAtomIdx &cgatomidx) const
{
    return this->operator[](cgatomidx);
}

/** Return the property for the atom at index 'cgatomidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& AtomProperty<T>::get(const CGAtomIdx &cgatomidx) const
{
    return this->operator[](cgatomidx);
}

/** Set the value of the property for the atom at index 'cgatomidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>& AtomProperty<T>::set(const CGAtomIdx &cgatomidx, const T &value)
{
    quint32 cgidx = cgatomidx.cutGroup().map(props.count());
    quint32 atomidx = cgatomidx.atom().map( props.at(cgidx).count() );

    props(cgidx, atomidx) = value;
    
    return *this;
}

/** Set the values for all atoms in the CutGroup at index 'cgidx'

    \throw SireError::incompatible_error
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>& AtomProperty<T>::set(CGIdx cgidx, const QVector<T> &values)
{
    props.update( cgidx.map(props.count()), values );
    
    return *this;
}

/** Return a const-reference to the PackedArray2D used to store
    all of the atom properties */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const PackedArray2D<T>& AtomProperty<T>::array() const
{
    return props;
}

/** Return a string representation of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString AtomProperty<T>::toString() const
{
    return QString("AtomProperty<%1>( %2 )")
                .arg( QMetaType::typeName( qMetaTypeId<T>() ) )
                .arg( Sire::toString(this->toVector()) );
}

/** Return a raw pointer to the array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array* AtomProperty<T>::data() const
{
    return props.constData();
}

/** Return a raw pointer to the array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array* AtomProperty<T>::constData() const
{
    return props.constData();
}

/** Return a raw pointer to the array of properties for 
    the atoms in the CutGroup at index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* AtomProperty<T>::data(CGIdx cgidx) const
{
    return this->at(cgidx).constData();
}

/** Return a raw pointer to the array of properties for 
    the atoms in the CutGroup at index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* AtomProperty<T>::constData(CGIdx cgidx) const
{
    return this->at(cgidx).constData();
}

/** Return the number of CutGroups in the molecule */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::size() const
{
    return props.count();
}

/** Return the number of CutGroups in the molecule */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::count() const
{
    return props.count();
}

/** Return the number of CutGroups in the molecule */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::nCutGroups() const
{
    return props.count();
}

/** Return the total number of atoms in the molecule */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::nAtoms() const
{
    return props.nValues();
}

/** Return whether or not this is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::isEmpty() const
{
    return props.isEmpty();
}

/** Return the number of atoms in the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::nAtoms(CGIdx cgidx) const
{
    return this->at(cgidx).count();
}

/** Return whether or not this property is compatible with the
    Molecule whose layout is described in 'molinfo'
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    if (props.nValues() != molinfo.nAtoms())
        return false;

    int ncg = molinfo.nCutGroups();

    if (ncg != props.count())
        return false;
        
    const typename PackedArray2D<T>::Array *props_array = props.constData();
    
    for (CGIdx i(0); i<ncg; ++i)
    {
        if (molinfo.nAtoms(i) != props_array[i].count())
            return false;
    }
    
    return true;
}

/** Return whether or not this property is compatible with the
    Molecule whose layout is described in 'molinfo'
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::isCompatibleWith(const MoleculeInfo &molinfo) const
{
    if (props.nValues() != molinfo.nAtoms())
        return false;

    int ncg = molinfo.nCutGroups();

    if (ncg != props.count())
        return false;
        
    const typename PackedArray2D<T>::Array *props_array = props.constData();
    
    for (CGIdx i(0); i<ncg; ++i)
    {
        if (molinfo.nAtoms(i) != props_array[i].count())
            return false;
    }
    
    return true;
}

/** Convert this atom property to an array of values. The values
    are written in CGAtomIdx order */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QVector<T> AtomProperty<T>::toVector() const
{
    if (this->nAtoms() == 0)
        return QVector<T>();
    
    QVector<T> ret( this->nAtoms() );
    
    SireBase::quickCopy<T>(ret.data(), props.constValueData(), this->nAtoms());

    return ret;
}

/** Convert the properties of the atoms selected in 'selection' to an 
    array of values. The values are written in CGAtomIdx order
    
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QVector<T> AtomProperty<T>::toVector(const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*this);
    
    if (selected_atoms.selectedAll())
        return this->toVector();
    
    else if (selected_atoms.selectedNone())
        return QVector<T>();
        
    else if (selected_atoms.selectedAllCutGroups())
    {
        QVector<T> vals( selected_atoms.nSelected() );
        T *value = vals.data();
        
        const int ncg = selected_atoms.nCutGroups();
        const typename PackedArray2D<T>::Array *props_array = props.constData();
        
        for (CGIdx i(0); i<ncg; ++i)
        {
            const T *group_props = props_array[i].constData();
        
            if (selected_atoms.selectedAll(i))
            {
                const int nats = props_array[i].nValues();
            
                SireBase::quickCopy<T>(value, group_props, nats);
                
                value += nats;
            }
            else
            {
                QList<SireID::Index> idxs = selected_atoms.selectedAtoms(i).toList();
                qSort(idxs);
                
                foreach (SireID::Index idx, idxs)
                {
                    *value = group_props[idx];
                    ++value;
                }
            }
        }
        
        return vals;
    }
    else
    {
        QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();
        qSort(cgidxs);
        
        QVector<T> vals( selected_atoms.nSelected() );
        T *value = vals.data();

        const typename PackedArray2D<T>::Array *props_array = props.constData();

        foreach (CGIdx i, cgidxs)
        {
            const T *group_props = props_array[i].constData();
        
            if (selected_atoms.selectedAll(i))
            {
                const int nats = props_array[i].nValues();
            
                SireBase::quickCopy<T>(value, group_props, nats);
                
                value += nats;
            }
            else
            {
                QList<SireID::Index> idxs = selected_atoms.selectedAtoms(i).toList();
                qSort(idxs);
                
                foreach (SireID::Index idx, idxs)
                {
                    *value = group_props[idx];
                    ++value;
                }
            }
        }
        
        return vals;
    }
}

/** Copy into this atom property set the values from 'values'. The values
    are copied in CGAtomIdx order, and there must be as many values
    as there are atoms
    
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomProperty<T>::copyFrom(const QVector<T> &values)
{
    if (values.count() != this->nAtoms())
        this->throwIncorrectNumberOfAtoms(values.count(), this->nAtoms());

    SireBase::quickCopy<T>(props.valueData(), values.constData(), values.count());
}

/** Copy into this atom property set the values from 'values', but only
    for the atoms selected in 'selection'. This copies the properties
    in in CGAtomIdx order, and there must be the same number of values
    as there are selected atoms
    
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomProperty<T>::copyFrom(const QVector<T> &values, 
                               const AtomSelection &selected_atoms)
{
    selected_atoms.assertCompatibleWith(*this);
    
    if (selected_atoms.selectedAll())
    {
        this->copyFrom(values);
        return;
    }

    if (values.count() != selected_atoms.nSelected())
        this->throwIncorrectNumberOfSelectedAtoms(values.count(),
                                                  selected_atoms.nSelected());
                                                  
    const T *values_array = values.constData();
    
    if (selected_atoms.selectedAllCutGroups())
    {
        const int ncg = selected_atoms.nCutGroups();
    
        for (CGIdx i(0); i<ncg; ++i)
        {
            T *group_props = props.data(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = props.nValues(i);
                SireBase::quickCopy<T>(group_props, values_array, nats);
                
                values_array += nats;
            }
            else
            {
                QList<SireID::Index> idxs = selected_atoms.selectedAtoms(i).toList();
                qSort(idxs);
                
                foreach (SireID::Index idx, idxs)
                {
                    group_props[idx] = *values_array;
                    ++values_array;
                }
            }
        }
    }
    else
    {
        QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();
        qSort(cgidxs);
        
        foreach (CGIdx i, cgidxs)
        {
            T *group_props = props.data(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = props.nValues(i);
                SireBase::quickCopy<T>(group_props, values_array, nats);
                
                values_array += nats;
            }
            else
            {
                QList<SireID::Index> idxs = selected_atoms.selectedAtoms(i).toList();
                qSort(idxs);
                
                foreach (SireID::Index idx, idxs)
                {
                    group_props[idx] = *values_array;
                    ++values_array;
                }
            }
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomSelection::assertCompatibleWith(const AtomProperty<T> &prop) const
{
    prop.assertCompatibleWith(this->info());
}

/** Match this property to the passed selection. This returns 
    the property only for the CutGroups that have been selected,
    and with default values for any atoms in those CutGroups that
    have not been selected. This is useful, e.g. for the forcefield
    classes, as this allows an AtomProperty<T> to be returned
    for only the atoms that are selected as part of the forcefield.
    
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T> AtomProperty<T>::matchToSelection(
                                      const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*this);

    if (selected_atoms.selectedAll())
        return *this;
    else if (selected_atoms.selectedNone())
        return AtomProperty<T>();
    
    else if (selected_atoms.selectedAllCutGroups())
    {
        PackedArray2D<T> new_props = props;
        
        int ncg = props.count();
        
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (not selected_atoms.selectedAll(i))
            {
                int nats = new_props.at(i).count();
                T *atom_props_array = new_props.data(i);

                const QSet<Index> &atoms = selected_atoms.selectedAtoms(i);
                
                for (Index j(0); j<nats; ++j)
                {
                    if (not atoms.contains(j))
                        atom_props_array[j] = T();
                }
            }
        }
        
        return AtomProperty<T>(new_props);
    }
    else
    {
        QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();
        
        QVector< QVector<T> > new_props = QVector< QVector<T> >(cgidxs.count());
        QVector<T> *new_props_array = new_props.data();
        
        const typename PackedArray2D<T>::Array *props_array = props.constData();
        
        int n = 0;
        
        foreach (CGIdx i, cgidxs)
        {
            if (selected_atoms.selectedAll(i))
            {
                new_props_array[n] = props_array[i].toQVector();
                ++n;
            }
            else
            {
                const QSet<Index> &atoms = selected_atoms.selectedAtoms(i);
                
                QVector<T> atom_props = props_array[i].toQVector();
                int nats = atom_props.count();
                
                T *atom_props_array = atom_props.data();
                
                for (Index j(0); j<nats; ++j)
                {
                    if (not atoms.contains(j))
                        atom_props_array[j] = T();
                }
                
                new_props_array[n] = atom_props;
                ++n;
            }
        }
        
        return AtomProperty<T>(new_props);
    }
}

/** Merge all of the atomic properties into a single array, with 
    the properties arranged in AtomIdx order */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropertyPtr AtomProperty<T>::merge(const MoleculeInfoData &moldata) const
{
    this->assertCompatibleWith(moldata);

    QVector<T> vals( moldata.nAtoms() );
    
    T *vals_array = vals.data();
    
    for (AtomIdx i(0); i<moldata.nAtoms(); ++i)
    {
        vals_array[i] = this->at( moldata.cgAtomIdx(i) );
    }
    
    return AtomProperty<T>(vals);
}

/** Divide the AtomProperty into beads according to the passed atom selections,
    and returning the properties in AtomIdx order within each bead
    
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropertyPtr AtomProperty<T>::divide(const QVector<AtomSelection> &beads) const
{
    if (beads.isEmpty())
        return PropertyPtr();

    const int nbeads = beads.count();
    const AtomSelection *beads_array = beads.constData();

    QVector< QVector<T> > bead_vals(nbeads);
    QVector<T> *bead_vals_array = bead_vals.data();
    
    for (int i=0; i<nbeads; ++i)
    {
        const AtomSelection &bead = beads_array[i];
        
        bead.assertCompatibleWith<T>(*this);
        
        QVector<T> vals( bead.nSelected() );
        T *vals_array = vals.data();
        
        if (bead.selectedAll())
        {
            for (AtomIdx j(0); j<bead.nSelected(); ++j)
            {
                vals_array[j] = this->at( bead.info().cgAtomIdx(j) );
            }
        }
        else
        {
            foreach (const AtomIdx &j, bead.selectedAtoms())
            {
                *vals_array = this->at( bead.info().cgAtomIdx(j) );
                ++vals_array;
            }
        }
        
        bead_vals_array[i] = vals;
    }
    
    return AtomProperty<T>(bead_vals);
}

/** Divide the properties into residues. This returns the values in 
    Residue/Index order

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropertyPtr AtomProperty<T>::divideByResidue(const MoleculeInfoData &molinfo) const
{
    this->assertCompatibleWith(molinfo);
    
    QVector< QVector<T> > res_vals( molinfo.nResidues() );
    QVector<T> *res_vals_array = res_vals.data();
    
    for (ResIdx i(0); i<molinfo.nResidues(); ++i)
    {
        const int nats = molinfo.nAtoms(i);
        
        QVector<T> vals(nats);
        T *vals_array = vals.data();
        
        for (int j=0; j<nats; ++j)
        {
            vals_array[j] = this->at( molinfo.cgAtomIdx(molinfo.getAtom(i,j)) );
        }
        
        res_vals_array[i] = vals;
    }
    
    return AtomProperty<T>(res_vals);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

/** Serialise this property to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMol::AtomProperty<T> &prop)
{
    //serialise the base class - this writes the header and version!
    ds << static_cast<const SireMol::AtomProp&>(prop);
    ds << prop.props;
    
    return ds;
}

/** Extract from an binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMol::AtomProperty<T> &prop)
{
    ds >> static_cast<SireMol::AtomProp&>(prop);
    ds >> prop.props;
        
    return ds;
}

///////
/////// Declare the basic types
///////

Q_DECLARE_METATYPE( SireMol::AtomStringProperty );
Q_DECLARE_METATYPE( SireMol::AtomIntProperty );
Q_DECLARE_METATYPE( SireMol::AtomFloatProperty );
Q_DECLARE_METATYPE( SireMol::AtomVariantProperty );

SIRE_EXPOSE_CLASS( SireMol::AtomProp )

SIRE_EXPOSE_ATOM_PROPERTY( QString, SireMol::AtomStringProperty )
SIRE_EXPOSE_ATOM_PROPERTY( qint64, SireMol::AtomIntProperty )
SIRE_EXPOSE_ATOM_PROPERTY( double, SireMol::AtomFloatProperty )
SIRE_EXPOSE_ATOM_PROPERTY( QVariant, SireMol::AtomVariantProperty )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMol::AtomProperty<QString>;
template class SireMol::AtomProperty<qint64>;
template class SireMol::AtomProperty<double>;
template class SireMol::AtomProperty<QVariant>;
#endif

SIRE_END_HEADER

#endif
