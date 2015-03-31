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

#ifndef SIREMOL_CGPROPERTY_HPP
#define SIREMOL_CGPROPERTY_HPP

#include <QVector>

#include "SireBase/qvariant_metatype.h"

#include "moleculeinfodata.h"
#include "molviewproperty.h"

#include "SireError/errors.h"

#include "tostring.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class CGProp;

template<class T>
class CGProperty;
}

QDataStream& operator<<(QDataStream&, const SireMol::CGProp&);
QDataStream& operator>>(QDataStream&, SireMol::CGProp&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireMol::CGProperty<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireMol::CGProperty<T>&);

namespace SireMol
{

////// Typedef the basic types
typedef CGProperty<QString>  CGStringProperty;
typedef CGProperty<qint64>   CGIntProperty;
typedef CGProperty<double>   CGFloatProperty;
typedef CGProperty<QVariant> CGVariantProperty;

/** Small class used to provide a common base for all CGProperty types */
class SIREMOL_EXPORT CGProp : public MolViewProperty
{
public:
    CGProp();
    CGProp(const CGProp &other);
    
    virtual ~CGProp();
    
    virtual bool canConvert(const QVariant &value) const=0;
    
    virtual void assignFrom(const CGProperty<QVariant> &values)=0;
    
    virtual CGProperty<QVariant> toVariant() const=0;
    
    virtual void assertCanConvert(const QVariant &value) const=0;
};

/** This is a property that can hold one value for each
    CutGroup in the molecule.
    
    mol.setProperty( "charge", CGCharges( [....] ) )
    mol.setProperty( "lj", CGLJs( [....] ) )

    cg.setProperty( "charge", 0.0 * mod_e )
    
    @author Christopher Woods
*/
template<class T>
class SIREMOL_EXPORT CGProperty 
    : public SireBase::ConcreteProperty<CGProperty<T>, CGProp>
{

friend QDataStream& ::operator<<<>(QDataStream&, const CGProperty<T>&);
friend QDataStream& ::operator>><>(QDataStream&, CGProperty<T>&);

public:
    CGProperty();

    CGProperty(const MoleculeInfoData &molinfo);
    
    CGProperty(const QVector<T> &values);
    
    CGProperty(const CGProperty<T> &other);
    
    ~CGProperty();
    
    CGProperty<T>& operator=(const CGProperty<T> &other);
    
    static const char* typeName();
    
    CGProperty<T>* clone() const;
    
    bool operator==(const CGProperty<T> &other) const;
    bool operator!=(const CGProperty<T> &other) const;

    const T& operator[](const CGIdx &cgidx) const;
    const T& at(const CGIdx &cgidx) const;
    const T& get(const CGIdx &cgidx) const;

    CGProperty<T>& set(CGIdx cgidx, const T &value);

    const T* data() const;
    const T* constData() const;

    QString toString() const;

    bool isEmpty() const;

    int size() const;
    int count() const;
    
    int nCutGroups() const;

    void assignFrom(const CGProperty<QVariant> &variant);

    const QVector<T>& array() const;
    
    CGProperty<QVariant> toVariant() const;
    
    static CGProperty<T> fromVariant(const CGProperty<QVariant> &variant);
    
    bool canConvert(const QVariant &value) const;
    
    void assertCanConvert(const QVariant &value) const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

private:
    /** The actual CutGroup property values */
    QVector<T> props;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
CGProperty<T>::CGProperty()
              : SireBase::ConcreteProperty<CGProperty<T>,CGProp>()
{}

/** Construct space for the values of the property for all of the 
    CutGroups in the molecule described by 'molinfo' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
CGProperty<T>::CGProperty(const MoleculeInfoData &molinfo)
              : SireBase::ConcreteProperty<CGProperty<T>,CGProp>()
{
    if (molinfo.nCutGroups() > 0)
    {
        props = QVector<T>(molinfo.nCutGroups());
        props.squeeze();
    }
}

/** Create CutGroup properties from the list of passed values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
CGProperty<T>::CGProperty(const QVector<T> &values)
              : SireBase::ConcreteProperty<CGProperty<T>,CGProp>()
{
    props = values;
    props.squeeze();
}

/** Assert that the variant can be converted to a value that can
    be held in this list of properties
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void CGProperty<T>::assertCanConvert(const QVariant &value) const
{
    if (not (value.isNull() or value.canConvert<T>()))
    {
        throw SireError::invalid_cast( QObject::tr(
            "Cannot convert an object of type %1 to an object "
            "of type %2, as required by a %3.")
                .arg(value.typeName()).arg( QMetaType::typeName(qMetaTypeId<T>()) )
                .arg(this->what()), CODELOC );
    }
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
CGProperty<T>::CGProperty(const CGProperty<T> &other)
              : SireBase::ConcreteProperty<CGProperty<T>,CGProp>(other),
                props(other.props)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
CGProperty<T>::~CGProperty()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
CGProperty<T>& CGProperty<T>::operator=(const CGProperty<T> &other)
{
    MolViewProperty::operator=(other);
    props = other.props;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool CGProperty<T>::operator==(const CGProperty<T> &other) const
{
    return props == other.props;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool CGProperty<T>::operator!=(const CGProperty<T> &other) const
{
    return props != other.props;
}

/** Return the property for the CutGroup at index 'cgidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& CGProperty<T>::operator[](const CGIdx &cgidx) const
{
    return props.constData()[cgidx.map(props.count())];
}
    
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* CGProperty<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< CGProperty<T> >() );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
CGProperty<T>* CGProperty<T>::clone() const
{
    return new CGProperty<T>(*this);
}

/** Return the underlying array holding the contents of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const QVector<T>& CGProperty<T>::array() const
{
    return props;
}

/** Return a string representation of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString CGProperty<T>::toString() const
{
    return QString("CGProperty<%1>( %2 )")
                .arg( QMetaType::typeName( qMetaTypeId<T>() ) )
                .arg( Sire::toString(this->array()) );
}

/** Return whether or not it is possible to convert the variant
    'value' so that it can be part of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool CGProperty<T>::canConvert(const QVariant &value) const
{
    return value.isNull() or value.canConvert<T>();
}

template<class T>
CGProperty<T> CGProperty<T>::fromVariant(const CGProperty<QVariant> &variant)
{
    CGProperty<T> array;
    array.assignFrom(variant);
    
    return array;
}

/** Assign the values of this property from the array of variants
    in 'values'
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void CGProperty<T>::assignFrom(const CGProperty<QVariant> &variant)
{
    if (variant.count() == 0)
    {
        props.clear();
        return;
    }
        
    int nvals = variant.count();
    const QVariant *variant_array = variant.constData();
    
    props = QVector<T>(nvals);
    props.squeeze();
    T *props_array = props.data();
    
    for (int i=0; i<nvals; ++i)
    {
        const QVariant &value = variant_array[i];
        CGProperty<T>::assertCanConvert(value);
        
        if (value.isNull())
            props_array[i] = T();
        else
            props_array[i] = value.value<T>();
    }
}

/** Convert the properties into an array of QVariants */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
CGProperty<QVariant> CGProperty<T>::toVariant() const
{
    if (props.isEmpty())
        return CGProperty<QVariant>();
        
    int nvals = props.count();
    const T *props_array = props.constData();
    
    QVector<QVariant> converted_vals(nvals);
    converted_vals.squeeze();
    QVariant *converted_vals_array = converted_vals.data();

    for (int i=0; i<nvals; ++i)
    {
        converted_vals_array[i].setValue<T>(props_array[i]);
    }
    
    return CGProperty<QVariant>(converted_vals);
}

/** Return the property for the CutGroup at index 'cgidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& CGProperty<T>::at(const CGIdx &cgidx) const
{
    return this->operator[](cgidx);
}

/** Return the property for the CutGroup at index 'cgidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& CGProperty<T>::get(const CGIdx &cgidx) const
{
    return this->operator[](cgidx);
}

/** Set the value of the property for the CutGroup at 
    index 'cgidx' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
CGProperty<T>& CGProperty<T>::set(CGIdx cgidx, const T &value)
{
    props.data()[cgidx.map(props.count())] = value;
    return *this;
}

/** Return a raw pointer to the array of property values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* CGProperty<T>::data() const
{
    return props.constData();
}

/** Return a raw pointer to the array of property values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* CGProperty<T>::constData() const
{
    return props.constData();
}

/** Return whether or not this property is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool CGProperty<T>::isEmpty() const
{
    return props.count() == 0;
}

/** Return the number of CutGroups */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int CGProperty<T>::size() const
{
    return props.count();
}

/** Return the number of CutGroups */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int CGProperty<T>::count() const
{
    return props.count();
}

/** Return the number of CutGroups */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int CGProperty<T>::nCutGroups() const
{
    return props.count();
}

/** Is this property compatible with the molecule that is represented
    by 'molinfo' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool CGProperty<T>::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return molinfo.nCutGroups() == this->nCutGroups();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

/** Serialise this property to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMol::CGProperty<T> &prop)
{
    //serialise the base class - this writes the header and version!
    ds << static_cast<const SireMol::CGProp&>(prop);
    ds << prop.props;
    
    return ds;
}

/** Extract from an binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMol::CGProperty<T> &prop)
{
    ds >> static_cast<SireMol::CGProp&>(prop);
    ds >> prop.props;
        
    return ds;
}

Q_DECLARE_METATYPE( SireMol::CGStringProperty );
Q_DECLARE_METATYPE( SireMol::CGIntProperty );
Q_DECLARE_METATYPE( SireMol::CGFloatProperty );
Q_DECLARE_METATYPE( SireMol::CGVariantProperty );

SIRE_EXPOSE_CLASS( SireMol::CGProp )

SIRE_EXPOSE_CUTGROUP_PROPERTY( QString, SireMol::CGStringProperty )
SIRE_EXPOSE_CUTGROUP_PROPERTY( qint64, SireMol::CGIntProperty )
SIRE_EXPOSE_CUTGROUP_PROPERTY( double, SireMol::CGFloatProperty )
SIRE_EXPOSE_CUTGROUP_PROPERTY( QVariant, SireMol::CGVariantProperty )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMol::CGProperty<QString>;
template class SireMol::CGProperty<qint64>;
template class SireMol::CGProperty<double>;
template class SireMol::CGProperty<QVariant>;
#endif

SIRE_END_HEADER

#endif
