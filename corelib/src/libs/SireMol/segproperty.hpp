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

#ifndef SIREMOL_SEGPROPERTY_HPP
#define SIREMOL_SEGPROPERTY_HPP

#include <QVector>

#include "SireBase/qvariant_metatype.h"
#include "SireBase/convert_property.hpp"

#include "moleculeinfodata.h"
#include "molviewproperty.h"

#include "SireError/errors.h"

#include "tostring.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class SegProp;

template<class T>
class SegProperty;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::SegProp&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::SegProp&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireMol::SegProperty<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireMol::SegProperty<T>&);

namespace SireMol
{

typedef SegProperty<QString>  SegStringProperty;
typedef SegProperty<qint64>   SegIntProperty;
typedef SegProperty<double>   SegFloatProperty;
typedef SegProperty<QVariant> SegVariantProperty;

/** Small class used to provide a common base for all SegProperty types */
class SIREMOL_EXPORT SegProp : public MolViewProperty
{
public:
    SegProp();
    SegProp(const SegProp &other);

    virtual ~SegProp();

    virtual bool canConvert(const QVariant &value) const=0;

    virtual void assignFrom(const SegProperty<QVariant> &values)=0;

    virtual QVariant getAsVariant(const SegIdx &segidx) const=0;
    virtual SireBase::PropertyPtr getAsProperty(const SegIdx &segidx) const=0;

    virtual SegProperty<QVariant> toVariant() const=0;

    virtual void assertCanConvert(const QVariant &value) const=0;
};

/** This is a property that can hold one value for each
    segment in the molecule.

    mol.setProperty( "charge", SegCharges( [....] ) )
    mol.setProperty( "lj", SegLJs( [....] ) )

    seg.setProperty( "charge", 0.0 * mod_e )

    @author Christopher Woods
*/
template<class T>
class SIREMOL_EXPORT SegProperty
    : public SireBase::ConcreteProperty<SegProperty<T>, SegProp>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const SegProperty<T>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, SegProperty<T>&);

public:
    SegProperty();

    SegProperty(const MoleculeInfoData &molinfo);

    SegProperty(const QVector<T> &values);

    SegProperty(const SegProperty<T> &other);

    ~SegProperty();

    SegProperty<T>& operator=(const SegProperty<T> &other);

    static const char* typeName();

    SegProperty<T>* clone() const;

    bool operator==(const SegProperty<T> &other) const;
    bool operator!=(const SegProperty<T> &other) const;

    const T& operator[](const SegIdx &segidx) const;
    const T& at(const SegIdx &segidx) const;
    const T& get(const SegIdx &segidx) const;

    QVariant getAsVariant(const SegIdx &idx) const;
    SireBase::PropertyPtr getAsProperty(const SegIdx &idx) const;

    SegProperty<T>& set(SegIdx segidx, const T &value);

    const T* data() const;
    const T* constData() const;

    bool isEmpty() const;

    int size() const;
    int count() const;

    int nSegments() const;

    QString toString() const;

    const QVector<T>& array() const;

    void assignFrom(const SegProperty<QVariant> &values);

    static SegProperty<T> fromVariant(const SegProperty<QVariant> &values);

    SegProperty<QVariant> toVariant() const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    bool canConvert(const QVariant &value) const;

    void assertCanConvert(const QVariant &value) const;

private:
    /** The actual segment property values */
    QVector<T> props;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
SegProperty<T>::SegProperty()
              : SireBase::ConcreteProperty<SegProperty<T>,SegProp>()
{}

/** Construct space for the values of the property for all of the
    segments in the molecule described by 'molinfo' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
SegProperty<T>::SegProperty(const MoleculeInfoData &molinfo)
              : SireBase::ConcreteProperty<SegProperty<T>,SegProp>()
{
    if (molinfo.nSegments() > 0)
    {
        props = QVector<T>(molinfo.nSegments());
        props.squeeze();
    }
}

/** Create segment properties from the list of passed values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
SegProperty<T>::SegProperty(const QVector<T> &values)
              : SireBase::ConcreteProperty<SegProperty<T>,SegProp>()
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
void SegProperty<T>::assertCanConvert(const QVariant &value) const
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
SegProperty<T>::SegProperty(const SegProperty<T> &other)
              : SireBase::ConcreteProperty<SegProperty<T>,SegProp>(other),
                props(other.props)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
SegProperty<T>::~SegProperty()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
SegProperty<T>& SegProperty<T>::operator=(const SegProperty<T> &other)
{
    MolViewProperty::operator=(other);
    props = other.props;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SegProperty<T>::operator==(const SegProperty<T> &other) const
{
    return props == other.props;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SegProperty<T>::operator!=(const SegProperty<T> &other) const
{
    return props != other.props;
}

/** Return the property for the segment at index 'segidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& SegProperty<T>::operator[](const SegIdx &segidx) const
{
    return props.constData()[segidx.map(props.count())];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* SegProperty<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< SegProperty<T> >() );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SegProperty<T>* SegProperty<T>::clone() const
{
    return new SegProperty<T>(*this);
}

/** Return the underlying array holding the contents of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const QVector<T>& SegProperty<T>::array() const
{
    return props;
}

/** Return a string representation of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString SegProperty<T>::toString() const
{
    return QString("SegProperty<%1>( %2 )")
                .arg( QMetaType::typeName( qMetaTypeId<T>() ) )
                .arg( Sire::toString(this->array()) );
}

/** Return whether or not it is possible to convert the variant
    'value' so that it can be part of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SegProperty<T>::canConvert(const QVariant &value) const
{
    return value.isNull() or value.canConvert<T>();
}

template<class T>
SegProperty<T> SegProperty<T>::fromVariant(const SegProperty<QVariant> &variant)
{
    SegProperty<T> array;
    array.assignFrom(variant);

    return array;
}

/** Assign the values of this property from the array of variants
    in 'values'

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SegProperty<T>::assignFrom(const SegProperty<QVariant> &variant)
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
        SegProperty<T>::assertCanConvert(value);

        if (value.isNull())
            props_array[i] = T();
        else
            props_array[i] = value.value<T>();
    }
}

/** Convert the properties into an array of QVariants */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
SegProperty<QVariant> SegProperty<T>::toVariant() const
{
    if (props.isEmpty())
        return SegProperty<QVariant>();

    int nvals = props.count();
    const T *props_array = props.constData();

    QVector<QVariant> converted_vals(nvals);
    converted_vals.squeeze();
    QVariant *converted_vals_array = converted_vals.data();

    for (int i=0; i<nvals; ++i)
    {
        converted_vals_array[i].setValue<T>(props_array[i]);
    }

    return SegProperty<QVariant>(converted_vals);
}

/** Return the property for the segment at index 'segidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& SegProperty<T>::at(const SegIdx &segidx) const
{
    return this->operator[](segidx);
}

/** Return the property for the segment at index 'segidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& SegProperty<T>::get(const SegIdx &segidx) const
{
    return this->operator[](segidx);
}

/** Return the value for the passed index, as
    a QVariant. This lets you get the value without knowing the
    actual type of this property

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QVariant SegProperty<T>::getAsVariant(const SegIdx &segidx) const
{
    const T &value = this->get(segidx);
    return QVariant::fromValue(value);
}

/** Return the value for this index as a
    Property. This lets you get the value without knowing the
    actual type of this property

   \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
SireBase::PropertyPtr SegProperty<T>::getAsProperty(
                                const SegIdx &segidx) const
{
    return SireBase::convert_property(this->get(segidx));
}

/** Set the value of the property for the segment at
    index 'segidx' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
SegProperty<T>& SegProperty<T>::set(SegIdx segidx, const T &value)
{
    props.data()[segidx.map(props.count())] = value;
    return *this;
}

/** Return a raw pointer to the array of property values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* SegProperty<T>::data() const
{
    return props.constData();
}

/** Return a raw pointer to the array of property values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* SegProperty<T>::constData() const
{
    return props.constData();
}

/** Return whether or not this property is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SegProperty<T>::isEmpty() const
{
    return props.count() == 0;
}

/** Return the number of segments */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SegProperty<T>::size() const
{
    return props.count();
}

/** Return the number of segments */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SegProperty<T>::count() const
{
    return props.count();
}

/** Return the number of segments */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SegProperty<T>::nSegments() const
{
    return props.count();
}

/** Is this property compatible with the molecule that is represented
    by 'molinfo' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SegProperty<T>::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return molinfo.nSegments() == this->nSegments();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

/** Serialise this property to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMol::SegProperty<T> &prop)
{
    //serialise the base class - this writes the header and version!
    ds << static_cast<const SireMol::SegProp&>(prop);
    ds << prop.props;

    return ds;
}

/** Extract from an binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMol::SegProperty<T> &prop)
{
    ds >> static_cast<SireMol::SegProp&>(prop);
    ds >> prop.props;

    return ds;
}

Q_DECLARE_METATYPE( SireMol::SegStringProperty );
Q_DECLARE_METATYPE( SireMol::SegIntProperty );
Q_DECLARE_METATYPE( SireMol::SegFloatProperty );
Q_DECLARE_METATYPE( SireMol::SegVariantProperty );

SIRE_EXPOSE_CLASS( SireMol::SegProp )

SIRE_EXPOSE_SEGMENT_PROPERTY( QString, SireMol::SegStringProperty )
SIRE_EXPOSE_SEGMENT_PROPERTY( qint64, SireMol::SegIntProperty )
SIRE_EXPOSE_SEGMENT_PROPERTY( double, SireMol::SegFloatProperty )
SIRE_EXPOSE_SEGMENT_PROPERTY( QVariant, SireMol::SegVariantProperty )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMol::SegProperty<QString>;
template class SireMol::SegProperty<qint64>;
template class SireMol::SegProperty<double>;
template class SireMol::SegProperty<QVariant>;
#endif

SIRE_END_HEADER

#endif

