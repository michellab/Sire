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

#ifndef SIREMOL_CHAINPROPERTY_HPP
#define SIREMOL_CHAINPROPERTY_HPP

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
class ChainProp;

template<class T>
class ChainProperty;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ChainProp&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ChainProp&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireMol::ChainProperty<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireMol::ChainProperty<T>&);

namespace SireMol
{

typedef ChainProperty<QString>  ChainStringProperty;
typedef ChainProperty<qint64>   ChainIntProperty;
typedef ChainProperty<double>   ChainFloatProperty;
typedef ChainProperty<QVariant> ChainVariantProperty;

/** Small class used to provide a common base for all ChainProperty types */
class SIREMOL_EXPORT ChainProp : public MolViewProperty
{
public:
    ChainProp();
    ChainProp(const ChainProp &other);

    virtual ~ChainProp();

    virtual bool canConvert(const QVariant &value) const=0;

    virtual void assignFrom(const ChainProperty<QVariant> &values)=0;

    virtual QVariant getAsVariant(const ChainIdx &chainidx) const=0;
    virtual SireBase::PropertyPtr getAsProperty(const ChainIdx &chainidx) const=0;

    virtual ChainProperty<QVariant> toVariant() const=0;

    virtual void assertCanConvert(const QVariant &value) const=0;
};

/** This is a property that can hold one value for each
    chain in the molecule.

    mol.setProperty( "charge", ChainCharges( [....] ) )
    mol.setProperty( "lj", ChainLJs( [....] ) )

    chain.setProperty( "charge", 0.0 * mod_e )

    @author Christopher Woods
*/
template<class T>
class SIREMOL_EXPORT ChainProperty
    : public SireBase::ConcreteProperty<ChainProperty<T>, ChainProp>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const ChainProperty<T>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, ChainProperty<T>&);

public:
    ChainProperty();

    ChainProperty(const MoleculeInfoData &molinfo);

    ChainProperty(const QVector<T> &values);

    ChainProperty(const ChainProperty<T> &other);

    ~ChainProperty();

    ChainProperty<T>& operator=(const ChainProperty<T> &other);

    static const char* typeName();

    ChainProperty<T>* clone() const;

    bool operator==(const ChainProperty<T> &other) const;
    bool operator!=(const ChainProperty<T> &other) const;

    const T& operator[](const ChainIdx &chainidx) const;
    const T& at(const ChainIdx &chainidx) const;
    const T& get(const ChainIdx &chainidx) const;

    QVariant getAsVariant(const ChainIdx &idx) const;
    SireBase::PropertyPtr getAsProperty(const ChainIdx &idx) const;

    ChainProperty<T>& set(ChainIdx chainidx, const T &value);

    const T* data() const;
    const T* constData() const;

    bool isEmpty() const;

    int size() const;
    int count() const;

    QString toString() const;

    const QVector<T>& array() const;

    int nChains() const;

    void assignFrom(const ChainProperty<QVariant> &values);

    static ChainProperty<T> fromVariant(const ChainProperty<QVariant> &values);

    ChainProperty<QVariant> toVariant() const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    bool canConvert(const QVariant &value) const;

    void assertCanConvert(const QVariant &value) const;

private:
    /** The actual chain property values */
    QVector<T> props;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
ChainProperty<T>::ChainProperty()
              : SireBase::ConcreteProperty<ChainProperty<T>,ChainProp>()
{}

/** Construct space for the values of the property for all of the
    chains in the molecule described by 'molinfo' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
ChainProperty<T>::ChainProperty(const MoleculeInfoData &molinfo)
              : SireBase::ConcreteProperty<ChainProperty<T>,ChainProp>()
{
    if (molinfo.nChains() > 0)
    {
        props = QVector<T>(molinfo.nChains());
        props.squeeze();
    }
}

/** Create chain properties from the list of passed values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
ChainProperty<T>::ChainProperty(const QVector<T> &values)
              : SireBase::ConcreteProperty<ChainProperty<T>,ChainProp>()
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
void ChainProperty<T>::assertCanConvert(const QVariant &value) const
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
ChainProperty<T>::ChainProperty(const ChainProperty<T> &other)
              : SireBase::ConcreteProperty<ChainProperty<T>,ChainProp>(other),
                props(other.props)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
ChainProperty<T>::~ChainProperty()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
ChainProperty<T>& ChainProperty<T>::operator=(const ChainProperty<T> &other)
{
    MolViewProperty::operator=(other);
    props = other.props;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool ChainProperty<T>::operator==(const ChainProperty<T> &other) const
{
    return props == other.props;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool ChainProperty<T>::operator!=(const ChainProperty<T> &other) const
{
    return props != other.props;
}

/** Return the property for the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& ChainProperty<T>::operator[](const ChainIdx &chainidx) const
{
    return props.constData()[chainidx.map(props.count())];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* ChainProperty<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< ChainProperty<T> >() );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
ChainProperty<T>* ChainProperty<T>::clone() const
{
    return new ChainProperty<T>(*this);
}

/** Return the underlying array holding the contents of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const QVector<T>& ChainProperty<T>::array() const
{
    return props;
}

/** Return a string representation of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString ChainProperty<T>::toString() const
{
    return QString("ChainProperty<%1>( %2 )")
                .arg( QMetaType::typeName( qMetaTypeId<T>() ) )
                .arg( Sire::toString(this->array()) );
}

/** Return whether or not it is possible to convert the variant
    'value' so that it can be part of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool ChainProperty<T>::canConvert(const QVariant &value) const
{
    return value.isNull() or value.canConvert<T>();
}

template<class T>
ChainProperty<T> ChainProperty<T>::fromVariant(const ChainProperty<QVariant> &variant)
{
    ChainProperty<T> array;
    array.assignFrom(variant);

    return array;
}

/** Assign the values of this property from the array of variants
    in 'values'

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void ChainProperty<T>::assignFrom(const ChainProperty<QVariant> &variant)
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
        ChainProperty<T>::assertCanConvert(value);

        if (value.isNull())
            props_array[i] = T();
        else
            props_array[i] = value.value<T>();
    }
}

/** Convert the properties into an array of QVariants */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
ChainProperty<QVariant> ChainProperty<T>::toVariant() const
{
    if (props.isEmpty())
        return ChainProperty<QVariant>();

    int nvals = props.count();
    const T *props_array = props.constData();

    QVector<QVariant> converted_vals(nvals);
    converted_vals.squeeze();
    QVariant *converted_vals_array = converted_vals.data();

    for (int i=0; i<nvals; ++i)
    {
        converted_vals_array[i].setValue<T>(props_array[i]);
    }

    return ChainProperty<QVariant>(converted_vals);
}

/** Return the property for the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& ChainProperty<T>::at(const ChainIdx &chainidx) const
{
    return this->operator[](chainidx);
}

/** Return the property for the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& ChainProperty<T>::get(const ChainIdx &chainidx) const
{
    return this->operator[](chainidx);
}

/** Return the value for the passed index, as
    a QVariant. This lets you get the value without knowing the
    actual type of this property

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QVariant ChainProperty<T>::getAsVariant(const ChainIdx &chainidx) const
{
    const T &value = this->get(chainidx);
    return QVariant::fromValue(value);
}

/** Return the value for this index as a
    Property. This lets you get the value without knowing the
    actual type of this property

   \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
SireBase::PropertyPtr ChainProperty<T>::getAsProperty(
                                const ChainIdx &chainidx) const
{
    return SireBase::convert_property(this->get(chainidx));
}

/** Set the value of the property for the chain at
    index 'chainidx' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
ChainProperty<T>& ChainProperty<T>::set(ChainIdx chainidx, const T &value)
{
    props.data()[chainidx.map(props.count())] = value;
    return *this;
}

/** Return a raw pointer to the array of property values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* ChainProperty<T>::data() const
{
    return props.constData();
}

/** Return a raw pointer to the array of property values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* ChainProperty<T>::constData() const
{
    return props.constData();
}

/** Return whether or not this property is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool ChainProperty<T>::isEmpty() const
{
    return props.count() == 0;
}

/** Return the number of chains */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int ChainProperty<T>::size() const
{
    return props.count();
}

/** Return the number of chains */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int ChainProperty<T>::count() const
{
    return props.count();
}

/** Return the number of chains */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int ChainProperty<T>::nChains() const
{
    return props.count();
}

/** Is this property compatible with the molecule that is represented
    by 'molinfo' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool ChainProperty<T>::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return molinfo.nChains() == this->nChains();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

/** Serialise this property to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMol::ChainProperty<T> &prop)
{
    //serialise the base class - this writes the header and version!
    ds << static_cast<const SireMol::ChainProp&>(prop);
    ds << prop.props;

    return ds;
}

/** Extract from an binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMol::ChainProperty<T> &prop)
{
    ds >> static_cast<SireMol::ChainProp&>(prop);
    ds >> prop.props;

    return ds;
}

Q_DECLARE_METATYPE( SireMol::ChainStringProperty );
Q_DECLARE_METATYPE( SireMol::ChainIntProperty );
Q_DECLARE_METATYPE( SireMol::ChainFloatProperty );
Q_DECLARE_METATYPE( SireMol::ChainVariantProperty );

SIRE_EXPOSE_CLASS( SireMol::ChainProp )

SIRE_EXPOSE_CHAIN_PROPERTY( QString, SireMol::ChainStringProperty )
SIRE_EXPOSE_CHAIN_PROPERTY( qint64, SireMol::ChainIntProperty )
SIRE_EXPOSE_CHAIN_PROPERTY( double, SireMol::ChainFloatProperty )
SIRE_EXPOSE_CHAIN_PROPERTY( QVariant, SireMol::ChainVariantProperty )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMol::ChainProperty<QString>;
template class SireMol::ChainProperty<qint64>;
template class SireMol::ChainProperty<double>;
template class SireMol::ChainProperty<QVariant>;
#endif

SIRE_END_HEADER

#endif

