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

#ifndef SIREMOL_BEADPROPERTY_HPP
#define SIREMOL_BEADPROPERTY_HPP

#include <QVector>

#include "SireBase/qvariant_metatype.h"

#include "moleculeinfodata.h"
#include "molviewproperty.h"
#include "beadidx.h"
#include "beading.h"

#include "SireError/errors.h"

#include "tostring.h"

#include "SireStream/shareddatastream.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class BeadProp;

template<class T>
class BeadProperty;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::BeadProp&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::BeadProp&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireMol::BeadProperty<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireMol::BeadProperty<T>&);

namespace SireMol
{

typedef BeadProperty<QString>  BeadStringProperty;
typedef BeadProperty<qint64>   BeadIntProperty;
typedef BeadProperty<double>   BeadFloatProperty;
typedef BeadProperty<QVariant> BeadVariantProperty;

/** Small class used to provide a common base for all BeadProperty types */
class SIREMOL_EXPORT BeadProp : public MolViewProperty
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const BeadProp&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, BeadProp&);

public:
    BeadProp();
    BeadProp(const Beading &beading);
    BeadProp(const BeadProp &other);
   
    virtual ~BeadProp();
    
    virtual bool canConvert(const QVariant &value) const=0;
    
    virtual void assignFrom(const BeadProperty<QVariant> &values)=0;
    
    virtual BeadProperty<QVariant> toVariant() const=0;
    
    virtual void assertCanConvert(const QVariant &value) const=0;

    const Beading& beading() const;

    void setBeading(const Beading &beading);

protected:
    BeadProp& operator=(const BeadProp &other);
    bool operator==(const BeadProp &other) const;
    bool operator!=(const BeadProp &other) const;

    int getNBeads(const MoleculeInfoData &molinfo) const;

private:
    /** The beading used to create the beads */
    BeadingPtr bdng;
};

/** This is a property that can hold one value for each
    bead in the molecule.
    
    mol.setProperty( "charge", BeadCharges( [....] ) )
    mol.setProperty( "lj", BeadLJs( [....] ) )

    bead.setProperty( "charge", 0.0 * mod_e )
    
    @author Christopher Woods
*/
template<class T>
class SIREMOL_EXPORT BeadProperty 
    : public SireBase::ConcreteProperty<BeadProperty<T>, BeadProp>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const BeadProperty<T>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, BeadProperty<T>&);

public:
    BeadProperty();

    BeadProperty(const MoleculeInfoData &molinfo,
                 const Beading &beading);

    BeadProperty(const QVector<T> &values,
                 const Beading &beading);
    
    BeadProperty(const BeadProperty<T> &other);
    
    ~BeadProperty();
    
    BeadProperty<T>& operator=(const BeadProperty<T> &other);
    
    static const char* typeName();
    
    BeadProperty<T>* clone() const;
    
    bool operator==(const BeadProperty<T> &other) const;
    bool operator!=(const BeadProperty<T> &other) const;

    const T& operator[](const BeadIdx &beadidx) const;
    const T& at(const BeadIdx &beadidx) const;
    const T& get(const BeadIdx &beadidx) const;

    BeadProperty<T>& set(BeadIdx beadidx, const T &value);

    const T* data() const;
    const T* constData() const;

    bool isEmpty() const;

    int size() const;
    int count() const;
    
    int nBeads() const;

    QString toString() const;
    
    const QVector<T>& array() const;

    void assignFrom(const BeadProperty<QVariant> &values);
    
    static BeadProperty<T> fromVariant(const BeadProperty<QVariant> &values);
    
    BeadProperty<QVariant> toVariant() const;
    
    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;
    
    bool canConvert(const QVariant &value) const;
    
    void assertCanConvert(const QVariant &value) const;

private:
    /** The actual bead property values */
    QVector<T> props;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
BeadProperty<T>::BeadProperty()
                : SireBase::ConcreteProperty<BeadProperty<T>,BeadProp>()
{}

/** Construct space for the values of the property for all of the 
    beads in the molecule described by 'molinfo' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
BeadProperty<T>::BeadProperty(const MoleculeInfoData &molinfo,
                              const Beading &beadin)
                : SireBase::ConcreteProperty<BeadProperty<T>,BeadProp>(beadin)
{
    const int nbeads = this->getNBeads(molinfo);

    if (nbeads > 0)
    {
        props = QVector<T>(nbeads);
        props.squeeze();
    }
}

/** Create bead properties from the list of passed values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
BeadProperty<T>::BeadProperty(const QVector<T> &values,
                              const Beading &beading)
                : SireBase::ConcreteProperty<BeadProperty<T>,BeadProp>(beading)
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
void BeadProperty<T>::assertCanConvert(const QVariant &value) const
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
BeadProperty<T>::BeadProperty(const BeadProperty<T> &other)
              : SireBase::ConcreteProperty<BeadProperty<T>,BeadProp>(other),
                props(other.props)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
BeadProperty<T>::~BeadProperty()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
BeadProperty<T>& BeadProperty<T>::operator=(const BeadProperty<T> &other)
{
    BeadProp::operator=(other);
    props = other.props;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool BeadProperty<T>::operator==(const BeadProperty<T> &other) const
{
    return props == other.props and BeadProp::operator==(other);
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool BeadProperty<T>::operator!=(const BeadProperty<T> &other) const
{
    return not BeadProperty<T>::operator==(other);
}

/** Return the property for the bead at index 'beadidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& BeadProperty<T>::operator[](const BeadIdx &beadidx) const
{
    return props.constData()[beadidx.map(props.count())];
}
    
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* BeadProperty<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< BeadProperty<T> >() );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
BeadProperty<T>* BeadProperty<T>::clone() const
{
    return new BeadProperty<T>(*this);
}

/** Return the underlying array holding the contents of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const QVector<T>& BeadProperty<T>::array() const
{
    return props;
}

/** Return a string representation of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString BeadProperty<T>::toString() const
{
    return QString("BeadProperty<%1>( %2, beading() == %3 )")
                .arg( QMetaType::typeName( qMetaTypeId<T>() ) )
                .arg( Sire::toString(this->array()) )
                .arg( this->beading().toString() );
}

/** Return whether or not it is possible to convert the variant
    'value' so that it can be part of this property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool BeadProperty<T>::canConvert(const QVariant &value) const
{
    return value.isNull() or value.canConvert<T>();
}

template<class T>
BeadProperty<T> BeadProperty<T>::fromVariant(const BeadProperty<QVariant> &variant)
{
    BeadProperty<T> array;
    array.setBeading(variant.beading());
    
    array.assignFrom(variant);
    
    return array;
}

/** Assign the values of this property from the array of variants
    in 'values'
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void BeadProperty<T>::assignFrom(const BeadProperty<QVariant> &variant)
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
        BeadProperty<T>::assertCanConvert(value);
        
        if (value.isNull())
            props_array[i] = T();
        else
            props_array[i] = value.value<T>();
    }

    this->setBeading(variant.beading());
}

/** Convert the properties into an array of QVariants */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
BeadProperty<QVariant> BeadProperty<T>::toVariant() const
{
    if (props.isEmpty())
        return BeadProperty<QVariant>();
        
    int nvals = props.count();
    const T *props_array = props.constData();
    
    QVector<QVariant> converted_vals(nvals);
    converted_vals.squeeze();
    QVariant *converted_vals_array = converted_vals.data();

    for (int i=0; i<nvals; ++i)
    {
        converted_vals_array[i].setValue<T>(props_array[i]);
    }
    
    return BeadProperty<QVariant>(converted_vals, this->beading());
}

/** Return the property for the bead at index 'beadidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& BeadProperty<T>::at(const BeadIdx &beadidx) const
{
    return this->operator[](beadidx);
}

/** Return the property for the bead at index 'beadidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& BeadProperty<T>::get(const BeadIdx &beadidx) const
{
    return this->operator[](beadidx);
}

/** Set the value of the property for the bead at 
    index 'beadidx' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
BeadProperty<T>& BeadProperty<T>::set(BeadIdx beadidx, const T &value)
{
    props.data()[beadidx.map(props.count())] = value;
    return *this;
}

/** Return a raw pointer to the array of property values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* BeadProperty<T>::data() const
{
    return props.constData();
}

/** Return a raw pointer to the array of property values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* BeadProperty<T>::constData() const
{
    return props.constData();
}

/** Return whether or not this property is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool BeadProperty<T>::isEmpty() const
{
    return props.count() == 0;
}

/** Return the number of beads */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int BeadProperty<T>::size() const
{
    return props.count();
}

/** Return the number of beads */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int BeadProperty<T>::count() const
{
    return props.count();
}

/** Return the number of beads */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int BeadProperty<T>::nBeads() const
{
    return props.count();
}

/** Is this property compatible with the molecule that is represented
    by 'molinfo' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool BeadProperty<T>::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return this->getNBeads(molinfo) == this->nBeads();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

/** Serialise this property to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMol::BeadProperty<T> &prop)
{
    //serialise the base class - this writes the header and version!
    ds << static_cast<const SireMol::BeadProp&>(prop);
    ds << prop.props;
    
    return ds;
}

/** Extract from an binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMol::BeadProperty<T> &prop)
{
    ds >> static_cast<SireMol::BeadProp&>(prop);
    ds >> prop.props;
        
    return ds;
}

Q_DECLARE_METATYPE( SireMol::BeadStringProperty );
Q_DECLARE_METATYPE( SireMol::BeadIntProperty );
Q_DECLARE_METATYPE( SireMol::BeadFloatProperty );
Q_DECLARE_METATYPE( SireMol::BeadVariantProperty );

SIRE_EXPOSE_CLASS( SireMol::BeadProp )

SIRE_EXPOSE_BEAD_PROPERTY( QString, SireMol::BeadStringProperty )
SIRE_EXPOSE_BEAD_PROPERTY( qint64, SireMol::BeadIntProperty )
SIRE_EXPOSE_BEAD_PROPERTY( double, SireMol::BeadFloatProperty )
SIRE_EXPOSE_BEAD_PROPERTY( QVariant, SireMol::BeadVariantProperty )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMol::BeadProperty<QString>;
template class SireMol::BeadProperty<qint64>;
template class SireMol::BeadProperty<double>;
template class SireMol::BeadProperty<QVariant>;
#endif

SIRE_END_HEADER

#endif
