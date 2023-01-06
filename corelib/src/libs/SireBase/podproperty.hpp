#ifndef SIREBASE_PODPROPERTY_HPP
#define SIREBASE_PODPROPERTY_HPP

#include "SireBase/property.h"
#include "SireBase/arrayproperty.hpp"

SIRE_BEGIN_HEADER

namespace SireBase
{
template<class T>
class PODProperty;

template<class T>
class PODArrayProperty;
}

template<class T>
QDataStream& operator<<(QDataStream&, const SireBase::PODProperty<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireBase::PODProperty<T>&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireBase::PODArrayProperty<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireBase::PODArrayProperty<T>&);

namespace SireBase
{

/** This class creates a Property wrapper around any plain old data (POD)
 *  type
 */
template<class T>
class PODProperty :
        public ConcreteProperty<PODProperty<T>, Property>,
        public T
{
public:
    PODProperty();
    PODProperty(const T &value);
    PODProperty(const PODProperty<T> &other);

    ~PODProperty();

    static const char* typeName();
    const char* what() const;

    QString toString() const;

    PODProperty<T>& operator=(const PODProperty<T> &other);

    bool operator==(const PODProperty<T> &other) const;
    bool operator!=(const PODProperty<T> &other) const;

    bool operator==(const T &other) const;
    bool operator!=(const T &other) const;

    PODProperty<T>* clone() const;
};

/** This class creates an ArrayProperty around plain old data (POD) types */
template<class T>
class PODArrayProperty :
    public ConcreteProperty<PODArrayProperty<T>, ArrayProperty<T>>
{
public:
    PODArrayProperty();
    PODArrayProperty(const T &value);
    PODArrayProperty(const QVector<T> &values);
    PODArrayProperty(const QList<T> &values);
    PODArrayProperty(const PODProperty<T> &value);
    PODArrayProperty(const QVector<PODProperty<T>> &values);
    PODArrayProperty(const QList<PODProperty<T>> &values);

    PODArrayProperty(const PODArrayProperty<T> &other);

    ~PODArrayProperty();

    static const char* typeName();
    const char* what() const;

    QString toString() const;

    PODArrayProperty<T>& operator=(const PODArrayProperty<T> &other);

    bool operator==(const PODArrayProperty<T> &other) const;
    bool operator!=(const PODArrayProperty<T> &other) const;

    PODArrayProperty<T>* clone() const;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODProperty<T>::PODProperty()
               : ConcreteProperty<PODProperty<T>, Property>(), T()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODProperty<T>::PODProperty(const T &value)
               : ConcreteProperty<PODProperty<T>, Property>(), T(value)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODProperty<T>::PODProperty(const PODProperty<T> &other)
               : ConcreteProperty<PODProperty<T>, Property>(other),
                 T(static_cast<const T&>(other))
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODProperty<T>::~PODProperty()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* PODProperty<T>::typeName()
{
    return QMetaType::typeName(qMetaTypeId<PODProperty<T>>());
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* PODProperty<T>::what() const
{
    return PODProperty<T>::typeName();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString PODProperty<T>::toString() const
{
    return T::toString();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODProperty<T>& PODProperty<T>::operator=(const PODProperty<T> &other)
{
    T::operator=(static_cast<const T&>(other));
    return *this;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PODProperty<T>::operator==(const PODProperty<T> &other) const
{
    return T::operator==(static_cast<const T&>(other));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PODProperty<T>::operator!=(const PODProperty<T> &other) const
{
    return T::operator!=(static_cast<const T&>(other));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PODProperty<T>::operator==(const T &other) const
{
    return T::operator==(other);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PODProperty<T>::operator!=(const T &other) const
{
    return T::operator!=(other);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODProperty<T>* PODProperty<T>::clone() const
{
    return new PODProperty<T>(*this);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>::PODArrayProperty()
                    : ConcreteProperty<PODArrayProperty<T>, ArrayProperty<T> >()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>::PODArrayProperty(const T &value)
                    : ConcreteProperty<PODArrayProperty<T>, ArrayProperty<T> >(value)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>::PODArrayProperty(const QVector<T> &values)
                    : ConcreteProperty<PODArrayProperty<T>, ArrayProperty<T> >(values)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>::PODArrayProperty(const QList<T> &values)
                    : ConcreteProperty<PODArrayProperty<T>, ArrayProperty<T> >(values)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>::PODArrayProperty(const PODProperty<T> &value)
                    : ConcreteProperty<PODArrayProperty<T>, ArrayProperty<T> >(value)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>::PODArrayProperty(const QVector<PODProperty<T>> &values)
                    : ConcreteProperty<PODArrayProperty<T>, ArrayProperty<T> >()
{
    QVector<T> vals;
    vals.reserve(values.count());

    for (const auto &val : values)
    {
        vals.append(val);
    }

    this->operator=(PODArrayProperty<T>(vals));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>::PODArrayProperty(const QList<PODProperty<T>> &values)
                    : ConcreteProperty<PODArrayProperty<T>, ArrayProperty<T> >()
{
    this->operator=(PODArrayProperty<T>(QVector<PODProperty<T>>(values.toVector())));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>::PODArrayProperty(const PODArrayProperty<T> &other)
                    : ConcreteProperty<PODArrayProperty<T>, ArrayProperty<T> >(other)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>::~PODArrayProperty()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* PODArrayProperty<T>::typeName()
{
    return QMetaType::typeName(qMetaTypeId<PODArrayProperty<T>>());
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* PODArrayProperty<T>::what() const
{
    return PODArrayProperty<T>::typeName();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString PODArrayProperty<T>::toString() const
{
    return ArrayProperty<T>::toString();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>& PODArrayProperty<T>::operator=(const PODArrayProperty<T> &other)
{
    ArrayProperty<T>::operator=(other);
    return *this;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PODArrayProperty<T>::operator==(const PODArrayProperty<T> &other) const
{
    return ArrayProperty<T>::operator==(other);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PODArrayProperty<T>::operator!=(const PODArrayProperty<T> &other) const
{
    return ArrayProperty<T>::operator!=(other);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
PODArrayProperty<T>* PODArrayProperty<T>::clone() const
{
    return new PODArrayProperty<T>(*this);
}

#endif

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds,
                        const SireBase::PODProperty<T> &p)
{
    ds << static_cast<const T&>(p);
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds,
                        SireBase::PODProperty<T> &p)
{
    ds >> static_cast<T&>(p);
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds,
                        const SireBase::PODArrayProperty<T> &p)
{
    ds << p.toVector();
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds,
                        SireBase::PODArrayProperty<T> &p)
{
    QVector<T> vals;
    ds >> vals;

    p = SireBase::PODArrayProperty<T>(vals);

    return ds;
}

#endif

SIRE_END_HEADER

#endif
