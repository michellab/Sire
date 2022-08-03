
#include "SireBase/generalunitproperty.h"

using namespace SireBase;
using namespace SireUnits::Dimension;

namespace SireBase
{
    PropertyPtr SIREBASE_EXPORT wrap(const GeneralUnit &unit)
    {
        return PropertyPtr( GeneralUnitProperty(unit) );
    }

    PropertyPtr SIREBASE_EXPORT wrap(const QVector<GeneralUnit> &units)
    {
        return PropertyPtr( GeneralUnitArrayProperty(units) );
    }

    PropertyPtr SIREBASE_EXPORT wrap(const QList<GeneralUnit> &units)
    {
        return PropertyPtr( GeneralUnitArrayProperty(units) );
    }
}

SIREBASE_EXPORT QDataStream& operator<<(QDataStream &ds, const SireBase::GeneralUnitProperty &p)
{
    ds << static_cast<const GeneralUnit&>(p);
    return ds;
}

SIREBASE_EXPORT QDataStream& operator>>(QDataStream &ds, SireBase::GeneralUnitProperty &p)
{
    ds >> static_cast<GeneralUnit&>(p);
    return ds;
}

SIREBASE_EXPORT QDataStream& operator<<(QDataStream &ds, const SireBase::GeneralUnitArrayProperty &p)
{
    ds << p.toVector();
    return ds;
}

SIREBASE_EXPORT QDataStream& operator>>(QDataStream &ds, SireBase::GeneralUnitArrayProperty &p)
{
    QVector<GeneralUnit> vals;
    ds >> vals;

    p = GeneralUnitArrayProperty(vals);

    return ds;
}

GeneralUnitProperty::GeneralUnitProperty()
                    : ConcreteProperty<GeneralUnitProperty, Property>(),
                      GeneralUnit()
{}

GeneralUnitProperty::GeneralUnitProperty(const SireUnits::Dimension::GeneralUnit &unit)
                    : ConcreteProperty<GeneralUnitProperty, Property>(),
                      GeneralUnit(unit)
{}

GeneralUnitProperty::GeneralUnitProperty(const GeneralUnitProperty &other)
                    : ConcreteProperty<GeneralUnitProperty, Property>(other),
                      GeneralUnit(static_cast<const GeneralUnit&>(other))
{}

GeneralUnitProperty::~GeneralUnitProperty()
{}

const char* GeneralUnitProperty::typeName()
{
    return "something";
}

const char* GeneralUnitProperty::what() const
{
    return GeneralUnitProperty::typeName();
}

QString GeneralUnitProperty::toString() const
{
    return GeneralUnit::toString();
}

GeneralUnitProperty* GeneralUnitProperty::clone() const
{
    return new GeneralUnitProperty(*this);
}

GeneralUnitProperty& GeneralUnitProperty::operator=(const GeneralUnitProperty &other)
{
    GeneralUnit::operator=(static_cast<const GeneralUnit&>(other));
    return *this;
}

bool GeneralUnitProperty::operator==(const GeneralUnitProperty &other) const
{
    return GeneralUnit::operator==(static_cast<const GeneralUnit&>(other));
}

bool GeneralUnitProperty::operator!=(const GeneralUnitProperty &other) const
{
    return GeneralUnit::operator!=(static_cast<const GeneralUnit&>(other));
}

GeneralUnitArrayProperty::GeneralUnitArrayProperty()
                         : ConcreteProperty<GeneralUnitArrayProperty, ArrayProperty<GeneralUnit>>()
{}

GeneralUnitArrayProperty::GeneralUnitArrayProperty(const QVector<SireUnits::Dimension::GeneralUnit> &units)
                         : ConcreteProperty<GeneralUnitArrayProperty, ArrayProperty<GeneralUnit>>(units)
{}

GeneralUnitArrayProperty::GeneralUnitArrayProperty(const QList<SireUnits::Dimension::GeneralUnit> &units)
                         : ConcreteProperty<GeneralUnitArrayProperty, ArrayProperty<GeneralUnit>>(units)
{}

GeneralUnitArrayProperty::GeneralUnitArrayProperty(const GeneralUnitArrayProperty &other)
                         : ConcreteProperty<GeneralUnitArrayProperty, ArrayProperty<GeneralUnit>>(other)
{}

GeneralUnitArrayProperty::~GeneralUnitArrayProperty()
{}

const char* GeneralUnitArrayProperty::typeName()
{
    return "something";
}

const char* GeneralUnitArrayProperty::what() const
{
    return GeneralUnitArrayProperty::typeName();
}

QString GeneralUnitArrayProperty::toString() const
{
    return ArrayProperty<GeneralUnit>::toString();
}

GeneralUnitArrayProperty* GeneralUnitArrayProperty::clone() const
{
    return new GeneralUnitArrayProperty(*this);
}

GeneralUnitArrayProperty& GeneralUnitArrayProperty::operator=(const GeneralUnitArrayProperty &other)
{
    ArrayProperty<GeneralUnit>::operator=(other);
    return *this;
}

bool GeneralUnitArrayProperty::operator==(const GeneralUnitArrayProperty &other) const
{
    return ArrayProperty<GeneralUnit>::operator==(other);
}

bool GeneralUnitArrayProperty::operator!=(const GeneralUnitArrayProperty &other) const
{
    return ArrayProperty<GeneralUnit>::operator!=(other);
}
