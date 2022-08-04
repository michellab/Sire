
#include "SireBase/generalunitproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireUnits::Dimension;

using namespace SireStream;

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

static const RegisterMetaType<GeneralUnitProperty> r_genprop;
static const RegisterMetaType<GeneralUnitArrayProperty> r_genarrayprop;

SIREBASE_EXPORT QDataStream& operator<<(QDataStream &ds, const SireBase::GeneralUnitProperty &p)
{
    writeHeader(ds, r_genprop, 1);
    ds << static_cast<const GeneralUnit&>(p);
    return ds;
}

SIREBASE_EXPORT QDataStream& operator>>(QDataStream &ds, SireBase::GeneralUnitProperty &p)
{
    VersionID v = readHeader(ds, r_genprop);

    if (v == 1)
        ds >> static_cast<GeneralUnit&>(p);
    else
        throw version_error(v, "1", r_genprop, CODELOC);

    return ds;
}

SIREBASE_EXPORT QDataStream& operator<<(QDataStream &ds, const SireBase::GeneralUnitArrayProperty &p)
{
    writeHeader(ds, r_genarrayprop, 1);

    SharedDataStream sds(ds);
    sds << p.toVector();

    return ds;
}

SIREBASE_EXPORT QDataStream& operator>>(QDataStream &ds, SireBase::GeneralUnitArrayProperty &p)
{
    VersionID v = readHeader(ds, r_genarrayprop);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        QVector<GeneralUnit> vals;
        sds >> vals;

        p = GeneralUnitArrayProperty(vals);
    }
    else
        throw version_error(v, "1", r_genarrayprop, CODELOC);

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
    return QMetaType::typeName(qMetaTypeId<GeneralUnitProperty>());
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
    return QMetaType::typeName(qMetaTypeId<GeneralUnitArrayProperty>());
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
