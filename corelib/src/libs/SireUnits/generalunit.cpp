
#include <QStringList>
#include <QMutex>
#include <QHash>

#include "generalunit.h"

#include "SireError/errors.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/temperature.h"

#include <QDebug>

using namespace SireUnits;
using namespace SireUnits::Dimension;

QDataStream& operator<<(QDataStream &ds, const SireUnits::Dimension::GeneralUnit &u)
{
    qint8 a = u.ANGLE();
    qint8 c = u.CHARGE();
    qint8 l = u.LENGTH();
    qint8 m = u.MASS();
    qint8 t1 = u.TEMPERATURE();
    qint8 t2 = u.TIME();
    qint8 q = u.QUANTITY();

    ds << a << c << l << m << t1 << t2 << q
       << static_cast<const Unit&>(u);

    return ds;
}

QDataStream& operator>>(QDataStream &ds, SireUnits::Dimension::GeneralUnit &u)
{
    qint8 a, c, l, m, t1, t2, q;

    ds >> a >> c >> l >> m >> t1 >> t2 >> q
       >> static_cast<Unit&>(u);

    u.Angle = a;
    u.Charge = c;
    u.Length = l;
    u.Mass = m;
    u.temperature = t1;
    u.Time = t2;
    u.Quantity = q;

    return ds;
}

namespace SireUnits
{
namespace Dimension
{
namespace detail
{

static QHash<QString,QString> typename_registry;
static QMutex registry_mutex;

static QString getKey(const GeneralUnit &unit)
{
    return QString("%1-%2-%3-%4-%5-%6-%7")
              .arg(unit.MASS()).arg(unit.LENGTH()).arg(unit.TIME())
              .arg(unit.CHARGE()).arg(unit.TEMPERATURE())
              .arg(unit.QUANTITY()).arg(unit.ANGLE());
}

void SIREUNITS_EXPORT registerTypeName(const GeneralUnit &unit, const char *typnam)
{
    QString key = getKey(unit);

    QMutexLocker lkr(&registry_mutex);
    if (not typename_registry.contains(key))
    {
        typename_registry.insert( key, QString(typnam) );
    }
}

static QString getTypeName(const GeneralUnit &unit)
{
    QString key = getKey(unit);

    QMutexLocker lkr(&registry_mutex);

    if (typename_registry.contains(key))
        return typename_registry.value(key);
    else
        return QString("SireUnits::Dimension::PhysUnit<%1,%2,%3,%4,%5,%6,%7>")
                 .arg(unit.MASS()).arg(unit.LENGTH()).arg(unit.TIME())
                 .arg(unit.CHARGE()).arg(unit.TEMPERATURE())
                 .arg(unit.QUANTITY()).arg(unit.ANGLE());
}

} // end of namespace detail
} // end of namespace Dimension
} // end of namespace SireUnits

GeneralUnit::GeneralUnit() : Unit(0)
{
    Mass = 0;
    Length = 0;
    Time = 0;
    Charge = 0;
    temperature = 0;
    Quantity = 0;
    Angle = 0;
}

GeneralUnit::GeneralUnit(double value) : Unit(value)
{
    Mass = 0;
    Length = 0;
    Time = 0;
    Charge = 0;
    temperature = 0;
    Quantity = 0;
    Angle = 0;
}

GeneralUnit::GeneralUnit(const GeneralUnit &other) : Unit(other)
{
    Mass = other.Mass;
    Length = other.Length;
    Time = other.Time;
    Charge = other.Charge;
    temperature = other.temperature;
    Quantity = other.Quantity;
    Angle = other.Angle;
}

GeneralUnit::~GeneralUnit()
{}

const char* GeneralUnit::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GeneralUnit>() );
}

const char* GeneralUnit::what() const
{
    return GeneralUnit::typeName();
}

/** Return the C++ type that this particular GeneralUnit corresponds to */
QString GeneralUnit::_to_cpp_type() const
{
    return detail::getTypeName(*this);
}

void GeneralUnit::assertCompatible(const GeneralUnit &other) const
{
    if (this->isZero() or other.isZero())
        return;

    if (Mass != other.Mass or
        Length != other.Length or
        Time != other.Time or
        Charge != other.Charge or
        temperature != other.temperature or
        Quantity != other.Quantity or
        Angle != other.Angle)
    {
        throw SireError::incompatible_error( QObject::tr(
            "Units for values %1 and %2 are incompatible.")
                .arg(this->toString()).arg(other.toString()));
    }
}

QString GeneralUnit::toString() const
{
    return SireUnits::Dimension::getUnitString(value(), Mass, Length,
                                               Time, Charge,
                                               temperature, Quantity,
                                               Angle);
}

double GeneralUnit::to(const GeneralUnit &units) const
{
    assertCompatible(units);
    return units.convertFromInternal(value());
}

double GeneralUnit::to(const TempBase &other) const
{
    //this must be a temperature!
    GeneralUnit general_temp;
    general_temp.temperature = 1;
    general_temp.setScale(other);

    this->assertCompatible(general_temp);

    return other.convertFromInternal(value()) / other.convertFromInternal();
}

bool GeneralUnit::isDimensionless() const
{
    return Mass == 0 and
           Length == 0 and
           Time == 0 and
           Charge == 0 and
           temperature == 0 and
           Quantity == 0 and
           Angle == 0;
}

bool GeneralUnit::isZero() const
{
    return (std::abs(this->value()) < std::numeric_limits<double>::epsilon());
}

int GeneralUnit::MASS() const
{
    return Mass;
}

int GeneralUnit::LENGTH() const
{
    return Length;
}

int GeneralUnit::TIME() const
{
    return Time;
}

int GeneralUnit::CHARGE() const
{
    return Charge;
}

int GeneralUnit::TEMPERATURE() const
{
    return temperature;
}

int GeneralUnit::QUANTITY() const
{
    return Quantity;
}

int GeneralUnit::ANGLE() const
{
    return Angle;
}

/** Return the units of this value */
GeneralUnit GeneralUnit::units() const
{
    GeneralUnit ret(*this);
    ret.setScale(1.0);
    return ret;
}

bool GeneralUnit::hasSameUnits(const GeneralUnit &other) const
{
    return Mass == other.Mass and
           Length == other.Length and
           Time == other.Time and
           Charge == other.Charge and
           temperature == other.temperature and
           Quantity == other.Quantity and
           Angle == other.Angle;
}

GeneralUnit& GeneralUnit::operator=(const GeneralUnit &other)
{
    setScale(other.value());

    Mass = other.MASS();
    Length = other.LENGTH();
    Time = other.TIME();
    Charge = other.CHARGE();
    temperature = other.TEMPERATURE();
    Quantity = other.QUANTITY();
    Angle = other.ANGLE();

    return *this;
}

bool GeneralUnit::operator==(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() == other.value();
}

bool GeneralUnit::operator!=(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() != other.value();
}

bool GeneralUnit::operator>(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() > other.value();
}

bool GeneralUnit::operator>=(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() >= other.value();
}

bool GeneralUnit::operator<(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() < other.value();
}

bool GeneralUnit::operator<=(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() <= other.value();
}

GeneralUnit GeneralUnit::operator-() const
{
    GeneralUnit ret = *this;
    ret.setScale( -value() );
    return ret;
}

GeneralUnit& GeneralUnit::operator+=(const GeneralUnit &other)
{
    if (this->isZero())
    {
        this->operator=(other);
        return *this;
    }

    assertCompatible(other);
    setScale(value() + other.value());

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit& GeneralUnit::operator-=(const GeneralUnit &other)
{
    if (this->isZero())
    {
        this->operator=(-other);
        return *this;
    }

    assertCompatible(other);
    setScale(value() - other.value());

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit GeneralUnit::operator+(const GeneralUnit &other) const
{
    GeneralUnit ret = *this;
    ret += other;
    return ret;
}

GeneralUnit GeneralUnit::operator-(const GeneralUnit &other) const
{
    GeneralUnit ret = *this;
    ret -= other;
    return ret;
}

GeneralUnit& GeneralUnit::operator+=(double val)
{
    assertCompatible(GeneralUnit(val));
    setScale(value() + val);

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit& GeneralUnit::operator-=(double val)
{
    assertCompatible(GeneralUnit(val));
    setScale(value() - val);

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit GeneralUnit::operator+(double val) const
{
    GeneralUnit ret = *this;
    ret += val;
    return ret;
}

GeneralUnit GeneralUnit::operator-(double val) const
{
    GeneralUnit ret = *this;
    ret -= val;
    return ret;
}

GeneralUnit GeneralUnit::operator*=(const GeneralUnit &other)
{
    setScale(value() * other.value());
    Mass += other.Mass;
    Length += other.Length;
    Time += other.Time;
    Charge += other.Charge;
    temperature += other.temperature;
    Quantity += other.Quantity;
    Angle += other.Angle;

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit GeneralUnit::operator/=(const GeneralUnit &other)
{
    setScale(value() / other.value());
    Mass -= other.Mass;
    Length -= other.Length;
    Time -= other.Time;
    Charge -= other.Charge;
    temperature -= other.temperature;
    Quantity -= other.Quantity;
    Angle -= other.Angle;

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit GeneralUnit::operator*(const GeneralUnit &other) const
{
    GeneralUnit ret = *this;
    ret *= other;

    return ret;
}

GeneralUnit GeneralUnit::operator/(const GeneralUnit &other) const
{
    GeneralUnit ret = *this;
    ret /= other;
    return ret;
}

GeneralUnit& GeneralUnit::operator*=(double val)
{
    setScale(value() * val);

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit& GeneralUnit::operator/=(double val)
{
    setScale(value() / val);

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit& GeneralUnit::operator*=(int val)
{
    setScale(value() * val);

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit& GeneralUnit::operator/=(int val)
{
    setScale(value() / val);

    if (this->isZero())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit GeneralUnit::operator*(double val) const
{
    GeneralUnit ret = *this;
    ret *= val;
    return ret;
}

GeneralUnit GeneralUnit::operator/(double val) const
{
    GeneralUnit ret = *this;
    ret /= val;
    return ret;
}

GeneralUnit GeneralUnit::operator*(int val) const
{
    GeneralUnit ret = *this;
    ret *= val;
    return ret;
}

GeneralUnit GeneralUnit::operator/(int val) const
{
    GeneralUnit ret = *this;
    ret /= val;
    return ret;
}

GeneralUnit GeneralUnit::invert() const
{
    GeneralUnit ret;

    ret.setScale( 1.0 / value() );

    ret.Mass = -Mass;
    ret.Length = -Length;
    ret.Time = -Time;
    ret.Charge = -Charge;
    ret.temperature = -temperature;
    ret.Quantity = -Quantity;
    ret.Angle = -Angle;

    return ret;
}
