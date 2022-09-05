#ifndef SIREUNITS_GENERALUNIT_H
#define SIREUNITS_GENERALUNIT_H

#include <QString>
#include <QHash>

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireUnits{ namespace Dimension {
class GeneralUnit;
}}

SIREUNITS_EXPORT QDataStream& operator<<(QDataStream&, const SireUnits::Dimension::GeneralUnit&);
SIREUNITS_EXPORT QDataStream& operator>>(QDataStream&, SireUnits::Dimension::GeneralUnit&);

namespace SireUnits
{

namespace Dimension
{

class GeneralUnit;
class TempBase;
class _Pvt_Kelvin;

namespace detail
{
void registerTypeName(const GeneralUnit &unit, const char *typnam);
}

class SIREUNITS_EXPORT GeneralUnit : public Unit
{

friend QDataStream& ::operator<<(QDataStream&, const GeneralUnit&);
friend QDataStream& ::operator>>(QDataStream&, GeneralUnit&);

public:
    GeneralUnit();

    explicit GeneralUnit(const TempBase &temperature);

    explicit GeneralUnit(double value);

    template<int M, int L, int T,
             int C, int t, int Q, int A>
    explicit GeneralUnit(const PhysUnit<M,L,T,C,t,Q,A> &unit) : Unit(unit)
    {
        Mass = M;
        Length = L;
        Time = T;
        Charge = C;
        temperature = t;
        Quantity = Q;
        Angle = A;
        detail::registerTypeName(*this, PhysUnit<M,L,T,C,t,Q,A>::typeName());

        if (this->isZero())
            this->operator=(GeneralUnit());
    }

    GeneralUnit(const GeneralUnit &other);

    ~GeneralUnit();

    static const char* typeName();
    const char* what() const;

    QString _to_cpp_type() const;

    QString unitString() const;

    GeneralUnit units() const;

    int MASS() const;
    int LENGTH() const;
    int TIME() const;
    int CHARGE() const;
    int TEMPERATURE() const;
    int QUANTITY() const;
    int ANGLE() const;

    bool hasSameUnits(const GeneralUnit &other) const;

    GeneralUnit& operator=(const GeneralUnit &other);

    bool operator==(const GeneralUnit &other) const;
    bool operator!=(const GeneralUnit &other) const;

    bool operator>(const GeneralUnit &other) const;
    bool operator>=(const GeneralUnit &other) const;
    bool operator<(const GeneralUnit &other) const;
    bool operator<=(const GeneralUnit &other) const;

    GeneralUnit operator-() const;

    GeneralUnit& operator+=(const GeneralUnit &other);

    GeneralUnit& operator-=(const GeneralUnit &other);

    GeneralUnit operator+(const GeneralUnit &other) const;

    GeneralUnit operator-(const GeneralUnit &other) const;

    GeneralUnit operator*=(const GeneralUnit &other);

    GeneralUnit operator/=(const GeneralUnit &other);

    GeneralUnit operator*(const GeneralUnit &other) const;

    GeneralUnit operator/(const GeneralUnit &other) const;

    GeneralUnit& operator+=(double val);
    GeneralUnit& operator-=(double val);

    GeneralUnit operator+(double val) const;
    GeneralUnit operator-(double val) const;

    GeneralUnit& operator*=(double val);
    GeneralUnit& operator/=(double val);

    GeneralUnit& operator*=(int val);
    GeneralUnit& operator/=(int val);

    GeneralUnit operator*(double val) const;
    GeneralUnit operator/(double val) const;

    GeneralUnit operator*(int val) const;
    GeneralUnit operator/(int val) const;

    GeneralUnit invert() const;

    double to(const TempBase &other) const;
    double to(const GeneralUnit &other) const;

    QString toString() const;

    template<class T>
    bool isUnit() const;

    template<class T>
    T toUnit() const;

    bool isDimensionless() const;

    bool isZero() const;

    void setAsDefault(const QString &unit_name) const;

    GeneralUnit getDefault() const;

    QHash<QString,GeneralUnit> components() const;

    GeneralUnit getComponent(const QString &component) const;

    void setComponent(const QString &component,
                      const GeneralUnit &value);

    void addComponent(const QString &component,
                      const GeneralUnit &value);

    void subtractComponent(const QString &component,
                           const GeneralUnit &value);

    template<int M, int L, int T,
             int C, int t, int Q, int A>
    void setComponent(const QString &component,
                      const PhysUnit<M,L,T,C,t,Q,A> &value)
    {
        this->setComponent(component, GeneralUnit(value));
    }

    template<int M, int L, int T,
             int C, int t, int Q, int A>
    void addComponent(const QString &component,
                      const PhysUnit<M,L,T,C,t,Q,A> &value)
    {
        this->addComponent(component, GeneralUnit(value));
    }

    template<int M, int L, int T,
             int C, int t, int Q, int A>
    void subtractComponent(const QString &component,
                      const PhysUnit<M,L,T,C,t,Q,A> &value)
    {
        this->subtractComponent(component, GeneralUnit(value));
    }

private:
    void assertCompatible(const GeneralUnit &other) const;

    int Mass, Length, Time, Charge, temperature, Quantity, Angle;

    QHash<QString,double> comps;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
inline bool GeneralUnit::isUnit() const
{
    return this->MASS() == T::MASS() and
           this->LENGTH() == T::LENGTH() and
           this->TIME() == T::TIME() and
           this->CHARGE() == T::CHARGE() and
           this->TEMPERATURE() == T::TEMPERATURE() and
           this->QUANTITY() == T::QUANTITY() and
           this->ANGLE() == T::ANGLE();
}

template<class T>
inline T GeneralUnit::toUnit() const
{
    this->assertCompatible(GeneralUnit(T(1.0)));
    return this->scaleFactor() * T(1.0);
}

inline GeneralUnit operator*(double val, const GeneralUnit &unit)
{
    return unit * val;
}

inline GeneralUnit operator*(int val, const GeneralUnit &unit)
{
    return unit * val;
}

inline GeneralUnit operator/(double val, const GeneralUnit &unit)
{
    return unit.invert() * val;
}

inline GeneralUnit operator/(int val, const GeneralUnit &unit)
{
    return unit.invert() * val;
}

template<int M, int L, int T,
         int C, int t, int Q, int A>
SIRE_OUTOFLINE_TEMPLATE
PhysUnit<M,L,T,C,t,Q,A>& PhysUnit<M,L,T,C,t,Q,A>::operator=(const GeneralUnit &other)
{
    this->operator=(other.toUnit<PhysUnit<M,L,T,C,t,Q,A>>());
    return *this;
}

template<int M, int L, int T,
         int C, int t, int Q, int A>
SIRE_OUTOFLINE_TEMPLATE
PhysUnit<M,L,T,C,t,Q,A>::PhysUnit(const GeneralUnit &other)
                        : Unit(0)
{
    this->operator=(other);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace Dimension

}

Q_DECLARE_METATYPE( SireUnits::Dimension::GeneralUnit );

SIRE_EXPOSE_CLASS( SireUnits::Dimension::GeneralUnit )

SIRE_END_HEADER

#endif
