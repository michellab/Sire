
#ifndef SIREMATHS_GENERALUNITPROPERTY_H
#define SIREMATHS_GENERALUNITPROPERTY_H

#include "SireBase/property.h"
#include "SireBase/arrayproperty.hpp"

#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class GeneralUnitProperty;
class GeneralUnitArrayProperty;
}

QDataStream& operator<<(QDataStream&, const SireBase::GeneralUnitProperty&);
QDataStream& operator>>(QDataStream&, SireBase::GeneralUnitProperty&);

QDataStream& operator<<(QDataStream&, const SireBase::GeneralUnitArrayProperty&);
QDataStream& operator>>(QDataStream&, SireBase::GeneralUnitArrayProperty&);

namespace SireBase
{

class SIREBASE_EXPORT GeneralUnitProperty
    : public ConcreteProperty<GeneralUnitProperty, Property>,
      public SireUnits::Dimension::GeneralUnit
{
public:
    GeneralUnitProperty();

    template<int M, int L, int T,
             int C, int t, int Q, int A>
    GeneralUnitProperty(const SireUnits::Dimension::PhysUnit<M,L,T,C,t,Q,A> &unit);

    GeneralUnitProperty(const SireUnits::Dimension::GeneralUnit &unit);

    GeneralUnitProperty(const GeneralUnitProperty &other);

    ~GeneralUnitProperty();

    static const char* typeName();
    const char *what() const;

    QString toString() const;

    GeneralUnitProperty* clone() const;

    GeneralUnitProperty& operator=(const GeneralUnitProperty &other);

    bool operator==(const GeneralUnitProperty &other) const;
    bool operator!=(const GeneralUnitProperty &other) const;

    virtual bool isAString() const;
    virtual bool isADouble() const;
    virtual bool isAnInteger() const;
    virtual bool isABoolean() const;
    virtual bool isAUnit() const;

    virtual QString asAString() const;
    virtual double asADouble() const;
    virtual int asAnInteger() const;
    virtual bool asABoolean() const;
    virtual SireUnits::Dimension::GeneralUnit asAUnit() const;
};

class SIREBASE_EXPORT GeneralUnitArrayProperty
    : public ConcreteProperty<GeneralUnitArrayProperty, ArrayProperty<SireUnits::Dimension::GeneralUnit> >
{
public:
    GeneralUnitArrayProperty();
    GeneralUnitArrayProperty(const QVector<SireUnits::Dimension::GeneralUnit> &units);
    GeneralUnitArrayProperty(const QList<SireUnits::Dimension::GeneralUnit> &units);
    GeneralUnitArrayProperty(const GeneralUnitArrayProperty &other);

    ~GeneralUnitArrayProperty();

    static const char* typeName();
    const char* what() const;

    QString toString() const;

    GeneralUnitArrayProperty* clone() const;

    GeneralUnitArrayProperty& operator=(const GeneralUnitArrayProperty &other);

    bool operator==(const GeneralUnitArrayProperty &other) const;
    bool operator!=(const GeneralUnitArrayProperty &other) const;
};

SireBase::PropertyPtr wrap(const SireUnits::Dimension::GeneralUnit &unit);
SireBase::PropertyPtr wrap(const QVector<SireUnits::Dimension::GeneralUnit> &units);
SireBase::PropertyPtr wrap(const QList<SireUnits::Dimension::GeneralUnit> &units);

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<int M, int L, int T,
         int C, int t, int Q, int A>
SIRE_OUTOFLINE_TEMPLATE
GeneralUnitProperty::GeneralUnitProperty(const SireUnits::Dimension::PhysUnit<M,L,T,C,t,Q,A> &unit)
                    : SireBase::ConcreteProperty<GeneralUnitProperty, SireBase::Property>(),
                      GeneralUnit(unit)
{}

#endif

}

Q_DECLARE_METATYPE(SireBase::GeneralUnitProperty);
Q_DECLARE_METATYPE(SireBase::GeneralUnitArrayProperty);

SIRE_EXPOSE_FUNCTION(SireBase::wrap)

SIRE_EXPOSE_CLASS(SireBase::GeneralUnitProperty)
SIRE_EXPOSE_CLASS(SireBase::GeneralUnitArrayProperty)

SIRE_END_HEADER

#endif
