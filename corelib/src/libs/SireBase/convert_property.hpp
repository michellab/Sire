#ifndef SIREBASE_CONVERT_PROPERTY_HPP
#define SIREBASE_CONVERT_PROPERTY_HPP

#include <boost/type_traits.hpp>

#include "property.h"
#include "variantproperty.h"
#include "generalunitproperty.h"

namespace SireBase
{

namespace detail
{
    template<int T>
    struct convert_property
    {
        template<class V>
        static SireBase::PropertyPtr convert(const V &value);
    };

    template<>
    struct convert_property<true>
    {
        template<class V>
        static SireBase::PropertyPtr convert(const V &value)
        {
            return SireBase::PropertyPtr(value);
        }
    };

    template<>
    struct convert_property<false>
    {
        template<class V>
        static SireBase::PropertyPtr convert(const V &value)
        {
            return SireBase::PropertyPtr(SireBase::VariantProperty(QVariant::fromValue(value)));
        }
    };

    template<int T>
    struct convert_unit_property
    {
        template<class V>
        static SireBase::PropertyPtr convert(const V &value);
    };

    template<>
    struct convert_unit_property<true>
    {
        template<class V>
        static SireBase::PropertyPtr convert(const V &value)
        {
            return SireBase::PropertyPtr(GeneralUnitProperty(value));
        }
    };

    template<>
    struct convert_unit_property<false>
    {
        template<class V>
        static SireBase::PropertyPtr convert(const V &value)
        {
            return detail::convert_property<boost::is_base_of<SireBase::Property, V>::value>::convert(value);
        }
    };
}

template<class T>
inline SireBase::PropertyPtr convert_property(const T &value)
{
    return detail::convert_unit_property<boost::is_base_of<SireUnits::Dimension::Unit, T>::value>::convert(value);
}

template<>
inline SireBase::PropertyPtr convert_property(const SireUnits::Dimension::GeneralUnit &value)
{
    return SireBase::PropertyPtr(GeneralUnitProperty(value));
}

}

#endif
