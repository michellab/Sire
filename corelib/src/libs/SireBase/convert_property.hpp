#ifndef SIREBASE_CONVERT_PROPERTY_HPP
#define SIREBASE_CONVERT_PROPERTY_HPP

#include <boost/type_traits.hpp>

#include "property.h"
#include "variantproperty.h"

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
}

template<class T>
SireBase::PropertyPtr convert_property(const T &value)
{
    return detail::convert_property<boost::is_base_of<SireBase::Property, T>::value>::convert(value);
}

}

#endif
