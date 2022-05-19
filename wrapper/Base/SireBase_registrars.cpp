//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireBase_registrars.h"

#include "Helpers/version_error_impl.h"

#include "cpuid.h"
#include "propertylist.h"
#include "majorminorversion.h"
#include "lengthproperty.h"
#include "stringproperty.h"
#include "ranges.h"
#include "stringmangler.h"
#include "properties.h"
#include "propertymap.h"
#include "property.h"
#include "booleanproperty.h"
#include "timeproperty.h"
#include "linktoproperty.h"
#include "numberproperty.h"
#include "variantproperty.h"

#include "Helpers/objectregistry.hpp"

void register_SireBase_objects()
{

    ObjectRegistry::registerConverterFor< SireBase::CPUID >();
    ObjectRegistry::registerConverterFor< SireBase::DoubleArrayProperty >();
    ObjectRegistry::registerConverterFor< SireBase::IntegerArrayProperty >();
    ObjectRegistry::registerConverterFor< SireBase::StringArrayProperty >();
    ObjectRegistry::registerConverterFor< SireBase::PropertyList >();
    ObjectRegistry::registerConverterFor< SireBase::MajorMinorVersion >();
    ObjectRegistry::registerConverterFor< SireBase::Version >();
    ObjectRegistry::registerConverterFor< SireBase::LengthProperty >();
    ObjectRegistry::registerConverterFor< SireBase::StringProperty >();
    ObjectRegistry::registerConverterFor< SireBase::SimpleRange >();
    ObjectRegistry::registerConverterFor< SireBase::NoMangling >();
    ObjectRegistry::registerConverterFor< SireBase::TrimString >();
    ObjectRegistry::registerConverterFor< SireBase::UpperCaseString >();
    ObjectRegistry::registerConverterFor< SireBase::LowerCaseString >();
    ObjectRegistry::registerConverterFor< SireBase::Properties >();
    ObjectRegistry::registerConverterFor< SireBase::PropertyName >();
    ObjectRegistry::registerConverterFor< SireBase::PropertyMap >();
    ObjectRegistry::registerConverterFor< SireBase::NullProperty >();
    ObjectRegistry::registerConverterFor< SireBase::BooleanProperty >();
    ObjectRegistry::registerConverterFor< SireBase::TimeProperty >();
    ObjectRegistry::registerConverterFor< SireBase::LinkToProperty >();
    ObjectRegistry::registerConverterFor< SireBase::NumberProperty >();
    ObjectRegistry::registerConverterFor< SireBase::VariantProperty >();

}

