#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireBase_properties.h"

#include "SireError/errors.h"
#include "SireError/getbacktrace.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "property.h"
#include <QDebug>
#include <QMutex>
#include "property.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "stringmangler.h"
#include <QMutex>
#include "stringmangler.h"
void register_SireBase_properties()
{
    register_property_container< SireBase::PropertyPtr, SireBase::Property >();
    register_property_container< SireBase::StringManglerPtr, SireBase::StringMangler >();
}
