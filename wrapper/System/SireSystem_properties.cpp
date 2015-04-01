#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireSystem_properties.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "system.h"
#include "systemmonitor.h"
#include <QMutex>
#include "systemmonitor.h"
#include "SireError/errors.h"
#include "SireMaths/maths.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireSystem/errors.h"
#include "constraint.h"
#include "delta.h"
#include "system.h"
#include <QDebug>
#include "constraint.h"
void register_SireSystem_properties()
{
    register_property_container< SireSystem::SysMonPtr, SireSystem::SystemMonitor >();
    register_property_container< SireSystem::ConstraintPtr, SireSystem::Constraint >();
}
