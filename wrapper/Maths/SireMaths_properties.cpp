#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireMaths_properties.h"

#include "SireError/errors.h"
#include "SireMaths/maths.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "accumulator.h"
#include "histogram.h"
#include <QDebug>
#include <QMutex>
#include <cmath>
#include "accumulator.h"
void register_SireMaths_properties()
{
    register_property_container< SireMaths::AccumulatorPtr, SireMaths::Accumulator >();
}
