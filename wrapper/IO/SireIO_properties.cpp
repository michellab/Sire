#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireIO_properties.h"

#include "SireError/errors.h"
#include "SireMol/cuttingfunction.h"
#include "SireMol/molecule.h"
#include "SireMol/molidx.h"
#include "SireMol/mover.hpp"
#include "SireStream/datastream.h"
#include "iobase.h"
#include <QDebug>
#include <QFile>
#include "iobase.h"
void register_SireIO_properties()
{
    register_property_container< SireIO::IOPtr, SireIO::IOBase >();
}
