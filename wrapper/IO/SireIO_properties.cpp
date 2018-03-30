#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireIO_properties.h"

#include "SireBase/booleanproperty.h"
#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireError/errors.h"
#include "SireFF/ffdetail.h"
#include "SireIO/errors.h"
#include "SireMM/mmdetail.h"
#include "SireMol/molecule.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireSystem/system.h"
#include "moleculeparser.h"
#include <QDebug>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QMutex>
#include <QTextStream>
#include "moleculeparser.h"
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
    register_property_container< SireIO::MoleculeParserPtr, SireIO::MoleculeParser >();
    register_property_container< SireIO::IOPtr, SireIO::IOBase >();
}
