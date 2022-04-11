#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "Squire_properties.h"

#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/molnum.h"
#include "SireMol/partialmolecule.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "qmchargecalculator.h"
#include "qmchargecalculator.h"
#include "SireError/errors.h"
#include "SireID/index.h"
#include "SireMaths/maths.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "gto.h"
#include "sgto.h"
#include "gto.h"
#include "SireError/errors.h"
#include "SireMol/molecule.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "latticecharges.h"
#include "qmprogram.h"
#include <QMutex>
#include "qmprogram.h"
void register_Squire_properties()
{
    register_property_container< Squire::QMChargeCalculatorPtr, Squire::QMChargeCalculator >();
    register_property_container< Squire::GTOPtr, Squire::GTO >();
    register_property_container< Squire::QMProgPtr, Squire::QMProgram >();
}
