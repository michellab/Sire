#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireMM_properties.h"

#include "SireFF/errors.h"
#include "SireMaths/maths.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireUnits/units.h"
#include "switchingfunction.h"
#include <QMutex>
#include <cmath>
#include <numeric>
#include "switchingfunction.h"
#include "SireBase/errors.h"
#include "SireBase/lengthproperty.h"
#include "SireBase/numberproperty.h"
#include "SireBase/properties.h"
#include "SireBase/stringproperty.h"
#include "SireError/errors.h"
#include "SireMaths/multidouble.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireVol/cartesian.h"
#include "SireVol/gridinfo.h"
#include "SireVol/periodicbox.h"
#include "cljboxes.h"
#include "cljfunction.h"
#include "switchingfunction.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tostring.h"
#include <QElapsedTimer>
#include "cljfunction.h"
#include "SireCAS/errors.h"
#include "SireCAS/expression.h"
#include "SireCAS/symbols.h"
#include "SireCAS/values.h"
#include "SireError/errors.h"
#include "SireFF/forcetable.h"
#include "SireMol/moleculedata.h"
#include "SireMol/molecules.h"
#include "SireMol/molid.h"
#include "SireMol/molnum.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "restraint.h"
#include "restraint.h"
void register_SireMM_properties()
{
    register_property_container< SireMM::SwitchFuncPtr, SireMM::SwitchingFunction >();
    register_property_container< SireMM::CLJFunctionPtr, SireMM::CLJFunction >();
    register_property_container< SireMM::RestraintPtr, SireMM::Restraint >();
    register_property_container< SireMM::Restraint3DPtr, SireMM::Restraint3D >();
}
