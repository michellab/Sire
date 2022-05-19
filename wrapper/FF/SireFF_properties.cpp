#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireFF_properties.h"

#include "SireBase/errors.h"
#include "SireError/errors.h"
#include "SireMol/errors.h"
#include "SireMol/mgnum.h"
#include "SireMol/molecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/mover.hpp"
#include "SireMol/partialmolecule.h"
#include "SireMol/viewsofmol.h"
#include "SireStream/datastream.h"
#include "forcefield.h"
#include <QDebug>
#include <QMutex>
#include "forcefield.h"
#include "SireError/errors.h"
#include "SireStream/datastream.h"
#include "probe.h"
#include "probe.h"
#include "SireMol/evaluator.h"
#include "SireMol/mgidx.h"
#include "SireMol/molecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/moleculegroups.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireVol/aabox.h"
#include "SireVol/errors.h"
#include "forcetable.h"
#include "point.h"
#include <boost/tuple/tuple.hpp>
#include "point.h"
void register_SireFF_properties()
{
    register_property_container< SireFF::FFPtr, SireFF::FF >();
    register_property_container< SireFF::ProbePtr, SireFF::Probe >();
    register_property_container< SireFF::PointPtr, SireFF::Point >();
}
