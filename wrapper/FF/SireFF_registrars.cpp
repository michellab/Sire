//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireFF_registrars.h"

#include "forcetable.h"
#include "ffname.h"
#include "patches.h"
#include "ffidx.h"
#include "potentialtable.h"
#include "ffmolgroup.h"
#include "ffcomponent.h"
#include "energytable.h"
#include "ffparameters.h"
#include "patch.h"
#include "fieldtable.h"
#include "forcefield.h"
#include "point.h"
#include "forcefields.h"
#include "probe.h"
#include "ffidentifier.h"

#include "Helpers/objectregistry.hpp"

void register_SireFF_objects()
{

    ObjectRegistry::registerConverterFor< SireFF::ForceTable >();
    ObjectRegistry::registerConverterFor< SireFF::MolForceTable >();
    ObjectRegistry::registerConverterFor< SireFF::FFName >();
    ObjectRegistry::registerConverterFor< SireFF::FFBead >();
    ObjectRegistry::registerConverterFor< SireFF::FFBeadChange >();
    ObjectRegistry::registerConverterFor< SireFF::Patches >();
    ObjectRegistry::registerConverterFor< SireFF::FFIdx >();
    ObjectRegistry::registerConverterFor< SireFF::PotentialTable >();
    ObjectRegistry::registerConverterFor< SireFF::GridPotentialTable >();
    ObjectRegistry::registerConverterFor< SireFF::MolPotentialTable >();
    ObjectRegistry::registerConverterFor< SireFF::FFMolGroup >();
    ObjectRegistry::registerConverterFor< SireFF::SingleComponent >();
    ObjectRegistry::registerConverterFor< SireFF::EnergyTable >();
    ObjectRegistry::registerConverterFor< SireFF::MolEnergyTable >();
    ObjectRegistry::registerConverterFor< SireFF::NullFFParameters >();
    ObjectRegistry::registerConverterFor< SireFF::NullFFParametersArray >();
    ObjectRegistry::registerConverterFor< SireFF::Patch >();
    ObjectRegistry::registerConverterFor< SireFF::FieldTable >();
    ObjectRegistry::registerConverterFor< SireFF::GridFieldTable >();
    ObjectRegistry::registerConverterFor< SireFF::MolFieldTable >();
    ObjectRegistry::registerConverterFor< SireFF::NullFF >();
    ObjectRegistry::registerConverterFor< SireFF::AtomPoint >();
    ObjectRegistry::registerConverterFor< SireFF::VectorPoint >();
    ObjectRegistry::registerConverterFor< SireFF::Center >();
    ObjectRegistry::registerConverterFor< SireFF::CenterOfGeometry >();
    ObjectRegistry::registerConverterFor< SireFF::CenterOfMass >();
    ObjectRegistry::registerConverterFor< SireFF::ForceFields >();
    ObjectRegistry::registerConverterFor< SireFF::NullProbe >();
    ObjectRegistry::registerConverterFor< SireID::Specify<SireFF::FFID> >();
    ObjectRegistry::registerConverterFor< SireID::IDAndSet<SireFF::FFID> >();
    ObjectRegistry::registerConverterFor< SireID::IDOrSet<SireFF::FFID> >();
    ObjectRegistry::registerConverterFor< SireFF::FFIdentifier >();

}

