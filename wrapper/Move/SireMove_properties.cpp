#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireMove_properties.h"

#include "SireCAS/symbol.h"
#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/moleculedata.h"
#include "SireMol/moleculeinfodata.h"
#include "SireMol/moleculeview.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireUnits/dimensions.h"
#include "SireUnits/temperature.h"
#include "SireUnits/units.h"
#include "velocitygenerator.h"
#include "velocitygenerator.h"
#include "SireError/errors.h"
#include "SireID/index.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireSystem/system.h"
#include "SireSystem/systemmonitor.h"
#include "moves.h"
#include "simstore.h"
#include "suprasystem.h"
#include "suprasystem.h"
#include "SireMaths/quaternion.h"
#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireSystem/system.h"
#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"
#include "SireVol/space.h"
#include "molinserter.h"
#include "molinserter.h"
#include "SireCAS/symbol.h"
#include "SireFF/forcefields.h"
#include "SireFF/forcetable.h"
#include "SireMaths/rangenerator.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/moleculeview.h"
#include "SireStream/datastream.h"
#include "SireSystem/system.h"
#include "ensemble.h"
#include "integrator.h"
#include "integratorworkspace.h"
#include "integrator.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "suprasubsystem.h"
#include "suprasubsystem.h"
#include "SireError/errors.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "suprasubmoves.h"
#include "suprasubsystem.h"
#include "suprasubmoves.h"
#include "SireError/errors.h"
#include "SireID/index.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "supramoves.h"
#include "suprasystem.h"
#include "supramoves.h"
#include "SireError/errors.h"
#include "SireMaths/rangenerator.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireSystem/constraints.h"
#include "SireSystem/errors.h"
#include "SireSystem/system.h"
#include "SireUnits/dimensions.h"
#include "SireUnits/temperature.h"
#include "SireUnits/units.h"
#include "ensemble.h"
#include "moves.h"
#include "simulation.h"
#include <QDebug>
#include <QElapsedTimer>
#include <QMutex>
#include "moves.h"
#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireSystem/system.h"
#include "sampler.h"
#include "uniformsampler.h"
#include "sampler.h"
#include "SireError/errors.h"
#include "SireMaths/rangenerator.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireSystem/system.h"
#include "SireUnits/temperature.h"
#include "SireUnits/units.h"
#include "ensemble.h"
#include "move.h"
#include <QMutex>
#include "move.h"
#include "SireError/errors.h"
#include "SireMaths/vector.h"
#include "SireMol/atom.h"
#include "SireMol/atomcoords.h"
#include "SireMol/evaluator.h"
#include "SireMol/moleculeview.h"
#include "SireMol/mover.hpp"
#include "SireMol/selector.hpp"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "getpoint.h"
#include <QDebug>
#include "getpoint.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "supramove.h"
#include "suprasystem.h"
#include "supramove.h"
#include "SireMol/molecule.h"
#include "SireMol/moleculegroup.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireSystem/system.h"
#include "moldeleter.h"
#include "uniformsampler.h"
#include "moldeleter.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "suprasubmove.h"
#include "suprasubsystem.h"
#include "suprasubmove.h"
void register_SireMove_properties()
{
    register_property_container< SireMove::VelGenPtr, SireMove::VelocityGenerator >();
    register_property_container< SireMove::SupraSystemPtr, SireMove::SupraSystem >();
    register_property_container< SireMove::MolInserterPtr, SireMove::MolInserter >();
    register_property_container< SireMove::IntegratorPtr, SireMove::Integrator >();
    register_property_container< SireMove::SupraSubSystemPtr, SireMove::SupraSubSystem >();
    register_property_container< SireMove::SupraSubMovesPtr, SireMove::SupraSubMoves >();
    register_property_container< SireMove::SupraMovesPtr, SireMove::SupraMoves >();
    register_property_container< SireMove::MovesPtr, SireMove::Moves >();
    register_property_container< SireMove::SamplerPtr, SireMove::Sampler >();
    register_property_container< SireMove::MovePtr, SireMove::Move >();
    register_property_container< SireMove::GetPointPtr, SireMove::GetPoint >();
    register_property_container< SireMove::SupraMovePtr, SireMove::SupraMove >();
    register_property_container< SireMove::MolDeleterPtr, SireMove::MolDeleter >();
    register_property_container< SireMove::SupraSubMovePtr, SireMove::SupraSubMove >();
}
