//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireMM_registrars.h"

#include "amberparams.h"
#include "cljboxes.h"
#include "intercoulombff.h"
#include "cljrffunction.h"
#include "threeatomfunctions.h"
#include "gridff2.h"
#include "cljcalculator.h"
#include "cljatoms.h"
#include "anglerestraint.h"
#include "cljgrid.h"
#include "switchingfunction.h"
#include "intersoftcljff.h"
#include "intergroupff.h"
#include "ljparameter.h"
#include "ljperturbation.h"
#include "intraff.h"
#include "softcljcomponent.h"
#include "cljdelta.h"
#include "cljextractor.h"
#include "fouratomfunctions.h"
#include "interff.h"
#include "intraljff.h"
#include "cljprobe.h"
#include "gromacsparams.h"
#include "ljpair.h"
#include "intrasoftcljff.h"
#include "internalparameters.h"
#include "clj14group.h"
#include "cljcomponent.h"
#include "intragroupff.h"
#include "cljshiftfunction.h"
#include "interljff.h"
#include "cljnbpairs.h"
#include "internalperturbation.h"
#include "distancerestraint.h"
#include "twoatomfunctions.h"
#include "restraintff.h"
#include "multicljcomponent.h"
#include "internalff.h"
#include "cljgroup.h"
#include "intracoulombff.h"
#include "intercljff.h"
#include "cljworkspace.h"
#include "dihedralrestraint.h"
#include "internalcomponent.h"
#include "gridff.h"
#include "intracljff.h"
#include "atomljs.h"
#include "mmdetail.h"
#include "cljfunction.h"
#include "cljparam.h"
#include "restraintcomponent.h"
#include "restraint.h"

#include "Helpers/objectregistry.hpp"

void register_SireMM_objects()
{

    ObjectRegistry::registerConverterFor< SireMM::AmberParams >();
    ObjectRegistry::registerConverterFor< SireMM::AmberBond >();
    ObjectRegistry::registerConverterFor< SireMM::AmberAngle >();
    ObjectRegistry::registerConverterFor< SireMM::AmberDihPart >();
    ObjectRegistry::registerConverterFor< SireMM::AmberDihedral >();
    ObjectRegistry::registerConverterFor< SireMM::AmberNB14 >();
    ObjectRegistry::registerConverterFor< SireMM::AmberNBDihPart >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBox >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBoxIndex >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBoxes >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBoxDistance >();
    ObjectRegistry::registerConverterFor< SireMM::InterCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::CLJRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJIntraRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftIntraRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::ThreeAtomFunctions >();
    ObjectRegistry::registerConverterFor< SireMM::GridFF2 >();
    ObjectRegistry::registerConverterFor< SireMM::CLJCalculator >();
    ObjectRegistry::registerConverterFor< SireMM::CLJAtom >();
    ObjectRegistry::registerConverterFor< SireMM::CLJAtoms >();
    ObjectRegistry::registerConverterFor< SireMM::AngleRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::CLJGrid >();
    ObjectRegistry::registerConverterFor< SireMM::NoCutoff >();
    ObjectRegistry::registerConverterFor< SireMM::HarmonicSwitchingFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CHARMMSwitchingFunction >();
    ObjectRegistry::registerConverterFor< SireMM::InterSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupFF >();
    ObjectRegistry::registerConverterFor< SireMM::LJParameter >();
    ObjectRegistry::registerConverterFor< SireMM::LJPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::IntraFF >();
    ObjectRegistry::registerConverterFor< SireMM::SoftCLJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::CLJDelta >();
    ObjectRegistry::registerConverterFor< SireMM::CLJExtractor >();
    ObjectRegistry::registerConverterFor< SireMM::FourAtomFunctions >();
    ObjectRegistry::registerConverterFor< SireMM::InterFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombProbe >();
    ObjectRegistry::registerConverterFor< SireMM::LJProbe >();
    ObjectRegistry::registerConverterFor< SireMM::CLJProbe >();
    ObjectRegistry::registerConverterFor< SireMM::GromacsAtomType >();
    ObjectRegistry::registerConverterFor< SireMM::GromacsBond >();
    ObjectRegistry::registerConverterFor< SireMM::GromacsAngle >();
    ObjectRegistry::registerConverterFor< SireMM::GromacsDihedral >();
    ObjectRegistry::registerConverterFor< SireMM::LJPair >();
    ObjectRegistry::registerConverterFor< SireMM::IntraSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InternalParameters >();
    ObjectRegistry::registerConverterFor< SireMM::InternalParameters3D >();
    ObjectRegistry::registerConverterFor< SireMM::CLJ14Group >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombComponent >();
    ObjectRegistry::registerConverterFor< SireMM::LJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::CLJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupFF >();
    ObjectRegistry::registerConverterFor< SireMM::CLJShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJIntraShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftIntraShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::InterLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombScaleFactor >();
    ObjectRegistry::registerConverterFor< SireMM::LJScaleFactor >();
    ObjectRegistry::registerConverterFor< SireMM::CLJScaleFactor >();
    ObjectRegistry::registerConverterFor< SireMM::CLJNBPairs >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombNBPairs >();
    ObjectRegistry::registerConverterFor< SireMM::LJNBPairs >();
    ObjectRegistry::registerConverterFor< SireMM::TwoAtomPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::ThreeAtomPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::FourAtomPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::DistanceRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::DoubleDistanceRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::TripleDistanceRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::HarmonicDistanceForceConstant >();
    ObjectRegistry::registerConverterFor< SireMM::TwoAtomFunctions >();
    ObjectRegistry::registerConverterFor< SireMM::RestraintFF >();
    ObjectRegistry::registerConverterFor< SireMM::MultiCLJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::InternalFF >();
    ObjectRegistry::registerConverterFor< SireMM::CLJGroup >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::CLJWorkspace >();
    ObjectRegistry::registerConverterFor< SireMM::DihedralRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::BondComponent >();
    ObjectRegistry::registerConverterFor< SireMM::AngleComponent >();
    ObjectRegistry::registerConverterFor< SireMM::DihedralComponent >();
    ObjectRegistry::registerConverterFor< SireMM::ImproperComponent >();
    ObjectRegistry::registerConverterFor< SireMM::UreyBradleyComponent >();
    ObjectRegistry::registerConverterFor< SireMM::StretchStretchComponent >();
    ObjectRegistry::registerConverterFor< SireMM::StretchBendComponent >();
    ObjectRegistry::registerConverterFor< SireMM::BendBendComponent >();
    ObjectRegistry::registerConverterFor< SireMM::StretchBendTorsionComponent >();
    ObjectRegistry::registerConverterFor< SireMM::Intra14CoulombComponent >();
    ObjectRegistry::registerConverterFor< SireMM::Intra14LJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::Intra14Component >();
    ObjectRegistry::registerConverterFor< SireMM::InternalComponent >();
    ObjectRegistry::registerConverterFor< SireMM::GridFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::AtomLJs >();
    ObjectRegistry::registerConverterFor< SireMM::MMDetail >();
    ObjectRegistry::registerConverterFor< SireMM::NullCLJFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJParams >();
    ObjectRegistry::registerConverterFor< SireMM::CLJParamsArray >();
    ObjectRegistry::registerConverterFor< SireMM::RestraintComponent >();
    ObjectRegistry::registerConverterFor< SireMM::NullRestraint >();

}

