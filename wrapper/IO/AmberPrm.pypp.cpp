// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AmberPrm.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/findexe.h"

#include "SireBase/parallel.h"

#include "SireBase/stringproperty.h"

#include "SireBase/tempdir.h"

#include "SireBase/unittest.h"

#include "SireCAS/trigfuncs.h"

#include "SireError/errors.h"

#include "SireIO/amberformat.h"

#include "SireIO/amberprm.h"

#include "SireIO/amberrst7.h"

#include "SireIO/errors.h"

#include "SireMM/amberparams.h"

#include "SireMM/atomljs.h"

#include "SireMM/cljnbpairs.h"

#include "SireMM/internalff.h"

#include "SireMM/ljparameter.h"

#include "SireMaths/maths.h"

#include "SireMol/amberparameters.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomcutting.h"

#include "SireMol/atomeditor.h"

#include "SireMol/atomelements.h"

#include "SireMol/atomidx.h"

#include "SireMol/atommasses.h"

#include "SireMol/atomvelocities.h"

#include "SireMol/cgatomidx.h"

#include "SireMol/connectivity.h"

#include "SireMol/element.h"

#include "SireMol/mgname.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireMol/reseditor.h"

#include "SireMol/residuecutting.h"

#include "SireMol/select.h"

#include "SireMol/selector.hpp"

#include "SireMol/trajectory.h"

#include "SireMove/flexibility.h"

#include "SireMove/internalmove.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/cartesian.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "amberprm.h"

#include <QDateTime>

#include <QDebug>

#include <QElapsedTimer>

#include <QFile>

#include <QHash>

#include <QRegularExpression>

#include <QSet>

#include <QTextStream>

#include <tuple>

#include "amberprm.h"

SireIO::AmberPrm __copy__(const SireIO::AmberPrm &other){ return SireIO::AmberPrm(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_AmberPrm_class(){

    { //::SireIO::AmberPrm
        typedef bp::class_< SireIO::AmberPrm, bp::bases< SireIO::MoleculeParser, SireBase::Property > > AmberPrm_exposer_t;
        AmberPrm_exposer_t AmberPrm_exposer = AmberPrm_exposer_t( "AmberPrm", "This class represents an Amber-format parameter file, currently\nsupporting top files produced from Amber7 until Amber16\n\nThe format of this file is described here;\n\nhttp:ambermd.orgformats.html\n\n(specifically the PARM parametertopology file specification)\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope AmberPrm_scope( AmberPrm_exposer );
        bp::enum_< SireIO::AmberPrm::FLAG_TYPE>("FLAG_TYPE")
            .value("UNKNOWN", SireIO::AmberPrm::UNKNOWN)
            .value("INTEGER", SireIO::AmberPrm::INTEGER)
            .value("FLOAT", SireIO::AmberPrm::FLOAT)
            .value("STRING", SireIO::AmberPrm::STRING)
            .export_values()
            ;
        AmberPrm_exposer.def( bp::init< QString const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("filename"), bp::arg("map")=SireBase::PropertyMap() ), "Construct by reading from the file called filename") );
        AmberPrm_exposer.def( bp::init< QStringList const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("lines"), bp::arg("map")=SireBase::PropertyMap() ), "Construct by reading from the contained in the passed\nset of lines") );
        AmberPrm_exposer.def( bp::init< SireSystem::System const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("system"), bp::arg("map")=SireBase::PropertyMap() ), "Construct by converting from the passed system, using the passed property\nmap to find the right properties") );
        AmberPrm_exposer.def( bp::init< SireIO::AmberPrm const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::AmberPrm::assertSane
        
            typedef void ( ::SireIO::AmberPrm::*assertSane_function_type)(  ) const;
            assertSane_function_type assertSane_function_value( &::SireIO::AmberPrm::assertSane );
            
            AmberPrm_exposer.def( 
                "assertSane"
                , assertSane_function_value
                , bp::release_gil_policy()
                , "Run through all of the data that has been read and perform a series\nof tests that will see if the prm7 data is sane. If any test fails,\nthen an exception will be thrown" );
        
        }
        { //::SireIO::AmberPrm::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::AmberPrm::*construct_function_type)( ::QString const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::AmberPrm::construct );
            
            AmberPrm_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("filename"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return this parser constructed from the passed filename" );
        
        }
        { //::SireIO::AmberPrm::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::AmberPrm::*construct_function_type)( ::QStringList const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::AmberPrm::construct );
            
            AmberPrm_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("lines"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return this parser constructed from the passed set of lines" );
        
        }
        { //::SireIO::AmberPrm::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::AmberPrm::*construct_function_type)( ::SireSystem::System const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::AmberPrm::construct );
            
            AmberPrm_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("system"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return this parser constructed from the passed SireSystem::System" );
        
        }
        { //::SireIO::AmberPrm::flagType
        
            typedef ::SireIO::AmberPrm::FLAG_TYPE ( ::SireIO::AmberPrm::*flagType_function_type)( ::QString const & ) const;
            flagType_function_type flagType_function_value( &::SireIO::AmberPrm::flagType );
            
            AmberPrm_exposer.def( 
                "flagType"
                , flagType_function_value
                , ( bp::arg("flag") )
                , bp::release_gil_policy()
                , "Return the flag type for the data associated with the passed flag.\nThis returns UNKNOWN if this is not known" );
        
        }
        { //::SireIO::AmberPrm::flags
        
            typedef ::QStringList ( ::SireIO::AmberPrm::*flags_function_type)(  ) const;
            flags_function_type flags_function_value( &::SireIO::AmberPrm::flags );
            
            AmberPrm_exposer.def( 
                "flags"
                , flags_function_value
                , bp::release_gil_policy()
                , "Return all of the flags that are held in this file" );
        
        }
        { //::SireIO::AmberPrm::floatData
        
            typedef ::QVector< double > ( ::SireIO::AmberPrm::*floatData_function_type)( ::QString const & ) const;
            floatData_function_type floatData_function_value( &::SireIO::AmberPrm::floatData );
            
            AmberPrm_exposer.def( 
                "floatData"
                , floatData_function_value
                , ( bp::arg("flag") )
                , bp::release_gil_policy()
                , "Return the float data for the passed flag. This returns an empty\nlist if there is no data associated with this flag. This raises\nan invalid_cast error if data exists, but it is the wrong type" );
        
        }
        { //::SireIO::AmberPrm::forcefield
        
            typedef ::SireMM::MMDetail ( ::SireIO::AmberPrm::*forcefield_function_type)(  ) const;
            forcefield_function_type forcefield_function_value( &::SireIO::AmberPrm::forcefield );
            
            AmberPrm_exposer.def( 
                "forcefield"
                , forcefield_function_value
                , bp::release_gil_policy()
                , "Return the forcefield for the molecules in this file" );
        
        }
        { //::SireIO::AmberPrm::formatDescription
        
            typedef ::QString ( ::SireIO::AmberPrm::*formatDescription_function_type)(  ) const;
            formatDescription_function_type formatDescription_function_value( &::SireIO::AmberPrm::formatDescription );
            
            AmberPrm_exposer.def( 
                "formatDescription"
                , formatDescription_function_value
                , bp::release_gil_policy()
                , "Return a description of the file format" );
        
        }
        { //::SireIO::AmberPrm::formatName
        
            typedef ::QString ( ::SireIO::AmberPrm::*formatName_function_type)(  ) const;
            formatName_function_type formatName_function_value( &::SireIO::AmberPrm::formatName );
            
            AmberPrm_exposer.def( 
                "formatName"
                , formatName_function_value
                , bp::release_gil_policy()
                , "Return the format name that is used to identify this file format within Sire" );
        
        }
        { //::SireIO::AmberPrm::formatSuffix
        
            typedef ::QStringList ( ::SireIO::AmberPrm::*formatSuffix_function_type)(  ) const;
            formatSuffix_function_type formatSuffix_function_value( &::SireIO::AmberPrm::formatSuffix );
            
            AmberPrm_exposer.def( 
                "formatSuffix"
                , formatSuffix_function_value
                , bp::release_gil_policy()
                , "Return the suffixes that AmberPrm files are normally associated with" );
        
        }
        { //::SireIO::AmberPrm::getMolecule
        
            typedef ::SireMol::Molecule ( ::SireIO::AmberPrm::*getMolecule_function_type)( int,::SireBase::PropertyMap const & ) const;
            getMolecule_function_type getMolecule_function_value( &::SireIO::AmberPrm::getMolecule );
            
            AmberPrm_exposer.def( 
                "getMolecule"
                , getMolecule_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the ith molecule that is described by this AmberPrm file. Note\nthat this molecule wont have any coordinate data, as this is not\nprovided in this file" );
        
        }
        { //::SireIO::AmberPrm::getMolecule
        
            typedef ::SireMol::Molecule ( ::SireIO::AmberPrm::*getMolecule_function_type)( int,::SireIO::AmberRst7 const &,::SireBase::PropertyMap const & ) const;
            getMolecule_function_type getMolecule_function_value( &::SireIO::AmberPrm::getMolecule );
            
            AmberPrm_exposer.def( 
                "getMolecule"
                , getMolecule_function_value
                , ( bp::arg("i"), bp::arg("rst"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the ith molecule that is described by this AmberPrm file, getting\nthe coordinate (and possibly velocity) data from the passed AmberRst file" );
        
        }
        { //::SireIO::AmberPrm::intData
        
            typedef ::QVector< long long > ( ::SireIO::AmberPrm::*intData_function_type)( ::QString const & ) const;
            intData_function_type intData_function_value( &::SireIO::AmberPrm::intData );
            
            AmberPrm_exposer.def( 
                "intData"
                , intData_function_value
                , ( bp::arg("flag") )
                , bp::release_gil_policy()
                , "Return the integer data for the passed flag. This returns an empty\nlist if there is no data associated with this flag. This raises\nan invalid_cast error if data exists, but it is the wrong type" );
        
        }
        { //::SireIO::AmberPrm::isTopology
        
            typedef bool ( ::SireIO::AmberPrm::*isTopology_function_type)(  ) const;
            isTopology_function_type isTopology_function_value( &::SireIO::AmberPrm::isTopology );
            
            AmberPrm_exposer.def( 
                "isTopology"
                , isTopology_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::AmberPrm::linesForFlag
        
            typedef ::QVector< QString > ( ::SireIO::AmberPrm::*linesForFlag_function_type)( ::QString const & ) const;
            linesForFlag_function_type linesForFlag_function_value( &::SireIO::AmberPrm::linesForFlag );
            
            AmberPrm_exposer.def( 
                "linesForFlag"
                , linesForFlag_function_value
                , ( bp::arg("flag") )
                , bp::release_gil_policy()
                , "Return the lines that correspond to the passed flag. This returns an\nempty list of there are no lines associated with the passed flag" );
        
        }
        { //::SireIO::AmberPrm::nAngles
        
            typedef int ( ::SireIO::AmberPrm::*nAngles_function_type)(  ) const;
            nAngles_function_type nAngles_function_value( &::SireIO::AmberPrm::nAngles );
            
            AmberPrm_exposer.def( 
                "nAngles"
                , nAngles_function_value
                , bp::release_gil_policy()
                , "Return the number of angles" );
        
        }
        { //::SireIO::AmberPrm::nAnglesNoHydrogen
        
            typedef int ( ::SireIO::AmberPrm::*nAnglesNoHydrogen_function_type)(  ) const;
            nAnglesNoHydrogen_function_type nAnglesNoHydrogen_function_value( &::SireIO::AmberPrm::nAnglesNoHydrogen );
            
            AmberPrm_exposer.def( 
                "nAnglesNoHydrogen"
                , nAnglesNoHydrogen_function_value
                , bp::release_gil_policy()
                , "Return the number of angles without hydrogen" );
        
        }
        { //::SireIO::AmberPrm::nAnglesWithHydrogen
        
            typedef int ( ::SireIO::AmberPrm::*nAnglesWithHydrogen_function_type)(  ) const;
            nAnglesWithHydrogen_function_type nAnglesWithHydrogen_function_value( &::SireIO::AmberPrm::nAnglesWithHydrogen );
            
            AmberPrm_exposer.def( 
                "nAnglesWithHydrogen"
                , nAnglesWithHydrogen_function_value
                , bp::release_gil_policy()
                , "Return the number of angles containing hydrogen" );
        
        }
        { //::SireIO::AmberPrm::nAtoms
        
            typedef int ( ::SireIO::AmberPrm::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::AmberPrm::nAtoms );
            
            AmberPrm_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the total number of atoms in the file" );
        
        }
        { //::SireIO::AmberPrm::nAtoms
        
            typedef int ( ::SireIO::AmberPrm::*nAtoms_function_type)( int ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::AmberPrm::nAtoms );
            
            AmberPrm_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("molidx") )
                , bp::release_gil_policy()
                , "Return the total number of atoms in the ith molecule in the file" );
        
        }
        { //::SireIO::AmberPrm::nBonds
        
            typedef int ( ::SireIO::AmberPrm::*nBonds_function_type)(  ) const;
            nBonds_function_type nBonds_function_value( &::SireIO::AmberPrm::nBonds );
            
            AmberPrm_exposer.def( 
                "nBonds"
                , nBonds_function_value
                , bp::release_gil_policy()
                , "Return the total number of bonds" );
        
        }
        { //::SireIO::AmberPrm::nBondsNoHydrogen
        
            typedef int ( ::SireIO::AmberPrm::*nBondsNoHydrogen_function_type)(  ) const;
            nBondsNoHydrogen_function_type nBondsNoHydrogen_function_value( &::SireIO::AmberPrm::nBondsNoHydrogen );
            
            AmberPrm_exposer.def( 
                "nBondsNoHydrogen"
                , nBondsNoHydrogen_function_value
                , bp::release_gil_policy()
                , "Return the number of bonds no containing hydrogen" );
        
        }
        { //::SireIO::AmberPrm::nBondsWithHydrogen
        
            typedef int ( ::SireIO::AmberPrm::*nBondsWithHydrogen_function_type)(  ) const;
            nBondsWithHydrogen_function_type nBondsWithHydrogen_function_value( &::SireIO::AmberPrm::nBondsWithHydrogen );
            
            AmberPrm_exposer.def( 
                "nBondsWithHydrogen"
                , nBondsWithHydrogen_function_value
                , bp::release_gil_policy()
                , "Return the number of bonds including hydrogen" );
        
        }
        { //::SireIO::AmberPrm::nDihedrals
        
            typedef int ( ::SireIO::AmberPrm::*nDihedrals_function_type)(  ) const;
            nDihedrals_function_type nDihedrals_function_value( &::SireIO::AmberPrm::nDihedrals );
            
            AmberPrm_exposer.def( 
                "nDihedrals"
                , nDihedrals_function_value
                , bp::release_gil_policy()
                , "Return the number of dihedrals" );
        
        }
        { //::SireIO::AmberPrm::nDihedralsNoHydrogen
        
            typedef int ( ::SireIO::AmberPrm::*nDihedralsNoHydrogen_function_type)(  ) const;
            nDihedralsNoHydrogen_function_type nDihedralsNoHydrogen_function_value( &::SireIO::AmberPrm::nDihedralsNoHydrogen );
            
            AmberPrm_exposer.def( 
                "nDihedralsNoHydrogen"
                , nDihedralsNoHydrogen_function_value
                , bp::release_gil_policy()
                , "Return the number of dihedrals without hydrogen" );
        
        }
        { //::SireIO::AmberPrm::nDihedralsWithHydrogen
        
            typedef int ( ::SireIO::AmberPrm::*nDihedralsWithHydrogen_function_type)(  ) const;
            nDihedralsWithHydrogen_function_type nDihedralsWithHydrogen_function_value( &::SireIO::AmberPrm::nDihedralsWithHydrogen );
            
            AmberPrm_exposer.def( 
                "nDihedralsWithHydrogen"
                , nDihedralsWithHydrogen_function_value
                , bp::release_gil_policy()
                , "Return the number of dihedrals containing hydrogen" );
        
        }
        { //::SireIO::AmberPrm::nExcluded
        
            typedef int ( ::SireIO::AmberPrm::*nExcluded_function_type)(  ) const;
            nExcluded_function_type nExcluded_function_value( &::SireIO::AmberPrm::nExcluded );
            
            AmberPrm_exposer.def( 
                "nExcluded"
                , nExcluded_function_value
                , bp::release_gil_policy()
                , "Return the number of excluded atoms" );
        
        }
        { //::SireIO::AmberPrm::nMolecules
        
            typedef int ( ::SireIO::AmberPrm::*nMolecules_function_type)(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireIO::AmberPrm::nMolecules );
            
            AmberPrm_exposer.def( 
                "nMolecules"
                , nMolecules_function_value
                , bp::release_gil_policy()
                , "Return the number of molecules in the file" );
        
        }
        { //::SireIO::AmberPrm::nResidues
        
            typedef int ( ::SireIO::AmberPrm::*nResidues_function_type)(  ) const;
            nResidues_function_type nResidues_function_value( &::SireIO::AmberPrm::nResidues );
            
            AmberPrm_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , bp::release_gil_policy()
                , "Return the number of residues" );
        
        }
        { //::SireIO::AmberPrm::nTypes
        
            typedef int ( ::SireIO::AmberPrm::*nTypes_function_type)(  ) const;
            nTypes_function_type nTypes_function_value( &::SireIO::AmberPrm::nTypes );
            
            AmberPrm_exposer.def( 
                "nTypes"
                , nTypes_function_value
                , bp::release_gil_policy()
                , "Return the number of distinct atom types" );
        
        }
        AmberPrm_exposer.def( bp::self != bp::self );
        { //::SireIO::AmberPrm::operator=
        
            typedef ::SireIO::AmberPrm & ( ::SireIO::AmberPrm::*assign_function_type)( ::SireIO::AmberPrm const & ) ;
            assign_function_type assign_function_value( &::SireIO::AmberPrm::operator= );
            
            AmberPrm_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AmberPrm_exposer.def( bp::self == bp::self );
        { //::SireIO::AmberPrm::parse
        
            typedef ::SireIO::AmberPrm ( *parse_function_type )( ::QString const &,::SireBase::PropertyMap const & );
            parse_function_type parse_function_value( &::SireIO::AmberPrm::parse );
            
            AmberPrm_exposer.def( 
                "parse"
                , parse_function_value
                , ( bp::arg("filename"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return an AmberPrm object parsed from the passed file" );
        
        }
        { //::SireIO::AmberPrm::stringData
        
            typedef ::QVector< QString > ( ::SireIO::AmberPrm::*stringData_function_type)( ::QString const & ) const;
            stringData_function_type stringData_function_value( &::SireIO::AmberPrm::stringData );
            
            AmberPrm_exposer.def( 
                "stringData"
                , stringData_function_value
                , ( bp::arg("flag") )
                , bp::release_gil_policy()
                , "Return the string data for the passed flag. This returns an empty\nlist if there is no data associated with this flag. This raises\nan invalid_cast error if data exists, but it is the wrong type" );
        
        }
        { //::SireIO::AmberPrm::title
        
            typedef ::QString ( ::SireIO::AmberPrm::*title_function_type)(  ) const;
            title_function_type title_function_value( &::SireIO::AmberPrm::title );
            
            AmberPrm_exposer.def( 
                "title"
                , title_function_value
                , bp::release_gil_policy()
                , "Return the title of the parameter file" );
        
        }
        { //::SireIO::AmberPrm::toString
        
            typedef ::QString ( ::SireIO::AmberPrm::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireIO::AmberPrm::toString );
            
            AmberPrm_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this object" );
        
        }
        { //::SireIO::AmberPrm::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::AmberPrm::typeName );
            
            AmberPrm_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::AmberPrm::what
        
            typedef char const * ( ::SireIO::AmberPrm::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireIO::AmberPrm::what );
            
            AmberPrm_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AmberPrm_exposer.staticmethod( "parse" );
        AmberPrm_exposer.staticmethod( "typeName" );
        AmberPrm_exposer.def( "__copy__", &__copy__);
        AmberPrm_exposer.def( "__deepcopy__", &__copy__);
        AmberPrm_exposer.def( "clone", &__copy__);
        AmberPrm_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::AmberPrm >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AmberPrm_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::AmberPrm >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AmberPrm_exposer.def_pickle(sire_pickle_suite< ::SireIO::AmberPrm >());
        AmberPrm_exposer.def( "__str__", &__str__< ::SireIO::AmberPrm > );
        AmberPrm_exposer.def( "__repr__", &__str__< ::SireIO::AmberPrm > );
    }

}
