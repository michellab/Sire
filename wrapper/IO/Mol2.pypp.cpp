// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Mol2.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/parallel.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireIO/mol2.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "mol2.h"

#include "mol2.h"

SireIO::Mol2 __copy__(const SireIO::Mol2 &other){ return SireIO::Mol2(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Mol2_class(){

    { //::SireIO::Mol2
        typedef bp::class_< SireIO::Mol2, bp::bases< SireIO::MoleculeParser, SireBase::Property > > Mol2_exposer_t;
        Mol2_exposer_t Mol2_exposer = Mol2_exposer_t( "Mol2", "This class holds a parser for reading and writing Tripos Mol2 files.\n\nAuthor: Lester Hedges\n", bp::init< >("Constructor") );
        bp::scope Mol2_scope( Mol2_exposer );
        Mol2_exposer.def( bp::init< QString const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("filename"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the file called filename. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        Mol2_exposer.def( bp::init< QStringList const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("lines"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the passed text lines. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        Mol2_exposer.def( bp::init< SireSystem::System const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("system"), bp::arg("map")=SireBase::PropertyMap() ), "Construct this parser by extracting all necessary information from the\npassed SireSystem::System, looking for the properties that are specified\nin the passed property map") );
        Mol2_exposer.def( bp::init< SireIO::Mol2 const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::Mol2::canFollow
        
            typedef bool ( ::SireIO::Mol2::*canFollow_function_type)(  ) const;
            canFollow_function_type canFollow_function_value( &::SireIO::Mol2::canFollow );
            
            Mol2_exposer.def( 
                "canFollow"
                , canFollow_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::Mol2::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::Mol2::*construct_function_type)( ::QString const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::Mol2::construct );
            
            Mol2_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("filename"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return the parser that has been constructed by reading in the passed\nfile using the passed properties" );
        
        }
        { //::SireIO::Mol2::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::Mol2::*construct_function_type)( ::QStringList const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::Mol2::construct );
            
            Mol2_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("lines"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return the parser that has been constructed by reading in the passed\ntext lines using the passed properties" );
        
        }
        { //::SireIO::Mol2::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::Mol2::*construct_function_type)( ::SireSystem::System const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::Mol2::construct );
            
            Mol2_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("system"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return the parser that has been constructed by extract all necessary\ndata from the passed SireSystem::System using the specified properties" );
        
        }
        { //::SireIO::Mol2::formatDescription
        
            typedef ::QString ( ::SireIO::Mol2::*formatDescription_function_type)(  ) const;
            formatDescription_function_type formatDescription_function_value( &::SireIO::Mol2::formatDescription );
            
            Mol2_exposer.def( 
                "formatDescription"
                , formatDescription_function_value
                , bp::release_gil_policy()
                , "Return a description of the file format" );
        
        }
        { //::SireIO::Mol2::formatName
        
            typedef ::QString ( ::SireIO::Mol2::*formatName_function_type)(  ) const;
            formatName_function_type formatName_function_value( &::SireIO::Mol2::formatName );
            
            Mol2_exposer.def( 
                "formatName"
                , formatName_function_value
                , bp::release_gil_policy()
                , "Return the format name that is used to identify this file format within Sire" );
        
        }
        { //::SireIO::Mol2::formatSuffix
        
            typedef ::QStringList ( ::SireIO::Mol2::*formatSuffix_function_type)(  ) const;
            formatSuffix_function_type formatSuffix_function_value( &::SireIO::Mol2::formatSuffix );
            
            Mol2_exposer.def( 
                "formatSuffix"
                , formatSuffix_function_value
                , bp::release_gil_policy()
                , "Return the suffixes that these files are normally associated with" );
        
        }
        { //::SireIO::Mol2::isLead
        
            typedef bool ( ::SireIO::Mol2::*isLead_function_type)(  ) const;
            isLead_function_type isLead_function_value( &::SireIO::Mol2::isLead );
            
            Mol2_exposer.def( 
                "isLead"
                , isLead_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::Mol2::nAtoms
        
            typedef int ( ::SireIO::Mol2::*nAtoms_function_type)( int ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::Mol2::nAtoms );
            
            Mol2_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return the number of atoms in a specific molecule." );
        
        }
        { //::SireIO::Mol2::nAtoms
        
            typedef int ( ::SireIO::Mol2::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::Mol2::nAtoms );
            
            Mol2_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the total number of atoms in all molecules." );
        
        }
        { //::SireIO::Mol2::nBonds
        
            typedef int ( ::SireIO::Mol2::*nBonds_function_type)( int ) const;
            nBonds_function_type nBonds_function_value( &::SireIO::Mol2::nBonds );
            
            Mol2_exposer.def( 
                "nBonds"
                , nBonds_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return the number of bonds in a specific molecule." );
        
        }
        { //::SireIO::Mol2::nBonds
        
            typedef int ( ::SireIO::Mol2::*nBonds_function_type)(  ) const;
            nBonds_function_type nBonds_function_value( &::SireIO::Mol2::nBonds );
            
            Mol2_exposer.def( 
                "nBonds"
                , nBonds_function_value
                , bp::release_gil_policy()
                , "Return the total number of bonds in all molecules." );
        
        }
        { //::SireIO::Mol2::nMolAtoms
        
            typedef ::QVector< int > ( ::SireIO::Mol2::*nMolAtoms_function_type)(  ) const;
            nMolAtoms_function_type nMolAtoms_function_value( &::SireIO::Mol2::nMolAtoms );
            
            Mol2_exposer.def( 
                "nMolAtoms"
                , nMolAtoms_function_value
                , bp::release_gil_policy()
                , "Return the number of atoms in each molecule." );
        
        }
        { //::SireIO::Mol2::nMolBonds
        
            typedef ::QVector< int > ( ::SireIO::Mol2::*nMolBonds_function_type)(  ) const;
            nMolBonds_function_type nMolBonds_function_value( &::SireIO::Mol2::nMolBonds );
            
            Mol2_exposer.def( 
                "nMolBonds"
                , nMolBonds_function_value
                , bp::release_gil_policy()
                , "Return the number of bonds in each molecule." );
        
        }
        { //::SireIO::Mol2::nMolSubstructures
        
            typedef ::QVector< int > ( ::SireIO::Mol2::*nMolSubstructures_function_type)(  ) const;
            nMolSubstructures_function_type nMolSubstructures_function_value( &::SireIO::Mol2::nMolSubstructures );
            
            Mol2_exposer.def( 
                "nMolSubstructures"
                , nMolSubstructures_function_value
                , bp::release_gil_policy()
                , "Return the number of substructures in each molecule." );
        
        }
        { //::SireIO::Mol2::nMolecules
        
            typedef int ( ::SireIO::Mol2::*nMolecules_function_type)(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireIO::Mol2::nMolecules );
            
            Mol2_exposer.def( 
                "nMolecules"
                , nMolecules_function_value
                , bp::release_gil_policy()
                , "Return the number of molecules in the system." );
        
        }
        { //::SireIO::Mol2::nSubstructures
        
            typedef int ( ::SireIO::Mol2::*nSubstructures_function_type)( int ) const;
            nSubstructures_function_type nSubstructures_function_value( &::SireIO::Mol2::nSubstructures );
            
            Mol2_exposer.def( 
                "nSubstructures"
                , nSubstructures_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return the number of substructures in a specific molecule." );
        
        }
        { //::SireIO::Mol2::nSubstructures
        
            typedef int ( ::SireIO::Mol2::*nSubstructures_function_type)(  ) const;
            nSubstructures_function_type nSubstructures_function_value( &::SireIO::Mol2::nSubstructures );
            
            Mol2_exposer.def( 
                "nSubstructures"
                , nSubstructures_function_value
                , bp::release_gil_policy()
                , "Return the total number of substructures in all molecules." );
        
        }
        Mol2_exposer.def( bp::self != bp::self );
        { //::SireIO::Mol2::operator=
        
            typedef ::SireIO::Mol2 & ( ::SireIO::Mol2::*assign_function_type)( ::SireIO::Mol2 const & ) ;
            assign_function_type assign_function_value( &::SireIO::Mol2::operator= );
            
            Mol2_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Mol2_exposer.def( bp::self == bp::self );
        { //::SireIO::Mol2::toLines
        
            typedef ::QVector< QString > ( ::SireIO::Mol2::*toLines_function_type)(  ) const;
            toLines_function_type toLines_function_value( &::SireIO::Mol2::toLines );
            
            Mol2_exposer.def( 
                "toLines"
                , toLines_function_value
                , bp::release_gil_policy()
                , "Convert the the parsed data to a collection of Mol2 record lines." );
        
        }
        { //::SireIO::Mol2::toString
        
            typedef ::QString ( ::SireIO::Mol2::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireIO::Mol2::toString );
            
            Mol2_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this parser" );
        
        }
        { //::SireIO::Mol2::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::Mol2::typeName );
            
            Mol2_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "Return the C++ name for this class" );
        
        }
        { //::SireIO::Mol2::what
        
            typedef char const * ( ::SireIO::Mol2::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireIO::Mol2::what );
            
            Mol2_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "Return the C++ name for this class" );
        
        }
        Mol2_exposer.staticmethod( "typeName" );
        Mol2_exposer.def( "__copy__", &__copy__);
        Mol2_exposer.def( "__deepcopy__", &__copy__);
        Mol2_exposer.def( "clone", &__copy__);
        Mol2_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::Mol2 >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Mol2_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::Mol2 >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Mol2_exposer.def_pickle(sire_pickle_suite< ::SireIO::Mol2 >());
        Mol2_exposer.def( "__str__", &__str__< ::SireIO::Mol2 > );
        Mol2_exposer.def( "__repr__", &__str__< ::SireIO::Mol2 > );
    }

}
