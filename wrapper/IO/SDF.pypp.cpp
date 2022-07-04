// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SDF.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/parallel.h"

#include "SireBase/propertylist.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireIO/sdf.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/atompropertylist.h"

#include "SireMol/atomradicals.h"

#include "SireMol/bondid.h"

#include "SireMol/bondtype.h"

#include "SireMol/connectivity.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireMol/radical.h"

#include "SireMol/stereoscopy.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "sdf.h"

#include "sire_version.h"

#include <QFile>

#include <QtMath>

#include "sdf.h"

SireIO::SDF __copy__(const SireIO::SDF &other){ return SireIO::SDF(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SDF_class(){

    { //::SireIO::SDF
        typedef bp::class_< SireIO::SDF, bp::bases< SireIO::MoleculeParser, SireBase::Property > > SDF_exposer_t;
        SDF_exposer_t SDF_exposer = SDF_exposer_t( "SDF", "This class holds a parser for reading and writing\nStructure Data File (SDF) molecular file formats\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope SDF_scope( SDF_exposer );
        SDF_exposer.def( bp::init< QString const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("filename"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the file called filename. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        SDF_exposer.def( bp::init< QStringList const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("lines"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the passed text lines. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        SDF_exposer.def( bp::init< SireSystem::System const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("system"), bp::arg("map")=SireBase::PropertyMap() ), "Construct this parser by extracting all necessary information from the\npassed SireSystem::System, looking for the properties that are specified\nin the passed property map") );
        SDF_exposer.def( bp::init< SireIO::SDF const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::SDF::canFollow
        
            typedef bool ( ::SireIO::SDF::*canFollow_function_type)(  ) const;
            canFollow_function_type canFollow_function_value( &::SireIO::SDF::canFollow );
            
            SDF_exposer.def( 
                "canFollow"
                , canFollow_function_value
                , "The SDF cannot follow another lead parsers." );
        
        }
        { //::SireIO::SDF::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::SDF::*construct_function_type)( ::QString const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::SDF::construct );
            
            SDF_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("filename"), bp::arg("map") )
                , "Return the parser that has been constructed by reading in the passed\nfile using the passed properties" );
        
        }
        { //::SireIO::SDF::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::SDF::*construct_function_type)( ::QStringList const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::SDF::construct );
            
            SDF_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("lines"), bp::arg("map") )
                , "Return the parser that has been constructed by reading in the passed\ntext lines using the passed properties" );
        
        }
        { //::SireIO::SDF::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::SDF::*construct_function_type)( ::SireSystem::System const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::SDF::construct );
            
            SDF_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("system"), bp::arg("map") )
                , "Return the parser that has been constructed by extract all necessary\ndata from the passed SireSystem::System using the specified properties" );
        
        }
        { //::SireIO::SDF::formatDescription
        
            typedef ::QString ( ::SireIO::SDF::*formatDescription_function_type)(  ) const;
            formatDescription_function_type formatDescription_function_value( &::SireIO::SDF::formatDescription );
            
            SDF_exposer.def( 
                "formatDescription"
                , formatDescription_function_value
                , "Return a description of the file format" );
        
        }
        { //::SireIO::SDF::formatName
        
            typedef ::QString ( ::SireIO::SDF::*formatName_function_type)(  ) const;
            formatName_function_type formatName_function_value( &::SireIO::SDF::formatName );
            
            SDF_exposer.def( 
                "formatName"
                , formatName_function_value
                , "Return the format name that is used to identify this file format within Sire" );
        
        }
        { //::SireIO::SDF::formatSuffix
        
            typedef ::QStringList ( ::SireIO::SDF::*formatSuffix_function_type)(  ) const;
            formatSuffix_function_type formatSuffix_function_value( &::SireIO::SDF::formatSuffix );
            
            SDF_exposer.def( 
                "formatSuffix"
                , formatSuffix_function_value
                , "Return the suffixes that these files are normally associated with" );
        
        }
        { //::SireIO::SDF::isLead
        
            typedef bool ( ::SireIO::SDF::*isLead_function_type)(  ) const;
            isLead_function_type isLead_function_value( &::SireIO::SDF::isLead );
            
            SDF_exposer.def( 
                "isLead"
                , isLead_function_value
                , "Return whether or not this is a lead parser. The lead parser is responsible\nfor starting the process of turning the parsed file into the System. There\nmust be one and one-only lead parser in a set of parsers creating a System" );
        
        }
        { //::SireIO::SDF::nAtoms
        
            typedef int ( ::SireIO::SDF::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::SDF::nAtoms );
            
            SDF_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , "Return the total number of atoms." );
        
        }
        { //::SireIO::SDF::nAtoms
        
            typedef int ( ::SireIO::SDF::*nAtoms_function_type)( int ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::SDF::nAtoms );
            
            SDF_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("i") )
                , "Return the number of atoms in molecule i." );
        
        }
        { //::SireIO::SDF::nMolecules
        
            typedef int ( ::SireIO::SDF::*nMolecules_function_type)(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireIO::SDF::nMolecules );
            
            SDF_exposer.def( 
                "nMolecules"
                , nMolecules_function_value
                , "Return the number of molecules loaded in this file" );
        
        }
        SDF_exposer.def( bp::self != bp::self );
        { //::SireIO::SDF::operator=
        
            typedef ::SireIO::SDF & ( ::SireIO::SDF::*assign_function_type)( ::SireIO::SDF const & ) ;
            assign_function_type assign_function_value( &::SireIO::SDF::operator= );
            
            SDF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        SDF_exposer.def( bp::self == bp::self );
        { //::SireIO::SDF::parseWarnings
        
            typedef ::QStringList ( ::SireIO::SDF::*parseWarnings_function_type)(  ) const;
            parseWarnings_function_type parseWarnings_function_value( &::SireIO::SDF::parseWarnings );
            
            SDF_exposer.def( 
                "parseWarnings"
                , parseWarnings_function_value
                , "Return any warnings raised when parsing this file" );
        
        }
        { //::SireIO::SDF::toString
        
            typedef ::QString ( ::SireIO::SDF::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireIO::SDF::toString );
            
            SDF_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this parser" );
        
        }
        { //::SireIO::SDF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::SDF::typeName );
            
            SDF_exposer.def( 
                "typeName"
                , typeName_function_value
                , "Return the C++ name for this class" );
        
        }
        { //::SireIO::SDF::what
        
            typedef char const * ( ::SireIO::SDF::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireIO::SDF::what );
            
            SDF_exposer.def( 
                "what"
                , what_function_value
                , "Return the C++ name for this class" );
        
        }
        SDF_exposer.staticmethod( "typeName" );
        SDF_exposer.def( "__copy__", &__copy__);
        SDF_exposer.def( "__deepcopy__", &__copy__);
        SDF_exposer.def( "clone", &__copy__);
        SDF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::SDF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SDF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::SDF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SDF_exposer.def_pickle(sire_pickle_suite< ::SireIO::SDF >());
        SDF_exposer.def( "__str__", &__str__< ::SireIO::SDF > );
        SDF_exposer.def( "__repr__", &__str__< ::SireIO::SDF > );
    }

}
