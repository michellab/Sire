// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "CGsWithAtoms.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "withatoms.h"

#include "withatoms.h"

SireMol::CGsWithAtoms __copy__(const SireMol::CGsWithAtoms &other){ return SireMol::CGsWithAtoms(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CGsWithAtoms_class(){

    { //::SireMol::CGsWithAtoms
        typedef bp::class_< SireMol::CGsWithAtoms, bp::bases< SireMol::CGID, SireID::ID > > CGsWithAtoms_exposer_t;
        CGsWithAtoms_exposer_t CGsWithAtoms_exposer = CGsWithAtoms_exposer_t( "CGsWithAtoms", "This ID class identifies CutGroups that contain atoms that\nmatch the passed AtomID\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope CGsWithAtoms_scope( CGsWithAtoms_exposer );
        CGsWithAtoms_exposer.def( bp::init< SireMol::AtomID const & >(( bp::arg("atomid") ), "Construct from the passed AtomID") );
        CGsWithAtoms_exposer.def( bp::init< SireMol::CGsWithAtoms const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::CGsWithAtoms::atomID
        
            typedef ::SireMol::AtomID const & ( ::SireMol::CGsWithAtoms::*atomID_function_type)(  ) const;
            atomID_function_type atomID_function_value( &::SireMol::CGsWithAtoms::atomID );
            
            CGsWithAtoms_exposer.def( 
                "atomID"
                , atomID_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the atom ID" );
        
        }
        { //::SireMol::CGsWithAtoms::hash
        
            typedef ::uint ( ::SireMol::CGsWithAtoms::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::CGsWithAtoms::hash );
            
            CGsWithAtoms_exposer.def( 
                "hash"
                , hash_function_value
                , "Return a hash of this identifier" );
        
        }
        { //::SireMol::CGsWithAtoms::isNull
        
            typedef bool ( ::SireMol::CGsWithAtoms::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::CGsWithAtoms::isNull );
            
            CGsWithAtoms_exposer.def( 
                "isNull"
                , isNull_function_value
                , "Is this selection null?" );
        
        }
        { //::SireMol::CGsWithAtoms::map
        
            typedef ::QList< SireMol::CGIdx > ( ::SireMol::CGsWithAtoms::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::CGsWithAtoms::map );
            
            CGsWithAtoms_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , "Map this ID to the list of indicies of CutGroups that match this ID\nThrow: SireMol::missing_cutgroup\nThrow: SireError::invalid_index\n" );
        
        }
        CGsWithAtoms_exposer.def( bp::self != bp::self );
        { //::SireMol::CGsWithAtoms::operator=
        
            typedef ::SireMol::CGsWithAtoms & ( ::SireMol::CGsWithAtoms::*assign_function_type)( ::SireMol::CGsWithAtoms const & ) ;
            assign_function_type assign_function_value( &::SireMol::CGsWithAtoms::operator= );
            
            CGsWithAtoms_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CGsWithAtoms_exposer.def( bp::self == bp::other< SireID::ID >() );
        CGsWithAtoms_exposer.def( bp::self == bp::self );
        { //::SireMol::CGsWithAtoms::toString
        
            typedef ::QString ( ::SireMol::CGsWithAtoms::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::CGsWithAtoms::toString );
            
            CGsWithAtoms_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representatio of this ID" );
        
        }
        { //::SireMol::CGsWithAtoms::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::CGsWithAtoms::typeName );
            
            CGsWithAtoms_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::CGsWithAtoms::what
        
            typedef char const * ( ::SireMol::CGsWithAtoms::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::CGsWithAtoms::what );
            
            CGsWithAtoms_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        CGsWithAtoms_exposer.staticmethod( "typeName" );
        CGsWithAtoms_exposer.def( "__copy__", &__copy__);
        CGsWithAtoms_exposer.def( "__deepcopy__", &__copy__);
        CGsWithAtoms_exposer.def( "clone", &__copy__);
        CGsWithAtoms_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::CGsWithAtoms >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGsWithAtoms_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::CGsWithAtoms >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGsWithAtoms_exposer.def( "__str__", &__str__< ::SireMol::CGsWithAtoms > );
        CGsWithAtoms_exposer.def( "__repr__", &__str__< ::SireMol::CGsWithAtoms > );
        CGsWithAtoms_exposer.def( "__hash__", &::SireMol::CGsWithAtoms::hash );
    }

}
