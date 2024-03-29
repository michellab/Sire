// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Values.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "symbol.h"

#include "values.h"

#include "values.h"

SireCAS::Values __copy__(const SireCAS::Values &other){ return SireCAS::Values(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Values_class(){

    { //::SireCAS::Values
        typedef bp::class_< SireCAS::Values > Values_exposer_t;
        Values_exposer_t Values_exposer = Values_exposer_t( "Values", "\nThis class holds a set of Symbols and their associated values. This is used\nwhen numerically evaluating an equation.\n\nAuthor: Christopher Woods\n", bp::init< >("Construct an empty set of values") );
        bp::scope Values_scope( Values_exposer );
        Values_exposer.def( bp::init< QList< SireCAS::SymbolValue > const & >(( bp::arg("values") ), "Construct from a list of values") );
        Values_exposer.def( bp::init< QHash< SireCAS::Symbol, double > const & >(( bp::arg("values") ), "Construct from a hash of values indexed by symbols") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const & >(( bp::arg("symval0") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const &, SireCAS::SymbolValue const & >(( bp::arg("symval0"), bp::arg("symval1") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const & >(( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const & >(( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const & >(( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const & >(( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const & >(( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5"), bp::arg("symval6") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const & >(( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5"), bp::arg("symval6"), bp::arg("symval7") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const & >(( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5"), bp::arg("symval6"), bp::arg("symval7"), bp::arg("symval8") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const &, SireCAS::SymbolValue const & >(( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5"), bp::arg("symval6"), bp::arg("symval7"), bp::arg("symval8"), bp::arg("symval9") ), "Construct from the passed values") );
        Values_exposer.def( bp::init< SireCAS::Values const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0"), bp::arg("symval1") )
                , bp::release_gil_policy()
                , "Add the passed values" );
        
        }
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2") )
                , bp::release_gil_policy()
                , "Add the passed values" );
        
        }
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3") )
                , bp::release_gil_policy()
                , "Add the passed values" );
        
        }
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4") )
                , bp::release_gil_policy()
                , "Add the passed values" );
        
        }
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5") )
                , bp::release_gil_policy()
                , "Add the passed values" );
        
        }
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5"), bp::arg("symval6") )
                , bp::release_gil_policy()
                , "Add the passed values" );
        
        }
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5"), bp::arg("symval6"), bp::arg("symval7") )
                , bp::release_gil_policy()
                , "Add the passed values" );
        
        }
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5"), bp::arg("symval6"), bp::arg("symval7"), bp::arg("symval8") )
                , bp::release_gil_policy()
                , "Add the passed values" );
        
        }
        { //::SireCAS::Values::add
        
            typedef void ( ::SireCAS::Values::*add_function_type)( ::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const &,::SireCAS::SymbolValue const & ) ;
            add_function_type add_function_value( &::SireCAS::Values::add );
            
            Values_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("symval0"), bp::arg("symval1"), bp::arg("symval2"), bp::arg("symval3"), bp::arg("symval4"), bp::arg("symval5"), bp::arg("symval6"), bp::arg("symval7"), bp::arg("symval8"), bp::arg("symval9") )
                , bp::release_gil_policy()
                , "Add the passed values" );
        
        }
        { //::SireCAS::Values::contains
        
            typedef bool ( ::SireCAS::Values::*contains_function_type)( ::SireCAS::Symbol const & ) const;
            contains_function_type contains_function_value( &::SireCAS::Values::contains );
            
            Values_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Values::count
        
            typedef int ( ::SireCAS::Values::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireCAS::Values::count );
            
            Values_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Values::isEmpty
        
            typedef bool ( ::SireCAS::Values::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireCAS::Values::isEmpty );
            
            Values_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Values::keys
        
            typedef ::QList< SireCAS::Symbol > ( ::SireCAS::Values::*keys_function_type)(  ) const;
            keys_function_type keys_function_value( &::SireCAS::Values::keys );
            
            Values_exposer.def( 
                "keys"
                , keys_function_value
                , bp::release_gil_policy()
                , "Return a list of the symbols that are present in this set" );
        
        }
        Values_exposer.def( bp::self != bp::self );
        { //::SireCAS::Values::operator()
        
            typedef double ( ::SireCAS::Values::*__call___function_type)( ::SireCAS::Symbol const & ) const;
            __call___function_type __call___function_value( &::SireCAS::Values::operator() );
            
            Values_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("sym") )
                , "" );
        
        }
        Values_exposer.def( bp::self == bp::self );
        { //::SireCAS::Values::operator[]
        
            typedef double ( ::SireCAS::Values::*__getitem___function_type)( ::SireCAS::Symbol const & ) const;
            __getitem___function_type __getitem___function_value( &::SireCAS::Values::operator[] );
            
            Values_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("sym") )
                , "" );
        
        }
        { //::SireCAS::Values::remove
        
            typedef void ( ::SireCAS::Values::*remove_function_type)( ::SireCAS::Symbol const & ) ;
            remove_function_type remove_function_value( &::SireCAS::Values::remove );
            
            Values_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Remove the value for the symbol symbol" );
        
        }
        { //::SireCAS::Values::remove
        
            typedef void ( ::SireCAS::Values::*remove_function_type)( ::SireCAS::SymbolID const & ) ;
            remove_function_type remove_function_value( &::SireCAS::Values::remove );
            
            Values_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("symbolid") )
                , bp::release_gil_policy()
                , "Remove the value for the symbol with ID symbolid" );
        
        }
        { //::SireCAS::Values::reserve
        
            typedef void ( ::SireCAS::Values::*reserve_function_type)( int ) ;
            reserve_function_type reserve_function_value( &::SireCAS::Values::reserve );
            
            Values_exposer.def( 
                "reserve"
                , reserve_function_value
                , ( bp::arg("n") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Values::set
        
            typedef void ( ::SireCAS::Values::*set_function_type)( ::SireCAS::Symbol const &,double ) ;
            set_function_type set_function_value( &::SireCAS::Values::set );
            
            Values_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("symbol"), bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Values::set
        
            typedef void ( ::SireCAS::Values::*set_function_type)( ::SireCAS::Values::const_iterator const & ) ;
            set_function_type set_function_value( &::SireCAS::Values::set );
            
            Values_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("it") )
                , bp::release_gil_policy()
                , "Set the value of the symbolvalue pair pointed to by the\niterator it in this set" );
        
        }
        { //::SireCAS::Values::symbols
        
            typedef ::QList< SireCAS::Symbol > ( ::SireCAS::Values::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::Values::symbols );
            
            Values_exposer.def( 
                "symbols"
                , symbols_function_value
                , bp::release_gil_policy()
                , "Return a list of the symbols that are present in this set" );
        
        }
        { //::SireCAS::Values::toString
        
            typedef ::QString ( ::SireCAS::Values::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireCAS::Values::toString );
            
            Values_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of these values" );
        
        }
        { //::SireCAS::Values::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::Values::typeName );
            
            Values_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Values::value
        
            typedef double ( ::SireCAS::Values::*value_function_type)( ::SireCAS::Symbol const & ) const;
            value_function_type value_function_value( &::SireCAS::Values::value );
            
            Values_exposer.def( 
                "value"
                , value_function_value
                , ( bp::arg("sym") )
                , bp::release_gil_policy()
                , "Return the value of the Symbol with ID id, or 0.0 if there is no such symbol" );
        
        }
        { //::SireCAS::Values::values
        
            typedef ::QHash< unsigned int, double > const & ( ::SireCAS::Values::*values_function_type)(  ) const;
            values_function_type values_function_value( &::SireCAS::Values::values );
            
            Values_exposer.def( 
                "values"
                , values_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireCAS::Values::what
        
            typedef char const * ( ::SireCAS::Values::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::Values::what );
            
            Values_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Values_exposer.staticmethod( "typeName" );
        Values_exposer.def( bp::self + bp::other< SireCAS::SymbolValue >() );
        Values_exposer.def( bp::self + bp::self );
        Values_exposer.def( "__copy__", &__copy__);
        Values_exposer.def( "__deepcopy__", &__copy__);
        Values_exposer.def( "clone", &__copy__);
        Values_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::Values >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Values_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::Values >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Values_exposer.def_pickle(sire_pickle_suite< ::SireCAS::Values >());
        Values_exposer.def( "__str__", &__str__< ::SireCAS::Values > );
        Values_exposer.def( "__repr__", &__str__< ::SireCAS::Values > );
        Values_exposer.def( "__len__", &__len_count< ::SireCAS::Values > );
    }

}
