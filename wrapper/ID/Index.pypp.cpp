// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Index.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "index.h"

#include <QDataStream>

#include "index.h"

SireID::Index __copy__(const SireID::Index &other){ return SireID::Index(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Index_class(){

    { //::SireID::Index
        typedef bp::class_< SireID::Index, bp::bases< SireID::IndexBase > > Index_exposer_t;
        Index_exposer_t Index_exposer = Index_exposer_t( "Index", "", bp::init< bp::optional< qint32 > >(( bp::arg("idx")=(::qint32)(SireID::IndexBase::null()) ), "") );
        bp::scope Index_scope( Index_exposer );
        Index_exposer.def( bp::init< SireID::Index const & >(( bp::arg("other") ), "") );
        { //::SireID::Index::null
        
            typedef ::SireID::Index ( *null_function_type )(  );
            null_function_type null_function_value( &::SireID::Index::null );
            
            Index_exposer.def( 
                "null"
                , null_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::Index::toString
        
            typedef ::QString ( ::SireID::Index::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::Index::toString );
            
            Index_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::Index::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::Index::typeName );
            
            Index_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Index_exposer.staticmethod( "null" );
        Index_exposer.staticmethod( "typeName" );
        Index_exposer.def( "__copy__", &__copy__);
        Index_exposer.def( "__deepcopy__", &__copy__);
        Index_exposer.def( "clone", &__copy__);
        Index_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::Index >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Index_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::Index >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Index_exposer.def_pickle(sire_pickle_suite< ::SireID::Index >());
        Index_exposer.def( "__str__", &__str__< ::SireID::Index > );
        Index_exposer.def( "__repr__", &__str__< ::SireID::Index > );
        Index_exposer.def( "__hash__", &::SireID::Index::hash );
    }

}
