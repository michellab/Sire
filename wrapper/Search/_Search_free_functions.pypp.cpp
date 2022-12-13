// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "_Search_free_functions.pypp.hpp"

namespace bp = boost::python;

#include "approx_equal.h"

#include <QDebug>

#include <QReadWriteLock>

#include "approx_equal.h"

#include "approx_equal.h"

#include <QDebug>

#include <QReadWriteLock>

#include "approx_equal.h"

#include "approx_equal.h"

#include <QDebug>

#include <QReadWriteLock>

#include "approx_equal.h"

#include "approx_equal.h"

#include <QDebug>

#include <QReadWriteLock>

#include "approx_equal.h"

#include "approx_equal.h"

#include <QDebug>

#include <QReadWriteLock>

#include "approx_equal.h"

#include "approx_equal.h"

#include <QDebug>

#include <QReadWriteLock>

#include "approx_equal.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

#include "approx_equal.h"

#include <QDebug>

#include <QReadWriteLock>

#include "approx_equal.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

#include "approx_equal.h"

#include <QDebug>

#include <QReadWriteLock>

#include "approx_equal.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

#include "SireMol/core.h"

#include "SireSearch/parser.h"

#include "helper_funcs.h"

#include <QReadWriteLock>

#include "helper_funcs.h"

void register_free_functions(){

    { //::SireSearch::approx_equal
    
        typedef bool ( *approx_equal_function_type )( double,double );
        approx_equal_function_type approx_equal_function_value( &::SireSearch::approx_equal );
        
        bp::def( 
            "approx_equal"
            , approx_equal_function_value
            , ( bp::arg("val0"), bp::arg("val1") )
            , "" );
    
    }

    { //::SireSearch::approx_greater
    
        typedef bool ( *approx_greater_function_type )( double,double );
        approx_greater_function_type approx_greater_function_value( &::SireSearch::approx_greater );
        
        bp::def( 
            "approx_greater"
            , approx_greater_function_value
            , ( bp::arg("val0"), bp::arg("val1") )
            , "" );
    
    }

    { //::SireSearch::approx_greater_equal
    
        typedef bool ( *approx_greater_equal_function_type )( double,double );
        approx_greater_equal_function_type approx_greater_equal_function_value( &::SireSearch::approx_greater_equal );
        
        bp::def( 
            "approx_greater_equal"
            , approx_greater_equal_function_value
            , ( bp::arg("val0"), bp::arg("val1") )
            , "" );
    
    }

    { //::SireSearch::approx_less
    
        typedef bool ( *approx_less_function_type )( double,double );
        approx_less_function_type approx_less_function_value( &::SireSearch::approx_less );
        
        bp::def( 
            "approx_less"
            , approx_less_function_value
            , ( bp::arg("val0"), bp::arg("val1") )
            , "" );
    
    }

    { //::SireSearch::approx_less_equal
    
        typedef bool ( *approx_less_equal_function_type )( double,double );
        approx_less_equal_function_type approx_less_equal_function_value( &::SireSearch::approx_less_equal );
        
        bp::def( 
            "approx_less_equal"
            , approx_less_equal_function_value
            , ( bp::arg("val0"), bp::arg("val1") )
            , "" );
    
    }

    { //::SireSearch::approx_not_equal
    
        typedef bool ( *approx_not_equal_function_type )( double,double );
        approx_not_equal_function_type approx_not_equal_function_value( &::SireSearch::approx_not_equal );
        
        bp::def( 
            "approx_not_equal"
            , approx_not_equal_function_value
            , ( bp::arg("val0"), bp::arg("val1") )
            , "" );
    
    }

    { //::SireSearch::delete_all_tokens
    
        typedef void ( *delete_all_tokens_function_type )(  );
        delete_all_tokens_function_type delete_all_tokens_function_value( &::SireSearch::delete_all_tokens );
        
        bp::def( 
            "delete_all_tokens"
            , delete_all_tokens_function_value
            , "" );
    
    }

    { //::SireSearch::delete_token
    
        typedef void ( *delete_token_function_type )( ::QString const & );
        delete_token_function_type delete_token_function_value( &::SireSearch::delete_token );
        
        bp::def( 
            "delete_token"
            , delete_token_function_value
            , ( bp::arg("token") )
            , "" );
    
    }

    { //::SireSearch::get_approx_epsilon
    
        typedef double ( *get_approx_epsilon_function_type )(  );
        get_approx_epsilon_function_type get_approx_epsilon_function_value( &::SireSearch::get_approx_epsilon );
        
        bp::def( 
            "get_approx_epsilon"
            , get_approx_epsilon_function_value
            , "" );
    
    }

    { //::SireSearch::get_min_protein_residues
    
        typedef int ( *get_min_protein_residues_function_type )(  );
        get_min_protein_residues_function_type get_min_protein_residues_function_value( &::SireSearch::get_min_protein_residues );
        
        bp::def( 
            "get_min_protein_residues"
            , get_min_protein_residues_function_value
            , "" );
    
    }

    { //::SireSearch::get_protein_residue_names
    
        typedef ::QSet< QString > ( *get_protein_residue_names_function_type )(  );
        get_protein_residue_names_function_type get_protein_residue_names_function_value( &::SireSearch::get_protein_residue_names );
        
        bp::def( 
            "get_protein_residue_names"
            , get_protein_residue_names_function_value
            , "" );
    
    }

    { //::SireSearch::get_token
    
        typedef ::QString ( *get_token_function_type )( ::QString const & );
        get_token_function_type get_token_function_value( &::SireSearch::get_token );
        
        bp::def( 
            "get_token"
            , get_token_function_value
            , ( bp::arg("token") )
            , "" );
    
    }

    { //::SireSearch::has_token
    
        typedef bool ( *has_token_function_type )( ::QString const & );
        has_token_function_type has_token_function_value( &::SireSearch::has_token );
        
        bp::def( 
            "has_token"
            , has_token_function_value
            , ( bp::arg("token") )
            , "" );
    
    }

    { //::SireSearch::install_search_parser
    
        typedef void ( *install_search_parser_function_type )(  );
        install_search_parser_function_type install_search_parser_function_value( &::SireSearch::install_search_parser );
        
        bp::def( 
            "install_search_parser"
            , install_search_parser_function_value
            , "" );
    
    }

    { //::SireSearch::set_approx_epsilon
    
        typedef void ( *set_approx_epsilon_function_type )( double );
        set_approx_epsilon_function_type set_approx_epsilon_function_value( &::SireSearch::set_approx_epsilon );
        
        bp::def( 
            "set_approx_epsilon"
            , set_approx_epsilon_function_value
            , ( bp::arg("eps") )
            , "" );
    
    }

    { //::SireSearch::set_min_protein_residues
    
        typedef void ( *set_min_protein_residues_function_type )( int );
        set_min_protein_residues_function_type set_min_protein_residues_function_value( &::SireSearch::set_min_protein_residues );
        
        bp::def( 
            "set_min_protein_residues"
            , set_min_protein_residues_function_value
            , ( bp::arg("nres") )
            , "" );
    
    }

    { //::SireSearch::set_protein_residue_names
    
        typedef void ( *set_protein_residue_names_function_type )( ::QSet< QString > const & );
        set_protein_residue_names_function_type set_protein_residue_names_function_value( &::SireSearch::set_protein_residue_names );
        
        bp::def( 
            "set_protein_residue_names"
            , set_protein_residue_names_function_value
            , ( bp::arg("names") )
            , "" );
    
    }

    { //::SireSearch::set_token
    
        typedef void ( *set_token_function_type )( ::QString const &,::QString const & );
        set_token_function_type set_token_function_value( &::SireSearch::set_token );
        
        bp::def( 
            "set_token"
            , set_token_function_value
            , ( bp::arg("token"), bp::arg("search") )
            , "" );
    
    }

}
