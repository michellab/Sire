// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomMatcher.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "atomidentifier.h"

#include "atomidx.h"

#include "atommatcher.h"

#include "atommatchers.h"

#include "atomname.h"

#include "atomselection.h"

#include "evaluator.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "mover.hpp"

#include "tostring.h"

#include "atommatcher.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_AtomMatcher_class(){

    { //::SireMol::AtomMatcher
        typedef bp::class_< SireMol::AtomMatcher, bp::bases< SireBase::Property >, boost::noncopyable > AtomMatcher_exposer_t;
        AtomMatcher_exposer_t AtomMatcher_exposer = AtomMatcher_exposer_t( "AtomMatcher", "Virtual base class of all of the functions used to match\natoms in one molecule layout with atoms in another layout\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope AtomMatcher_scope( AtomMatcher_exposer );
        { //::SireMol::AtomMatcher::add
        
            typedef ::SireMol::AtomMultiMatcher ( ::SireMol::AtomMatcher::*add_function_type)( ::SireMol::AtomMatcher const & ) const;
            add_function_type add_function_value( &::SireMol::AtomMatcher::add );
            
            AtomMatcher_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Return the matcher that matches using this matcher, and then other (in that order)" );
        
        }
        { //::SireMol::AtomMatcher::changesOrder
        
            typedef bool ( ::SireMol::AtomMatcher::*changesOrder_function_type)( ::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const & ) const;
            changesOrder_function_type changesOrder_function_value( &::SireMol::AtomMatcher::changesOrder );
            
            AtomMatcher_exposer.def( 
                "changesOrder"
                , changesOrder_function_value
                , ( bp::arg("molinfo0"), bp::arg("molinfo1") )
                , bp::release_gil_policy()
                , "Return whether or not this match changes the order of number of atoms" );
        
        }
        { //::SireMol::AtomMatcher::changesOrder
        
            typedef bool ( ::SireMol::AtomMatcher::*changesOrder_function_type)( ::SireMol::MoleculeView const &,::SireMol::MoleculeView const & ) const;
            changesOrder_function_type changesOrder_function_value( &::SireMol::AtomMatcher::changesOrder );
            
            AtomMatcher_exposer.def( 
                "changesOrder"
                , changesOrder_function_value
                , ( bp::arg("molview0"), bp::arg("molview1") )
                , bp::release_gil_policy()
                , "Return whether or not this match changes the order of number of atoms" );
        
        }
        { //::SireMol::AtomMatcher::changesOrder
        
            typedef bool ( ::SireMol::AtomMatcher::*changesOrder_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            changesOrder_function_type changesOrder_function_value( &::SireMol::AtomMatcher::changesOrder );
            
            AtomMatcher_exposer.def( 
                "changesOrder"
                , changesOrder_function_value
                , ( bp::arg("molview0"), bp::arg("map0"), bp::arg("molview1"), bp::arg("map1") )
                , bp::release_gil_policy()
                , "Return whether or not this match changes the order or number of viewed atoms" );
        
        }
        { //::SireMol::AtomMatcher::changesOrder
        
            typedef bool ( ::SireMol::AtomMatcher::*changesOrder_function_type)( ::SireMol::MoleculeView const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            changesOrder_function_type changesOrder_function_value( &::SireMol::AtomMatcher::changesOrder );
            
            AtomMatcher_exposer.def( 
                "changesOrder"
                , changesOrder_function_value
                , ( bp::arg("molview0"), bp::arg("molview1"), bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomMatcher::isNull
        
            typedef bool ( ::SireMol::AtomMatcher::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::AtomMatcher::isNull );
            
            AtomMatcher_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "Return whether or not this matcher is null (cannot be used for matching)" );
        
        }
        { //::SireMol::AtomMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomMatcher::*match_function_type)( ::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const & ) const;
            match_function_type match_function_value( &::SireMol::AtomMatcher::match );
            
            AtomMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molinfo0"), bp::arg("molinfo1") )
                , bp::release_gil_policy()
                , "Match atoms based only on the data in the MoleculeInfoData of the molecules." );
        
        }
        { //::SireMol::AtomMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomMatcher::*match_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            match_function_type match_function_value( &::SireMol::AtomMatcher::match );
            
            AtomMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molview0"), bp::arg("map0"), bp::arg("molview1"), bp::arg("map1") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomMatcher::*match_function_type)( ::SireMol::MoleculeView const &,::SireMol::MoleculeView const & ) const;
            match_function_type match_function_value( &::SireMol::AtomMatcher::match );
            
            AtomMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molview0"), bp::arg("molview1") )
                , bp::release_gil_policy()
                , "Match atoms based only on the data in the MoleculeInfoData of the molecules." );
        
        }
        { //::SireMol::AtomMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomMatcher::*match_function_type)( ::SireMol::MoleculeView const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            match_function_type match_function_value( &::SireMol::AtomMatcher::match );
            
            AtomMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molview0"), bp::arg("molview1"), bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomMatcher::null
        
            typedef ::SireMol::AtomMultiMatcher const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::AtomMatcher::null );
            
            AtomMatcher_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        AtomMatcher_exposer.def( bp::self + bp::self );
        { //::SireMol::AtomMatcher::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomMatcher::typeName );
            
            AtomMatcher_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomMatcher_exposer.staticmethod( "null" );
        AtomMatcher_exposer.staticmethod( "typeName" );
        AtomMatcher_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomMatcher_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomMatcher_exposer.def_pickle(sire_pickle_suite< ::SireMol::AtomMatcher >());
        AtomMatcher_exposer.def( "__str__", &__str__< ::SireMol::AtomMatcher > );
        AtomMatcher_exposer.def( "__repr__", &__str__< ::SireMol::AtomMatcher > );
    }

}
