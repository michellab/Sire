// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "BondHunter.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/coordgroup.h"

#include "atom.h"

#include "atomcoords.h"

#include "atomelements.h"

#include "atomselection.h"

#include "bondhunter.h"

#include "connectivity.h"

#include "molecule.h"

#include "moleculedata.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "mover.hpp"

#include "selector.hpp"

#include <QDebug>

#include <QMutex>

#include "bondhunter.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_BondHunter_class(){

    { //::SireMol::BondHunter
        typedef bp::class_< SireMol::BondHunter, bp::bases< SireBase::Property >, boost::noncopyable > BondHunter_exposer_t;
        BondHunter_exposer_t BondHunter_exposer = BondHunter_exposer_t( "BondHunter", "Base class of all functions used to hunt for bonds in a molecule\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope BondHunter_scope( BondHunter_exposer );
        { //::SireMol::BondHunter::null
        
            typedef ::SireMol::NullBondHunter const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::BondHunter::null );
            
            BondHunter_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::BondHunter::operator()
        
            typedef ::SireMol::Connectivity ( ::SireMol::BondHunter::*__call___function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::BondHunter::operator() );
            
            BondHunter_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the connectivity of the molecule viewed in molview\nusing this function to hunt for all of the bonded atoms.\nThis only searches for the bonds between atoms that are\npart of this view" );
        
        }
        { //::SireMol::BondHunter::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::BondHunter::typeName );
            
            BondHunter_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        BondHunter_exposer.staticmethod( "null" );
        BondHunter_exposer.staticmethod( "typeName" );
        BondHunter_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::BondHunter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BondHunter_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::BondHunter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BondHunter_exposer.def_pickle(sire_pickle_suite< ::SireMol::BondHunter >());
        BondHunter_exposer.def( "__str__", &__str__< ::SireMol::BondHunter > );
        BondHunter_exposer.def( "__repr__", &__str__< ::SireMol::BondHunter > );
    }

}
