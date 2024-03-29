// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "GetCOMPoint.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireMol/atom.h"

#include "SireMol/atomcoords.h"

#include "SireMol/evaluator.h"

#include "SireMol/moleculeview.h"

#include "SireMol/mover.hpp"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "getpoint.h"

#include <QDebug>

#include "getpoint.h"

SireMove::GetCOMPoint __copy__(const SireMove::GetCOMPoint &other){ return SireMove::GetCOMPoint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_GetCOMPoint_class(){

    { //::SireMove::GetCOMPoint
        typedef bp::class_< SireMove::GetCOMPoint, bp::bases< SireMove::GetPoint, SireBase::Property > > GetCOMPoint_exposer_t;
        GetCOMPoint_exposer_t GetCOMPoint_exposer = GetCOMPoint_exposer_t( "GetCOMPoint", "This function returns the center of mass (COG) of the\natoms in the passed view of the molecule\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope GetCOMPoint_scope( GetCOMPoint_exposer );
        GetCOMPoint_exposer.def( bp::init< SireMol::AtomID const & >(( bp::arg("atomid") ), "Construct to get the COM of the atoms in the molecule that match\nthe passed AtomID") );
        GetCOMPoint_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("atomid0"), bp::arg("atomid1") ), "Construct to get the COM of the atoms in the molecule that\nmatch either of the two passed AtomIDs") );
        GetCOMPoint_exposer.def( bp::init< QList< SireMol::AtomIdentifier > const & >(( bp::arg("atomids") ), "Construct to get the COM of the atoms in the molecule that\nmatch any of the passed AtomIDs") );
        GetCOMPoint_exposer.def( bp::init< SireMove::GetCOMPoint const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::GetCOMPoint::atomID
        
            typedef ::SireMol::AtomID const & ( ::SireMove::GetCOMPoint::*atomID_function_type)(  ) const;
            atomID_function_type atomID_function_value( &::SireMove::GetCOMPoint::atomID );
            
            GetCOMPoint_exposer.def( 
                "atomID"
                , atomID_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the AtomID(s) used to limit the search for the point" );
        
        }
        { //::SireMove::GetCOMPoint::getPoint
        
            typedef ::SireMaths::Vector ( ::SireMove::GetCOMPoint::*getPoint_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            getPoint_function_type getPoint_function_value( &::SireMove::GetCOMPoint::getPoint );
            
            GetCOMPoint_exposer.def( 
                "getPoint"
                , getPoint_function_value
                , ( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        GetCOMPoint_exposer.def( bp::self != bp::self );
        { //::SireMove::GetCOMPoint::operator=
        
            typedef ::SireMove::GetCOMPoint & ( ::SireMove::GetCOMPoint::*assign_function_type)( ::SireMove::GetCOMPoint const & ) ;
            assign_function_type assign_function_value( &::SireMove::GetCOMPoint::operator= );
            
            GetCOMPoint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        GetCOMPoint_exposer.def( bp::self == bp::self );
        { //::SireMove::GetCOMPoint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::GetCOMPoint::typeName );
            
            GetCOMPoint_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        GetCOMPoint_exposer.staticmethod( "typeName" );
        GetCOMPoint_exposer.def( "__copy__", &__copy__);
        GetCOMPoint_exposer.def( "__deepcopy__", &__copy__);
        GetCOMPoint_exposer.def( "clone", &__copy__);
        GetCOMPoint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::GetCOMPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GetCOMPoint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::GetCOMPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GetCOMPoint_exposer.def_pickle(sire_pickle_suite< ::SireMove::GetCOMPoint >());
        GetCOMPoint_exposer.def( "__str__", &__str__< ::SireMove::GetCOMPoint > );
        GetCOMPoint_exposer.def( "__repr__", &__str__< ::SireMove::GetCOMPoint > );
    }

}
