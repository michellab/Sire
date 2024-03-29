// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PS_GTO.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/array2d.hpp"

#include "SireError/errors.h"

#include "SireMaths/boys.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "pgto.h"

#include "pointcharge.h"

#include "pointdipole.h"

#include "pgto.h"

Squire::PS_GTO __copy__(const Squire::PS_GTO &other){ return Squire::PS_GTO(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_PS_GTO_class(){

    { //::Squire::PS_GTO
        typedef bp::class_< Squire::PS_GTO, bp::bases< Squire::GTOPair, Squire::ShellPair, SireBase::Property > > PS_GTO_exposer_t;
        PS_GTO_exposer_t PS_GTO_exposer = PS_GTO_exposer_t( "PS_GTO", "This is a combined S-P GTO shell pair", bp::init< >("Constructor") );
        bp::scope PS_GTO_scope( PS_GTO_exposer );
        PS_GTO_exposer.def( bp::init< SireMaths::Vector const &, Squire::S_GTO const &, SireMaths::Vector const &, Squire::P_GTO const & >(( bp::arg("A"), bp::arg("a"), bp::arg("B"), bp::arg("b") ), "Construct combining orbital a at position A with orbital b at\nposition B") );
        PS_GTO_exposer.def( bp::init< SireMaths::Vector const &, Squire::P_GTO const &, SireMaths::Vector const &, Squire::S_GTO const & >(( bp::arg("A"), bp::arg("a"), bp::arg("B"), bp::arg("b") ), "Construct combining orbital a at position A with orbital b at\nposition B") );
        PS_GTO_exposer.def( bp::init< Squire::PS_GTO const & >(( bp::arg("other") ), "Copy constructor") );
        { //::Squire::PS_GTO::P_minus_A
        
            typedef ::SireMaths::Vector const & ( ::Squire::PS_GTO::*P_minus_A_function_type)(  ) const;
            P_minus_A_function_type P_minus_A_function_value( &::Squire::PS_GTO::P_minus_A );
            
            PS_GTO_exposer.def( 
                "P_minus_A"
                , P_minus_A_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::Squire::PS_GTO::P_minus_B
        
            typedef ::SireMaths::Vector const & ( ::Squire::PS_GTO::*P_minus_B_function_type)(  ) const;
            P_minus_B_function_type P_minus_B_function_value( &::Squire::PS_GTO::P_minus_B );
            
            PS_GTO_exposer.def( 
                "P_minus_B"
                , P_minus_B_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::Squire::PS_GTO::Q_minus_C
        
            typedef ::SireMaths::Vector const & ( ::Squire::PS_GTO::*Q_minus_C_function_type)(  ) const;
            Q_minus_C_function_type Q_minus_C_function_value( &::Squire::PS_GTO::Q_minus_C );
            
            PS_GTO_exposer.def( 
                "Q_minus_C"
                , Q_minus_C_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::Squire::PS_GTO::Q_minus_D
        
            typedef ::SireMaths::Vector const & ( ::Squire::PS_GTO::*Q_minus_D_function_type)(  ) const;
            Q_minus_D_function_type Q_minus_D_function_value( &::Squire::PS_GTO::Q_minus_D );
            
            PS_GTO_exposer.def( 
                "Q_minus_D"
                , Q_minus_D_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::Squire::PS_GTO::angularMomentum0
        
            typedef int ( ::Squire::PS_GTO::*angularMomentum0_function_type)(  ) const;
            angularMomentum0_function_type angularMomentum0_function_value( &::Squire::PS_GTO::angularMomentum0 );
            
            PS_GTO_exposer.def( 
                "angularMomentum0"
                , angularMomentum0_function_value
                , bp::release_gil_policy()
                , "Return the angular momentum of the first GTO shell in this pair" );
        
        }
        { //::Squire::PS_GTO::angularMomentum1
        
            typedef int ( ::Squire::PS_GTO::*angularMomentum1_function_type)(  ) const;
            angularMomentum1_function_type angularMomentum1_function_value( &::Squire::PS_GTO::angularMomentum1 );
            
            PS_GTO_exposer.def( 
                "angularMomentum1"
                , angularMomentum1_function_value
                , bp::release_gil_policy()
                , "Return the angular momentum of the second GTO shell in this pair" );
        
        }
        { //::Squire::PS_GTO::nOrbitals0
        
            typedef int ( ::Squire::PS_GTO::*nOrbitals0_function_type)(  ) const;
            nOrbitals0_function_type nOrbitals0_function_value( &::Squire::PS_GTO::nOrbitals0 );
            
            PS_GTO_exposer.def( 
                "nOrbitals0"
                , nOrbitals0_function_value
                , bp::release_gil_policy()
                , "Return the number of orbitals in the first GTO shell in this pair" );
        
        }
        { //::Squire::PS_GTO::nOrbitals1
        
            typedef int ( ::Squire::PS_GTO::*nOrbitals1_function_type)(  ) const;
            nOrbitals1_function_type nOrbitals1_function_value( &::Squire::PS_GTO::nOrbitals1 );
            
            PS_GTO_exposer.def( 
                "nOrbitals1"
                , nOrbitals1_function_value
                , bp::release_gil_policy()
                , "Return the number of orbitals in the second GTO shell in this pair" );
        
        }
        PS_GTO_exposer.def( bp::self != bp::self );
        { //::Squire::PS_GTO::operator=
        
            typedef ::Squire::PS_GTO & ( ::Squire::PS_GTO::*assign_function_type)( ::Squire::PS_GTO const & ) ;
            assign_function_type assign_function_value( &::Squire::PS_GTO::operator= );
            
            PS_GTO_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PS_GTO_exposer.def( bp::self == bp::self );
        { //::Squire::PS_GTO::scale
        
            typedef double ( ::Squire::PS_GTO::*scale_function_type)(  ) const;
            scale_function_type scale_function_value( &::Squire::PS_GTO::scale );
            
            PS_GTO_exposer.def( 
                "scale"
                , scale_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::Squire::PS_GTO::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::PS_GTO::typeName );
            
            PS_GTO_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        PS_GTO_exposer.staticmethod( "typeName" );
        PS_GTO_exposer.def( "__copy__", &__copy__);
        PS_GTO_exposer.def( "__deepcopy__", &__copy__);
        PS_GTO_exposer.def( "clone", &__copy__);
        PS_GTO_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::PS_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PS_GTO_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::PS_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PS_GTO_exposer.def_pickle(sire_pickle_suite< ::Squire::PS_GTO >());
        PS_GTO_exposer.def( "__str__", &__str__< ::Squire::PS_GTO > );
        PS_GTO_exposer.def( "__repr__", &__str__< ::Squire::PS_GTO > );
    }

}
