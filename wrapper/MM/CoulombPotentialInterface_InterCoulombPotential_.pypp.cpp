// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "CoulombPotentialInterface_InterCoulombPotential_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intercoulombff.h"

#include "intercoulombff.h"

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential>&){ return "SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential>";}

void register_CoulombPotentialInterface_InterCoulombPotential__class(){

    { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >
        typedef bp::class_< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >, boost::noncopyable > CoulombPotentialInterface_InterCoulombPotential__exposer_t;
        CoulombPotentialInterface_InterCoulombPotential__exposer_t CoulombPotentialInterface_InterCoulombPotential__exposer = CoulombPotentialInterface_InterCoulombPotential__exposer_t( "CoulombPotentialInterface_InterCoulombPotential_", "", bp::no_init );
        bp::scope CoulombPotentialInterface_InterCoulombPotential__scope( CoulombPotentialInterface_InterCoulombPotential__exposer );
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::containsProperty
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef bool ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::containsProperty );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::parameters
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef ::SireMM::ChargeParameterName3D ( *parameters_function_type )(  );
            parameters_function_type parameters_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::parameters );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "parameters"
                , parameters_function_value
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::properties
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::properties );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::property
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::property );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::setProperty
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef bool ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::setProperty );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::setShiftElectrostatics
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef bool ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*setShiftElectrostatics_function_type)( bool ) ;
            setShiftElectrostatics_function_type setShiftElectrostatics_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::setShiftElectrostatics );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "setShiftElectrostatics"
                , setShiftElectrostatics_function_value
                , ( bp::arg("switchelectro") )
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::setSpace
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef bool ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::setSpace );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("new_space") )
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::setSwitchingFunction
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef bool ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*setSwitchingFunction_function_type)( ::SireMM::SwitchingFunction const & ) ;
            setSwitchingFunction_function_type setSwitchingFunction_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::setSwitchingFunction );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "setSwitchingFunction"
                , setSwitchingFunction_function_value
                , ( bp::arg("new_switchfunc") )
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::shiftElectrostatics
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef bool ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*shiftElectrostatics_function_type)(  ) const;
            shiftElectrostatics_function_type shiftElectrostatics_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::shiftElectrostatics );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "shiftElectrostatics"
                , shiftElectrostatics_function_value
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::space
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef ::SireVol::Space const & ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*space_function_type)(  ) const;
            space_function_type space_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::space );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::switchingFunction
        
            typedef SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > exported_class_t;
            typedef ::SireMM::SwitchingFunction const & ( ::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::*switchingFunction_function_type)(  ) const;
            switchingFunction_function_type switchingFunction_function_value( &::SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential >::switchingFunction );
            
            CoulombPotentialInterface_InterCoulombPotential__exposer.def( 
                "switchingFunction"
                , switchingFunction_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        CoulombPotentialInterface_InterCoulombPotential__exposer.staticmethod( "parameters" );
        CoulombPotentialInterface_InterCoulombPotential__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CoulombPotentialInterface_InterCoulombPotential__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CoulombPotentialInterface_InterCoulombPotential__exposer.def( "__str__", &pvt_get_name);
        CoulombPotentialInterface_InterCoulombPotential__exposer.def( "__repr__", &pvt_get_name);
    }

}
