// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "CLJPotentialInterface_InterCLJPotential_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intercljff.h"

#include "intercljff.h"

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::CLJPotentialInterface<SireMM::InterCLJPotential>&){ return "SireMM::CLJPotentialInterface<SireMM::InterCLJPotential>";}

#include "Helpers/release_gil_policy.hpp"

void register_CLJPotentialInterface_InterCLJPotential__class(){

    { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >
        typedef bp::class_< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >, boost::noncopyable > CLJPotentialInterface_InterCLJPotential__exposer_t;
        CLJPotentialInterface_InterCLJPotential__exposer_t CLJPotentialInterface_InterCLJPotential__exposer = CLJPotentialInterface_InterCLJPotential__exposer_t( "CLJPotentialInterface_InterCLJPotential_", "", bp::no_init );
        bp::scope CLJPotentialInterface_InterCLJPotential__scope( CLJPotentialInterface_InterCLJPotential__exposer );
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::combiningRules
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef ::QString const & ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*combiningRules_function_type)(  ) const;
            combiningRules_function_type combiningRules_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::combiningRules );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "combiningRules"
                , combiningRules_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::containsProperty
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::containsProperty );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::parameters
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef ::SireMM::InterCLJPotential::ParameterNames ( *parameters_function_type )(  );
            parameters_function_type parameters_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::parameters );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "parameters"
                , parameters_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::properties
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::properties );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::property
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::property );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::reactionFieldDielectric
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef double ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*reactionFieldDielectric_function_type)(  ) const;
            reactionFieldDielectric_function_type reactionFieldDielectric_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::reactionFieldDielectric );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "reactionFieldDielectric"
                , reactionFieldDielectric_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setCombiningRules
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*setCombiningRules_function_type)( ::QString const & ) ;
            setCombiningRules_function_type setCombiningRules_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setCombiningRules );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "setCombiningRules"
                , setCombiningRules_function_value
                , ( bp::arg("combiningrules") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setProperty
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setProperty );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setReactionFieldDielectric
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*setReactionFieldDielectric_function_type)( double ) ;
            setReactionFieldDielectric_function_type setReactionFieldDielectric_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setReactionFieldDielectric );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "setReactionFieldDielectric"
                , setReactionFieldDielectric_function_value
                , ( bp::arg("dielectric") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setShiftElectrostatics
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*setShiftElectrostatics_function_type)( bool ) ;
            setShiftElectrostatics_function_type setShiftElectrostatics_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setShiftElectrostatics );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "setShiftElectrostatics"
                , setShiftElectrostatics_function_value
                , ( bp::arg("switchelectro") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setSpace
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setSpace );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("new_space") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setSwitchingFunction
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*setSwitchingFunction_function_type)( ::SireMM::SwitchingFunction const & ) ;
            setSwitchingFunction_function_type setSwitchingFunction_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setSwitchingFunction );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "setSwitchingFunction"
                , setSwitchingFunction_function_value
                , ( bp::arg("new_switchfunc") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setUseAtomisticCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*setUseAtomisticCutoff_function_type)( bool ) ;
            setUseAtomisticCutoff_function_type setUseAtomisticCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setUseAtomisticCutoff );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "setUseAtomisticCutoff"
                , setUseAtomisticCutoff_function_value
                , ( bp::arg("switchatomistic") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setUseGroupCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*setUseGroupCutoff_function_type)( bool ) ;
            setUseGroupCutoff_function_type setUseGroupCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setUseGroupCutoff );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "setUseGroupCutoff"
                , setUseGroupCutoff_function_value
                , ( bp::arg("switchgroup") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setUseReactionField
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*setUseReactionField_function_type)( bool ) ;
            setUseReactionField_function_type setUseReactionField_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::setUseReactionField );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "setUseReactionField"
                , setUseReactionField_function_value
                , ( bp::arg("switchrf") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::shiftElectrostatics
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*shiftElectrostatics_function_type)(  ) const;
            shiftElectrostatics_function_type shiftElectrostatics_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::shiftElectrostatics );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "shiftElectrostatics"
                , shiftElectrostatics_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::space
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef ::SireVol::Space const & ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*space_function_type)(  ) const;
            space_function_type space_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::space );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::switchingFunction
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef ::SireMM::SwitchingFunction const & ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*switchingFunction_function_type)(  ) const;
            switchingFunction_function_type switchingFunction_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::switchingFunction );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "switchingFunction"
                , switchingFunction_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::useAtomisticCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*useAtomisticCutoff_function_type)(  ) const;
            useAtomisticCutoff_function_type useAtomisticCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::useAtomisticCutoff );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "useAtomisticCutoff"
                , useAtomisticCutoff_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::useGroupCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*useGroupCutoff_function_type)(  ) const;
            useGroupCutoff_function_type useGroupCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::useGroupCutoff );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "useGroupCutoff"
                , useGroupCutoff_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::useReactionField
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::*useReactionField_function_type)(  ) const;
            useReactionField_function_type useReactionField_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterCLJPotential >::useReactionField );
            
            CLJPotentialInterface_InterCLJPotential__exposer.def( 
                "useReactionField"
                , useReactionField_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CLJPotentialInterface_InterCLJPotential__exposer.staticmethod( "parameters" );
        CLJPotentialInterface_InterCLJPotential__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJPotentialInterface_InterCLJPotential__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJPotentialInterface_InterCLJPotential__exposer.def_pickle(sire_pickle_suite< ::SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >());
        CLJPotentialInterface_InterCLJPotential__exposer.def( "__str__", &pvt_get_name);
        CLJPotentialInterface_InterCLJPotential__exposer.def( "__repr__", &pvt_get_name);
    }

}
