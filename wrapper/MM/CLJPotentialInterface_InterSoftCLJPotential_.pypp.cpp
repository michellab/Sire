// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "CLJPotentialInterface_InterSoftCLJPotential_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intersoftcljff.h"

#include "intersoftcljff.h"

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential>&){ return "SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential>";}

#include "Helpers/release_gil_policy.hpp"

void register_CLJPotentialInterface_InterSoftCLJPotential__class(){

    { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >
        typedef bp::class_< SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >, boost::noncopyable > CLJPotentialInterface_InterSoftCLJPotential__exposer_t;
        CLJPotentialInterface_InterSoftCLJPotential__exposer_t CLJPotentialInterface_InterSoftCLJPotential__exposer = CLJPotentialInterface_InterSoftCLJPotential__exposer_t( "CLJPotentialInterface_InterSoftCLJPotential_", "", bp::no_init );
        bp::scope CLJPotentialInterface_InterSoftCLJPotential__scope( CLJPotentialInterface_InterSoftCLJPotential__exposer );
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::combiningRules
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef ::QString const & ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*combiningRules_function_type)(  ) const;
            combiningRules_function_type combiningRules_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::combiningRules );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "combiningRules"
                , combiningRules_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::containsProperty
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::containsProperty );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::parameters
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef ::SireMM::InterSoftCLJPotential::ParameterNames ( *parameters_function_type )(  );
            parameters_function_type parameters_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::parameters );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "parameters"
                , parameters_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::properties
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::properties );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::property
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::property );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::reactionFieldDielectric
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef double ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*reactionFieldDielectric_function_type)(  ) const;
            reactionFieldDielectric_function_type reactionFieldDielectric_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::reactionFieldDielectric );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "reactionFieldDielectric"
                , reactionFieldDielectric_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setCombiningRules
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setCombiningRules_function_type)( ::QString const & ) ;
            setCombiningRules_function_type setCombiningRules_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setCombiningRules );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setCombiningRules"
                , setCombiningRules_function_value
                , ( bp::arg("combiningrules") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setProperty
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setProperty );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setReactionFieldDielectric
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setReactionFieldDielectric_function_type)( double ) ;
            setReactionFieldDielectric_function_type setReactionFieldDielectric_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setReactionFieldDielectric );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setReactionFieldDielectric"
                , setReactionFieldDielectric_function_value
                , ( bp::arg("dielectric") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setShiftElectrostatics
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setShiftElectrostatics_function_type)( bool ) ;
            setShiftElectrostatics_function_type setShiftElectrostatics_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setShiftElectrostatics );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setShiftElectrostatics"
                , setShiftElectrostatics_function_value
                , ( bp::arg("switchelectro") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setSpace
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setSpace );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("new_space") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setSwitchingFunction
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setSwitchingFunction_function_type)( ::SireMM::SwitchingFunction const & ) ;
            setSwitchingFunction_function_type setSwitchingFunction_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setSwitchingFunction );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setSwitchingFunction"
                , setSwitchingFunction_function_value
                , ( bp::arg("new_switchfunc") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setUseAtomisticCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setUseAtomisticCutoff_function_type)( bool ) ;
            setUseAtomisticCutoff_function_type setUseAtomisticCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setUseAtomisticCutoff );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setUseAtomisticCutoff"
                , setUseAtomisticCutoff_function_value
                , ( bp::arg("switchatomistic") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setUseGroupCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setUseGroupCutoff_function_type)( bool ) ;
            setUseGroupCutoff_function_type setUseGroupCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setUseGroupCutoff );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setUseGroupCutoff"
                , setUseGroupCutoff_function_value
                , ( bp::arg("switchgroup") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setUseReactionField
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setUseReactionField_function_type)( bool ) ;
            setUseReactionField_function_type setUseReactionField_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::setUseReactionField );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setUseReactionField"
                , setUseReactionField_function_value
                , ( bp::arg("switchrf") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::shiftElectrostatics
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*shiftElectrostatics_function_type)(  ) const;
            shiftElectrostatics_function_type shiftElectrostatics_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::shiftElectrostatics );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "shiftElectrostatics"
                , shiftElectrostatics_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::space
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef ::SireVol::Space const & ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*space_function_type)(  ) const;
            space_function_type space_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::space );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::switchingFunction
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef ::SireMM::SwitchingFunction const & ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*switchingFunction_function_type)(  ) const;
            switchingFunction_function_type switchingFunction_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::switchingFunction );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "switchingFunction"
                , switchingFunction_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::useAtomisticCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*useAtomisticCutoff_function_type)(  ) const;
            useAtomisticCutoff_function_type useAtomisticCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::useAtomisticCutoff );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "useAtomisticCutoff"
                , useAtomisticCutoff_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::useGroupCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*useGroupCutoff_function_type)(  ) const;
            useGroupCutoff_function_type useGroupCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::useGroupCutoff );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "useGroupCutoff"
                , useGroupCutoff_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::useReactionField
        
            typedef SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::*useReactionField_function_type)(  ) const;
            useReactionField_function_type useReactionField_function_value( &::SireMM::CLJPotentialInterface< SireMM::InterSoftCLJPotential >::useReactionField );
            
            CLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "useReactionField"
                , useReactionField_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CLJPotentialInterface_InterSoftCLJPotential__exposer.staticmethod( "parameters" );
        CLJPotentialInterface_InterSoftCLJPotential__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJPotentialInterface_InterSoftCLJPotential__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJPotentialInterface_InterSoftCLJPotential__exposer.def_pickle(sire_pickle_suite< ::SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential> >());
        CLJPotentialInterface_InterSoftCLJPotential__exposer.def( "__str__", &pvt_get_name);
        CLJPotentialInterface_InterSoftCLJPotential__exposer.def( "__repr__", &pvt_get_name);
    }

}
