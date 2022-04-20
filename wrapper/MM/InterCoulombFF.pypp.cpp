// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "InterCoulombFF.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intercoulombff.h"

#include "intercoulombff.h"

SireFF::Inter2B3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> > __copy__(const SireFF::Inter2B3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> > &other){ return SireFF::Inter2B3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_InterCoulombFF_class(){

    { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >
        typedef bp::class_< SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >, bp::bases< SireFF::FF3D, SireFF::Inter2BFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >, SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential>, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > InterCoulombFF_exposer_t;
        InterCoulombFF_exposer_t InterCoulombFF_exposer = InterCoulombFF_exposer_t( "InterCoulombFF", "", bp::init< >("") );
        bp::scope InterCoulombFF_scope( InterCoulombFF_exposer );
        InterCoulombFF_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        InterCoulombFF_exposer.def( bp::init< SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*energy_function_type)(  ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::energy );
            
            InterCoulombFF_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*energy_function_type)( ::SireCAS::Symbol const & ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::energy );
            
            InterCoulombFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("component") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::energy );
            
            InterCoulombFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("scale_energy")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,::SireCAS::Symbol const &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::energy );
            
            InterCoulombFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("symbol"), bp::arg("scale_energy")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*field_function_type)( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::field );
            
            InterCoulombFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::field );
            
            InterCoulombFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::field );
            
            InterCoulombFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::field );
            
            InterCoulombFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::force
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*force_function_type)( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::force );
            
            InterCoulombFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::force
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*force_function_type)( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::force );
            
            InterCoulombFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        InterCoulombFF_exposer.def( bp::self != bp::self );
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::operator=
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > & ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*assign_function_type)( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::operator= );
            
            InterCoulombFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        InterCoulombFF_exposer.def( bp::self == bp::self );
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::packCoordinates
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*packCoordinates_function_type)(  ) ;
            packCoordinates_function_type packCoordinates_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::packCoordinates );
            
            InterCoulombFF_exposer.def( 
                "packCoordinates"
                , packCoordinates_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::potential );
            
            InterCoulombFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::potential );
            
            InterCoulombFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::potential );
            
            InterCoulombFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::potential );
            
            InterCoulombFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::typeName
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::typeName );
            
            InterCoulombFF_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::what
        
            typedef SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Inter2B3DFF< SireMM::CoulombPotentialInterface< SireMM::InterCoulombPotential > >::what );
            
            InterCoulombFF_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        InterCoulombFF_exposer.staticmethod( "typeName" );
        InterCoulombFF_exposer.def( "__copy__", &__copy__);
        InterCoulombFF_exposer.def( "__deepcopy__", &__copy__);
        InterCoulombFF_exposer.def( "clone", &__copy__);
        InterCoulombFF_exposer.def( "__str__", &__str__< ::SireFF::Inter2B3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> > > );
        InterCoulombFF_exposer.def( "__repr__", &__str__< ::SireFF::Inter2B3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> > > );
        InterCoulombFF_exposer.def( "__len__", &__len_count< ::SireFF::Inter2B3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> > > );
    }

}
