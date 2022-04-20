// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "FF3D.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/symbol.h"

#include "ff3d.h"

#include "ff3d.h"

const char* pvt_get_name(const SireFF::FF3D&){ return "SireFF::FF3D";}

#include "Helpers/release_gil_policy.hpp"

void register_FF3D_class(){

    { //::SireFF::FF3D
        typedef bp::class_< SireFF::FF3D, boost::noncopyable > FF3D_exposer_t;
        FF3D_exposer_t FF3D_exposer = FF3D_exposer_t( "FF3D", "This class provides the virtual interface for the 3D\nforcefields. These are forcefields that use 3D coordinates\nfor the atoms, and therefore you can calculate 3D forces\non the atoms. This class provides a virtual interface,\nand should be multiply inherited with FF to be used.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope FF3D_scope( FF3D_exposer );
        { //::SireFF::FF3D::energy
        
            typedef void ( ::SireFF::FF3D::*energy_function_type)( ::SireFF::EnergyTable &,double ) ;
            energy_function_type energy_function_value( &::SireFF::FF3D::energy );
            
            FF3D_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("scale_energy")=1 )
                , bp::release_gil_policy()
                , "Calculate all of the energies of the\nmolecules in the forcetable forcetable due to the\nmolecules interactions in this forcefield" );
        
        }
        { //::SireFF::FF3D::energy
        
            typedef void ( ::SireFF::FF3D::*energy_function_type)( ::SireFF::EnergyTable &,::SireCAS::Symbol const &,double ) ;
            energy_function_type energy_function_value( &::SireFF::FF3D::energy );
            
            FF3D_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("component"), bp::arg("scale_energy")=1 )
                , bp::release_gil_policy()
                , "Calculate all of the energies acting on all of the\nmolecules in the energytable energytable due to the\nspecified component of the molecules interactions in\nthis forcefield" );
        
        }
        { //::SireFF::FF3D::field
        
            typedef void ( ::SireFF::FF3D::*field_function_type)( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireFF::FF3D::field );
            
            FF3D_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 )
                , bp::release_gil_policy()
                , "Calculate the fields acting at all of the points\nin fieldtable due to the molecules in this forcefield" );
        
        }
        { //::SireFF::FF3D::field
        
            typedef void ( ::SireFF::FF3D::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireFF::FF3D::field );
            
            FF3D_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 )
                , bp::release_gil_policy()
                , "Calculate the fields acting at all of the points\nin fieldtable due to the specified component component\nfrom the molecules in this forcefield" );
        
        }
        { //::SireFF::FF3D::field
        
            typedef void ( ::SireFF::FF3D::*field_function_type)( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::FF3D::field );
            
            FF3D_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , bp::release_gil_policy()
                , "Calculate the fields acting at all of the points\nin fieldtable due to the molecules in this forcefield" );
        
        }
        { //::SireFF::FF3D::field
        
            typedef void ( ::SireFF::FF3D::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::FF3D::field );
            
            FF3D_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , bp::release_gil_policy()
                , "Calculate the fields acting at all of the points\nin fieldtable due to the specified component component\nfrom the molecules in this forcefield" );
        
        }
        { //::SireFF::FF3D::force
        
            typedef void ( ::SireFF::FF3D::*force_function_type)( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireFF::FF3D::force );
            
            FF3D_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , bp::release_gil_policy()
                , "Calculate all of the forces acting on all of the\nmolecules in the forcetable forcetable due to the\nmolecules interactions in this forcefield" );
        
        }
        { //::SireFF::FF3D::force
        
            typedef void ( ::SireFF::FF3D::*force_function_type)( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireFF::FF3D::force );
            
            FF3D_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("component"), bp::arg("scale_force")=1 )
                , bp::release_gil_policy()
                , "Calculate all of the forces acting on all of the\nmolecules in the forcetable forcetable due to the\nspecified component of the molecules interactions in\nthis forcefield" );
        
        }
        { //::SireFF::FF3D::potential
        
            typedef void ( ::SireFF::FF3D::*potential_function_type)( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireFF::FF3D::potential );
            
            FF3D_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 )
                , bp::release_gil_policy()
                , "Calculate the potential acting at the points in potentialtable\ndue to the molecules in this forcefield" );
        
        }
        { //::SireFF::FF3D::potential
        
            typedef void ( ::SireFF::FF3D::*potential_function_type)( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::FF3D::potential );
            
            FF3D_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , bp::release_gil_policy()
                , "Calculate the potential acting at the points in potentialtable\ndue to the molecules in this forcefield" );
        
        }
        { //::SireFF::FF3D::potential
        
            typedef void ( ::SireFF::FF3D::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::FF3D::potential );
            
            FF3D_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 )
                , bp::release_gil_policy()
                , "Calculate the potential acting at the points in potentialtable\ndue to the specified component component from the molecules in\nthis table" );
        
        }
        { //::SireFF::FF3D::potential
        
            typedef void ( ::SireFF::FF3D::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::FF3D::potential );
            
            FF3D_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , bp::release_gil_policy()
                , "Calculate the potential acting at the points in potentialtable\ndue to the specified component component from the molecules in\nthis table" );
        
        }
        { //::SireFF::FF3D::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::FF3D::typeName );
            
            FF3D_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        FF3D_exposer.staticmethod( "typeName" );
        FF3D_exposer.def( "__str__", &pvt_get_name);
        FF3D_exposer.def( "__repr__", &pvt_get_name);
    }

}
