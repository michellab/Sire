// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "AnglePerturbation.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/identities.h"

#include "SireCAS/values.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "geometryperturbation.h"

#include "molecule.h"

#include "moleditor.h"

#include "mover.hpp"

#include "geometryperturbation.h"

SireMol::AnglePerturbation __copy__(const SireMol::AnglePerturbation &other){ return SireMol::AnglePerturbation(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_AnglePerturbation_class(){

    { //::SireMol::AnglePerturbation
        typedef bp::class_< SireMol::AnglePerturbation, bp::bases< SireMol::GeometryPerturbation, SireMol::Perturbation, SireBase::Property > > AnglePerturbation_exposer_t;
        AnglePerturbation_exposer_t AnglePerturbation_exposer = AnglePerturbation_exposer_t( "AnglePerturbation", "This perturbation moves an angle between two sizes.\n\nThis uses the anchors property to anchor parts\nof the molecule, the weight function property\nto weight the motion of the parts of the molecule,\nand the coordinates property to get the coordinates\nto move\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope AnglePerturbation_scope( AnglePerturbation_exposer );
        AnglePerturbation_exposer.def( bp::init< SireMol::AngleID const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("angle"), bp::arg("start"), bp::arg("end"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to perturb the angle angle from start to end") );
        AnglePerturbation_exposer.def( bp::init< SireMol::AngleID const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("angle"), bp::arg("start"), bp::arg("end"), bp::arg("mapping_function"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to perturb the angle angle from start to end\nusing the passed mapping function") );
        AnglePerturbation_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("start"), bp::arg("end"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to perturb the angle between atoms atom0, atom1 and atom2\nfrom start to end") );
        AnglePerturbation_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("start"), bp::arg("end"), bp::arg("mapping_function"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to perturb the angle between atoms atom0, atom1 and atom2\nfrom start to end using the passed mapping function") );
        AnglePerturbation_exposer.def( bp::init< SireMol::AnglePerturbation const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::AnglePerturbation::angle
        
            typedef ::SireMol::AngleID const & ( ::SireMol::AnglePerturbation::*angle_function_type)(  ) const;
            angle_function_type angle_function_value( &::SireMol::AnglePerturbation::angle );
            
            AnglePerturbation_exposer.def( 
                "angle"
                , angle_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the ID that identifies that angle that will be perturbed" );
        
        }
        { //::SireMol::AnglePerturbation::end
        
            typedef ::SireUnits::Dimension::Angle const & ( ::SireMol::AnglePerturbation::*end_function_type)(  ) const;
            end_function_type end_function_value( &::SireMol::AnglePerturbation::end );
            
            AnglePerturbation_exposer.def( 
                "end"
                , end_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the end length of the angle" );
        
        }
        AnglePerturbation_exposer.def( bp::self != bp::self );
        { //::SireMol::AnglePerturbation::operator=
        
            typedef ::SireMol::AnglePerturbation & ( ::SireMol::AnglePerturbation::*assign_function_type)( ::SireMol::AnglePerturbation const & ) ;
            assign_function_type assign_function_value( &::SireMol::AnglePerturbation::operator= );
            
            AnglePerturbation_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AnglePerturbation_exposer.def( bp::self == bp::self );
        { //::SireMol::AnglePerturbation::start
        
            typedef ::SireUnits::Dimension::Angle const & ( ::SireMol::AnglePerturbation::*start_function_type)(  ) const;
            start_function_type start_function_value( &::SireMol::AnglePerturbation::start );
            
            AnglePerturbation_exposer.def( 
                "start"
                , start_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the start length of the angle" );
        
        }
        { //::SireMol::AnglePerturbation::toString
        
            typedef ::QString ( ::SireMol::AnglePerturbation::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AnglePerturbation::toString );
            
            AnglePerturbation_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AnglePerturbation::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AnglePerturbation::typeName );
            
            AnglePerturbation_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AnglePerturbation::wouldChange
        
            typedef bool ( ::SireMol::AnglePerturbation::*wouldChange_function_type)( ::SireMol::Molecule const &,::SireCAS::Values const & ) const;
            wouldChange_function_type wouldChange_function_value( &::SireMol::AnglePerturbation::wouldChange );
            
            AnglePerturbation_exposer.def( 
                "wouldChange"
                , wouldChange_function_value
                , ( bp::arg("molecule"), bp::arg("values") )
                , bp::release_gil_policy()
                , "Return whether or not this perturbation with the passed values would\nchange the molecule molecule" );
        
        }
        AnglePerturbation_exposer.staticmethod( "typeName" );
        AnglePerturbation_exposer.def( "__copy__", &__copy__);
        AnglePerturbation_exposer.def( "__deepcopy__", &__copy__);
        AnglePerturbation_exposer.def( "clone", &__copy__);
        AnglePerturbation_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AnglePerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AnglePerturbation_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AnglePerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AnglePerturbation_exposer.def_pickle(sire_pickle_suite< ::SireMol::AnglePerturbation >());
        AnglePerturbation_exposer.def( "__str__", &__str__< ::SireMol::AnglePerturbation > );
        AnglePerturbation_exposer.def( "__repr__", &__str__< ::SireMol::AnglePerturbation > );
    }

}
