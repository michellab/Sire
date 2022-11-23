// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SphereArrayProperty.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/sphereproperty.h"

#include "sphereproperty.h"

#include "sphereproperty.h"

SireBase::PODArrayProperty<SireMaths::Sphere> __copy__(const SireBase::PODArrayProperty<SireMaths::Sphere> &other){ return SireBase::PODArrayProperty<SireMaths::Sphere>(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_SphereArrayProperty_class(){

    { //::SireBase::PODArrayProperty< SireMaths::Sphere >
        typedef bp::class_< SireBase::PODArrayProperty< SireMaths::Sphere >, bp::bases< SireBase::Property > > SphereArrayProperty_exposer_t;
        SphereArrayProperty_exposer_t SphereArrayProperty_exposer = SphereArrayProperty_exposer_t( "SphereArrayProperty", "", bp::init< >("") );
        bp::scope SphereArrayProperty_scope( SphereArrayProperty_exposer );
        SphereArrayProperty_exposer.def( bp::init< SireMaths::Sphere const & >(( bp::arg("value") ), "") );
        SphereArrayProperty_exposer.def( bp::init< QVector< SireMaths::Sphere > const & >(( bp::arg("values") ), "") );
        SphereArrayProperty_exposer.def( bp::init< QList< SireMaths::Sphere > const & >(( bp::arg("values") ), "") );
        SphereArrayProperty_exposer.def( bp::init< SireBase::PODProperty< SireMaths::Sphere > const & >(( bp::arg("value") ), "") );
        SphereArrayProperty_exposer.def( bp::init< QVector< SireBase::PODProperty< SireMaths::Sphere > > const & >(( bp::arg("values") ), "") );
        SphereArrayProperty_exposer.def( bp::init< QList< SireBase::PODProperty< SireMaths::Sphere > > const & >(( bp::arg("values") ), "") );
        SphereArrayProperty_exposer.def( bp::init< SireBase::PODArrayProperty< SireMaths::Sphere > const & >(( bp::arg("other") ), "") );
        SphereArrayProperty_exposer.def( bp::self != bp::self );
        { //::SireBase::PODArrayProperty< SireMaths::Sphere >::operator=
        
            typedef SireBase::PODArrayProperty< SireMaths::Sphere > exported_class_t;
            typedef ::SireBase::PODArrayProperty< SireMaths::Sphere > & ( ::SireBase::PODArrayProperty< SireMaths::Sphere >::*assign_function_type)( ::SireBase::PODArrayProperty< SireMaths::Sphere > const & ) ;
            assign_function_type assign_function_value( &::SireBase::PODArrayProperty< SireMaths::Sphere >::operator= );
            
            SphereArrayProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        SphereArrayProperty_exposer.def( bp::self == bp::self );
        { //::SireBase::PODArrayProperty< SireMaths::Sphere >::toString
        
            typedef SireBase::PODArrayProperty< SireMaths::Sphere > exported_class_t;
            typedef ::QString ( ::SireBase::PODArrayProperty< SireMaths::Sphere >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::PODArrayProperty< SireMaths::Sphere >::toString );
            
            SphereArrayProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::PODArrayProperty< SireMaths::Sphere >::typeName
        
            typedef SireBase::PODArrayProperty< SireMaths::Sphere > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::PODArrayProperty< SireMaths::Sphere >::typeName );
            
            SphereArrayProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::PODArrayProperty< SireMaths::Sphere >::what
        
            typedef SireBase::PODArrayProperty< SireMaths::Sphere > exported_class_t;
            typedef char const * ( ::SireBase::PODArrayProperty< SireMaths::Sphere >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireBase::PODArrayProperty< SireMaths::Sphere >::what );
            
            SphereArrayProperty_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SphereArrayProperty_exposer.staticmethod( "typeName" );
        SphereArrayProperty_exposer.def( "__copy__", &__copy__);
        SphereArrayProperty_exposer.def( "__deepcopy__", &__copy__);
        SphereArrayProperty_exposer.def( "clone", &__copy__);
        SphereArrayProperty_exposer.def( "__str__", &__str__< ::SireBase::PODArrayProperty<SireMaths::Sphere> > );
        SphereArrayProperty_exposer.def( "__repr__", &__str__< ::SireBase::PODArrayProperty<SireMaths::Sphere> > );
        SphereArrayProperty_exposer.def( "__len__", &__len_size< ::SireBase::PODArrayProperty<SireMaths::Sphere> > );
        SphereArrayProperty_exposer.def( "__getitem__", &::SireBase::PODArrayProperty<SireMaths::Sphere>::getitem );
    }

}