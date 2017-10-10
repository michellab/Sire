// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "IntraGroupLJFFBase.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intraljff.h"

#include "intraljff.h"

SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > __copy__(const SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > &other){ return SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_IntraGroupLJFFBase_class(){

    { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >
        typedef bp::class_< SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >, bp::bases< SireMM::LJPotentialInterface<SireMM::IntraLJPotential>, SireFF::G2FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > IntraGroupLJFFBase_exposer_t;
        IntraGroupLJFFBase_exposer_t IntraGroupLJFFBase_exposer = IntraGroupLJFFBase_exposer_t( "IntraGroupLJFFBase", "", bp::init< >("") );
        bp::scope IntraGroupLJFFBase_scope( IntraGroupLJFFBase_exposer );
        IntraGroupLJFFBase_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        IntraGroupLJFFBase_exposer.def( bp::init< SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::components
        
            typedef SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef ::SireMM::LJComponent const & ( ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*components_function_type)(  ) const;
            components_function_type components_function_value( &::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::components );
            
            IntraGroupLJFFBase_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::containsProperty
        
            typedef SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::containsProperty );
            
            IntraGroupLJFFBase_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::mustNowRecalculateFromScratch
        
            typedef SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::mustNowRecalculateFromScratch );
            
            IntraGroupLJFFBase_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value
                , "" );
        
        }
        IntraGroupLJFFBase_exposer.def( bp::self != bp::self );
        { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::operator=
        
            typedef SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > & ( ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*assign_function_type)( ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::operator= );
            
            IntraGroupLJFFBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IntraGroupLJFFBase_exposer.def( bp::self == bp::self );
        { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::properties
        
            typedef SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::properties );
            
            IntraGroupLJFFBase_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::property
        
            typedef SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::property );
            
            IntraGroupLJFFBase_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::setProperty
        
            typedef SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::setProperty );
            
            IntraGroupLJFFBase_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") )
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::typeName
        
            typedef SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::typeName );
            
            IntraGroupLJFFBase_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::what
        
            typedef SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Intra2B2GFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::what );
            
            IntraGroupLJFFBase_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        IntraGroupLJFFBase_exposer.staticmethod( "typeName" );
        IntraGroupLJFFBase_exposer.def( "__copy__", &__copy__);
        IntraGroupLJFFBase_exposer.def( "__deepcopy__", &__copy__);
        IntraGroupLJFFBase_exposer.def( "clone", &__copy__);
        IntraGroupLJFFBase_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraGroupLJFFBase_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraGroupLJFFBase_exposer.def( "__str__", &__str__< ::SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > > );
        IntraGroupLJFFBase_exposer.def( "__repr__", &__str__< ::SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > > );
        IntraGroupLJFFBase_exposer.def( "__len__", &__len_count< ::SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > > );
    }

}
