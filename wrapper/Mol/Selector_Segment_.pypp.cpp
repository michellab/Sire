// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Selector_Segment_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "evaluator.h"

#include "groupatomids.h"

#include "molecule.h"

#include "mover.hpp"

#include "mover_metaid.h"

#include "segeditor.h"

#include "segment.h"

#include "selector.hpp"

#include "segment.h"

SireMol::Selector<SireMol::Segment> __copy__(const SireMol::Selector<SireMol::Segment> &other){ return SireMol::Selector<SireMol::Segment>(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Selector_Segment__class(){

    { //::SireMol::Selector< SireMol::Segment >
        typedef bp::class_< SireMol::Selector< SireMol::Segment >, bp::bases< SireMol::MoleculeView, SireBase::Property > > Selector_Segment__exposer_t;
        Selector_Segment__exposer_t Selector_Segment__exposer = Selector_Segment__exposer_t( "Selector_Segment_", "", bp::init< >("") );
        bp::scope Selector_Segment__scope( Selector_Segment__exposer );
        Selector_Segment__exposer.def( bp::init< SireMol::Segment const & >(( bp::arg("view") ), "") );
        Selector_Segment__exposer.def( bp::init< SireMol::MoleculeData const & >(( bp::arg("moldata") ), "") );
        Selector_Segment__exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomSelection const & >(( bp::arg("moldata"), bp::arg("selected_atoms") ), "") );
        Selector_Segment__exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::Segment::ID const & >(( bp::arg("moldata"), bp::arg("viewid") ), "") );
        Selector_Segment__exposer.def( bp::init< SireMol::MoleculeData const &, QList< SireMol::SegIdx > const & >(( bp::arg("moldata"), bp::arg("idxs") ), "") );
        Selector_Segment__exposer.def( bp::init< SireMol::Selector< SireMol::Segment > const & >(( bp::arg("other") ), "") );
        { //::SireMol::Selector< SireMol::Segment >::add
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*add_function_type)( ::SireMol::Selector< SireMol::Segment > const & ) const;
            add_function_type add_function_value( &::SireMol::Selector< SireMol::Segment >::add );
            
            Selector_Segment__exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::add
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*add_function_type)( ::SireMol::Segment const & ) const;
            add_function_type add_function_value( &::SireMol::Selector< SireMol::Segment >::add );
            
            Selector_Segment__exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::add
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*add_function_type)( ::SireMol::Segment::ID const & ) const;
            add_function_type add_function_value( &::SireMol::Selector< SireMol::Segment >::add );
            
            Selector_Segment__exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::contains
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*contains_function_type)( ::SireMol::Selector< SireMol::Segment > const & ) const;
            contains_function_type contains_function_value( &::SireMol::Selector< SireMol::Segment >::contains );
            
            Selector_Segment__exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::contains
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*contains_function_type)( ::SireMol::Segment const & ) const;
            contains_function_type contains_function_value( &::SireMol::Selector< SireMol::Segment >::contains );
            
            Selector_Segment__exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::contains
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*contains_function_type)( ::SireMol::Segment::ID const & ) const;
            contains_function_type contains_function_value( &::SireMol::Selector< SireMol::Segment >::contains );
            
            Selector_Segment__exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::evaluate
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Evaluator ( ::SireMol::Selector< SireMol::Segment >::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Selector< SireMol::Segment >::evaluate );
            
            Selector_Segment__exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::evaluate
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Evaluator ( ::SireMol::Selector< SireMol::Segment >::*evaluate_function_type)( int ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Selector< SireMol::Segment >::evaluate );
            
            Selector_Segment__exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::evaluate
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Evaluator ( ::SireMol::Selector< SireMol::Segment >::*evaluate_function_type)( int,int ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Selector< SireMol::Segment >::evaluate );
            
            Selector_Segment__exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::hasMetadata
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Selector< SireMol::Segment >::hasMetadata );
            
            Selector_Segment__exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::hasMetadata
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Selector< SireMol::Segment >::hasMetadata );
            
            Selector_Segment__exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::hasProperty
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::Selector< SireMol::Segment >::hasProperty );
            
            Selector_Segment__exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::index
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Segment::Index ( ::SireMol::Selector< SireMol::Segment >::*index_function_type)( int ) const;
            index_function_type index_function_value( &::SireMol::Selector< SireMol::Segment >::index );
            
            Selector_Segment__exposer.def( 
                "index"
                , index_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::intersection
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*intersection_function_type)( ::SireMol::Selector< SireMol::Segment > const & ) const;
            intersection_function_type intersection_function_value( &::SireMol::Selector< SireMol::Segment >::intersection );
            
            Selector_Segment__exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::intersection
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*intersection_function_type)( ::SireMol::Segment const & ) const;
            intersection_function_type intersection_function_value( &::SireMol::Selector< SireMol::Segment >::intersection );
            
            Selector_Segment__exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::intersection
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*intersection_function_type)( ::SireMol::Segment::ID const & ) const;
            intersection_function_type intersection_function_value( &::SireMol::Selector< SireMol::Segment >::intersection );
            
            Selector_Segment__exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::intersects
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*intersects_function_type)( ::SireMol::Selector< SireMol::Segment > const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Selector< SireMol::Segment >::intersects );
            
            Selector_Segment__exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::intersects
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*intersects_function_type)( ::SireMol::Segment const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Selector< SireMol::Segment >::intersects );
            
            Selector_Segment__exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::intersects
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*intersects_function_type)( ::SireMol::Segment::ID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Selector< SireMol::Segment >::intersects );
            
            Selector_Segment__exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::invert
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMol::Selector< SireMol::Segment >::invert );
            
            Selector_Segment__exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::isEmpty
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Selector< SireMol::Segment >::isEmpty );
            
            Selector_Segment__exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::metadataKeys
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::QStringList ( ::SireMol::Selector< SireMol::Segment >::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Selector< SireMol::Segment >::metadataKeys );
            
            Selector_Segment__exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::metadataKeys
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::QStringList ( ::SireMol::Selector< SireMol::Segment >::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Selector< SireMol::Segment >::metadataKeys );
            
            Selector_Segment__exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::move
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Segment > > ( ::SireMol::Selector< SireMol::Segment >::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::Selector< SireMol::Segment >::move );
            
            Selector_Segment__exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::move
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Segment > > ( ::SireMol::Selector< SireMol::Segment >::*move_function_type)( int ) const;
            move_function_type move_function_value( &::SireMol::Selector< SireMol::Segment >::move );
            
            Selector_Segment__exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::move
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Segment > > ( ::SireMol::Selector< SireMol::Segment >::*move_function_type)( int,int ) const;
            move_function_type move_function_value( &::SireMol::Selector< SireMol::Segment >::move );
            
            Selector_Segment__exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::nViews
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef int ( ::SireMol::Selector< SireMol::Segment >::*nViews_function_type)(  ) const;
            nViews_function_type nViews_function_value( &::SireMol::Selector< SireMol::Segment >::nViews );
            
            Selector_Segment__exposer.def( 
                "nViews"
                , nViews_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Selector_Segment__exposer.def( bp::self != bp::self );
        { //::SireMol::Selector< SireMol::Segment >::operator()
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Segment ( ::SireMol::Selector< SireMol::Segment >::*__call___function_type)( int ) const;
            __call___function_type __call___function_value( &::SireMol::Selector< SireMol::Segment >::operator() );
            
            Selector_Segment__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator()
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireMol::Selector< SireMol::Segment >::operator() );
            
            Selector_Segment__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , "" );
        
        }
        Selector_Segment__exposer.def( bp::self + bp::self );
        Selector_Segment__exposer.def( bp::self + bp::other< SireMol::SegID >() );
        Selector_Segment__exposer.def( bp::self + bp::other< SireMol::Segment >() );
        Selector_Segment__exposer.def( bp::self - bp::self );
        Selector_Segment__exposer.def( bp::self - bp::other< SireMol::SegID >() );
        Selector_Segment__exposer.def( bp::self - bp::other< SireMol::Segment >() );
        { //::SireMol::Selector< SireMol::Segment >::operator=
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > & ( ::SireMol::Selector< SireMol::Segment >::*assign_function_type)( ::SireMol::Selector< SireMol::Segment > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Selector< SireMol::Segment >::operator= );
            
            Selector_Segment__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator=
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > & ( ::SireMol::Selector< SireMol::Segment >::*assign_function_type)( ::SireMol::Segment const & ) ;
            assign_function_type assign_function_value( &::SireMol::Selector< SireMol::Segment >::operator= );
            
            Selector_Segment__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("view") )
                , bp::return_self< >()
                , "" );
        
        }
        Selector_Segment__exposer.def( bp::self == bp::self );
        { //::SireMol::Selector< SireMol::Segment >::operator[]
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Segment >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Segment >::operator[] );
            
            Selector_Segment__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator[]
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Segment >::*__getitem___function_type)( ::QString const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Segment >::operator[] );
            
            Selector_Segment__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator[]
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Segment >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Segment >::operator[] );
            
            Selector_Segment__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator[]
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Segment >::*__getitem___function_type)( ::SireMol::AtomID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Segment >::operator[] );
            
            Selector_Segment__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("atomid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator[]
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Segment >::*__getitem___function_type)( ::SireMol::ResID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Segment >::operator[] );
            
            Selector_Segment__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("resid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator[]
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Segment >::*__getitem___function_type)( ::SireMol::CGID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Segment >::operator[] );
            
            Selector_Segment__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator[]
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Segment >::*__getitem___function_type)( ::SireMol::ChainID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Segment >::operator[] );
            
            Selector_Segment__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator[]
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Segment >::*__getitem___function_type)( ::SireMol::SegID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Segment >::operator[] );
            
            Selector_Segment__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("segid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::operator[]
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Segment >::*__getitem___function_type)( ::SireID::Index const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Segment >::operator[] );
            
            Selector_Segment__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idx") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::propertyKeys
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::QStringList ( ::SireMol::Selector< SireMol::Segment >::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::Selector< SireMol::Segment >::propertyKeys );
            
            Selector_Segment__exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::selectedAll
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Segment >::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::Selector< SireMol::Segment >::selectedAll );
            
            Selector_Segment__exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::selection
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::AtomSelection ( ::SireMol::Selector< SireMol::Segment >::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::Selector< SireMol::Segment >::selection );
            
            Selector_Segment__exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::selection
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::AtomSelection ( ::SireMol::Selector< SireMol::Segment >::*selection_function_type)( int ) const;
            selection_function_type selection_function_value( &::SireMol::Selector< SireMol::Segment >::selection );
            
            Selector_Segment__exposer.def( 
                "selection"
                , selection_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::selection
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::AtomSelection ( ::SireMol::Selector< SireMol::Segment >::*selection_function_type)( int,int ) const;
            selection_function_type selection_function_value( &::SireMol::Selector< SireMol::Segment >::selection );
            
            Selector_Segment__exposer.def( 
                "selection"
                , selection_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::selector
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*selector_function_type)(  ) const;
            selector_function_type selector_function_value( &::SireMol::Selector< SireMol::Segment >::selector );
            
            Selector_Segment__exposer.def( 
                "selector"
                , selector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::selector
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*selector_function_type)( int ) const;
            selector_function_type selector_function_value( &::SireMol::Selector< SireMol::Segment >::selector );
            
            Selector_Segment__exposer.def( 
                "selector"
                , selector_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::selector
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*selector_function_type)( int,int ) const;
            selector_function_type selector_function_value( &::SireMol::Selector< SireMol::Segment >::selector );
            
            Selector_Segment__exposer.def( 
                "selector"
                , selector_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::subtract
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*subtract_function_type)( ::SireMol::Selector< SireMol::Segment > const & ) const;
            subtract_function_type subtract_function_value( &::SireMol::Selector< SireMol::Segment >::subtract );
            
            Selector_Segment__exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::subtract
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*subtract_function_type)( ::SireMol::Segment const & ) const;
            subtract_function_type subtract_function_value( &::SireMol::Selector< SireMol::Segment >::subtract );
            
            Selector_Segment__exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::subtract
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Selector< SireMol::Segment >::*subtract_function_type)( ::SireMol::Segment::ID const & ) const;
            subtract_function_type subtract_function_value( &::SireMol::Selector< SireMol::Segment >::subtract );
            
            Selector_Segment__exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::toString
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef ::QString ( ::SireMol::Selector< SireMol::Segment >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Selector< SireMol::Segment >::toString );
            
            Selector_Segment__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Segment >::typeName
        
            typedef SireMol::Selector< SireMol::Segment > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Selector< SireMol::Segment >::typeName );
            
            Selector_Segment__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Selector_Segment__exposer.staticmethod( "typeName" );
        Selector_Segment__exposer.def( "__copy__", &__copy__);
        Selector_Segment__exposer.def( "__deepcopy__", &__copy__);
        Selector_Segment__exposer.def( "clone", &__copy__);
        Selector_Segment__exposer.def( "__str__", &__str__< ::SireMol::Selector<SireMol::Segment> > );
        Selector_Segment__exposer.def( "__repr__", &__str__< ::SireMol::Selector<SireMol::Segment> > );
        Selector_Segment__exposer.def( "__len__", &__len_size< ::SireMol::Selector<SireMol::Segment> > );
    }

}
