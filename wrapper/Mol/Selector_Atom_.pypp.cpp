// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Selector_Atom_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "atomcharges.h"

#include "atomeditor.h"

#include "atomproperty.hpp"

#include "chain.h"

#include "cutgroup.h"

#include "evaluator.h"

#include "molecule.h"

#include "mover.hpp"

#include "mover_metaid.h"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include <QDebug>

#include "atom.h"

SireMol::Selector<SireMol::Atom> __copy__(const SireMol::Selector<SireMol::Atom> &other){ return SireMol::Selector<SireMol::Atom>(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Selector_Atom__class(){

    { //::SireMol::Selector< SireMol::Atom >
        typedef bp::class_< SireMol::Selector< SireMol::Atom >, bp::bases< SireMol::MoleculeView, SireBase::Property > > Selector_Atom__exposer_t;
        Selector_Atom__exposer_t Selector_Atom__exposer = Selector_Atom__exposer_t( "Selector_Atom_", "", bp::init< >("") );
        bp::scope Selector_Atom__scope( Selector_Atom__exposer );
        Selector_Atom__exposer.def( bp::init< SireMol::Atom const & >(( bp::arg("view") ), "") );
        Selector_Atom__exposer.def( bp::init< SireMol::MoleculeData const & >(( bp::arg("moldata") ), "") );
        Selector_Atom__exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomSelection const & >(( bp::arg("moldata"), bp::arg("selected_atoms") ), "") );
        Selector_Atom__exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::Atom::ID const & >(( bp::arg("moldata"), bp::arg("viewid") ), "") );
        Selector_Atom__exposer.def( bp::init< SireMol::MoleculeData const &, QList< SireMol::AtomIdx > const & >(( bp::arg("moldata"), bp::arg("idxs") ), "") );
        Selector_Atom__exposer.def( bp::init< SireMol::Selector< SireMol::Atom > const & >(( bp::arg("other") ), "") );
        { //::SireMol::Selector< SireMol::Atom >::IDs
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::Selector< SireMol::Atom >::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireMol::Selector< SireMol::Atom >::IDs );
            
            Selector_Atom__exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::add
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*add_function_type)( ::SireMol::Selector< SireMol::Atom > const & ) const;
            add_function_type add_function_value( &::SireMol::Selector< SireMol::Atom >::add );
            
            Selector_Atom__exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::add
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*add_function_type)( ::SireMol::Atom const & ) const;
            add_function_type add_function_value( &::SireMol::Selector< SireMol::Atom >::add );
            
            Selector_Atom__exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::add
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*add_function_type)( ::SireMol::Atom::ID const & ) const;
            add_function_type add_function_value( &::SireMol::Selector< SireMol::Atom >::add );
            
            Selector_Atom__exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::contains
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*contains_function_type)( ::SireMol::Selector< SireMol::Atom > const & ) const;
            contains_function_type contains_function_value( &::SireMol::Selector< SireMol::Atom >::contains );
            
            Selector_Atom__exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::contains
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*contains_function_type)( ::SireMol::Atom const & ) const;
            contains_function_type contains_function_value( &::SireMol::Selector< SireMol::Atom >::contains );
            
            Selector_Atom__exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::contains
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*contains_function_type)( ::SireMol::Atom::ID const & ) const;
            contains_function_type contains_function_value( &::SireMol::Selector< SireMol::Atom >::contains );
            
            Selector_Atom__exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::evaluate
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Evaluator ( ::SireMol::Selector< SireMol::Atom >::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Selector< SireMol::Atom >::evaluate );
            
            Selector_Atom__exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::evaluate
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Evaluator ( ::SireMol::Selector< SireMol::Atom >::*evaluate_function_type)( int ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Selector< SireMol::Atom >::evaluate );
            
            Selector_Atom__exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::evaluate
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Evaluator ( ::SireMol::Selector< SireMol::Atom >::*evaluate_function_type)( int,int ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Selector< SireMol::Atom >::evaluate );
            
            Selector_Atom__exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::hasMetadata
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Selector< SireMol::Atom >::hasMetadata );
            
            Selector_Atom__exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::hasMetadata
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Selector< SireMol::Atom >::hasMetadata );
            
            Selector_Atom__exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::hasProperty
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::Selector< SireMol::Atom >::hasProperty );
            
            Selector_Atom__exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::index
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Atom::Index ( ::SireMol::Selector< SireMol::Atom >::*index_function_type)( int ) const;
            index_function_type index_function_value( &::SireMol::Selector< SireMol::Atom >::index );
            
            Selector_Atom__exposer.def( 
                "index"
                , index_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::indexes
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::Selector< SireMol::Atom >::*indexes_function_type)(  ) const;
            indexes_function_type indexes_function_value( &::SireMol::Selector< SireMol::Atom >::indexes );
            
            Selector_Atom__exposer.def( 
                "indexes"
                , indexes_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::intersection
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*intersection_function_type)( ::SireMol::Selector< SireMol::Atom > const & ) const;
            intersection_function_type intersection_function_value( &::SireMol::Selector< SireMol::Atom >::intersection );
            
            Selector_Atom__exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::intersection
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*intersection_function_type)( ::SireMol::Atom const & ) const;
            intersection_function_type intersection_function_value( &::SireMol::Selector< SireMol::Atom >::intersection );
            
            Selector_Atom__exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::intersection
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*intersection_function_type)( ::SireMol::Atom::ID const & ) const;
            intersection_function_type intersection_function_value( &::SireMol::Selector< SireMol::Atom >::intersection );
            
            Selector_Atom__exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::intersects
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*intersects_function_type)( ::SireMol::Selector< SireMol::Atom > const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Selector< SireMol::Atom >::intersects );
            
            Selector_Atom__exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::intersects
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*intersects_function_type)( ::SireMol::Atom const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Selector< SireMol::Atom >::intersects );
            
            Selector_Atom__exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::intersects
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*intersects_function_type)( ::SireMol::Atom::ID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Selector< SireMol::Atom >::intersects );
            
            Selector_Atom__exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::invert
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMol::Selector< SireMol::Atom >::invert );
            
            Selector_Atom__exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::isEmpty
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Selector< SireMol::Atom >::isEmpty );
            
            Selector_Atom__exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::metadataKeys
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::QStringList ( ::SireMol::Selector< SireMol::Atom >::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Selector< SireMol::Atom >::metadataKeys );
            
            Selector_Atom__exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::metadataKeys
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::QStringList ( ::SireMol::Selector< SireMol::Atom >::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Selector< SireMol::Atom >::metadataKeys );
            
            Selector_Atom__exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::move
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Atom > > ( ::SireMol::Selector< SireMol::Atom >::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::Selector< SireMol::Atom >::move );
            
            Selector_Atom__exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::move
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Atom > > ( ::SireMol::Selector< SireMol::Atom >::*move_function_type)( int ) const;
            move_function_type move_function_value( &::SireMol::Selector< SireMol::Atom >::move );
            
            Selector_Atom__exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::move
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Atom > > ( ::SireMol::Selector< SireMol::Atom >::*move_function_type)( int,int ) const;
            move_function_type move_function_value( &::SireMol::Selector< SireMol::Atom >::move );
            
            Selector_Atom__exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::nViews
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef int ( ::SireMol::Selector< SireMol::Atom >::*nViews_function_type)(  ) const;
            nViews_function_type nViews_function_value( &::SireMol::Selector< SireMol::Atom >::nViews );
            
            Selector_Atom__exposer.def( 
                "nViews"
                , nViews_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::names
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::QList< SireMol::AtomName > ( ::SireMol::Selector< SireMol::Atom >::*names_function_type)(  ) const;
            names_function_type names_function_value( &::SireMol::Selector< SireMol::Atom >::names );
            
            Selector_Atom__exposer.def( 
                "names"
                , names_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::numbers
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::QList< SireMol::AtomNum > ( ::SireMol::Selector< SireMol::Atom >::*numbers_function_type)(  ) const;
            numbers_function_type numbers_function_value( &::SireMol::Selector< SireMol::Atom >::numbers );
            
            Selector_Atom__exposer.def( 
                "numbers"
                , numbers_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Selector_Atom__exposer.def( bp::self != bp::self );
        { //::SireMol::Selector< SireMol::Atom >::operator()
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Atom ( ::SireMol::Selector< SireMol::Atom >::*__call___function_type)( int ) const;
            __call___function_type __call___function_value( &::SireMol::Selector< SireMol::Atom >::operator() );
            
            Selector_Atom__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator()
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireMol::Selector< SireMol::Atom >::operator() );
            
            Selector_Atom__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , "" );
        
        }
        Selector_Atom__exposer.def( bp::self + bp::self );
        Selector_Atom__exposer.def( bp::self + bp::other< SireMol::AtomID >() );
        Selector_Atom__exposer.def( bp::self + bp::other< SireMol::Atom >() );
        Selector_Atom__exposer.def( bp::self - bp::self );
        Selector_Atom__exposer.def( bp::self - bp::other< SireMol::AtomID >() );
        Selector_Atom__exposer.def( bp::self - bp::other< SireMol::Atom >() );
        { //::SireMol::Selector< SireMol::Atom >::operator=
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > & ( ::SireMol::Selector< SireMol::Atom >::*assign_function_type)( ::SireMol::Selector< SireMol::Atom > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Selector< SireMol::Atom >::operator= );
            
            Selector_Atom__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator=
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > & ( ::SireMol::Selector< SireMol::Atom >::*assign_function_type)( ::SireMol::Atom const & ) ;
            assign_function_type assign_function_value( &::SireMol::Selector< SireMol::Atom >::operator= );
            
            Selector_Atom__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("view") )
                , bp::return_self< >()
                , "" );
        
        }
        Selector_Atom__exposer.def( bp::self == bp::self );
        { //::SireMol::Selector< SireMol::Atom >::operator[]
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Atom >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Atom >::operator[] );
            
            Selector_Atom__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator[]
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Atom >::*__getitem___function_type)( ::QString const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Atom >::operator[] );
            
            Selector_Atom__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator[]
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Atom >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Atom >::operator[] );
            
            Selector_Atom__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator[]
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Atom >::*__getitem___function_type)( ::SireMol::AtomID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Atom >::operator[] );
            
            Selector_Atom__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("atomid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator[]
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Atom >::*__getitem___function_type)( ::SireMol::ResID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Atom >::operator[] );
            
            Selector_Atom__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("resid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator[]
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Atom >::*__getitem___function_type)( ::SireMol::CGID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Atom >::operator[] );
            
            Selector_Atom__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator[]
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Atom >::*__getitem___function_type)( ::SireMol::ChainID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Atom >::operator[] );
            
            Selector_Atom__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator[]
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Atom >::*__getitem___function_type)( ::SireMol::SegID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Atom >::operator[] );
            
            Selector_Atom__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("segid") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::operator[]
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Selector< SireMol::Atom >::*__getitem___function_type)( ::SireID::Index const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Atom >::operator[] );
            
            Selector_Atom__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idx") )
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::propertyKeys
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::QStringList ( ::SireMol::Selector< SireMol::Atom >::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::Selector< SireMol::Atom >::propertyKeys );
            
            Selector_Atom__exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::selectedAll
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Atom >::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::Selector< SireMol::Atom >::selectedAll );
            
            Selector_Atom__exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::selection
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::AtomSelection ( ::SireMol::Selector< SireMol::Atom >::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::Selector< SireMol::Atom >::selection );
            
            Selector_Atom__exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::selection
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::AtomSelection ( ::SireMol::Selector< SireMol::Atom >::*selection_function_type)( int ) const;
            selection_function_type selection_function_value( &::SireMol::Selector< SireMol::Atom >::selection );
            
            Selector_Atom__exposer.def( 
                "selection"
                , selection_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::selection
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::AtomSelection ( ::SireMol::Selector< SireMol::Atom >::*selection_function_type)( int,int ) const;
            selection_function_type selection_function_value( &::SireMol::Selector< SireMol::Atom >::selection );
            
            Selector_Atom__exposer.def( 
                "selection"
                , selection_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::selector
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*selector_function_type)(  ) const;
            selector_function_type selector_function_value( &::SireMol::Selector< SireMol::Atom >::selector );
            
            Selector_Atom__exposer.def( 
                "selector"
                , selector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::selector
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*selector_function_type)( int ) const;
            selector_function_type selector_function_value( &::SireMol::Selector< SireMol::Atom >::selector );
            
            Selector_Atom__exposer.def( 
                "selector"
                , selector_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::selector
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*selector_function_type)( int,int ) const;
            selector_function_type selector_function_value( &::SireMol::Selector< SireMol::Atom >::selector );
            
            Selector_Atom__exposer.def( 
                "selector"
                , selector_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::subtract
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*subtract_function_type)( ::SireMol::Selector< SireMol::Atom > const & ) const;
            subtract_function_type subtract_function_value( &::SireMol::Selector< SireMol::Atom >::subtract );
            
            Selector_Atom__exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::subtract
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*subtract_function_type)( ::SireMol::Atom const & ) const;
            subtract_function_type subtract_function_value( &::SireMol::Selector< SireMol::Atom >::subtract );
            
            Selector_Atom__exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("view") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::subtract
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::Selector< SireMol::Atom >::*subtract_function_type)( ::SireMol::Atom::ID const & ) const;
            subtract_function_type subtract_function_value( &::SireMol::Selector< SireMol::Atom >::subtract );
            
            Selector_Atom__exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::toString
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef ::QString ( ::SireMol::Selector< SireMol::Atom >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Selector< SireMol::Atom >::toString );
            
            Selector_Atom__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Selector< SireMol::Atom >::typeName
        
            typedef SireMol::Selector< SireMol::Atom > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Selector< SireMol::Atom >::typeName );
            
            Selector_Atom__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Selector_Atom__exposer.staticmethod( "typeName" );
        Selector_Atom__exposer.def( "__copy__", &__copy__);
        Selector_Atom__exposer.def( "__deepcopy__", &__copy__);
        Selector_Atom__exposer.def( "clone", &__copy__);
        Selector_Atom__exposer.def( "__str__", &__str__< ::SireMol::Selector<SireMol::Atom> > );
        Selector_Atom__exposer.def( "__repr__", &__str__< ::SireMol::Selector<SireMol::Atom> > );
        Selector_Atom__exposer.def( "__len__", &__len_size< ::SireMol::Selector<SireMol::Atom> > );
    }

}
