// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CGEditorBase.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "atom.h"

#include "atomeditor.h"

#include "cgeditor.h"

#include "chain.h"

#include "chaineditor.h"

#include "cutgroup.h"

#include "groupatomids.h"

#include "molecule.h"

#include "moleditor.h"

#include "mover.hpp"

#include "reseditor.h"

#include "residue.h"

#include "segeditor.h"

#include "segment.h"

#include "selector.hpp"

#include "cgeditor.h"

#include "cgproperty.hpp"

SireMol::CGEditorBase& set_Metadata_SireMol_CGStringProperty_function1(
                                  SireMol::CGEditorBase &molview,
                                   const QString &metakey, const QString &p)
                                   { return molview.setMetadata< QString >(metakey, p); }

SireMol::CGEditorBase& set_Metadata_SireMol_CGStringProperty_function2(
                                  SireMol::CGEditorBase &molview,
                                   const QString &key, const QString &metakey, const QString &p)
                                   { return molview.setMetadata< QString >(key, metakey, p); }

SireMol::CGEditorBase& set_Metadata_SireMol_CGIntProperty_function1(
                                  SireMol::CGEditorBase &molview,
                                   const QString &metakey, const qint64 &p)
                                   { return molview.setMetadata< qint64 >(metakey, p); }

SireMol::CGEditorBase& set_Metadata_SireMol_CGIntProperty_function2(
                                  SireMol::CGEditorBase &molview,
                                   const QString &key, const QString &metakey, const qint64 &p)
                                   { return molview.setMetadata< qint64 >(key, metakey, p); }

SireMol::CGEditorBase& set_Metadata_SireMol_CGFloatProperty_function1(
                                  SireMol::CGEditorBase &molview,
                                   const QString &metakey, const double &p)
                                   { return molview.setMetadata< double >(metakey, p); }

SireMol::CGEditorBase& set_Metadata_SireMol_CGFloatProperty_function2(
                                  SireMol::CGEditorBase &molview,
                                   const QString &key, const QString &metakey, const double &p)
                                   { return molview.setMetadata< double >(key, metakey, p); }

SireMol::CGEditorBase& set_Metadata_SireMol_CGVariantProperty_function1(
                                  SireMol::CGEditorBase &molview,
                                   const QString &metakey, const QVariant &p)
                                   { return molview.setMetadata< QVariant >(metakey, p); }

SireMol::CGEditorBase& set_Metadata_SireMol_CGVariantProperty_function2(
                                  SireMol::CGEditorBase &molview,
                                   const QString &key, const QString &metakey, const QVariant &p)
                                   { return molview.setMetadata< QVariant >(key, metakey, p); }

SireMol::CGEditorBase& set_Metadata_SireMol_CGPropertyProperty_function1(
                                  SireMol::CGEditorBase &molview,
                                   const QString &metakey, const SireBase::PropertyPtr &p)
                                   { return molview.setMetadata< SireBase::PropertyPtr >(metakey, p); }

SireMol::CGEditorBase& set_Metadata_SireMol_CGPropertyProperty_function2(
                                  SireMol::CGEditorBase &molview,
                                   const QString &key, const QString &metakey, const SireBase::PropertyPtr &p)
                                   { return molview.setMetadata< SireBase::PropertyPtr >(key, metakey, p); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_CGEditorBase_class(){

    { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >
        typedef bp::class_< SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >, bp::bases< SireMol::CutGroup, SireMol::MoleculeView, SireBase::Property >, boost::noncopyable > CGEditorBase_exposer_t;
        CGEditorBase_exposer_t CGEditorBase_exposer = CGEditorBase_exposer_t( "CGEditorBase", "", bp::no_init );
        bp::scope CGEditorBase_scope( CGEditorBase_exposer );
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::atom
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*atom_function_type)(  ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::atom );
            
            CGEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::atom
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*atom_function_type)( int,::SireBase::PropertyMap const & ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::atom );
            
            CGEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::atom
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*atom_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::atom );
            
            CGEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::atom
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*atom_function_type)( ::SireMol::AtomID const &,::SireBase::PropertyMap const & ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::atom );
            
            CGEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("atomid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::chain
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*chain_function_type)(  ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::chain );
            
            CGEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::chain
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*chain_function_type)( int,::SireBase::PropertyMap const & ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::chain );
            
            CGEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::chain
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*chain_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::chain );
            
            CGEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::chain
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*chain_function_type)( ::SireMol::ChainID const &,::SireBase::PropertyMap const & ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::chain );
            
            CGEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("chainid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::cutGroup
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*cutGroup_function_type)(  ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::cutGroup );
            
            CGEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::cutGroup
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*cutGroup_function_type)( int,::SireBase::PropertyMap const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::cutGroup );
            
            CGEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::cutGroup
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*cutGroup_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::cutGroup );
            
            CGEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::cutGroup
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*cutGroup_function_type)( ::SireMol::CGID const &,::SireBase::PropertyMap const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::cutGroup );
            
            CGEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("cgid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::molecule
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::MolEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*molecule_function_type)(  ) ;
            molecule_function_type molecule_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::molecule );
            
            CGEditorBase_exposer.def( 
                "molecule"
                , molecule_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator=
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > & ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*assign_function_type)( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator= );
            
            CGEditorBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator=
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > & ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*assign_function_type)( ::SireMol::CutGroup const & ) ;
            assign_function_type assign_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator= );
            
            CGEditorBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[]
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*__getitem___function_type)( int ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[] );
            
            CGEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[]
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*__getitem___function_type)( ::QString const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[] );
            
            CGEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[]
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*__getitem___function_type)( ::SireMol::AtomID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[] );
            
            CGEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("atomid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[]
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*__getitem___function_type)( ::SireMol::ResID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[] );
            
            CGEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("resid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[]
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*__getitem___function_type)( ::SireMol::CGID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[] );
            
            CGEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[]
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*__getitem___function_type)( ::SireMol::ChainID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[] );
            
            CGEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[]
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*__getitem___function_type)( ::SireMol::SegID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[] );
            
            CGEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("segid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[]
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*__getitem___function_type)( ::SireID::Index const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::operator[] );
            
            CGEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idx") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::removeMetadata
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::CGEditor & ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*removeMetadata_function_type)( ::SireBase::PropertyName const & ) ;
            removeMetadata_function_type removeMetadata_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::removeMetadata );
            
            CGEditorBase_exposer.def( 
                "removeMetadata"
                , removeMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::removeMetadata
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::CGEditor & ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*removeMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) ;
            removeMetadata_function_type removeMetadata_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::removeMetadata );
            
            CGEditorBase_exposer.def( 
                "removeMetadata"
                , removeMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::removeProperty
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::CGEditor & ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*removeProperty_function_type)( ::SireBase::PropertyName const & ) ;
            removeProperty_function_type removeProperty_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::removeProperty );
            
            CGEditorBase_exposer.def( 
                "removeProperty"
                , removeProperty_function_value
                , ( bp::arg("key") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::residue
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*residue_function_type)(  ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::residue );
            
            CGEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::residue
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*residue_function_type)( int,::SireBase::PropertyMap const & ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::residue );
            
            CGEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::residue
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*residue_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::residue );
            
            CGEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::residue
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*residue_function_type)( ::SireMol::ResID const &,::SireBase::PropertyMap const & ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::residue );
            
            CGEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("resid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::segment
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*segment_function_type)(  ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::segment );
            
            CGEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::segment
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*segment_function_type)( int,::SireBase::PropertyMap const & ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::segment );
            
            CGEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::segment
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*segment_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::segment );
            
            CGEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::segment
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*segment_function_type)( ::SireMol::SegID const &,::SireBase::PropertyMap const & ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::segment );
            
            CGEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("segid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*select_function_type)( ::SireMol::AtomID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select );
            
            CGEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("atomid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*select_function_type)( ::SireMol::CGID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select );
            
            CGEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("cgid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*select_function_type)( ::SireMol::ResID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select );
            
            CGEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("resid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*select_function_type)( ::SireMol::ChainID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select );
            
            CGEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("chainid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select
        
            typedef SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::*select_function_type)( ::SireMol::SegID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::CGEditor, SireMol::CutGroup >::select );
            
            CGEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("segid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        CGEditorBase_exposer.def( "_set_property_QString",
                                           &SireMol::CGEditorBase::setProperty< QString >, bp::return_self< >() );
        CGEditorBase_exposer.def( "_set_metadata_QString", &set_Metadata_SireMol_CGStringProperty_function1, bp::return_self< >());
        CGEditorBase_exposer.def( "_set_metadata_QString", &set_Metadata_SireMol_CGStringProperty_function2, bp::return_self< >());
        CGEditorBase_exposer.def( "_set_property_qint64",
                                           &SireMol::CGEditorBase::setProperty< qint64 >, bp::return_self< >() );
        CGEditorBase_exposer.def( "_set_metadata_qint64", &set_Metadata_SireMol_CGIntProperty_function1, bp::return_self< >());
        CGEditorBase_exposer.def( "_set_metadata_qint64", &set_Metadata_SireMol_CGIntProperty_function2, bp::return_self< >());
        CGEditorBase_exposer.def( "_set_property_double",
                                           &SireMol::CGEditorBase::setProperty< double >, bp::return_self< >() );
        CGEditorBase_exposer.def( "_set_metadata_double", &set_Metadata_SireMol_CGFloatProperty_function1, bp::return_self< >());
        CGEditorBase_exposer.def( "_set_metadata_double", &set_Metadata_SireMol_CGFloatProperty_function2, bp::return_self< >());
        CGEditorBase_exposer.def( "_set_property_QVariant",
                                           &SireMol::CGEditorBase::setProperty< QVariant >, bp::return_self< >() );
        CGEditorBase_exposer.def( "_set_metadata_QVariant", &set_Metadata_SireMol_CGVariantProperty_function1, bp::return_self< >());
        CGEditorBase_exposer.def( "_set_metadata_QVariant", &set_Metadata_SireMol_CGVariantProperty_function2, bp::return_self< >());
        CGEditorBase_exposer.def( "_set_property_SireBase_PropertyPtr",
                                           &SireMol::CGEditorBase::setProperty< SireBase::PropertyPtr >, bp::return_self< >() );
        CGEditorBase_exposer.def( "_set_metadata_SireBase_PropertyPtr", &set_Metadata_SireMol_CGPropertyProperty_function1, bp::return_self< >());
        CGEditorBase_exposer.def( "_set_metadata_SireBase_PropertyPtr", &set_Metadata_SireMol_CGPropertyProperty_function2, bp::return_self< >());
        CGEditorBase_exposer.def( "__str__", &__str__< ::SireMol::Editor<SireMol::CGEditor, SireMol::CutGroup> > );
        CGEditorBase_exposer.def( "__repr__", &__str__< ::SireMol::Editor<SireMol::CGEditor, SireMol::CutGroup> > );
        CGEditorBase_exposer.def( "__len__", &__len_size< ::SireMol::Editor<SireMol::CGEditor, SireMol::CutGroup> > );
    }

}
