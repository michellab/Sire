// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ResEditorBase.pypp.hpp"

namespace bp = boost::python;

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

#include "reseditor.h"

#include "resproperty.hpp"

SireMol::ResEditorBase& set_Metadata_SireMol_ResStringProperty_function1(
                                  SireMol::ResEditorBase &molview,
                                   const QString &metakey, const QString &p)
                                   { return molview.setMetadata< QString >(metakey, p); }

SireMol::ResEditorBase& set_Metadata_SireMol_ResStringProperty_function2(
                                  SireMol::ResEditorBase &molview,
                                   const QString &key, const QString &metakey, const QString &p)
                                   { return molview.setMetadata< QString >(key, metakey, p); }

SireMol::ResEditorBase& set_Metadata_SireMol_ResIntProperty_function1(
                                  SireMol::ResEditorBase &molview,
                                   const QString &metakey, const qint64 &p)
                                   { return molview.setMetadata< qint64 >(metakey, p); }

SireMol::ResEditorBase& set_Metadata_SireMol_ResIntProperty_function2(
                                  SireMol::ResEditorBase &molview,
                                   const QString &key, const QString &metakey, const qint64 &p)
                                   { return molview.setMetadata< qint64 >(key, metakey, p); }

SireMol::ResEditorBase& set_Metadata_SireMol_ResFloatProperty_function1(
                                  SireMol::ResEditorBase &molview,
                                   const QString &metakey, const double &p)
                                   { return molview.setMetadata< double >(metakey, p); }

SireMol::ResEditorBase& set_Metadata_SireMol_ResFloatProperty_function2(
                                  SireMol::ResEditorBase &molview,
                                   const QString &key, const QString &metakey, const double &p)
                                   { return molview.setMetadata< double >(key, metakey, p); }

SireMol::ResEditorBase& set_Metadata_SireMol_ResVariantProperty_function1(
                                  SireMol::ResEditorBase &molview,
                                   const QString &metakey, const QVariant &p)
                                   { return molview.setMetadata< QVariant >(metakey, p); }

SireMol::ResEditorBase& set_Metadata_SireMol_ResVariantProperty_function2(
                                  SireMol::ResEditorBase &molview,
                                   const QString &key, const QString &metakey, const QVariant &p)
                                   { return molview.setMetadata< QVariant >(key, metakey, p); }

SireMol::ResEditorBase& set_Metadata_SireMol_ResPropertyProperty_function1(
                                  SireMol::ResEditorBase &molview,
                                   const QString &metakey, const SireBase::PropertyPtr &p)
                                   { return molview.setMetadata< SireBase::PropertyPtr >(metakey, p); }

SireMol::ResEditorBase& set_Metadata_SireMol_ResPropertyProperty_function2(
                                  SireMol::ResEditorBase &molview,
                                   const QString &key, const QString &metakey, const SireBase::PropertyPtr &p)
                                   { return molview.setMetadata< SireBase::PropertyPtr >(key, metakey, p); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_ResEditorBase_class(){

    { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >
        typedef bp::class_< SireMol::Editor< SireMol::ResEditor, SireMol::Residue >, bp::bases< SireMol::Residue, SireMol::MoleculeView, SireBase::Property >, boost::noncopyable > ResEditorBase_exposer_t;
        ResEditorBase_exposer_t ResEditorBase_exposer = ResEditorBase_exposer_t( "ResEditorBase", "", bp::no_init );
        bp::scope ResEditorBase_scope( ResEditorBase_exposer );
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::atom
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*atom_function_type)(  ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::atom );
            
            ResEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::atom
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*atom_function_type)( int,::SireBase::PropertyMap const & ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::atom );
            
            ResEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::atom
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*atom_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::atom );
            
            ResEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::atom
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*atom_function_type)( ::SireMol::AtomID const &,::SireBase::PropertyMap const & ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::atom );
            
            ResEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("atomid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::chain
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*chain_function_type)(  ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::chain );
            
            ResEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::chain
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*chain_function_type)( int,::SireBase::PropertyMap const & ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::chain );
            
            ResEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::chain
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*chain_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::chain );
            
            ResEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::chain
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*chain_function_type)( ::SireMol::ChainID const &,::SireBase::PropertyMap const & ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::chain );
            
            ResEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("chainid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::cutGroup
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*cutGroup_function_type)(  ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::cutGroup );
            
            ResEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::cutGroup
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*cutGroup_function_type)( int,::SireBase::PropertyMap const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::cutGroup );
            
            ResEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::cutGroup
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*cutGroup_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::cutGroup );
            
            ResEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::cutGroup
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*cutGroup_function_type)( ::SireMol::CGID const &,::SireBase::PropertyMap const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::cutGroup );
            
            ResEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("cgid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::molecule
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::MolEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*molecule_function_type)(  ) ;
            molecule_function_type molecule_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::molecule );
            
            ResEditorBase_exposer.def( 
                "molecule"
                , molecule_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator=
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue > & ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*assign_function_type)( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator= );
            
            ResEditorBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator=
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue > & ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*assign_function_type)( ::SireMol::Residue const & ) ;
            assign_function_type assign_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator= );
            
            ResEditorBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[]
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*__getitem___function_type)( int ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[] );
            
            ResEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[]
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*__getitem___function_type)( ::QString const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[] );
            
            ResEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[]
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*__getitem___function_type)( ::SireMol::AtomID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[] );
            
            ResEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("atomid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[]
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*__getitem___function_type)( ::SireMol::ResID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[] );
            
            ResEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("resid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[]
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*__getitem___function_type)( ::SireMol::CGID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[] );
            
            ResEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[]
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*__getitem___function_type)( ::SireMol::ChainID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[] );
            
            ResEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[]
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*__getitem___function_type)( ::SireMol::SegID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[] );
            
            ResEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("segid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[]
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*__getitem___function_type)( ::SireID::Index const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::operator[] );
            
            ResEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idx") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::removeMetadata
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ResEditor & ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*removeMetadata_function_type)( ::SireBase::PropertyName const & ) ;
            removeMetadata_function_type removeMetadata_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::removeMetadata );
            
            ResEditorBase_exposer.def( 
                "removeMetadata"
                , removeMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::removeMetadata
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ResEditor & ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*removeMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) ;
            removeMetadata_function_type removeMetadata_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::removeMetadata );
            
            ResEditorBase_exposer.def( 
                "removeMetadata"
                , removeMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::removeProperty
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ResEditor & ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*removeProperty_function_type)( ::SireBase::PropertyName const & ) ;
            removeProperty_function_type removeProperty_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::removeProperty );
            
            ResEditorBase_exposer.def( 
                "removeProperty"
                , removeProperty_function_value
                , ( bp::arg("key") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::residue
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*residue_function_type)(  ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::residue );
            
            ResEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::residue
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*residue_function_type)( int,::SireBase::PropertyMap const & ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::residue );
            
            ResEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::residue
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*residue_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::residue );
            
            ResEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::residue
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*residue_function_type)( ::SireMol::ResID const &,::SireBase::PropertyMap const & ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::residue );
            
            ResEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("resid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::segment
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*segment_function_type)(  ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::segment );
            
            ResEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::segment
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*segment_function_type)( int,::SireBase::PropertyMap const & ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::segment );
            
            ResEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::segment
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*segment_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::segment );
            
            ResEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::segment
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*segment_function_type)( ::SireMol::SegID const &,::SireBase::PropertyMap const & ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::segment );
            
            ResEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("segid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*select_function_type)( ::SireMol::AtomID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select );
            
            ResEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("atomid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*select_function_type)( ::SireMol::CGID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select );
            
            ResEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("cgid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*select_function_type)( ::SireMol::ResID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select );
            
            ResEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("resid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*select_function_type)( ::SireMol::ChainID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select );
            
            ResEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("chainid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select
        
            typedef SireMol::Editor< SireMol::ResEditor, SireMol::Residue > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::*select_function_type)( ::SireMol::SegID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::ResEditor, SireMol::Residue >::select );
            
            ResEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("segid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        ResEditorBase_exposer.def( "_set_property_QString",
                                           &SireMol::ResEditorBase::setProperty< QString >, bp::return_self< >() );
        ResEditorBase_exposer.def( "_set_metadata_QString", &set_Metadata_SireMol_ResStringProperty_function1, bp::return_self< >());
        ResEditorBase_exposer.def( "_set_metadata_QString", &set_Metadata_SireMol_ResStringProperty_function2, bp::return_self< >());
        ResEditorBase_exposer.def( "_set_property_qint64",
                                           &SireMol::ResEditorBase::setProperty< qint64 >, bp::return_self< >() );
        ResEditorBase_exposer.def( "_set_metadata_qint64", &set_Metadata_SireMol_ResIntProperty_function1, bp::return_self< >());
        ResEditorBase_exposer.def( "_set_metadata_qint64", &set_Metadata_SireMol_ResIntProperty_function2, bp::return_self< >());
        ResEditorBase_exposer.def( "_set_property_double",
                                           &SireMol::ResEditorBase::setProperty< double >, bp::return_self< >() );
        ResEditorBase_exposer.def( "_set_metadata_double", &set_Metadata_SireMol_ResFloatProperty_function1, bp::return_self< >());
        ResEditorBase_exposer.def( "_set_metadata_double", &set_Metadata_SireMol_ResFloatProperty_function2, bp::return_self< >());
        ResEditorBase_exposer.def( "_set_property_QVariant",
                                           &SireMol::ResEditorBase::setProperty< QVariant >, bp::return_self< >() );
        ResEditorBase_exposer.def( "_set_metadata_QVariant", &set_Metadata_SireMol_ResVariantProperty_function1, bp::return_self< >());
        ResEditorBase_exposer.def( "_set_metadata_QVariant", &set_Metadata_SireMol_ResVariantProperty_function2, bp::return_self< >());
        ResEditorBase_exposer.def( "_set_property_SireBase_PropertyPtr",
                                           &SireMol::ResEditorBase::setProperty< SireBase::PropertyPtr >, bp::return_self< >() );
        ResEditorBase_exposer.def( "_set_metadata_SireBase_PropertyPtr", &set_Metadata_SireMol_ResPropertyProperty_function1, bp::return_self< >());
        ResEditorBase_exposer.def( "_set_metadata_SireBase_PropertyPtr", &set_Metadata_SireMol_ResPropertyProperty_function2, bp::return_self< >());
        ResEditorBase_exposer.def( "__str__", &__str__< ::SireMol::Editor<SireMol::ResEditor, SireMol::Residue> > );
        ResEditorBase_exposer.def( "__repr__", &__str__< ::SireMol::Editor<SireMol::ResEditor, SireMol::Residue> > );
        ResEditorBase_exposer.def( "__len__", &__len_size< ::SireMol::Editor<SireMol::ResEditor, SireMol::Residue> > );
    }

}
