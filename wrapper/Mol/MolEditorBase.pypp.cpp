// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MolEditorBase.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "atomeditor.h"

#include "cgeditor.h"

#include "chain.h"

#include "chaineditor.h"

#include "cutgroup.h"

#include "moleditor.h"

#include "mover.hpp"

#include "reseditor.h"

#include "residue.h"

#include "segeditor.h"

#include "segment.h"

#include "selector.hpp"

#include "moleditor.h"

SireMol::MolEditorBase& set_Metadata_function1(
                              SireMol::MolEditorBase &molview,
                              const QString &metakey, const SireBase::Property &p)
                              { return molview.setMetadata<SireBase::Property>(metakey, p); }

SireMol::MolEditorBase& set_Metadata_function2(
                              SireMol::MolEditorBase &molview,
                              const QString &key, const QString &metakey,
                              const SireBase::Property &p)
                              { return molview.setMetadata<SireBase::Property>(key, metakey, p); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_MolEditorBase_class(){

    { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >
        typedef bp::class_< SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >, bp::bases< SireMol::Molecule, SireMol::MoleculeView, SireBase::Property >, boost::noncopyable > MolEditorBase_exposer_t;
        MolEditorBase_exposer_t MolEditorBase_exposer = MolEditorBase_exposer_t( "MolEditorBase", "", bp::no_init );
        bp::scope MolEditorBase_scope( MolEditorBase_exposer );
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::atom
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*atom_function_type)(  ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::atom );
            
            MolEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::atom
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*atom_function_type)( int,::SireBase::PropertyMap const & ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::atom );
            
            MolEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::atom
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*atom_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::atom );
            
            MolEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::atom
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*atom_function_type)( ::SireMol::AtomID const &,::SireBase::PropertyMap const & ) ;
            atom_function_type atom_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::atom );
            
            MolEditorBase_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("atomid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::chain
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*chain_function_type)(  ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::chain );
            
            MolEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::chain
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*chain_function_type)( int,::SireBase::PropertyMap const & ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::chain );
            
            MolEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::chain
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*chain_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::chain );
            
            MolEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::chain
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*chain_function_type)( ::SireMol::ChainID const &,::SireBase::PropertyMap const & ) ;
            chain_function_type chain_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::chain );
            
            MolEditorBase_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("chainid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::cutGroup
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*cutGroup_function_type)(  ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::cutGroup );
            
            MolEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::cutGroup
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*cutGroup_function_type)( int,::SireBase::PropertyMap const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::cutGroup );
            
            MolEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::cutGroup
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*cutGroup_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::cutGroup );
            
            MolEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::cutGroup
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*cutGroup_function_type)( ::SireMol::CGID const &,::SireBase::PropertyMap const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::cutGroup );
            
            MolEditorBase_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("cgid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::molecule
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*molecule_function_type)(  ) ;
            molecule_function_type molecule_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::molecule );
            
            MolEditorBase_exposer.def( 
                "molecule"
                , molecule_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator=
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > & ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*assign_function_type)( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator= );
            
            MolEditorBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator=
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > & ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*assign_function_type)( ::SireMol::Molecule const & ) ;
            assign_function_type assign_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator= );
            
            MolEditorBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[]
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*__getitem___function_type)( int ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[] );
            
            MolEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[]
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*__getitem___function_type)( ::QString const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[] );
            
            MolEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[]
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*__getitem___function_type)( ::SireMol::AtomID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[] );
            
            MolEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("atomid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[]
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*__getitem___function_type)( ::SireMol::ResID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[] );
            
            MolEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("resid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[]
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*__getitem___function_type)( ::SireMol::CGID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[] );
            
            MolEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[]
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*__getitem___function_type)( ::SireMol::ChainID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[] );
            
            MolEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[]
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*__getitem___function_type)( ::SireMol::SegID const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[] );
            
            MolEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("segid") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[]
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolViewPtr ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*__getitem___function_type)( ::SireID::Index const & ) ;
            __getitem___function_type __getitem___function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::operator[] );
            
            MolEditorBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idx") )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::removeMetadata
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolEditor & ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*removeMetadata_function_type)( ::SireBase::PropertyName const & ) ;
            removeMetadata_function_type removeMetadata_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::removeMetadata );
            
            MolEditorBase_exposer.def( 
                "removeMetadata"
                , removeMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::removeMetadata
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolEditor & ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*removeMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) ;
            removeMetadata_function_type removeMetadata_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::removeMetadata );
            
            MolEditorBase_exposer.def( 
                "removeMetadata"
                , removeMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::removeProperty
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::MolEditor & ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*removeProperty_function_type)( ::SireBase::PropertyName const & ) ;
            removeProperty_function_type removeProperty_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::removeProperty );
            
            MolEditorBase_exposer.def( 
                "removeProperty"
                , removeProperty_function_value
                , ( bp::arg("key") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::residue
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*residue_function_type)(  ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::residue );
            
            MolEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::residue
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*residue_function_type)( int,::SireBase::PropertyMap const & ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::residue );
            
            MolEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::residue
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*residue_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::residue );
            
            MolEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::residue
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*residue_function_type)( ::SireMol::ResID const &,::SireBase::PropertyMap const & ) ;
            residue_function_type residue_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::residue );
            
            MolEditorBase_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("resid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::segment
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*segment_function_type)(  ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::segment );
            
            MolEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::segment
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*segment_function_type)( int,::SireBase::PropertyMap const & ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::segment );
            
            MolEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("i"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::segment
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*segment_function_type)( ::QString const &,::SireBase::PropertyMap const & ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::segment );
            
            MolEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("name"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::segment
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*segment_function_type)( ::SireMol::SegID const &,::SireBase::PropertyMap const & ) ;
            segment_function_type segment_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::segment );
            
            MolEditorBase_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("segid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::AtomEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*select_function_type)( ::SireMol::AtomID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select );
            
            MolEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("atomid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::CGEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*select_function_type)( ::SireMol::CGID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select );
            
            MolEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("cgid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ResEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*select_function_type)( ::SireMol::ResID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select );
            
            MolEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("resid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::ChainEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*select_function_type)( ::SireMol::ChainID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select );
            
            MolEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("chainid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select
        
            typedef SireMol::Editor< SireMol::MolEditor, SireMol::Molecule > exported_class_t;
            typedef ::SireMol::SegEditor ( ::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::*select_function_type)( ::SireMol::SegID const &,::SireBase::PropertyMap const & ) ;
            select_function_type select_function_value( &::SireMol::Editor< SireMol::MolEditor, SireMol::Molecule >::select );
            
            MolEditorBase_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("segid"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        MolEditorBase_exposer.def( "setProperty",
                                       &SireMol::MolEditorBase::setProperty<SireBase::Property>, bp::return_self< >() );
        MolEditorBase_exposer.def( "setMetadata", &set_Metadata_function1, bp::return_self< >());
        MolEditorBase_exposer.def( "setMetadata", &set_Metadata_function2, bp::return_self< >());
        MolEditorBase_exposer.def( "__str__", &__str__< ::SireMol::Editor<SireMol::MolEditor, SireMol::Molecule> > );
        MolEditorBase_exposer.def( "__repr__", &__str__< ::SireMol::Editor<SireMol::MolEditor, SireMol::Molecule> > );
        MolEditorBase_exposer.def( "__len__", &__len_size< ::SireMol::Editor<SireMol::MolEditor, SireMol::Molecule> > );
    }

}
