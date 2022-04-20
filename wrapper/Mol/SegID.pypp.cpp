// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SegID.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "atom.h"

#include "chain.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "groupatomids.h"

#include "groupgroupids.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "molecules.h"

#include "molinfo.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "residue.h"

#include "segid.h"

#include "segidentifier.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "segid.h"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_SegID_class(){

    { //::SireMol::SegID
        typedef bp::class_< SireMol::SegID, bp::bases< SireID::ID >, boost::noncopyable > SegID_exposer_t;
        SegID_exposer_t SegID_exposer = SegID_exposer_t( "SegID", "This is the base class of all identifiers that are used\nto identify a Segment within a Molecule\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope SegID_scope( SegID_exposer );
        { //::SireMol::SegID::any
        
            typedef ::SireID::MatchAll< SireMol::SegID > ( *any_function_type )(  );
            any_function_type any_function_value( &::SireMol::SegID::any );
            
            SegID_exposer.def( 
                "any"
                , any_function_value
                , bp::release_gil_policy()
                , "Return a match for anything" );
        
        }
        { //::SireMol::SegID::atom
        
            typedef ::SireMol::AtomsIn< SireMol::SegID > ( ::SireMol::SegID::*atom_function_type)( int ) const;
            atom_function_type atom_function_value( &::SireMol::SegID::atom );
            
            SegID_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return a specific atom in the matching residues" );
        
        }
        { //::SireMol::SegID::atoms
        
            typedef ::SireMol::AtomsIn< SireMol::SegID > ( ::SireMol::SegID::*atoms_function_type)(  ) const;
            atoms_function_type atoms_function_value( &::SireMol::SegID::atoms );
            
            SegID_exposer.def( 
                "atoms"
                , atoms_function_value
                , bp::release_gil_policy()
                , "Return the atoms in the matching residues" );
        
        }
        { //::SireMol::SegID::atoms
        
            typedef ::SireMol::AtomsIn< SireMol::SegID > ( ::SireMol::SegID::*atoms_function_type)( int,int ) const;
            atoms_function_type atoms_function_value( &::SireMol::SegID::atoms );
            
            SegID_exposer.def( 
                "atoms"
                , atoms_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "Return a range of atoms in the matching residues" );
        
        }
        { //::SireMol::SegID::fromString
        
            typedef ::SireMol::SegIdentifier ( *fromString_function_type )( ::QString const & );
            fromString_function_type fromString_function_value( &::SireMol::SegID::fromString );
            
            SegID_exposer.def( 
                "fromString"
                , fromString_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "Return an AtomID constructed from the passed string" );
        
        }
        { //::SireMol::SegID::inverse
        
            typedef ::SireID::InvertMatch< SireMol::SegID > ( ::SireMol::SegID::*inverse_function_type)(  ) const;
            inverse_function_type inverse_function_value( &::SireMol::SegID::inverse );
            
            SegID_exposer.def( 
                "inverse"
                , inverse_function_value
                , bp::release_gil_policy()
                , "Return the inverse of this match" );
        
        }
        { //::SireMol::SegID::invert
        
            typedef ::SireID::InvertMatch< SireMol::SegID > ( ::SireMol::SegID::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMol::SegID::invert );
            
            SegID_exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "Return the inverse of this match" );
        
        }
        { //::SireMol::SegID::map
        
            typedef ::QList< SireMol::SegIdx > ( ::SireMol::SegID::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::SegID::map );
            
            SegID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "Map this ID back to the indicies of the segments in the molecule,\nusing the passed MoleculeInfo to do the mapping" );
        
        }
        { //::SireMol::SegID::map
        
            typedef ::QList< SireMol::SegIdx > ( ::SireMol::SegID::*map_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            map_function_type map_function_value( &::SireMol::SegID::map );
            
            SegID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Map this SegID to the atoms in the passed molecule view\nThrow: SireMol::missing_segment\nThrow: SireError::invalid_index\n" );
        
        }
        SegID_exposer.def( !bp::self );
        SegID_exposer.def( bp::self & bp::self );
        SegID_exposer.def( bp::self & bp::other< SireMol::AtomID >() );
        SegID_exposer.def( bp::self & bp::other< SireMol::CGID >() );
        SegID_exposer.def( bp::self & bp::other< SireMol::ResID >() );
        SegID_exposer.def( bp::self & bp::other< SireMol::ChainID >() );
        { //::SireMol::SegID::operator()
        
            typedef ::SireID::Specify< SireMol::SegID > ( ::SireMol::SegID::*__call___function_type)( ::SireBase::Range const & ) const;
            __call___function_type __call___function_value( &::SireMol::SegID::operator() );
            
            SegID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("range") )
                , "" );
        
        }
        { //::SireMol::SegID::operator()
        
            typedef ::SireID::Specify< SireMol::SegID > ( ::SireMol::SegID::*__call___function_type)( ::qint64 ) const;
            __call___function_type __call___function_value( &::SireMol::SegID::operator() );
            
            SegID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::SegID::operator()
        
            typedef ::SireID::Specify< SireMol::SegID > ( ::SireMol::SegID::*__call___function_type)( ::qint64,::qint64 ) const;
            __call___function_type __call___function_value( &::SireMol::SegID::operator() );
            
            SegID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("start"), bp::arg("end") )
                , "" );
        
        }
        { //::SireMol::SegID::operator()
        
            typedef ::SireID::Specify< SireMol::SegID > ( ::SireMol::SegID::*__call___function_type)( ::qint64,::qint64,::qint64 ) const;
            __call___function_type __call___function_value( &::SireMol::SegID::operator() );
            
            SegID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("start"), bp::arg("end"), bp::arg("increment") )
                , "" );
        
        }
        SegID_exposer.def( bp::self * bp::self );
        SegID_exposer.def( bp::self * bp::other< SireMol::AtomID >() );
        SegID_exposer.def( bp::self + bp::self );
        SegID_exposer.def( bp::self + bp::other< SireMol::AtomID >() );
        SegID_exposer.def( bp::self + bp::other< SireMol::CGID >() );
        SegID_exposer.def( bp::self + bp::other< SireMol::ResID >() );
        SegID_exposer.def( bp::self + bp::other< SireMol::ChainID >() );
        SegID_exposer.def( bp::self - bp::self );
        SegID_exposer.def( bp::self - bp::other< SireMol::AtomID >() );
        SegID_exposer.def( bp::self - bp::other< SireMol::CGID >() );
        SegID_exposer.def( bp::self - bp::other< SireMol::ResID >() );
        SegID_exposer.def( bp::self - bp::other< SireMol::ChainID >() );
        SegID_exposer.def( -bp::self );
        { //::SireMol::SegID::operator[]
        
            typedef ::SireID::Specify< SireMol::SegID > ( ::SireMol::SegID::*__getitem___function_type)( ::qint64 ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::SegID::operator[] );
            
            SegID_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::SegID::operator[]
        
            typedef ::SireID::Specify< SireMol::SegID > ( ::SireMol::SegID::*__getitem___function_type)( ::SireBase::Range const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::SegID::operator[] );
            
            SegID_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("range") )
                , "" );
        
        }
        SegID_exposer.def( bp::self | bp::self );
        SegID_exposer.def( bp::self | bp::other< SireMol::AtomID >() );
        { //::SireMol::SegID::selectAllFrom
        
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::SegID::*selectAllFrom_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            selectAllFrom_function_type selectAllFrom_function_value( &::SireMol::SegID::selectAllFrom );
            
            SegID_exposer.def( 
                "selectAllFrom"
                , selectAllFrom_function_value
                , ( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Select all the atoms from the passed view that match this ID\nThrow: SireMol::missing_segment\nThrow: SireError::invalid_index\nThrow: SireMol::duplicate_segment\n" );
        
        }
        { //::SireMol::SegID::selectAllFrom
        
            typedef ::QHash< SireMol::MolNum, SireMol::Selector< SireMol::Segment > > ( ::SireMol::SegID::*selectAllFrom_function_type)( ::SireMol::Molecules const &,::SireBase::PropertyMap const & ) const;
            selectAllFrom_function_type selectAllFrom_function_value( &::SireMol::SegID::selectAllFrom );
            
            SegID_exposer.def( 
                "selectAllFrom"
                , selectAllFrom_function_value
                , ( bp::arg("molecules"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Return all of the atoms from the molecules that match\nthis ID\nThrow: SireMol::missing_segment\n" );
        
        }
        { //::SireMol::SegID::selectAllFrom
        
            typedef ::QHash< SireMol::MolNum, SireMol::Selector< SireMol::Segment > > ( ::SireMol::SegID::*selectAllFrom_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) const;
            selectAllFrom_function_type selectAllFrom_function_value( &::SireMol::SegID::selectAllFrom );
            
            SegID_exposer.def( 
                "selectAllFrom"
                , selectAllFrom_function_value
                , ( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Return the atoms from the molecule group molgroup that match\nthis ID\nThrow: SireMol::missing_segment\n" );
        
        }
        { //::SireMol::SegID::selectAllFrom
        
            typedef ::QHash< SireMol::MolNum, SireMol::Selector< SireMol::Segment > > ( ::SireMol::SegID::*selectAllFrom_function_type)( ::SireMol::MolGroupsBase const &,::SireBase::PropertyMap const & ) const;
            selectAllFrom_function_type selectAllFrom_function_value( &::SireMol::SegID::selectAllFrom );
            
            SegID_exposer.def( 
                "selectAllFrom"
                , selectAllFrom_function_value
                , ( bp::arg("molgroups"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Return the set of atoms that match this ID in the molecule groups\nset molgroups\nThrow: SireMol::missing_segment\n" );
        
        }
        { //::SireMol::SegID::selectFrom
        
            typedef ::SireMol::Segment ( ::SireMol::SegID::*selectFrom_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            selectFrom_function_type selectFrom_function_value( &::SireMol::SegID::selectFrom );
            
            SegID_exposer.def( 
                "selectFrom"
                , selectFrom_function_value
                , ( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Select the atom from the passed view that matches this ID\nThrow: SireMol::missing_segment\nThrow: SireError::invalid_index\nThrow: SireMol::duplicate_segment\n" );
        
        }
        { //::SireMol::SegID::selectFrom
        
            typedef ::SireMol::Segment ( ::SireMol::SegID::*selectFrom_function_type)( ::SireMol::Molecules const &,::SireBase::PropertyMap const & ) const;
            selectFrom_function_type selectFrom_function_value( &::SireMol::SegID::selectFrom );
            
            SegID_exposer.def( 
                "selectFrom"
                , selectFrom_function_value
                , ( bp::arg("molecules"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Return the atom from the molecules molecules that matches\nthis ID\nThrow: SireMol::missing_segment\nThrow: SireMol::duplicate_segment\n" );
        
        }
        { //::SireMol::SegID::selectFrom
        
            typedef ::SireMol::Segment ( ::SireMol::SegID::*selectFrom_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) const;
            selectFrom_function_type selectFrom_function_value( &::SireMol::SegID::selectFrom );
            
            SegID_exposer.def( 
                "selectFrom"
                , selectFrom_function_value
                , ( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Return the atom from the molecule group molgroup that matches\nthis ID\nThrow: SireMol::missing_segment\nThrow: SireMol::duplicate_segment\n" );
        
        }
        { //::SireMol::SegID::selectFrom
        
            typedef ::SireMol::Segment ( ::SireMol::SegID::*selectFrom_function_type)( ::SireMol::MolGroupsBase const &,::SireBase::PropertyMap const & ) const;
            selectFrom_function_type selectFrom_function_value( &::SireMol::SegID::selectFrom );
            
            SegID_exposer.def( 
                "selectFrom"
                , selectFrom_function_value
                , ( bp::arg("molgroups"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Return the atom from the molecule groups molgroups that matches\nthis ID\nThrow: SireMol::missing_segment\nThrow: SireMol::duplicate_segment\n" );
        
        }
        { //::SireMol::SegID::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::SegID::typeName );
            
            SegID_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SegID_exposer.staticmethod( "any" );
        SegID_exposer.staticmethod( "fromString" );
        SegID_exposer.staticmethod( "typeName" );
        SegID_exposer.def( "__str__", &__str__< ::SireMol::SegID > );
        SegID_exposer.def( "__repr__", &__str__< ::SireMol::SegID > );
        SegID_exposer.def( "__hash__", &::SireMol::SegID::hash );
    }

}
