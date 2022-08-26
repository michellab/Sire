// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Trajectory.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/generalunitproperty.h"

#include "SireBase/slice.h"

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMol/core.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "trajectory.h"

#include "trajectory.h"

SireMol::Trajectory __copy__(const SireMol::Trajectory &other){ return SireMol::Trajectory(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Trajectory_class(){

    { //::SireMol::Trajectory
        typedef bp::class_< SireMol::Trajectory, bp::bases< SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property > > Trajectory_exposer_t;
        Trajectory_exposer_t Trajectory_exposer = Trajectory_exposer_t( "Trajectory", "This is a molecular property that holds the handle to the\ntrajectory data for that molecule. In addition to the\nhandle, this also holds the index of the first atom\nin the underlying trajectory data (trajectory data is a\nvector of coordinates in atom index order for each molecule)\n", bp::init< >("") );
        bp::scope Trajectory_scope( Trajectory_exposer );
        Trajectory_exposer.def( bp::init< SireMol::TrajectoryData const & >(( bp::arg("trajectory") ), "") );
        Trajectory_exposer.def( bp::init< QList< SireBase::SharedPolyPointer< SireMol::TrajectoryData > > const & >(( bp::arg("trajectories") ), "") );
        Trajectory_exposer.def( bp::init< SireMol::TrajectoryData const &, int, int >(( bp::arg("trajectory"), bp::arg("start_atom"), bp::arg("natoms") ), "") );
        Trajectory_exposer.def( bp::init< QList< SireBase::SharedPolyPointer< SireMol::TrajectoryData > > const &, int, int >(( bp::arg("trajectories"), bp::arg("start_atom"), bp::arg("natoms") ), "") );
        Trajectory_exposer.def( bp::init< SireMol::Trajectory const & >(( bp::arg("other") ), "") );
        { //::SireMol::Trajectory::appendFrame
        
            typedef void ( ::SireMol::Trajectory::*appendFrame_function_type)( ::SireMol::Frame const & ) ;
            appendFrame_function_type appendFrame_function_value( &::SireMol::Trajectory::appendFrame );
            
            Trajectory_exposer.def( 
                "appendFrame"
                , appendFrame_function_value
                , ( bp::arg("frame") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::count
        
            typedef int ( ::SireMol::Trajectory::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::Trajectory::count );
            
            Trajectory_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::deleteFrame
        
            typedef void ( ::SireMol::Trajectory::*deleteFrame_function_type)( int ) ;
            deleteFrame_function_type deleteFrame_function_value( &::SireMol::Trajectory::deleteFrame );
            
            Trajectory_exposer.def( 
                "deleteFrame"
                , deleteFrame_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::getFrame
        
            typedef ::SireMol::Frame ( ::SireMol::Trajectory::*getFrame_function_type)( int ) const;
            getFrame_function_type getFrame_function_value( &::SireMol::Trajectory::getFrame );
            
            Trajectory_exposer.def( 
                "getFrame"
                , getFrame_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::insertFrame
        
            typedef void ( ::SireMol::Trajectory::*insertFrame_function_type)( int,::SireMol::Frame const & ) ;
            insertFrame_function_type insertFrame_function_value( &::SireMol::Trajectory::insertFrame );
            
            Trajectory_exposer.def( 
                "insertFrame"
                , insertFrame_function_value
                , ( bp::arg("i"), bp::arg("frame") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::isCompatibleWith
        
            typedef bool ( ::SireMol::Trajectory::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::Trajectory::isCompatibleWith );
            
            Trajectory_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::isEmpty
        
            typedef bool ( ::SireMol::Trajectory::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Trajectory::isEmpty );
            
            Trajectory_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::nAtoms
        
            typedef int ( ::SireMol::Trajectory::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Trajectory::nAtoms );
            
            Trajectory_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::nFrames
        
            typedef int ( ::SireMol::Trajectory::*nFrames_function_type)(  ) const;
            nFrames_function_type nFrames_function_value( &::SireMol::Trajectory::nFrames );
            
            Trajectory_exposer.def( 
                "nFrames"
                , nFrames_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Trajectory_exposer.def( bp::self != bp::self );
        { //::SireMol::Trajectory::operator=
        
            typedef ::SireMol::Trajectory & ( ::SireMol::Trajectory::*assign_function_type)( ::SireMol::Trajectory const & ) ;
            assign_function_type assign_function_value( &::SireMol::Trajectory::operator= );
            
            Trajectory_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Trajectory_exposer.def( bp::self == bp::self );
        { //::SireMol::Trajectory::operator[]
        
            typedef ::SireMol::Frame ( ::SireMol::Trajectory::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Trajectory::operator[] );
            
            Trajectory_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::Trajectory::operator[]
        
            typedef ::QList< SireMol::Frame > ( ::SireMol::Trajectory::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Trajectory::operator[] );
            
            Trajectory_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::Trajectory::operator[]
        
            typedef ::QList< SireMol::Frame > ( ::SireMol::Trajectory::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Trajectory::operator[] );
            
            Trajectory_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::Trajectory::setFrame
        
            typedef void ( ::SireMol::Trajectory::*setFrame_function_type)( int,::SireMol::Frame const & ) ;
            setFrame_function_type setFrame_function_value( &::SireMol::Trajectory::setFrame );
            
            Trajectory_exposer.def( 
                "setFrame"
                , setFrame_function_value
                , ( bp::arg("i"), bp::arg("frame") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::size
        
            typedef int ( ::SireMol::Trajectory::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::Trajectory::size );
            
            Trajectory_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::toString
        
            typedef ::QString ( ::SireMol::Trajectory::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Trajectory::toString );
            
            Trajectory_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Trajectory::typeName );
            
            Trajectory_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::what
        
            typedef char const * ( ::SireMol::Trajectory::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::Trajectory::what );
            
            Trajectory_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Trajectory_exposer.staticmethod( "typeName" );
        Trajectory_exposer.def( "__copy__", &__copy__);
        Trajectory_exposer.def( "__deepcopy__", &__copy__);
        Trajectory_exposer.def( "clone", &__copy__);
        Trajectory_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Trajectory >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Trajectory_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Trajectory >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Trajectory_exposer.def_pickle(sire_pickle_suite< ::SireMol::Trajectory >());
        Trajectory_exposer.def( "__str__", &__str__< ::SireMol::Trajectory > );
        Trajectory_exposer.def( "__repr__", &__str__< ::SireMol::Trajectory > );
        Trajectory_exposer.def( "__len__", &__len_size< ::SireMol::Trajectory > );
    }

}
