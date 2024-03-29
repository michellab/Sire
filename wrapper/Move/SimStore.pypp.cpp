// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "SimStore.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/errors.h"

#include "SireStream/shareddatastream.h"

#include "simstore.h"

#include <QDebug>

#include <QDir>

#include <QFileInfo>

#include "simstore.h"

SireMove::SimStore __copy__(const SireMove::SimStore &other){ return SireMove::SimStore(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMove::SimStore&){ return "SireMove::SimStore";}

#include "Helpers/release_gil_policy.hpp"

void register_SimStore_class(){

    { //::SireMove::SimStore
        typedef bp::class_< SireMove::SimStore > SimStore_exposer_t;
        SimStore_exposer_t SimStore_exposer = SimStore_exposer_t( "SimStore", "This is a simple class that provides a place to store the primary\ninformation for a simulation, namely the system being simulated and\nthe moves to be applied to the system. This information can be stored\ndirectly, for rapid access, or it can be all compressed down into\na binary array, so to save memory or diskspace.\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope SimStore_scope( SimStore_exposer );
        SimStore_exposer.def( bp::init< SireSystem::System const &, SireMove::Moves const &, bp::optional< bool > >(( bp::arg("system"), bp::arg("moves"), bp::arg("compress")=(bool)(false) ), "Construct from the passed system and moves - optionally specify whether\nor not to compress the data now") );
        SimStore_exposer.def( bp::init< SireMove::SimStore const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::SimStore::isPacked
        
            typedef bool ( ::SireMove::SimStore::*isPacked_function_type)(  ) const;
            isPacked_function_type isPacked_function_value( &::SireMove::SimStore::isPacked );
            
            SimStore_exposer.def( 
                "isPacked"
                , isPacked_function_value
                , bp::release_gil_policy()
                , "Return whether or not the data is packed into memory or\nis packed onto disk" );
        
        }
        { //::SireMove::SimStore::isPackedToDisk
        
            typedef bool ( ::SireMove::SimStore::*isPackedToDisk_function_type)(  ) const;
            isPackedToDisk_function_type isPackedToDisk_function_value( &::SireMove::SimStore::isPackedToDisk );
            
            SimStore_exposer.def( 
                "isPackedToDisk"
                , isPackedToDisk_function_value
                , bp::release_gil_policy()
                , "Return whether or not the data is packed to disk" );
        
        }
        { //::SireMove::SimStore::isPackedToMemory
        
            typedef bool ( ::SireMove::SimStore::*isPackedToMemory_function_type)(  ) const;
            isPackedToMemory_function_type isPackedToMemory_function_value( &::SireMove::SimStore::isPackedToMemory );
            
            SimStore_exposer.def( 
                "isPackedToMemory"
                , isPackedToMemory_function_value
                , bp::release_gil_policy()
                , "Return whether or not the data is packed into a compressed\nbinary array" );
        
        }
        { //::SireMove::SimStore::moves
        
            typedef ::SireMove::Moves const & ( ::SireMove::SimStore::*moves_function_type)(  ) const;
            moves_function_type moves_function_value( &::SireMove::SimStore::moves );
            
            SimStore_exposer.def( 
                "moves"
                , moves_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return a copy of the moves being stored" );
        
        }
        SimStore_exposer.def( bp::self != bp::self );
        { //::SireMove::SimStore::operator=
        
            typedef ::SireMove::SimStore & ( ::SireMove::SimStore::*assign_function_type)( ::SireMove::SimStore const & ) ;
            assign_function_type assign_function_value( &::SireMove::SimStore::operator= );
            
            SimStore_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        SimStore_exposer.def( bp::self == bp::self );
        { //::SireMove::SimStore::pack
        
            typedef void ( ::SireMove::SimStore::*pack_function_type)(  ) ;
            pack_function_type pack_function_value( &::SireMove::SimStore::pack );
            
            SimStore_exposer.def( 
                "pack"
                , pack_function_value
                , bp::release_gil_policy()
                , "Pack the system - this packs the system to the same state\nas the last time this function was called (so if it was previously\npacked to disk, then this will pack it to disk). This allows\nyou to call SimStore::unpack(), knowing that SimStore::pack()\nwill restore the packing state" );
        
        }
        { //::SireMove::SimStore::packToDisk
        
            typedef void ( ::SireMove::SimStore::*packToDisk_function_type)(  ) ;
            packToDisk_function_type packToDisk_function_value( &::SireMove::SimStore::packToDisk );
            
            SimStore_exposer.def( 
                "packToDisk"
                , packToDisk_function_value
                , bp::release_gil_policy()
                , "Pack the system and moves to disk - this places the data in\na temporary file in QDir::tempPath()" );
        
        }
        { //::SireMove::SimStore::packToDisk
        
            typedef void ( ::SireMove::SimStore::*packToDisk_function_type)( ::QString const & ) ;
            packToDisk_function_type packToDisk_function_value( &::SireMove::SimStore::packToDisk );
            
            SimStore_exposer.def( 
                "packToDisk"
                , packToDisk_function_value
                , ( bp::arg("tempdir") )
                , bp::release_gil_policy()
                , "Pack the system and moves to disk - this places the data\ninto a temporary file in tempdir" );
        
        }
        { //::SireMove::SimStore::packToMemory
        
            typedef void ( ::SireMove::SimStore::*packToMemory_function_type)(  ) ;
            packToMemory_function_type packToMemory_function_value( &::SireMove::SimStore::packToMemory );
            
            SimStore_exposer.def( 
                "packToMemory"
                , packToMemory_function_value
                , bp::release_gil_policy()
                , "Pack the system and moves to memory - this will compress\nthem to a compressed binary array if they are not packed,\nwill do nothing if they are already packed to memory,\nor will move them from disk to memory if they are already\npacked to disk" );
        
        }
        { //::SireMove::SimStore::setMoves
        
            typedef void ( ::SireMove::SimStore::*setMoves_function_type)( ::SireMove::Moves const & ) ;
            setMoves_function_type setMoves_function_value( &::SireMove::SimStore::setMoves );
            
            SimStore_exposer.def( 
                "setMoves"
                , setMoves_function_value
                , ( bp::arg("moves") )
                , bp::release_gil_policy()
                , "Set the moves to be stored" );
        
        }
        { //::SireMove::SimStore::setSystem
        
            typedef void ( ::SireMove::SimStore::*setSystem_function_type)( ::SireSystem::System const & ) ;
            setSystem_function_type setSystem_function_value( &::SireMove::SimStore::setSystem );
            
            SimStore_exposer.def( 
                "setSystem"
                , setSystem_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Set the system to be stored" );
        
        }
        { //::SireMove::SimStore::setSystemAndMoves
        
            typedef void ( ::SireMove::SimStore::*setSystemAndMoves_function_type)( ::SireSystem::System const &,::SireMove::Moves const & ) ;
            setSystemAndMoves_function_type setSystemAndMoves_function_value( &::SireMove::SimStore::setSystemAndMoves );
            
            SimStore_exposer.def( 
                "setSystemAndMoves"
                , setSystemAndMoves_function_value
                , ( bp::arg("system"), bp::arg("moves") )
                , bp::release_gil_policy()
                , "Set both the system and moves to be stored" );
        
        }
        { //::SireMove::SimStore::system
        
            typedef ::SireSystem::System const & ( ::SireMove::SimStore::*system_function_type)(  ) const;
            system_function_type system_function_value( &::SireMove::SimStore::system );
            
            SimStore_exposer.def( 
                "system"
                , system_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return a copy of the system being stored" );
        
        }
        { //::SireMove::SimStore::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::SimStore::typeName );
            
            SimStore_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::SimStore::unpack
        
            typedef void ( ::SireMove::SimStore::*unpack_function_type)(  ) ;
            unpack_function_type unpack_function_value( &::SireMove::SimStore::unpack );
            
            SimStore_exposer.def( 
                "unpack"
                , unpack_function_value
                , bp::release_gil_policy()
                , "Unpack the system and move from the compressed binary array" );
        
        }
        { //::SireMove::SimStore::what
        
            typedef char const * ( ::SireMove::SimStore::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMove::SimStore::what );
            
            SimStore_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SimStore_exposer.staticmethod( "typeName" );
        SimStore_exposer.def( "__copy__", &__copy__);
        SimStore_exposer.def( "__deepcopy__", &__copy__);
        SimStore_exposer.def( "clone", &__copy__);
        SimStore_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::SimStore >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SimStore_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::SimStore >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SimStore_exposer.def_pickle(sire_pickle_suite< ::SireMove::SimStore >());
        SimStore_exposer.def( "__str__", &pvt_get_name);
        SimStore_exposer.def( "__repr__", &pvt_get_name);
    }

}
