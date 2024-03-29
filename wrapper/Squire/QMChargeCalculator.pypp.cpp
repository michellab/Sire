// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "QMChargeCalculator.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/molecule.h"

#include "SireMol/molecules.h"

#include "SireMol/molnum.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "qmchargecalculator.h"

#include "qmchargecalculator.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_QMChargeCalculator_class(){

    { //::Squire::QMChargeCalculator
        typedef bp::class_< Squire::QMChargeCalculator, bp::bases< SireBase::Property >, boost::noncopyable > QMChargeCalculator_exposer_t;
        QMChargeCalculator_exposer_t QMChargeCalculator_exposer = QMChargeCalculator_exposer_t( "QMChargeCalculator", "This is the base class of all functions which are used\nto calculate atomic partial charges of molecules from\nan underlying QM Hamiltonian\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope QMChargeCalculator_scope( QMChargeCalculator_exposer );
        { //::Squire::QMChargeCalculator::calculate
        
            typedef ::SireMol::AtomCharges ( ::Squire::QMChargeCalculator::*calculate_function_type)( ::SireMol::PartialMolecule const &,::SireBase::PropertyMap const & ) const;
            calculate_function_type calculate_function_value( &::Squire::QMChargeCalculator::calculate );
            
            QMChargeCalculator_exposer.def( 
                "calculate"
                , calculate_function_value
                , ( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the partial charges on the molecule molecule" );
        
        }
        { //::Squire::QMChargeCalculator::mayChangeCharges
        
            typedef bool ( ::Squire::QMChargeCalculator::*mayChangeCharges_function_type)( ::SireMol::PartialMolecule const &,::SireMol::PartialMolecule const &,::SireBase::PropertyMap const & ) const;
            mayChangeCharges_function_type mayChangeCharges_function_value( &::Squire::QMChargeCalculator::mayChangeCharges );
            
            QMChargeCalculator_exposer.def( 
                "mayChangeCharges"
                , mayChangeCharges_function_value
                , ( bp::arg("oldmol"), bp::arg("newmol"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::Squire::QMChargeCalculator::null
        
            typedef ::Squire::NullQMChargeCalculator const & ( *null_function_type )(  );
            null_function_type null_function_value( &::Squire::QMChargeCalculator::null );
            
            QMChargeCalculator_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::Squire::QMChargeCalculator::operator()
        
            typedef ::SireMol::AtomCharges ( ::Squire::QMChargeCalculator::*__call___function_type)( ::SireMol::PartialMolecule const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::Squire::QMChargeCalculator::operator() );
            
            QMChargeCalculator_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::Squire::QMChargeCalculator::scaleFactor
        
            typedef double ( ::Squire::QMChargeCalculator::*scaleFactor_function_type)(  ) const;
            scaleFactor_function_type scaleFactor_function_value( &::Squire::QMChargeCalculator::scaleFactor );
            
            QMChargeCalculator_exposer.def( 
                "scaleFactor"
                , scaleFactor_function_value
                , bp::release_gil_policy()
                , "Return the scale factor for the charges" );
        
        }
        { //::Squire::QMChargeCalculator::setScaleFactor
        
            typedef void ( ::Squire::QMChargeCalculator::*setScaleFactor_function_type)( double ) ;
            setScaleFactor_function_type setScaleFactor_function_value( &::Squire::QMChargeCalculator::setScaleFactor );
            
            QMChargeCalculator_exposer.def( 
                "setScaleFactor"
                , setScaleFactor_function_value
                , ( bp::arg("sclfactor") )
                , bp::release_gil_policy()
                , "Set the scale factor for the charges" );
        
        }
        { //::Squire::QMChargeCalculator::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::QMChargeCalculator::typeName );
            
            QMChargeCalculator_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        QMChargeCalculator_exposer.staticmethod( "null" );
        QMChargeCalculator_exposer.staticmethod( "typeName" );
        QMChargeCalculator_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::QMChargeCalculator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QMChargeCalculator_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::QMChargeCalculator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QMChargeCalculator_exposer.def_pickle(sire_pickle_suite< ::Squire::QMChargeCalculator >());
        QMChargeCalculator_exposer.def( "__str__", &__str__< ::Squire::QMChargeCalculator > );
        QMChargeCalculator_exposer.def( "__repr__", &__str__< ::Squire::QMChargeCalculator > );
    }

}
