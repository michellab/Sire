// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GromacsAngle.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/conditional.h"

#include "SireCAS/exp.h"

#include "SireCAS/sum.h"

#include "SireCAS/trigfuncs.h"

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "amberparams.h"

#include "gromacsparams.h"

#include "gromacsparams.h"

SireMM::GromacsAngle __copy__(const SireMM::GromacsAngle &other){ return SireMM::GromacsAngle(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_GromacsAngle_class(){

    { //::SireMM::GromacsAngle
        typedef bp::class_< SireMM::GromacsAngle > GromacsAngle_exposer_t;
        GromacsAngle_exposer_t GromacsAngle_exposer = GromacsAngle_exposer_t( "GromacsAngle", "This class holds all of the information about a Gromacs Angle\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope GromacsAngle_scope( GromacsAngle_exposer );
        GromacsAngle_exposer.def( bp::init< int >(( bp::arg("function_type") ), "Construct an angle that is of the specified type, but the parameters have yet\nto be resolved. This is because Gromacs can indicate the required type of\nfunction in the molecule specification, without providing the parameters") );
        GromacsAngle_exposer.def( bp::init< int, double, bp::optional< double, double, double, double, double > >(( bp::arg("function_type"), bp::arg("k0"), bp::arg("k1")=0, bp::arg("k2")=0, bp::arg("k3")=0, bp::arg("k4")=0, bp::arg("k5")=0 ), "Construct an angle of the specified function type with specified parameters\n(the order should be the same as in the Gromacs Manual, table 5.5)") );
        GromacsAngle_exposer.def( bp::init< int, QList< double > const & >(( bp::arg("function_type"), bp::arg("params") ), "Construct an angle of the specified function type by interpreting the parameter\ndata from the passed list of parameter values. These should be in the\nsame order as in the Gromacs Manual, table 5.5") );
        GromacsAngle_exposer.def( bp::init< SireCAS::Expression const &, SireCAS::Symbol const & >(( bp::arg("angle"), bp::arg("theta") ), "Construct from the passed angle, using theta as the symbol for the theta value") );
        GromacsAngle_exposer.def( bp::init< SireMM::GromacsAngle const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::GromacsAngle::assertResolved
        
            typedef void ( ::SireMM::GromacsAngle::*assertResolved_function_type)(  ) const;
            assertResolved_function_type assertResolved_function_value( &::SireMM::GromacsAngle::assertResolved );
            
            GromacsAngle_exposer.def( 
                "assertResolved"
                , assertResolved_function_value
                , bp::release_gil_policy()
                , "Assert that the parameters for this angle have been resolved" );
        
        }
        { //::SireMM::GromacsAngle::at
        
            typedef double ( ::SireMM::GromacsAngle::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMM::GromacsAngle::at );
            
            GromacsAngle_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return the ith parameter for this angle" );
        
        }
        { //::SireMM::GromacsAngle::count
        
            typedef int ( ::SireMM::GromacsAngle::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMM::GromacsAngle::count );
            
            GromacsAngle_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of parameters associated with this angle type" );
        
        }
        { //::SireMM::GromacsAngle::functionType
        
            typedef int ( ::SireMM::GromacsAngle::*functionType_function_type)(  ) const;
            functionType_function_type functionType_function_value( &::SireMM::GromacsAngle::functionType );
            
            GromacsAngle_exposer.def( 
                "functionType"
                , functionType_function_value
                , bp::release_gil_policy()
                , "Return the Gromacs ID number for the function type for this angle. See table\n5.5 in the Gromacs manual for information" );
        
        }
        { //::SireMM::GromacsAngle::functionTypeString
        
            typedef ::QString ( ::SireMM::GromacsAngle::*functionTypeString_function_type)(  ) const;
            functionTypeString_function_type functionTypeString_function_value( &::SireMM::GromacsAngle::functionTypeString );
            
            GromacsAngle_exposer.def( 
                "functionTypeString"
                , functionTypeString_function_value
                , bp::release_gil_policy()
                , "Return the string description of the function type for this angle" );
        
        }
        { //::SireMM::GromacsAngle::hash
        
            typedef ::uint ( ::SireMM::GromacsAngle::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMM::GromacsAngle::hash );
            
            GromacsAngle_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "Return a hash for this bond" );
        
        }
        { //::SireMM::GromacsAngle::isBondAngleCrossTerm
        
            typedef bool ( ::SireMM::GromacsAngle::*isBondAngleCrossTerm_function_type)(  ) const;
            isBondAngleCrossTerm_function_type isBondAngleCrossTerm_function_value( &::SireMM::GromacsAngle::isBondAngleCrossTerm );
            
            GromacsAngle_exposer.def( 
                "isBondAngleCrossTerm"
                , isBondAngleCrossTerm_function_value
                , bp::release_gil_policy()
                , "Return whether or not this angle is really a mix of bond and angle terms" );
        
        }
        { //::SireMM::GromacsAngle::isBondBondCrossTerm
        
            typedef bool ( ::SireMM::GromacsAngle::*isBondBondCrossTerm_function_type)(  ) const;
            isBondBondCrossTerm_function_type isBondBondCrossTerm_function_value( &::SireMM::GromacsAngle::isBondBondCrossTerm );
            
            GromacsAngle_exposer.def( 
                "isBondBondCrossTerm"
                , isBondBondCrossTerm_function_value
                , bp::release_gil_policy()
                , "Return whether or not this angle is really a mix of multiple bond terms" );
        
        }
        { //::SireMM::GromacsAngle::isHarmonic
        
            typedef bool ( ::SireMM::GromacsAngle::*isHarmonic_function_type)(  ) const;
            isHarmonic_function_type isHarmonic_function_value( &::SireMM::GromacsAngle::isHarmonic );
            
            GromacsAngle_exposer.def( 
                "isHarmonic"
                , isHarmonic_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is a harmonic angle" );
        
        }
        { //::SireMM::GromacsAngle::isResolved
        
            typedef bool ( ::SireMM::GromacsAngle::*isResolved_function_type)(  ) const;
            isResolved_function_type isResolved_function_value( &::SireMM::GromacsAngle::isResolved );
            
            GromacsAngle_exposer.def( 
                "isResolved"
                , isResolved_function_value
                , bp::release_gil_policy()
                , "Return whether or not the parameters for this angle are resolved" );
        
        }
        { //::SireMM::GromacsAngle::isSimple
        
            typedef bool ( ::SireMM::GromacsAngle::*isSimple_function_type)(  ) const;
            isSimple_function_type isSimple_function_value( &::SireMM::GromacsAngle::isSimple );
            
            GromacsAngle_exposer.def( 
                "isSimple"
                , isSimple_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is a simple angle function, based only on the\nsize of the angle" );
        
        }
        { //::SireMM::GromacsAngle::needsResolving
        
            typedef bool ( ::SireMM::GromacsAngle::*needsResolving_function_type)(  ) const;
            needsResolving_function_type needsResolving_function_value( &::SireMM::GromacsAngle::needsResolving );
            
            GromacsAngle_exposer.def( 
                "needsResolving"
                , needsResolving_function_value
                , bp::release_gil_policy()
                , "Return whether or not this parameter needs resolving" );
        
        }
        GromacsAngle_exposer.def( bp::self != bp::self );
        GromacsAngle_exposer.def( bp::self < bp::self );
        GromacsAngle_exposer.def( bp::self <= bp::self );
        { //::SireMM::GromacsAngle::operator=
        
            typedef ::SireMM::GromacsAngle & ( ::SireMM::GromacsAngle::*assign_function_type)( ::SireMM::GromacsAngle const & ) ;
            assign_function_type assign_function_value( &::SireMM::GromacsAngle::operator= );
            
            GromacsAngle_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        GromacsAngle_exposer.def( bp::self == bp::self );
        GromacsAngle_exposer.def( bp::self > bp::self );
        GromacsAngle_exposer.def( bp::self >= bp::self );
        { //::SireMM::GromacsAngle::operator[]
        
            typedef double ( ::SireMM::GromacsAngle::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::GromacsAngle::operator[] );
            
            GromacsAngle_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMM::GromacsAngle::parameters
        
            typedef ::QList< double > ( ::SireMM::GromacsAngle::*parameters_function_type)(  ) const;
            parameters_function_type parameters_function_value( &::SireMM::GromacsAngle::parameters );
            
            GromacsAngle_exposer.def( 
                "parameters"
                , parameters_function_value
                , bp::release_gil_policy()
                , "Return all of the parameters for this angle" );
        
        }
        { //::SireMM::GromacsAngle::size
        
            typedef int ( ::SireMM::GromacsAngle::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::GromacsAngle::size );
            
            GromacsAngle_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of parameters associated with this angle type" );
        
        }
        { //::SireMM::GromacsAngle::toAngleTerm
        
            typedef ::SireMM::GromacsAngle ( ::SireMM::GromacsAngle::*toAngleTerm_function_type)(  ) const;
            toAngleTerm_function_type toAngleTerm_function_value( &::SireMM::GromacsAngle::toAngleTerm );
            
            GromacsAngle_exposer.def( 
                "toAngleTerm"
                , toAngleTerm_function_value
                , bp::release_gil_policy()
                , "Return only the angle term part of this angle" );
        
        }
        { //::SireMM::GromacsAngle::toBondAngleExpression
        
            typedef ::SireCAS::Expression ( ::SireMM::GromacsAngle::*toBondAngleExpression_function_type)( ::SireCAS::Symbol const &,::SireCAS::Symbol const & ) const;
            toBondAngleExpression_function_type toBondAngleExpression_function_value( &::SireMM::GromacsAngle::toBondAngleExpression );
            
            GromacsAngle_exposer.def( 
                "toBondAngleExpression"
                , toBondAngleExpression_function_value
                , ( bp::arg("r"), bp::arg("theta") )
                , bp::release_gil_policy()
                , "Return this function converted to a SireCAS::Expression using the passed symbol\nto represent the bond length (r02) and angle size (t012)" );
        
        }
        { //::SireMM::GromacsAngle::toBondBondExpression
        
            typedef ::SireCAS::Expression ( ::SireMM::GromacsAngle::*toBondBondExpression_function_type)( ::SireCAS::Symbol const &,::SireCAS::Symbol const &,::SireCAS::Symbol const & ) const;
            toBondBondExpression_function_type toBondBondExpression_function_value( &::SireMM::GromacsAngle::toBondBondExpression );
            
            GromacsAngle_exposer.def( 
                "toBondBondExpression"
                , toBondBondExpression_function_value
                , ( bp::arg("r01"), bp::arg("r12"), bp::arg("r02") )
                , bp::release_gil_policy()
                , "Return this function converted to a SireCAS::Expression using the passed symbol\nto represent the bond lengths r01, r12 and r02" );
        
        }
        { //::SireMM::GromacsAngle::toBondTerm
        
            typedef ::SireMM::GromacsBond ( ::SireMM::GromacsAngle::*toBondTerm_function_type)(  ) const;
            toBondTerm_function_type toBondTerm_function_value( &::SireMM::GromacsAngle::toBondTerm );
            
            GromacsAngle_exposer.def( 
                "toBondTerm"
                , toBondTerm_function_value
                , bp::release_gil_policy()
                , "Return only the bond term part of this angle" );
        
        }
        { //::SireMM::GromacsAngle::toExpression
        
            typedef ::SireCAS::Expression ( ::SireMM::GromacsAngle::*toExpression_function_type)( ::SireCAS::Symbol const & ) const;
            toExpression_function_type toExpression_function_value( &::SireMM::GromacsAngle::toExpression );
            
            GromacsAngle_exposer.def( 
                "toExpression"
                , toExpression_function_value
                , ( bp::arg("theta") )
                , bp::release_gil_policy()
                , "Return this function converted to a SireCAS::Expression using the passed symbol\nto represent the angle size" );
        
        }
        { //::SireMM::GromacsAngle::toString
        
            typedef ::QString ( ::SireMM::GromacsAngle::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::GromacsAngle::toString );
            
            GromacsAngle_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this angle" );
        
        }
        { //::SireMM::GromacsAngle::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::GromacsAngle::typeName );
            
            GromacsAngle_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::GromacsAngle::what
        
            typedef char const * ( ::SireMM::GromacsAngle::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::GromacsAngle::what );
            
            GromacsAngle_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        GromacsAngle_exposer.staticmethod( "typeName" );
        GromacsAngle_exposer.def( "__copy__", &__copy__);
        GromacsAngle_exposer.def( "__deepcopy__", &__copy__);
        GromacsAngle_exposer.def( "clone", &__copy__);
        GromacsAngle_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::GromacsAngle >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GromacsAngle_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::GromacsAngle >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GromacsAngle_exposer.def_pickle(sire_pickle_suite< ::SireMM::GromacsAngle >());
        GromacsAngle_exposer.def( "__str__", &__str__< ::SireMM::GromacsAngle > );
        GromacsAngle_exposer.def( "__repr__", &__str__< ::SireMM::GromacsAngle > );
        GromacsAngle_exposer.def( "__len__", &__len_size< ::SireMM::GromacsAngle > );
        GromacsAngle_exposer.def( "__hash__", &::SireMM::GromacsAngle::hash );
    }

}
