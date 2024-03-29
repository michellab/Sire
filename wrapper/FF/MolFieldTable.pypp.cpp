// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MolFieldTable.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMol/atomselection.h"

#include "SireMol/errors.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleculeview.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "fieldtable.h"

#include "fieldtable.h"

SireFF::MolFieldTable __copy__(const SireFF::MolFieldTable &other){ return SireFF::MolFieldTable(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_MolFieldTable_class(){

    { //::SireFF::MolFieldTable
        typedef bp::class_< SireFF::MolFieldTable > MolFieldTable_exposer_t;
        MolFieldTable_exposer_t MolFieldTable_exposer = MolFieldTable_exposer_t( "MolFieldTable", "This class holds the field at the points of all of the atoms of\nselected CutGroups in a molecule. The MolFieldTable is used\nto accumulate all of the fields acting on these atoms during\na field evaluation, and also to control which fields are\nevaluated (as only the fields on atoms in selected CutGroups\nare evaluated). This allows you to provide some control over\nthe calculation, e.g. only placing a few protein residues into\nthe field table, thereby preventing the fields on all atoms\nin a protein from being evaluated if they arent actually\nnecessary.\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope MolFieldTable_scope( MolFieldTable_exposer );
        MolFieldTable_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molview") ), "Construct to hold the field acting at the points of all of the atoms\nof all of the cutgroups viewed in molview") );
        MolFieldTable_exposer.def( bp::init< SireFF::MolFieldTable const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireFF::MolFieldTable::add
        
            typedef bool ( ::SireFF::MolFieldTable::*add_function_type)( ::SireMol::CGAtomIdx const &,::SireMaths::Vector const & ) ;
            add_function_type add_function_value( &::SireFF::MolFieldTable::add );
            
            MolFieldTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("cgatomidx"), bp::arg("field") )
                , bp::release_gil_policy()
                , "Add the field field onto this table - this returns whether or not the\natom is in this table\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireFF::MolFieldTable::add
        
            typedef bool ( ::SireFF::MolFieldTable::*add_function_type)( ::SireMol::AtomSelection const &,::SireMaths::Vector const & ) ;
            add_function_type add_function_value( &::SireFF::MolFieldTable::add );
            
            MolFieldTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("selected_atoms"), bp::arg("field") )
                , bp::release_gil_policy()
                , "Add the field field onto this table for all of the atoms\nin selected_atoms - this returns whether\nor not any selected atoms are in this table\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireFF::MolFieldTable::add
        
            typedef void ( ::SireFF::MolFieldTable::*add_function_type)( ::SireFF::MolFieldTable const & ) ;
            add_function_type add_function_value( &::SireFF::MolFieldTable::add );
            
            MolFieldTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Add the fields contained in other onto this field table. This will only\nadd the fields for CutGroups that are in both tables" );
        
        }
        { //::SireFF::MolFieldTable::add
        
            typedef void ( ::SireFF::MolFieldTable::*add_function_type)( ::SireMaths::Vector const & ) ;
            add_function_type add_function_value( &::SireFF::MolFieldTable::add );
            
            MolFieldTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("field") )
                , bp::release_gil_policy()
                , "Add the field field onto all of the atom points in this table" );
        
        }
        { //::SireFF::MolFieldTable::divide
        
            typedef void ( ::SireFF::MolFieldTable::*divide_function_type)( double ) ;
            divide_function_type divide_function_value( &::SireFF::MolFieldTable::divide );
            
            MolFieldTable_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Divide the field at all atom points by value" );
        
        }
        { //::SireFF::MolFieldTable::initialise
        
            typedef void ( ::SireFF::MolFieldTable::*initialise_function_type)(  ) ;
            initialise_function_type initialise_function_value( &::SireFF::MolFieldTable::initialise );
            
            MolFieldTable_exposer.def( 
                "initialise"
                , initialise_function_value
                , bp::release_gil_policy()
                , "Initialise this table - this clears all of the fields, resetting them to zero" );
        
        }
        { //::SireFF::MolFieldTable::map
        
            typedef int ( ::SireFF::MolFieldTable::*map_function_type)( ::SireMol::CGIdx ) const;
            map_function_type map_function_value( &::SireFF::MolFieldTable::map );
            
            MolFieldTable_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("cgidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::MolFieldTable::molNum
        
            typedef ::SireMol::MolNum ( ::SireFF::MolFieldTable::*molNum_function_type)(  ) const;
            molNum_function_type molNum_function_value( &::SireFF::MolFieldTable::molNum );
            
            MolFieldTable_exposer.def( 
                "molNum"
                , molNum_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::MolFieldTable::molUID
        
            typedef ::QUuid const & ( ::SireFF::MolFieldTable::*molUID_function_type)(  ) const;
            molUID_function_type molUID_function_value( &::SireFF::MolFieldTable::molUID );
            
            MolFieldTable_exposer.def( 
                "molUID"
                , molUID_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::MolFieldTable::multiply
        
            typedef void ( ::SireFF::MolFieldTable::*multiply_function_type)( double ) ;
            multiply_function_type multiply_function_value( &::SireFF::MolFieldTable::multiply );
            
            MolFieldTable_exposer.def( 
                "multiply"
                , multiply_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Multiply the field at all atom points by value" );
        
        }
        { //::SireFF::MolFieldTable::nCutGroups
        
            typedef int ( ::SireFF::MolFieldTable::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireFF::MolFieldTable::nCutGroups );
            
            MolFieldTable_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::MolFieldTable::nSelectedCutGroups
        
            typedef int ( ::SireFF::MolFieldTable::*nSelectedCutGroups_function_type)(  ) const;
            nSelectedCutGroups_function_type nSelectedCutGroups_function_value( &::SireFF::MolFieldTable::nSelectedCutGroups );
            
            MolFieldTable_exposer.def( 
                "nSelectedCutGroups"
                , nSelectedCutGroups_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        MolFieldTable_exposer.def( bp::self != bp::self );
        MolFieldTable_exposer.def( bp::self * bp::other< double >() );
        MolFieldTable_exposer.def( bp::self + bp::self );
        MolFieldTable_exposer.def( bp::self + bp::other< SireMaths::Vector >() );
        MolFieldTable_exposer.def( bp::self - bp::self );
        MolFieldTable_exposer.def( bp::self - bp::other< SireMaths::Vector >() );
        MolFieldTable_exposer.def( -bp::self );
        MolFieldTable_exposer.def( bp::self / bp::other< double >() );
        { //::SireFF::MolFieldTable::operator=
        
            typedef ::SireFF::MolFieldTable & ( ::SireFF::MolFieldTable::*assign_function_type)( ::SireFF::MolFieldTable const & ) ;
            assign_function_type assign_function_value( &::SireFF::MolFieldTable::operator= );
            
            MolFieldTable_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireFF::MolFieldTable::operator=
        
            typedef ::SireFF::MolFieldTable & ( ::SireFF::MolFieldTable::*assign_function_type)( ::SireMaths::Vector const & ) ;
            assign_function_type assign_function_value( &::SireFF::MolFieldTable::operator= );
            
            MolFieldTable_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("field") )
                , bp::return_self< >()
                , "" );
        
        }
        MolFieldTable_exposer.def( bp::self == bp::self );
        { //::SireFF::MolFieldTable::selected
        
            typedef bool ( ::SireFF::MolFieldTable::*selected_function_type)( ::SireMol::CGIdx ) const;
            selected_function_type selected_function_value( &::SireFF::MolFieldTable::selected );
            
            MolFieldTable_exposer.def( 
                "selected"
                , selected_function_value
                , ( bp::arg("cgidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::MolFieldTable::selectedAll
        
            typedef bool ( ::SireFF::MolFieldTable::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireFF::MolFieldTable::selectedAll );
            
            MolFieldTable_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::MolFieldTable::setAll
        
            typedef void ( ::SireFF::MolFieldTable::*setAll_function_type)( ::SireMaths::Vector const & ) ;
            setAll_function_type setAll_function_value( &::SireFF::MolFieldTable::setAll );
            
            MolFieldTable_exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Set all of the fields at the atom points equal to field" );
        
        }
        { //::SireFF::MolFieldTable::subtract
        
            typedef bool ( ::SireFF::MolFieldTable::*subtract_function_type)( ::SireMol::CGAtomIdx const &,::SireMaths::Vector const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::MolFieldTable::subtract );
            
            MolFieldTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("cgatomidx"), bp::arg("field") )
                , bp::release_gil_policy()
                , "Subtract the field field from this table - this returns whether or not the\natom is in this table\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireFF::MolFieldTable::subtract
        
            typedef bool ( ::SireFF::MolFieldTable::*subtract_function_type)( ::SireMol::AtomSelection const &,::SireMaths::Vector const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::MolFieldTable::subtract );
            
            MolFieldTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("selected_atoms"), bp::arg("field") )
                , bp::release_gil_policy()
                , "Subtract the field field from this table for all of the atoms\nin selected_atoms - this returns whether\nor not any selected atoms are in this table\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireFF::MolFieldTable::subtract
        
            typedef void ( ::SireFF::MolFieldTable::*subtract_function_type)( ::SireFF::MolFieldTable const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::MolFieldTable::subtract );
            
            MolFieldTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Subtract the fields contained in other from this field table. This will only\nsubtract the fields for CutGroups that are in both tables" );
        
        }
        { //::SireFF::MolFieldTable::subtract
        
            typedef void ( ::SireFF::MolFieldTable::*subtract_function_type)( ::SireMaths::Vector const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::MolFieldTable::subtract );
            
            MolFieldTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("field") )
                , bp::release_gil_policy()
                , "Subtract the field field from all of the atom points in this table" );
        
        }
        { //::SireFF::MolFieldTable::toVector
        
            typedef ::QVector< SireMaths::Vector > ( ::SireFF::MolFieldTable::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireFF::MolFieldTable::toVector );
            
            MolFieldTable_exposer.def( 
                "toVector"
                , toVector_function_value
                , bp::release_gil_policy()
                , "Return all of the fields in this table in a single array" );
        
        }
        { //::SireFF::MolFieldTable::toVector
        
            typedef ::QVector< SireMaths::Vector > ( ::SireFF::MolFieldTable::*toVector_function_type)( ::SireMol::AtomSelection const & ) const;
            toVector_function_type toVector_function_value( &::SireFF::MolFieldTable::toVector );
            
            MolFieldTable_exposer.def( 
                "toVector"
                , toVector_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "Return an array of all of the fields at the location of\nthe atoms selected in selection\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireFF::MolFieldTable::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::MolFieldTable::typeName );
            
            MolFieldTable_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::MolFieldTable::what
        
            typedef char const * ( ::SireFF::MolFieldTable::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::MolFieldTable::what );
            
            MolFieldTable_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        MolFieldTable_exposer.staticmethod( "typeName" );
        MolFieldTable_exposer.def( "__copy__", &__copy__);
        MolFieldTable_exposer.def( "__deepcopy__", &__copy__);
        MolFieldTable_exposer.def( "clone", &__copy__);
        MolFieldTable_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::MolFieldTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolFieldTable_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::MolFieldTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolFieldTable_exposer.def_pickle(sire_pickle_suite< ::SireFF::MolFieldTable >());
        MolFieldTable_exposer.def( "__str__", &__str__< ::SireFF::MolFieldTable > );
        MolFieldTable_exposer.def( "__repr__", &__str__< ::SireFF::MolFieldTable > );
        MolFieldTable_exposer.def( "__len__", &__len_size< ::SireFF::MolFieldTable > );
    }

}
