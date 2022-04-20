// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "PartialMolecule.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "chain.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "evaluator.h"

#include "molecule.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include <QDebug>

#include "partialmolecule.h"

SireMol::PartialMolecule __copy__(const SireMol::PartialMolecule &other){ return SireMol::PartialMolecule(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_PartialMolecule_class(){

    { //::SireMol::PartialMolecule
        typedef bp::class_< SireMol::PartialMolecule, bp::bases< SireMol::MoleculeView, SireBase::Property > > PartialMolecule_exposer_t;
        PartialMolecule_exposer_t PartialMolecule_exposer = PartialMolecule_exposer_t( "PartialMolecule", "This class provides a view to an arbitrary part of a molecule\n(which can range from just a single atom all the way through to\nthe entire molecule). As such, this class can be used to\nrepresent Molecule, Residue and Atom, as well as everything\nin between\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope PartialMolecule_scope( PartialMolecule_exposer );
        PartialMolecule_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molecule") ), "Construct from the passed view") );
        PartialMolecule_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomSelection const & >(( bp::arg("moldata"), bp::arg("atoms") ), "Construct from the selected atoms of the passed molecule whose\ndata is in moldata") );
        PartialMolecule_exposer.def( bp::init< SireMol::PartialMolecule const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::PartialMolecule::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::PartialMolecule::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::PartialMolecule::evaluate );
            
            PartialMolecule_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "Return an evaluator that can evaluate properties\nover all of the atoms in this view" );
        
        }
        { //::SireMol::PartialMolecule::extract
        
            typedef ::SireMol::PartialMolecule ( ::SireMol::PartialMolecule::*extract_function_type)(  ) const;
            extract_function_type extract_function_value( &::SireMol::PartialMolecule::extract );
            
            PartialMolecule_exposer.def( 
                "extract"
                , extract_function_value
                , bp::release_gil_policy()
                , "Extract a copy of this PartialMolecule which contains only the currently\nselected atoms. This allows the used to pull out parts of a larger molecule,\ne.g. if they want to have only selected residues in a protein and do not\nwant to have to store or manipulate the larger protein molecule" );
        
        }
        { //::SireMol::PartialMolecule::hasMetadata
        
            typedef bool ( ::SireMol::PartialMolecule::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::PartialMolecule::hasMetadata );
            
            PartialMolecule_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether or not this molecule has some metadata\nat the metakey metakey" );
        
        }
        { //::SireMol::PartialMolecule::hasMetadata
        
            typedef bool ( ::SireMol::PartialMolecule::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::PartialMolecule::hasMetadata );
            
            PartialMolecule_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether or not the property at key key has\nsome metadata at metakey metakey\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::PartialMolecule::hasProperty
        
            typedef bool ( ::SireMol::PartialMolecule::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::PartialMolecule::hasProperty );
            
            PartialMolecule_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return whether or not this molecule has a property at key key" );
        
        }
        { //::SireMol::PartialMolecule::isEmpty
        
            typedef bool ( ::SireMol::PartialMolecule::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::PartialMolecule::isEmpty );
            
            PartialMolecule_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMol::PartialMolecule::metadata
        
            typedef ::SireBase::Property const & ( ::SireMol::PartialMolecule::*metadata_function_type)( ::SireBase::PropertyName const & ) const;
            metadata_function_type metadata_function_value( &::SireMol::PartialMolecule::metadata );
            
            PartialMolecule_exposer.def( 
                "metadata"
                , metadata_function_value
                , ( bp::arg("metakey") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the metadata at metakey metakey. Note that this returns the\nmetadata for the molecule - no attempt is made to mask this\nmetadata to match the current selection\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::PartialMolecule::metadata
        
            typedef ::SireBase::Property const & ( ::SireMol::PartialMolecule::*metadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            metadata_function_type metadata_function_value( &::SireMol::PartialMolecule::metadata );
            
            PartialMolecule_exposer.def( 
                "metadata"
                , metadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the metadata at the metakey metakey for the property\nat key key. Note that this returns the\nmetadata for the molecule - no attempt is made to mask this\nmetadata to match the current selection\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::PartialMolecule::metadataKeys
        
            typedef ::QStringList ( ::SireMol::PartialMolecule::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::PartialMolecule::metadataKeys );
            
            PartialMolecule_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "Return the keys of all of the metadata contained directly by\nthis molecule" );
        
        }
        { //::SireMol::PartialMolecule::metadataKeys
        
            typedef ::QStringList ( ::SireMol::PartialMolecule::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::PartialMolecule::metadataKeys );
            
            PartialMolecule_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the keys of all metadata for the property at key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::PartialMolecule::move
        
            typedef ::SireMol::Mover< SireMol::PartialMolecule > ( ::SireMol::PartialMolecule::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::PartialMolecule::move );
            
            PartialMolecule_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "Return a mover that can move all of the atoms in this view" );
        
        }
        { //::SireMol::PartialMolecule::nAtoms
        
            typedef int ( ::SireMol::PartialMolecule::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::PartialMolecule::nAtoms );
            
            PartialMolecule_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the number of atoms in this view" );
        
        }
        { //::SireMol::PartialMolecule::nChains
        
            typedef int ( ::SireMol::PartialMolecule::*nChains_function_type)(  ) const;
            nChains_function_type nChains_function_value( &::SireMol::PartialMolecule::nChains );
            
            PartialMolecule_exposer.def( 
                "nChains"
                , nChains_function_value
                , bp::release_gil_policy()
                , "Return the number of chains in this view" );
        
        }
        { //::SireMol::PartialMolecule::nCutGroups
        
            typedef int ( ::SireMol::PartialMolecule::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::PartialMolecule::nCutGroups );
            
            PartialMolecule_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , bp::release_gil_policy()
                , "Return the number of CutGroups in this view" );
        
        }
        { //::SireMol::PartialMolecule::nResidues
        
            typedef int ( ::SireMol::PartialMolecule::*nResidues_function_type)(  ) const;
            nResidues_function_type nResidues_function_value( &::SireMol::PartialMolecule::nResidues );
            
            PartialMolecule_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , bp::release_gil_policy()
                , "Return the number of residues in this view" );
        
        }
        { //::SireMol::PartialMolecule::nSegments
        
            typedef int ( ::SireMol::PartialMolecule::*nSegments_function_type)(  ) const;
            nSegments_function_type nSegments_function_value( &::SireMol::PartialMolecule::nSegments );
            
            PartialMolecule_exposer.def( 
                "nSegments"
                , nSegments_function_value
                , bp::release_gil_policy()
                , "Return the number of segments in this view" );
        
        }
        { //::SireMol::PartialMolecule::name
        
            typedef ::SireMol::MolName const & ( ::SireMol::PartialMolecule::*name_function_type)(  ) const;
            name_function_type name_function_value( &::SireMol::PartialMolecule::name );
            
            PartialMolecule_exposer.def( 
                "name"
                , name_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the name of this molecule" );
        
        }
        { //::SireMol::PartialMolecule::number
        
            typedef ::SireMol::MolNum ( ::SireMol::PartialMolecule::*number_function_type)(  ) const;
            number_function_type number_function_value( &::SireMol::PartialMolecule::number );
            
            PartialMolecule_exposer.def( 
                "number"
                , number_function_value
                , bp::release_gil_policy()
                , "Return the identifying number of this molecule" );
        
        }
        PartialMolecule_exposer.def( bp::self != bp::self );
        { //::SireMol::PartialMolecule::operator=
        
            typedef ::SireMol::PartialMolecule & ( ::SireMol::PartialMolecule::*assign_function_type)( ::SireMol::MoleculeView const & ) ;
            assign_function_type assign_function_value( &::SireMol::PartialMolecule::operator= );
            
            PartialMolecule_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::PartialMolecule::operator=
        
            typedef ::SireMol::PartialMolecule & ( ::SireMol::PartialMolecule::*assign_function_type)( ::SireMol::PartialMolecule const & ) ;
            assign_function_type assign_function_value( &::SireMol::PartialMolecule::operator= );
            
            PartialMolecule_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PartialMolecule_exposer.def( bp::self == bp::self );
        { //::SireMol::PartialMolecule::property
        
            typedef ::SireBase::Property const & ( ::SireMol::PartialMolecule::*property_function_type)( ::SireBase::PropertyName const & ) const;
            property_function_type property_function_value( &::SireMol::PartialMolecule::property );
            
            PartialMolecule_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the property at key key. Note that this returns the\nproperty for the molecule - no attempt is made to mask this\nproperty to match the current selection\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::PartialMolecule::propertyKeys
        
            typedef ::QStringList ( ::SireMol::PartialMolecule::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::PartialMolecule::propertyKeys );
            
            PartialMolecule_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "Return the keys of all of the properties contained in this molecule" );
        
        }
        { //::SireMol::PartialMolecule::selectedAll
        
            typedef bool ( ::SireMol::PartialMolecule::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::PartialMolecule::selectedAll );
            
            PartialMolecule_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "Return whether or not this contains the entire molecule" );
        
        }
        { //::SireMol::PartialMolecule::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::PartialMolecule::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::PartialMolecule::selection );
            
            PartialMolecule_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "Return the atoms that are part of this view" );
        
        }
        { //::SireMol::PartialMolecule::toString
        
            typedef ::QString ( ::SireMol::PartialMolecule::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::PartialMolecule::toString );
            
            PartialMolecule_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this molecule" );
        
        }
        { //::SireMol::PartialMolecule::toUnit
        
            typedef ::SireMol::MolViewPtr ( ::SireMol::PartialMolecule::*toUnit_function_type)(  ) const;
            toUnit_function_type toUnit_function_value( &::SireMol::PartialMolecule::toUnit );
            
            PartialMolecule_exposer.def( 
                "toUnit"
                , toUnit_function_value
                , bp::release_gil_policy()
                , "Return a copy of this PartialMolecule that has been reduced to its unit\ntype, i.e. if this is a single Atom, this returns the Atom, if this is a single\nresidue, this returns the Residue etc." );
        
        }
        { //::SireMol::PartialMolecule::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::PartialMolecule::typeName );
            
            PartialMolecule_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::PartialMolecule::version
        
            typedef ::quint64 ( ::SireMol::PartialMolecule::*version_function_type)(  ) const;
            version_function_type version_function_value( &::SireMol::PartialMolecule::version );
            
            PartialMolecule_exposer.def( 
                "version"
                , version_function_value
                , bp::release_gil_policy()
                , "Return the version number of this molecule - all molecules\nwith the same ID number and version number must be identical" );
        
        }
        { //::SireMol::PartialMolecule::version
        
            typedef ::quint64 ( ::SireMol::PartialMolecule::*version_function_type)( ::SireBase::PropertyName const & ) const;
            version_function_type version_function_value( &::SireMol::PartialMolecule::version );
            
            PartialMolecule_exposer.def( 
                "version"
                , version_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the version number of the property at key key.\nAll molecules with the same ID number and same property version\nnumber must have the same value of this property\n(although this says nothing about any metadata associated\nwith this property)\nThrow: SireBase::missing_property\n" );
        
        }
        PartialMolecule_exposer.staticmethod( "typeName" );
        PartialMolecule_exposer.def( "__copy__", &__copy__);
        PartialMolecule_exposer.def( "__deepcopy__", &__copy__);
        PartialMolecule_exposer.def( "clone", &__copy__);
        PartialMolecule_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::PartialMolecule >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PartialMolecule_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::PartialMolecule >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PartialMolecule_exposer.def_pickle(sire_pickle_suite< ::SireMol::PartialMolecule >());
        PartialMolecule_exposer.def( "__str__", &__str__< ::SireMol::PartialMolecule > );
        PartialMolecule_exposer.def( "__repr__", &__str__< ::SireMol::PartialMolecule > );
        PartialMolecule_exposer.def( "__len__", &__len_size< ::SireMol::PartialMolecule > );
    }

}
