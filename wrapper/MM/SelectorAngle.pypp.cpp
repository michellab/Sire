// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SelectorAngle.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/slice.h"

#include "SireCAS/symbol.h"

#include "SireCAS/values.h"

#include "SireID/index.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "selectorangle.h"

#include "threeatomfunctions.h"

#include <QDebug>

#include "selectorangle.h"

SireMM::SelectorAngle __copy__(const SireMM::SelectorAngle &other){ return SireMM::SelectorAngle(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_SelectorAngle_class(){

    { //::SireMM::SelectorAngle
        typedef bp::class_< SireMM::SelectorAngle, bp::bases< SireMol::MoleculeView, SireBase::Property > > SelectorAngle_exposer_t;
        SelectorAngle_exposer_t SelectorAngle_exposer = SelectorAngle_exposer_t( "SelectorAngle", "This provides a Selector<T>-style interface for multiple angles", bp::init< >("") );
        bp::scope SelectorAngle_scope( SelectorAngle_exposer );
        SelectorAngle_exposer.def( bp::init< SireMM::Angle const & >(( bp::arg("angle") ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeData const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomID const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::AngleID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("angle"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeView const &, QList< SireMol::AngleID > const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("angles"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::AtomID const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::Selector< SireMol::Atom > const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atoms"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::Selector< SireMol::Atom > const &, SireMol::Selector< SireMol::Atom > const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atoms0"), bp::arg("atoms1"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMol::Selector< SireMol::Atom > const &, SireMol::Selector< SireMol::Atom > const &, SireMol::Selector< SireMol::Atom > const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atoms0"), bp::arg("atoms1"), bp::arg("atoms2"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorAngle_exposer.def( bp::init< SireMM::SelectorAngle const & >(( bp::arg("other") ), "") );
        { //::SireMM::SelectorAngle::IDs
        
            typedef ::QList< SireMol::AngleID > ( ::SireMM::SelectorAngle::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireMM::SelectorAngle::IDs );
            
            SelectorAngle_exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::add
        
            typedef ::SireMM::SelectorAngle ( ::SireMM::SelectorAngle::*add_function_type)( ::SireMM::Angle const & ) const;
            add_function_type add_function_value( &::SireMM::SelectorAngle::add );
            
            SelectorAngle_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("Angle") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::add
        
            typedef ::SireMM::SelectorAngle ( ::SireMM::SelectorAngle::*add_function_type)( ::SireMM::SelectorAngle const & ) const;
            add_function_type add_function_value( &::SireMM::SelectorAngle::add );
            
            SelectorAngle_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::count
        
            typedef int ( ::SireMM::SelectorAngle::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMM::SelectorAngle::count );
            
            SelectorAngle_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::energies
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > ( ::SireMM::SelectorAngle::*energies_function_type)(  ) const;
            energies_function_type energies_function_value( &::SireMM::SelectorAngle::energies );
            
            SelectorAngle_exposer.def( 
                "energies"
                , energies_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::energies
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > ( ::SireMM::SelectorAngle::*energies_function_type)( ::SireBase::PropertyMap const & ) const;
            energies_function_type energies_function_value( &::SireMM::SelectorAngle::energies );
            
            SelectorAngle_exposer.def( 
                "energies"
                , energies_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::energy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMM::SelectorAngle::*energy_function_type)(  ) const;
            energy_function_type energy_function_value( &::SireMM::SelectorAngle::energy );
            
            SelectorAngle_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::energy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMM::SelectorAngle::*energy_function_type)( ::SireBase::PropertyMap const & ) const;
            energy_function_type energy_function_value( &::SireMM::SelectorAngle::energy );
            
            SelectorAngle_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMM::SelectorAngle::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMM::SelectorAngle::evaluate );
            
            SelectorAngle_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::hasMetadata
        
            typedef bool ( ::SireMM::SelectorAngle::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMM::SelectorAngle::hasMetadata );
            
            SelectorAngle_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::hasMetadata
        
            typedef bool ( ::SireMM::SelectorAngle::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMM::SelectorAngle::hasMetadata );
            
            SelectorAngle_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::hasProperty
        
            typedef bool ( ::SireMM::SelectorAngle::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMM::SelectorAngle::hasProperty );
            
            SelectorAngle_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::intersection
        
            typedef ::SireMM::SelectorAngle ( ::SireMM::SelectorAngle::*intersection_function_type)( ::SireMM::SelectorAngle const & ) const;
            intersection_function_type intersection_function_value( &::SireMM::SelectorAngle::intersection );
            
            SelectorAngle_exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::invert
        
            typedef ::SireMM::SelectorAngle ( ::SireMM::SelectorAngle::*invert_function_type)( ::SireBase::PropertyMap const & ) const;
            invert_function_type invert_function_value( &::SireMM::SelectorAngle::invert );
            
            SelectorAngle_exposer.def( 
                "invert"
                , invert_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::invert
        
            typedef ::SireMM::SelectorAngle ( ::SireMM::SelectorAngle::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMM::SelectorAngle::invert );
            
            SelectorAngle_exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::isEmpty
        
            typedef bool ( ::SireMM::SelectorAngle::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::SelectorAngle::isEmpty );
            
            SelectorAngle_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::measures
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 0, 0, 0, 0, 0, 0, 1 > > ( ::SireMM::SelectorAngle::*measures_function_type)(  ) const;
            measures_function_type measures_function_value( &::SireMM::SelectorAngle::measures );
            
            SelectorAngle_exposer.def( 
                "measures"
                , measures_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::measures
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 0, 0, 0, 0, 0, 0, 1 > > ( ::SireMM::SelectorAngle::*measures_function_type)( ::SireBase::PropertyMap const & ) const;
            measures_function_type measures_function_value( &::SireMM::SelectorAngle::measures );
            
            SelectorAngle_exposer.def( 
                "measures"
                , measures_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::metadataKeys
        
            typedef ::QStringList ( ::SireMM::SelectorAngle::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMM::SelectorAngle::metadataKeys );
            
            SelectorAngle_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::metadataKeys
        
            typedef ::QStringList ( ::SireMM::SelectorAngle::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMM::SelectorAngle::metadataKeys );
            
            SelectorAngle_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::move
        
            typedef ::SireMol::Mover< SireMM::SelectorAngle > ( ::SireMM::SelectorAngle::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMM::SelectorAngle::move );
            
            SelectorAngle_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::nViews
        
            typedef int ( ::SireMM::SelectorAngle::*nViews_function_type)(  ) const;
            nViews_function_type nViews_function_value( &::SireMM::SelectorAngle::nViews );
            
            SelectorAngle_exposer.def( 
                "nViews"
                , nViews_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SelectorAngle_exposer.def( bp::self != bp::self );
        { //::SireMM::SelectorAngle::operator()
        
            typedef ::SireMM::Angle ( ::SireMM::SelectorAngle::*__call___function_type)( int ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorAngle::operator() );
            
            SelectorAngle_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMM::SelectorAngle::operator()
        
            typedef ::SireMM::SelectorAngle ( ::SireMM::SelectorAngle::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorAngle::operator() );
            
            SelectorAngle_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , "" );
        
        }
        { //::SireMM::SelectorAngle::operator()
        
            typedef ::SireMM::SelectorAngle ( ::SireMM::SelectorAngle::*__call___function_type)( ::SireBase::Slice const & ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorAngle::operator() );
            
            SelectorAngle_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMM::SelectorAngle::operator()
        
            typedef ::SireMM::SelectorAngle ( ::SireMM::SelectorAngle::*__call___function_type)( ::QList< long long > const & ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorAngle::operator() );
            
            SelectorAngle_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMM::SelectorAngle::operator()
        
            typedef ::SireMM::SelectorAngle ( ::SireMM::SelectorAngle::*__call___function_type)( ::SireMol::AngleID const & ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorAngle::operator() );
            
            SelectorAngle_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("angle") )
                , "" );
        
        }
        { //::SireMM::SelectorAngle::operator=
        
            typedef ::SireMM::SelectorAngle & ( ::SireMM::SelectorAngle::*assign_function_type)( ::SireMM::SelectorAngle const & ) ;
            assign_function_type assign_function_value( &::SireMM::SelectorAngle::operator= );
            
            SelectorAngle_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("angle") )
                , bp::return_self< >()
                , "" );
        
        }
        SelectorAngle_exposer.def( bp::self == bp::self );
        { //::SireMM::SelectorAngle::operator[]
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorAngle::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::SelectorAngle::operator[] );
            
            SelectorAngle_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMM::SelectorAngle::operator[]
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorAngle::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::SelectorAngle::operator[] );
            
            SelectorAngle_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMM::SelectorAngle::operator[]
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorAngle::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::SelectorAngle::operator[] );
            
            SelectorAngle_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMM::SelectorAngle::operator[]
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorAngle::*__getitem___function_type)( ::SireMol::AngleID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::SelectorAngle::operator[] );
            
            SelectorAngle_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("angle") )
                , "" );
        
        }
        { //::SireMM::SelectorAngle::potentials
        
            typedef ::QList< SireCAS::Expression > ( ::SireMM::SelectorAngle::*potentials_function_type)(  ) const;
            potentials_function_type potentials_function_value( &::SireMM::SelectorAngle::potentials );
            
            SelectorAngle_exposer.def( 
                "potentials"
                , potentials_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::potentials
        
            typedef ::QList< SireCAS::Expression > ( ::SireMM::SelectorAngle::*potentials_function_type)( ::SireBase::PropertyMap const & ) const;
            potentials_function_type potentials_function_value( &::SireMM::SelectorAngle::potentials );
            
            SelectorAngle_exposer.def( 
                "potentials"
                , potentials_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::properties
        
            typedef ::QList< SireBase::Properties > ( ::SireMM::SelectorAngle::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::SelectorAngle::properties );
            
            SelectorAngle_exposer.def( 
                "properties"
                , properties_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::property
        
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireMM::SelectorAngle::*property_function_type)( ::SireBase::PropertyName const & ) const;
            property_function_type property_function_value( &::SireMM::SelectorAngle::property );
            
            SelectorAngle_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::property
        
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireMM::SelectorAngle::*property_function_type)( ::SireBase::PropertyName const &,::SireBase::Property const & ) const;
            property_function_type property_function_value( &::SireMM::SelectorAngle::property );
            
            SelectorAngle_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key"), bp::arg("default_value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::propertyKeys
        
            typedef ::QStringList ( ::SireMM::SelectorAngle::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMM::SelectorAngle::propertyKeys );
            
            SelectorAngle_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::selectedAll
        
            typedef bool ( ::SireMM::SelectorAngle::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMM::SelectorAngle::selectedAll );
            
            SelectorAngle_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMM::SelectorAngle::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMM::SelectorAngle::selection );
            
            SelectorAngle_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::size
        
            typedef int ( ::SireMM::SelectorAngle::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::SelectorAngle::size );
            
            SelectorAngle_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::sizes
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 0, 0, 0, 0, 0, 0, 1 > > ( ::SireMM::SelectorAngle::*sizes_function_type)(  ) const;
            sizes_function_type sizes_function_value( &::SireMM::SelectorAngle::sizes );
            
            SelectorAngle_exposer.def( 
                "sizes"
                , sizes_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::sizes
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 0, 0, 0, 0, 0, 0, 1 > > ( ::SireMM::SelectorAngle::*sizes_function_type)( ::SireBase::PropertyMap const & ) const;
            sizes_function_type sizes_function_value( &::SireMM::SelectorAngle::sizes );
            
            SelectorAngle_exposer.def( 
                "sizes"
                , sizes_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::toList
        
            typedef ::QList< SireBase::PropPtr< SireMol::MoleculeView > > ( ::SireMM::SelectorAngle::*toList_function_type)(  ) const;
            toList_function_type toList_function_value( &::SireMM::SelectorAngle::toList );
            
            SelectorAngle_exposer.def( 
                "toList"
                , toList_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::toSelector
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorAngle::*toSelector_function_type)(  ) const;
            toSelector_function_type toSelector_function_value( &::SireMM::SelectorAngle::toSelector );
            
            SelectorAngle_exposer.def( 
                "toSelector"
                , toSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::toString
        
            typedef ::QString ( ::SireMM::SelectorAngle::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::SelectorAngle::toString );
            
            SelectorAngle_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::SelectorAngle::typeName );
            
            SelectorAngle_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorAngle::what
        
            typedef char const * ( ::SireMM::SelectorAngle::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::SelectorAngle::what );
            
            SelectorAngle_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SelectorAngle_exposer.staticmethod( "typeName" );
        SelectorAngle_exposer.def( "__copy__", &__copy__);
        SelectorAngle_exposer.def( "__deepcopy__", &__copy__);
        SelectorAngle_exposer.def( "clone", &__copy__);
        SelectorAngle_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::SelectorAngle >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SelectorAngle_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::SelectorAngle >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SelectorAngle_exposer.def_pickle(sire_pickle_suite< ::SireMM::SelectorAngle >());
        SelectorAngle_exposer.def( "__str__", &__str__< ::SireMM::SelectorAngle > );
        SelectorAngle_exposer.def( "__repr__", &__str__< ::SireMM::SelectorAngle > );
        SelectorAngle_exposer.def( "__len__", &__len_size< ::SireMM::SelectorAngle > );
    }

}
