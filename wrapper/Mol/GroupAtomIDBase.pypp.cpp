// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GroupAtomIDBase.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "atomidentifier.h"

#include "cgidentifier.h"

#include "chainidentifier.h"

#include "groupatomids.h"

#include "residentifier.h"

#include "segidentifier.h"

#include "groupatomids.h"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_GroupAtomIDBase_class(){

    { //::SireMol::GroupAtomIDBase
        typedef bp::class_< SireMol::GroupAtomIDBase, bp::bases< SireMol::AtomID, SireID::ID >, boost::noncopyable > GroupAtomIDBase_exposer_t;
        GroupAtomIDBase_exposer_t GroupAtomIDBase_exposer = GroupAtomIDBase_exposer_t( "GroupAtomIDBase", "This is the base class of GroupAtomID, used to abstract\ntemplate-independent parts away from the template code.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope GroupAtomIDBase_scope( GroupAtomIDBase_exposer );
        GroupAtomIDBase_exposer.def( "__str__", &__str__< ::SireMol::GroupAtomIDBase > );
        GroupAtomIDBase_exposer.def( "__repr__", &__str__< ::SireMol::GroupAtomIDBase > );
        GroupAtomIDBase_exposer.def( "__hash__", &::SireMol::GroupAtomIDBase::hash );
    }

}
