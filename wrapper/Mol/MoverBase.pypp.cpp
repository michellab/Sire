// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MoverBase.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/align.h"

#include "SireMaths/axisset.h"

#include "SireMaths/matrix.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/rotate.h"

#include "SireMaths/vectorproperty.h"

#include "SireMol/errors.h"

#include "SireUnits/units.h"

#include "SireVol/coordgroup.h"

#include "SireVol/space.h"

#include "angleid.h"

#include "atomcoords.h"

#include "atommatcher.h"

#include "atommatchers.h"

#include "bondid.h"

#include "connectivity.h"

#include "dihedralid.h"

#include "improperid.h"

#include "mover.h"

#include "tostring.h"

#include "weightfunction.h"

#include "mover.h"

SireMol::MoverBase __copy__(const SireMol::MoverBase &other){ return SireMol::MoverBase(other); }

const char* pvt_get_name(const SireMol::MoverBase&){ return "SireMol::MoverBase";}

#include "Helpers/release_gil_policy.hpp"

void register_MoverBase_class(){

    { //::SireMol::MoverBase
        typedef bp::class_< SireMol::MoverBase > MoverBase_exposer_t;
        MoverBase_exposer_t MoverBase_exposer = MoverBase_exposer_t( "MoverBase", "This class provides the template-independent part\nof Mover<T>. This class is not designed to be used\non its own\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope MoverBase_scope( MoverBase_exposer );
        MoverBase_exposer.def( bp::init< SireMol::MoverBase const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::MoverBase::operator=
        
            typedef ::SireMol::MoverBase & ( ::SireMol::MoverBase::*assign_function_type)( ::SireMol::MoverBase const & ) ;
            assign_function_type assign_function_value( &::SireMol::MoverBase::operator= );
            
            MoverBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        MoverBase_exposer.def( "__copy__", &__copy__);
        MoverBase_exposer.def( "__deepcopy__", &__copy__);
        MoverBase_exposer.def( "clone", &__copy__);
        MoverBase_exposer.def( "__str__", &pvt_get_name);
        MoverBase_exposer.def( "__repr__", &pvt_get_name);
    }

}
