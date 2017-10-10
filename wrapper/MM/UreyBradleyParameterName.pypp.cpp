// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "UreyBradleyParameterName.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/property.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireFF/detail/atomiccoords3d.h"

#include "SireFF/errors.h"

#include "SireMaths/line.h"

#include "SireMaths/torsion.h"

#include "SireMaths/triangle.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/units.h"

#include "internalff.h"

#include "tostring.h"

#include <QDebug>

#include <cstdio>

#include "internalff.h"

SireMM::UreyBradleyParameterName __copy__(const SireMM::UreyBradleyParameterName &other){ return SireMM::UreyBradleyParameterName(other); }

const char* pvt_get_name(const SireMM::UreyBradleyParameterName&){ return "SireMM::UreyBradleyParameterName";}

void register_UreyBradleyParameterName_class(){

    { //::SireMM::UreyBradleyParameterName
        typedef bp::class_< SireMM::UreyBradleyParameterName > UreyBradleyParameterName_exposer_t;
        UreyBradleyParameterName_exposer_t UreyBradleyParameterName_exposer = UreyBradleyParameterName_exposer_t( "UreyBradleyParameterName", "", bp::init< >("") );
        bp::scope UreyBradleyParameterName_scope( UreyBradleyParameterName_exposer );
        { //::SireMM::UreyBradleyParameterName::ureyBradley
        
            typedef ::SireBase::PropertyName const & ( ::SireMM::UreyBradleyParameterName::*ureyBradley_function_type)(  ) const;
            ureyBradley_function_type ureyBradley_function_value( &::SireMM::UreyBradleyParameterName::ureyBradley );
            
            UreyBradleyParameterName_exposer.def( 
                "ureyBradley"
                , ureyBradley_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        UreyBradleyParameterName_exposer.def( "__copy__", &__copy__);
        UreyBradleyParameterName_exposer.def( "__deepcopy__", &__copy__);
        UreyBradleyParameterName_exposer.def( "clone", &__copy__);
        UreyBradleyParameterName_exposer.def( "__str__", &pvt_get_name);
        UreyBradleyParameterName_exposer.def( "__repr__", &pvt_get_name);
    }

}
