// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "StretchBendParameterName.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/property.h"

#include "SireBase/propertylist.h"

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

SireMM::StretchBendParameterName __copy__(const SireMM::StretchBendParameterName &other){ return SireMM::StretchBendParameterName(other); }

const char* pvt_get_name(const SireMM::StretchBendParameterName&){ return "SireMM::StretchBendParameterName";}

#include "Helpers/release_gil_policy.hpp"

void register_StretchBendParameterName_class(){

    { //::SireMM::StretchBendParameterName
        typedef bp::class_< SireMM::StretchBendParameterName > StretchBendParameterName_exposer_t;
        StretchBendParameterName_exposer_t StretchBendParameterName_exposer = StretchBendParameterName_exposer_t( "StretchBendParameterName", "This class provides the default name of the\nproperty that contains the stretch-bend parameters", bp::init< >("") );
        bp::scope StretchBendParameterName_scope( StretchBendParameterName_exposer );
        { //::SireMM::StretchBendParameterName::stretchBend
        
            typedef ::SireBase::PropertyName const & ( ::SireMM::StretchBendParameterName::*stretchBend_function_type)(  ) const;
            stretchBend_function_type stretchBend_function_value( &::SireMM::StretchBendParameterName::stretchBend );
            
            StretchBendParameterName_exposer.def( 
                "stretchBend"
                , stretchBend_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        StretchBendParameterName_exposer.def( "__copy__", &__copy__);
        StretchBendParameterName_exposer.def( "__deepcopy__", &__copy__);
        StretchBendParameterName_exposer.def( "clone", &__copy__);
        StretchBendParameterName_exposer.def( "__str__", &pvt_get_name);
        StretchBendParameterName_exposer.def( "__repr__", &pvt_get_name);
    }

}
