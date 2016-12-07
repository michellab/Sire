// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "TinkerParameters.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireMaths/vector.h"

#include "SireMol/element.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/molecules.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "tinker.h"

#include <QFile>

#include <QRegExp>

#include <QString>

#include <QStringList>

#include <QTextStream>

#include "tinker.h"

SireIO::TinkerParameters __copy__(const SireIO::TinkerParameters &other){ return SireIO::TinkerParameters(other); }

const char* pvt_get_name(const SireIO::TinkerParameters&){ return "SireIO::TinkerParameters";}

void register_TinkerParameters_class(){

    { //::SireIO::TinkerParameters
        typedef bp::class_< SireIO::TinkerParameters, bp::bases< SireIO::IOParametersBase > > TinkerParameters_exposer_t;
        TinkerParameters_exposer_t TinkerParameters_exposer = TinkerParameters_exposer_t( "TinkerParameters", "This class holds all of the sources and default values of the\nproperties and parameters used by the Tinker readerwriter\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope TinkerParameters_scope( TinkerParameters_exposer );
        TinkerParameters_exposer.def( "__copy__", &__copy__);
        TinkerParameters_exposer.def( "__deepcopy__", &__copy__);
        TinkerParameters_exposer.def( "clone", &__copy__);
        TinkerParameters_exposer.def( "__str__", &pvt_get_name);
        TinkerParameters_exposer.def( "__repr__", &pvt_get_name);
    }

}
