#include "SireStream/version_error.h"

namespace SireError
{

const char* version_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<version_error>() );
}

const char* version_error::what() const throw()
{
    return version_error::typeName();
}

static RegisterMetaType<version_error> r_version_error;

}
