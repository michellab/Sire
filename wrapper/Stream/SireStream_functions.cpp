
#include <boost/python.hpp>

#include "Helpers/objectregistry.hpp"

using namespace boost::python;

void register_SireStream_functions()
{
    {
        typedef object (*load_function_type)(const QByteArray&);
        load_function_type load_function_value( &ObjectRegistry::load );

        def( "load", load_function_value );
    }    
    {
        typedef object (*load_function_type)(const QString&);
        load_function_type load_function_value( &ObjectRegistry::load );
        def( "load", load_function_value );
    }
    {
        typedef QByteArray (*save_function_type)(const object&);
        save_function_type save_function_value( &ObjectRegistry::save );
        def( "save", save_function_value );
    }
    {
        typedef void (*save_function_type)(const object&, const QString&);
        save_function_type save_function_value( &ObjectRegistry::save );
        def( "save", save_function_value );
    }
}
