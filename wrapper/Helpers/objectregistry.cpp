
#include "objectregistry.hpp"

#include <QMutex>

#include <cstdlib>

#include "SireStream/streamdata.hpp"

#include "SireError/errors.h"

using boost::python::object;
using boost::python::extract;

using boost::tuples::tuple;
using boost::shared_ptr;

ObjectRegistry::ObjectRegistry()
{}

ObjectRegistry::~ObjectRegistry()
{}

void ObjectRegistry::throwExtractionError(const object &obj, 
                                          const QString &type_name) const
{
    throw SireError::program_bug( QObject::tr(
       "There was an error trying to extract a Python object of type %1 using "
       "a C++ type of %2.").arg( extract<const char*>(obj.attr("__class__").attr("__name__"))() )
                           .arg(type_name), CODELOC );
}

object ObjectRegistry::getObjects( 
                const QList< boost::tuple<shared_ptr<void>,QString> > &objects )
{
    if (objects.isEmpty())
    {
        //nothing of interest here!
        return object();
    }
    else if (objects.count() == 1)
    {
        const boost::tuple<shared_ptr<void>,QString> &obj = objects[0];
    
        if (obj.get<0>().get() == 0 or obj.get<1>().isEmpty())
            //again, nothing
            return object();
            
        return getConverter(obj.get<1>()).convertFromVoid(obj.get<0>().get());
    }
    else
    {
        //place the objects into a list
        boost::python::list l;
        
        for (int i=0; i<objects.count(); ++i)
        {
            const tuple<shared_ptr<void>,QString> &obj = objects[i];

            if (obj.get<0>().get() == 0 or obj.get<1>().isEmpty())
                //nothing
                l.append( object() );
            else
                l.append( 
                     getConverter(obj.get<1>()).convertFromVoid(obj.get<0>().get()) );
        }
        
        //return the tuple version of this list
        return boost::python::tuple( l );
    }
}

object ObjectRegistry::load(const QByteArray &data)
{
    return ObjectRegistry::getObjects( SireStream::load(data) );
}

object ObjectRegistry::load(const QString &filename)
{
    return ObjectRegistry::getObjects( SireStream::load(filename) );
}

namespace bp = boost::python;

boost::tuple<shared_ptr<void>,QString> 
ObjectRegistry::getObjectFromPython(const object &obj)
{
    object result = obj.attr("what")();

    extract<QString> test_result(result);

    if (not test_result.check())
    {
        throw SireError::invalid_arg( QObject::tr(
           "You cannot save this object to binary, as while it has a \".what()\" "
           "member function, it does not return a string, as expected. "
           "Ask the programmer of this class to provide this functionality."),
                CODELOC );
    }

    QString type_name = test_result();

    return boost::tuple<shared_ptr<void>,QString>( getConverter(type_name).getObject(obj),
                                                   type_name );
}

QList< boost::tuple<shared_ptr<void>,QString> > 
ObjectRegistry::getObjectsFromPython(const object &obj)
{
    PyObject *obj_ptr = obj.ptr();

    QList< boost::tuple<shared_ptr<void>,QString> > objects_to_save;

    if ( PyTuple_Check(obj_ptr) )
    {
        bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );
        
        int n = PyTuple_Size(obj_ptr);
        
        for (int i=0; i<n; ++i)
        {
            objects_to_save.append( getObjectFromPython(t[i]) );
        }
    }
    else if ( PyList_Check(obj_ptr) )
    {
        bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );
        
        int n = PyList_Size(obj_ptr);
        
        for (int i=0; i<n; ++i)
        {
            objects_to_save.append( getObjectFromPython(l[i]) );
        }
    }
    else
    {
        objects_to_save.append( getObjectFromPython(obj) );
    }

    return objects_to_save;
}

QByteArray ObjectRegistry::save(const object &obj)
{
    QList< boost::tuple<shared_ptr<void>,QString> > 
                                    objects_t = getObjectsFromPython(obj);
    
    return SireStream::detail::streamDataSave(objects_t);
}

void ObjectRegistry::save(const object &obj, const QString &filename)
{
    QList< boost::tuple<shared_ptr<void>,QString> > 
                                    objects_t = getObjectsFromPython(obj);
    
    SireStream::detail::streamDataSave(objects_t, filename);
}

Q_GLOBAL_STATIC( QMutex, registryMutex );
 
typedef QHash< QString,ObjectRegistry* > ObjectRegistryType;

Q_GLOBAL_STATIC( ObjectRegistryType, objectRegistry );

static ObjectRegistryType& getRegistry()
{ 
    ObjectRegistryType *registry = objectRegistry();
    
    if (registry == 0)
        std::abort();

    return *registry;
}

void ObjectRegistry::registerConverter(const char *type_name,
                                       const ObjectRegistry &converter)
{
    QMutexLocker lkr( registryMutex() );
    getRegistry().insert( type_name, converter.clone() );   
}

const ObjectRegistry& ObjectRegistry::getConverter(const QString &type_name)
{
    QMutexLocker lkr( registryMutex() );

    if (not getRegistry().contains(type_name))
    {
        throw SireError::unknown_type( QObject::tr(
           "The type %1 is not known in the registry. Make sure that you load the "
           "Python module that contains this type, and, if that doesn't work, then "
           "ask a programmer to see if they could add support.").arg(type_name),
               CODELOC );
    }

    const ObjectRegistry *registry = *(getRegistry().constFind(type_name));
    
    if (registry == 0)
        throw SireError::program_bug( QObject::tr(
                "How did the registry become null for %1?").arg(type_name), 
                    CODELOC );

    return *registry;
}

