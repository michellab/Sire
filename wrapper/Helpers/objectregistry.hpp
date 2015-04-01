#ifndef PYTHON2_HELPERS_OBJECTREGISTRY_HPP
#define PYTHON2_HELPERS_OBJECTREGISTRY_HPP

#include <Python.h>
#include <boost/python.hpp>

#include "SireStream/streamdata.hpp"

#include <QHash>
#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

class SIRE_EXPORT ObjectRegistry
{
public:
    virtual ~ObjectRegistry();

    static boost::python::object load(const QByteArray &data);
    static boost::python::object load(const QString &filename);

    static QByteArray save(const boost::python::object &object);
    static void save(const boost::python::object &object, const QString &filename);

    virtual ObjectRegistry* clone() const=0;
    virtual const char* what() const=0;

    template<class T>
    static void registerConverterFor();

protected:
    virtual boost::python::object convertFromVoid(const void * const ptr) const=0;

    virtual boost::shared_ptr<void> getObject(
                            const boost::python::object &object) const=0;

    static void registerConverter(const char *type_name,
                                  const ObjectRegistry &converter);

    static const ObjectRegistry& getConverter(const QString &type_name);

    static boost::python::object getObjects( 
                const QList< boost::tuple<boost::shared_ptr<void>,QString> > &objects );

    static boost::tuple<boost::shared_ptr<void>,QString> getObjectFromPython(
                                        const boost::python::object &obj);

    static QList< boost::tuple<boost::shared_ptr<void>,QString> > getObjectsFromPython(
                                        const boost::python::object &obj);
    
    void throwExtractionError(const boost::python::object &obj, 
                              const QString &type_name) const;

    ObjectRegistry();
};

namespace detail
{

template<class T>
const char* getString(T value)
{
    return value;
}

template<>
inline const char* getString<const char*>(const char *value)
{
    return value;
}

template<>
inline const char* getString<QString>(QString value)
{
    return value.toUtf8().constData();
}

template<class T>
class SIRE_EXPORT ObjectRegistryT : public ObjectRegistry
{

friend class ObjectRegistry;

public:
    ~ObjectRegistryT()
    {}

    ObjectRegistry* clone() const;

    const char* what() const;

protected:
    boost::python::object convertFromVoid(const void * const ptr) const
    {
        if (ptr == 0)
            return boost::python::object();

        const T * const t_ptr = static_cast<const T * const>(ptr);

        return boost::python::object( *t_ptr );
    }

    boost::shared_ptr<void> getObject(const boost::python::object &obj) const
    {
        boost::python::extract<const T&> t_object(obj);

        if (not t_object.check())
            this->throwExtractionError(obj, T::typeName());

        return boost::shared_ptr<T>( new T(t_object()), SireStream::detail::void_deleter(
                                                                 qMetaTypeId<T>()) );
    }

private:
    ObjectRegistryT() : ObjectRegistry()
    {}
};
    
template<class T>
SIRE_OUTOFLINE_TEMPLATE
ObjectRegistry* ObjectRegistryT<T>::clone() const
{
    return new ObjectRegistryT<T>();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* ObjectRegistryT<T>::what() const
{
    return getString( T::typeName() );
}

}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void ObjectRegistry::registerConverterFor()
{
    ObjectRegistry::registerConverter( detail::getString(T::typeName()), detail::ObjectRegistryT<T>() );
}

SIRE_END_HEADER

#endif

