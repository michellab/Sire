
#include "libfoo.h"

#include "SireStream/shareddatastream.h"

using namespace SireStream::detail;

FooBase::FooBase()
{}

FooBase::FooBase(const FooBase&)
{}

FooBase::~FooBase()
{}

FooBase Q_DECL_EXPORT *makeFoo_Int()
{
    return new Foo<int>(42);
}

FooBase Q_DECL_EXPORT *makeFoo_Double()
{
    return new Foo<double>(365.25);
}

bool Q_DECL_EXPORT testFoo_Int(const FooBase &foo)
{
    return foo.isA< Foo<int> >();
}

bool Q_DECL_EXPORT testFoo_Double(const FooBase &foo)
{
    return foo.isA< Foo<double> >();
}

bool Q_DECL_EXPORT testFoo_Holder_Double(const SharedDataHolder &holder)
{
    try
    {
        return holder.sharedData< Foo<double> >().isA< Foo<double> >();
    }
    catch(const SireError::exception &e)
    {
        std::cout << qPrintable( e.toString() );
        return false;
    }
}
