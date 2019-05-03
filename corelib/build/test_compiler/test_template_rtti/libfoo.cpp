
#include "libfoo.h"

#include "SireStream/shareddatastream.h"

using namespace SireStream::detail;

FooBase::FooBase()
{}

FooBase::FooBase(const FooBase&)
{}

FooBase::~FooBase()
{}

FooBase *makeFoo_Int()
{
    return new Foo<int>(42);
}

FooBase *makeFoo_Double()
{
    return new Foo<double>(365.25);
}

bool testFoo_Int(const FooBase &foo)
{
    return foo.isA< Foo<int> >();
}

bool testFoo_Double(const FooBase &foo)
{
    return foo.isA< Foo<double> >();
}

bool testFoo_Holder_Double(const SharedDataHolder &holder)
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
