
#include "libfoo.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <iostream>
#include <memory>

using std::cout;
using std::endl;
using std::auto_ptr;

int main(int argc, const char **argv)
{

    Foo<int> my_foo_int(21);
    Foo<double> my_foo_double(999.9);
    
    std::auto_ptr<FooBase> libfoo_int( makeFoo_Int() );
    std::auto_ptr<FooBase> libfoo_double( makeFoo_Double() );
    
    my_foo_int.helloWorld();
    my_foo_double.helloWorld();
    
    libfoo_int->helloWorld();
    libfoo_double->helloWorld();
    
    if (not my_foo_int.isA< Foo<int> >())
    {
        cout << "my_foo_int is not a Foo<int>!?!?" << endl;
        return -1;
    }
    
    if (not my_foo_double.isA< Foo<double> >())
    {
        cout << "my_foo_double is not a Foo<double>!?!?" << endl;
        return -1;
    }
    
    if (not testFoo_Int( *libfoo_int ))
    {
        cout << "libfoo_int is not a Foo<int>!?!?" << endl;
        return -1;
    }
    
    if (not testFoo_Double( *libfoo_double ))
    {
        cout << "libfoo_double is not a Foo<double>!?!?" << endl;
        return -1;
    }
    
    if (not testFoo_Int(my_foo_int))
    {
        cout << "my_foo_int is not a Foo<int>. "
             << "RTTI across shared library boundaries isn't working!" << endl;
             
        return -1;
    }
    
    if (not testFoo_Double(my_foo_double))
    {
        cout << "my_foo_double is not a Foo<double>. "
             << "RTTI across shared library boundaries isn't working!" << endl;
             
        return -1;
    }
    
    if (not libfoo_int->isA< Foo<int> >())
    {
        cout << "libfoo_int is not a Foo<int>. "
             << "RTTI across shared library boundaries isn't working!" << endl;
             
        return -1;
    }
    
    if (not libfoo_double->isA< Foo<double> >())
    {
        cout << "libfoo_int is not a Foo<int>. "
             << "RTTI across shared library boundaries isn't working!" << endl;
             
        return -1;
    }

    SireStream::detail::SharedDataHolderT< Foo<double> > holder_double = 
            SireStream::detail::SharedDataHolderT< Foo<double> >(
                    libfoo_double->asA< Foo<double> >() );
    
    testFoo_Holder_Double( holder_double );
    
    cout << "Everything is ok. RTTI across shared library boundaries "
            "is working. This is good, as the Sire property system "
            "requires that this works!" << endl;
    
    return 0;
}
