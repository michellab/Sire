
#ifndef TEST_TEMPLATE_RTTI_LIBFOO_H
#define TEST_TEMPLATE_RTTI_LIBFOO_H

#include <iostream>

#include <qglobal.h>

namespace SireStream
{
namespace detail
{
class SharedDataHolder;
}
}

class Q_DECL_EXPORT FooBase
{
public:
    FooBase();
    FooBase(const FooBase &other);
    
    virtual ~FooBase();
    
    virtual FooBase* clone() const=0;
    
    virtual void helloWorld() const=0;
    
    template<class T>
    bool isA() const
    {
        return dynamic_cast<const T*>(this) != 0;
    }
    
    template<class T>
    const T& asA() const
    {
        const T *this_as_T = dynamic_cast<const T*>(this);
        
        //assert(this_as_T);
        
        return *this_as_T;
    }
};

Q_DECL_EXPORT FooBase* makeFoo_Int();
Q_DECL_EXPORT FooBase* makeFoo_Double();

Q_DECL_EXPORT bool testFoo_Int(const FooBase &foo);
Q_DECL_EXPORT bool testFoo_Double(const FooBase &foo);

Q_DECL_EXPORT bool testFoo_Holder_Double(const SireStream::detail::SharedDataHolder &holder);
Q_DECL_EXPORT bool testFoo_Holder_Int(const SireStream::detail::SharedDataHolder &holder);

template<class T>
class Q_DECL_EXPORT Foo : public FooBase
{
public:
    Foo();
    Foo(const T &value);
    
    Foo(const Foo<T> &other);
    
    ~Foo();
    
    Foo<T>* clone() const;
    
    void helloWorld() const;

private:
    T value;
};

template<class T>
Foo<T>::Foo() : FooBase()
{}

template<class T>
Foo<T>::Foo(const T &val) : FooBase(), value(val)
{}

template<class T>
Foo<T>::Foo(const Foo<T> &other) : FooBase(other), value(other.value)
{}

template<class T>
Foo<T>::~Foo()
{}

template<class T>
Foo<T>* Foo<T>::clone() const
{
    return new Foo<T>(*this);
}

template<class T>
void Foo<T>::helloWorld() const
{
    std::cout << "Hello World! My value is " << value << std::endl;
}

#endif
