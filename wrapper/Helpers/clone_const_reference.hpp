#ifndef CLONE_CONST_REFERENCE_HPP
#define CLONE_CONST_REFERENCE_HPP

/** This return policy copies the returned pointer or reference using the
    object's .clone() function, and it then manages the new cloned object.
    This is a useful alternative to copy_const_reference if you are
    working with a polymorphic hierarchy of classes that span across
    shared library boundaries, as you can be sure that the object
    will be correctly copied into python

    @author Christopher Woods
*/

#include <boost/python.hpp>
#include <boost/python/detail/indirect_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/python/to_python_indirect.hpp>

namespace boost
{
namespace python
{

namespace detail
{

struct make_clone_reference_holder
{
    template <class T>
    static PyObject* execute(T* p)
    {
        if (p == 0)
            return python::detail::none();
        else
        {
            try
            {
                //try the first route, which will downcast to the derived type
                return detail::make_owning_holder::execute<T>(p->clone());
            }
            catch(...)
            {
                PyErr_Clear();

                //failed as the derived type is not wrapped. Instead,
                //return a python wrapper to the base type
                boost::python::object obj( p->clone() );
                return boost::python::incref( obj.ptr() );
            }
        }
    }
};

  template <class R>
  struct clone_const_reference_expects_a_const_reference_return_type
# if defined(__GNUC__) && __GNUC__ >= 3 || defined(__EDG__)
  {}
# endif
  ;

} // end of namespace detail

struct clone_const_reference
{
    template<class T>
    struct apply
    {
        typedef typename mpl::if_c<
            indirect_traits::is_reference_to_const<T>::value
          , to_python_indirect<T,detail::make_clone_reference_holder>
          , detail::clone_const_reference_expects_a_const_reference_return_type<T>
        >::type type;
    };
};

} // end of namespace python
} // end of namespace boost

#endif
