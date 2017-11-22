// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ArrayProperty_Vector_.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/vectorproperty.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "vectorproperty.h"

#include "vectorproperty.h"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_ArrayProperty_Vector__class(){

    { //::SireBase::ArrayProperty< SireMaths::Vector >
        typedef bp::class_< SireBase::ArrayProperty< SireMaths::Vector >, bp::bases< SireBase::Property >, boost::noncopyable > ArrayProperty_Vector__exposer_t;
        ArrayProperty_Vector__exposer_t ArrayProperty_Vector__exposer = ArrayProperty_Vector__exposer_t( "ArrayProperty_Vector_", "", bp::no_init );
        bp::scope ArrayProperty_Vector__scope( ArrayProperty_Vector__exposer );
        { //::SireBase::ArrayProperty< SireMaths::Vector >::append
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*append_function_type)( ::SireMaths::Vector ) ;
            append_function_type append_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::append );
            
            ArrayProperty_Vector__exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::append
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*append_function_type)( ::SireBase::ArrayProperty< SireMaths::Vector > const & ) ;
            append_function_type append_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::append );
            
            ArrayProperty_Vector__exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("values") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::array
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::QVector< SireMaths::Vector > ( ::SireBase::ArrayProperty< SireMaths::Vector >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::array );
            
            ArrayProperty_Vector__exposer.def( 
                "array"
                , array_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::at
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::SireMaths::Vector ( ::SireBase::ArrayProperty< SireMaths::Vector >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::at );
            
            ArrayProperty_Vector__exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::clear
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*clear_function_type)(  ) ;
            clear_function_type clear_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::clear );
            
            ArrayProperty_Vector__exposer.def( 
                "clear"
                , clear_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::count
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef int ( ::SireBase::ArrayProperty< SireMaths::Vector >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::count );
            
            ArrayProperty_Vector__exposer.def( 
                "count"
                , count_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::empty
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef bool ( ::SireBase::ArrayProperty< SireMaths::Vector >::*empty_function_type)(  ) const;
            empty_function_type empty_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::empty );
            
            ArrayProperty_Vector__exposer.def( 
                "empty"
                , empty_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::getitem
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::SireMaths::Vector ( ::SireBase::ArrayProperty< SireMaths::Vector >::*getitem_function_type)( int ) const;
            getitem_function_type getitem_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::getitem );
            
            ArrayProperty_Vector__exposer.def( 
                "getitem"
                , getitem_function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::insert
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*insert_function_type)( int,::SireMaths::Vector ) ;
            insert_function_type insert_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::insert );
            
            ArrayProperty_Vector__exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("i"), bp::arg("value") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::isEmpty
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef bool ( ::SireBase::ArrayProperty< SireMaths::Vector >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::isEmpty );
            
            ArrayProperty_Vector__exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::move
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*move_function_type)( int,int ) ;
            move_function_type move_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::move );
            
            ArrayProperty_Vector__exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("from"), bp::arg("to") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::operator[]
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::SireMaths::Vector ( ::SireBase::ArrayProperty< SireMaths::Vector >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::operator[] );
            
            ArrayProperty_Vector__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::pop_back
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*pop_back_function_type)(  ) ;
            pop_back_function_type pop_back_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::pop_back );
            
            ArrayProperty_Vector__exposer.def( 
                "pop_back"
                , pop_back_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::pop_front
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*pop_front_function_type)(  ) ;
            pop_front_function_type pop_front_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::pop_front );
            
            ArrayProperty_Vector__exposer.def( 
                "pop_front"
                , pop_front_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::prepend
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*prepend_function_type)( ::SireMaths::Vector ) ;
            prepend_function_type prepend_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::prepend );
            
            ArrayProperty_Vector__exposer.def( 
                "prepend"
                , prepend_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::push_back
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*push_back_function_type)( ::SireMaths::Vector ) ;
            push_back_function_type push_back_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::push_back );
            
            ArrayProperty_Vector__exposer.def( 
                "push_back"
                , push_back_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::push_front
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*push_front_function_type)( ::SireMaths::Vector ) ;
            push_front_function_type push_front_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::push_front );
            
            ArrayProperty_Vector__exposer.def( 
                "push_front"
                , push_front_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::removeAt
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*removeAt_function_type)( int ) ;
            removeAt_function_type removeAt_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::removeAt );
            
            ArrayProperty_Vector__exposer.def( 
                "removeAt"
                , removeAt_function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::removeFirst
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*removeFirst_function_type)(  ) ;
            removeFirst_function_type removeFirst_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::removeFirst );
            
            ArrayProperty_Vector__exposer.def( 
                "removeFirst"
                , removeFirst_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::removeLast
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*removeLast_function_type)(  ) ;
            removeLast_function_type removeLast_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::removeLast );
            
            ArrayProperty_Vector__exposer.def( 
                "removeLast"
                , removeLast_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::replace
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*replace_function_type)( int,::SireMaths::Vector ) ;
            replace_function_type replace_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::replace );
            
            ArrayProperty_Vector__exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("i"), bp::arg("value") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::size
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef int ( ::SireBase::ArrayProperty< SireMaths::Vector >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::size );
            
            ArrayProperty_Vector__exposer.def( 
                "size"
                , size_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::swap
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*swap_function_type)( ::SireBase::ArrayProperty< SireMaths::Vector > & ) ;
            swap_function_type swap_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::swap );
            
            ArrayProperty_Vector__exposer.def( 
                "swap"
                , swap_function_value
                , ( bp::arg("other") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::swap
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::ArrayProperty< SireMaths::Vector >::*swap_function_type)( int,int ) ;
            swap_function_type swap_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::swap );
            
            ArrayProperty_Vector__exposer.def( 
                "swap"
                , swap_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::takeAt
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::SireMaths::Vector ( ::SireBase::ArrayProperty< SireMaths::Vector >::*takeAt_function_type)( int ) ;
            takeAt_function_type takeAt_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::takeAt );
            
            ArrayProperty_Vector__exposer.def( 
                "takeAt"
                , takeAt_function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::takeFirst
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::SireMaths::Vector ( ::SireBase::ArrayProperty< SireMaths::Vector >::*takeFirst_function_type)(  ) ;
            takeFirst_function_type takeFirst_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::takeFirst );
            
            ArrayProperty_Vector__exposer.def( 
                "takeFirst"
                , takeFirst_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::takeLast
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::SireMaths::Vector ( ::SireBase::ArrayProperty< SireMaths::Vector >::*takeLast_function_type)(  ) ;
            takeLast_function_type takeLast_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::takeLast );
            
            ArrayProperty_Vector__exposer.def( 
                "takeLast"
                , takeLast_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::toList
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::QList< SireMaths::Vector > ( ::SireBase::ArrayProperty< SireMaths::Vector >::*toList_function_type)(  ) const;
            toList_function_type toList_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::toList );
            
            ArrayProperty_Vector__exposer.def( 
                "toList"
                , toList_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::toString
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::QString ( ::SireBase::ArrayProperty< SireMaths::Vector >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::toString );
            
            ArrayProperty_Vector__exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::toVector
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::QVector< SireMaths::Vector > ( ::SireBase::ArrayProperty< SireMaths::Vector >::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::toVector );
            
            ArrayProperty_Vector__exposer.def( 
                "toVector"
                , toVector_function_value
                , "" );
        
        }
        { //::SireBase::ArrayProperty< SireMaths::Vector >::value
        
            typedef SireBase::ArrayProperty< SireMaths::Vector > exported_class_t;
            typedef ::QVector< SireMaths::Vector > ( ::SireBase::ArrayProperty< SireMaths::Vector >::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireBase::ArrayProperty< SireMaths::Vector >::value );
            
            ArrayProperty_Vector__exposer.def( 
                "value"
                , value_function_value
                , "" );
        
        }
        ArrayProperty_Vector__exposer.def( "__str__", &__str__< ::SireBase::ArrayProperty<SireMaths::Vector> > );
        ArrayProperty_Vector__exposer.def( "__repr__", &__str__< ::SireBase::ArrayProperty<SireMaths::Vector> > );
        ArrayProperty_Vector__exposer.def( "__len__", &__len_size< ::SireBase::ArrayProperty<SireMaths::Vector> > );
        ArrayProperty_Vector__exposer.def( "__getitem__", &::SireBase::ArrayProperty<SireMaths::Vector>::getitem );
    }

}