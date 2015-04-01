/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#include <Python.h>

#include <boost/python.hpp>
#include <QChar>
#include <QString>
#include <QRegExp>
#include <wchar.h>

#include "Helpers/str.hpp"

using namespace boost::python;

/** This function converts a QChar to a python (8 bit) str object
*/
struct qchar_to_python_string
{
    static PyObject* convert(const QChar &cpp_string)
    {
        //get Qt to encode the string as UTF8, which python can then decode
        // - this is probably not as efficient as it could be, but Sire
        //   is not a text processing app ;-)
        QByteArray utf8 = QString(cpp_string).toUtf8();

        if (utf8.isEmpty())
            return PyUnicode_DecodeUTF8(0, 0, "strict");
        else
        {
            return PyUnicode_DecodeUTF8(utf8.constData(), utf8.count(), "strict");
        }
    }
};

/** This function convert a python string or unicode to a QChar */
void QChar_from_python_str_or_unicode( PyObject* obj_ptr,
                                       converter::rvalue_from_python_stage1_data* data )
{
    // First, convert the object to a python unicode object
    PyObject* unicode_object;

    if (PyUnicode_Check(obj_ptr))
    {
        unicode_object = obj_ptr;
        Py_INCREF(unicode_object);
    }
    else if (PyBytes_Check(obj_ptr))
    {
        // Coerce the `str' to `unicode' using UTF-8. UnicodeDecodeError
        // is thrown if the string doesn't make sense in that encoding.
        unicode_object = PyUnicode_FromEncodedObject(obj_ptr, "utf-8", "strict"); // new reference

        if (unicode_object == 0)
            boost::python::throw_error_already_set();
    }
    else
    {
        //we should never get here!!!
        return;
    }

    //now convert the Python unicode string to UTF-8
    PyObject *utf8 = PyUnicode_AsUTF8String( unicode_object );

    //we now don't need to the unicode string anymore
    Py_DECREF(unicode_object);

    if (utf8 == 0)
        boost::python::throw_error_already_set();

    //get the raw buffer from the string
    char *utf8_string;
    Py_ssize_t utf8_nchars;
    int ok = PyBytes_AsStringAndSize(utf8, &utf8_string, &utf8_nchars);

    if (ok != 0)
    {
        Py_DECREF(utf8);
        boost::python::throw_error_already_set();
    }

    //use Qt to convert the UTF8 string to QString unicode
    QString unicode_qstring;

    if (utf8_nchars > 0)
        unicode_qstring = QString::fromUtf8(utf8_string, utf8_nchars);

    Py_DECREF(utf8);

    QChar unicode_qchar;
    
    if (unicode_qstring.count() > 0)
       unicode_qchar = unicode_qstring[0];

    void* storage = ((converter::rvalue_from_python_storage<QChar>*) data)->storage.bytes;
    new (storage) QChar(unicode_qchar);
    data->convertible = storage;
}

/** The actual struct used to control the conversion */
struct QChar_from_python
{
    QChar_from_python()
    {
        converter::registry::push_back(  &convertible,
                                         &construct,
                                         type_id<QChar>() );
    }

    /** Can the python object pointed to by 'obj_ptr' be converted
        to a QChar?
    */
    static void* convertible(PyObject* obj_ptr)
    {
        if ( PyBytes_Check(obj_ptr) )
        {
             if (PyBytes_GET_SIZE(obj_ptr) == 1){ return obj_ptr; }
             else{ return 0; }
        }
        else if ( PyUnicode_Check(obj_ptr) )
        {
             if (PyUnicode_GET_SIZE(obj_ptr) == 1){ return obj_ptr; }
             else{ return 0; }
        }
        else
            return 0;
    }

    /** Perform the actual conversion */
    static void construct(  PyObject* obj_ptr,
                            converter::rvalue_from_python_stage1_data* data)
    {
        //use python-qt conversion function
        QChar_from_python_str_or_unicode(obj_ptr, data);
    }
};

void autoconvert_QChar()
{
    //code to get a QString from a python object
    QChar_from_python();
    to_python_converter< QChar, qchar_to_python_string >();
}
