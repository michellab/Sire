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
#include <QString>
#include <QRegExp>
#include <wchar.h>

#include "Helpers/str.hpp"
#include "Helpers/release_gil_policy.hpp"

using namespace boost::python;

/** This function converts a QString to a python unicode object

*/
struct qstring_to_python_string
{
    static PyObject* convert(const QString &cpp_string)
    {
        auto raii = boost::python::release_gil_policy::acquire_gil();

        //get Qt to encode the string as UTF8, which python can then decode
        // - this is probably not as efficient as it could be, but Sire
        //   is not a text processing app ;-)
        QByteArray utf8 = cpp_string.toUtf8();

        if (utf8.isEmpty())
            return PyUnicode_DecodeUTF8(0, 0, "strict");
        else
        {
            return PyUnicode_DecodeUTF8(utf8.constData(), utf8.count(), "strict");
        }
    }
};

/** This is the QString __str__ function */
PyObject* QString__str__(QString const& s)
{
    PyObject *unicode_object = qstring_to_python_string::convert(s);

    PyObject *str_object = PyUnicode_AsUTF8String(unicode_object);

    Py_DECREF(unicode_object);

    if (str_object == 0)
    {
        Py_DECREF(str_object);
        boost::python::throw_error_already_set();
    }

    return str_object;
}

/** This function convert a python string or unicode to a QString */
void QString_from_python_str_or_unicode( PyObject* obj_ptr,
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

    void* storage = ((converter::rvalue_from_python_storage<QString>*) data)->storage.bytes;
    new (storage) QString(unicode_qstring);
    data->convertible = storage;
}

/** The actual struct used to control the conversion */
struct QString_from_python
{
    QString_from_python()
    {
        converter::registry::push_back(  &convertible,
                                         &construct,
                                         type_id<QString>() );
    }

    /** Can the python object pointed to by 'obj_ptr' be converted
        to a QString?
    */
    static void* convertible(PyObject* obj_ptr)
    {
        if ( PyBytes_Check(obj_ptr) or
             PyUnicode_Check(obj_ptr) )
                return obj_ptr;
        else
            return 0;
    }

    /** Perform the actual conversion */
    static void construct(  PyObject* obj_ptr,
                            converter::rvalue_from_python_stage1_data* data)
    {
        //use python-qt conversion function
        auto raii = boost::python::release_gil_policy::acquire_gil();
        QString_from_python_str_or_unicode(obj_ptr, data);
    }
};

/** This exports the QString wrapper. I use the python-qt conversion to convert
    from python to QString, then a full wrapping of QString to convert from
    QString to python. This allows me to have a fully-featured QString type
    that is visible from python.

    @author Christopher Woods
*/
void autoconvert_QString()
{
    QString_from_python();
    to_python_converter< QString, qstring_to_python_string >();
}
