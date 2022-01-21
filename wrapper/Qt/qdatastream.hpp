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

#ifndef PYWRAP_SIREQT_QDATASTREAM_HPP
#define PYWRAP_SIREQT_QDATASTREAM_HPP

#include <Python.h>

#include <QDataStream>
#include <QByteArray>
#include <QDebug>

#include <boost/python.hpp>

#include "sireglobal.h"

namespace bp = boost::python;

SIRE_BEGIN_HEADER

/** Expose the QDataStream serialisation function */
template<class T>
QDataStream& __rlshift__QDataStream(const T &value, QDataStream &ds)
{
    ds << value;
    return ds;
}

/** Expose the QDataStream deserialisation function */
template<class T>
QDataStream& __rrshift__QDataStream(T &value, QDataStream &ds)
{
    ds >> value;
    return ds;
}

/** Also allow pickling */
template<class T>
bp::dict __getstate__base64(const T &value)
{
    QByteArray b;
    QDataStream ds(&b, QIODevice::WriteOnly);
    ds << value;

    bp::dict d;
    d["sire_pickle_version"] = 1;
    d["sire_pickle_data"] = bp::str(b.toBase64().constData());

    return d;
}

/** Also allow pickling */
template<class T>
void __setstate__base64(T &value, const bp::dict &d)
{
    qDebug() << "NEED TO WRITE THIS!";
}


SIRE_END_HEADER

#endif
