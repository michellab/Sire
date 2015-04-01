/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) <year>  <name of author>
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

#include <QVector>
#include <QSet>
#include <QStringList>
#include <QUuid>

#include <boost/tuple/tuple.hpp>

#include "Helpers/convertlist.hpp"
#include "Helpers/convertdict.hpp"
#include "Helpers/convertset.hpp"
#include "Base/convertpackedarray.hpp"

#include "sireqt_headers.h"

#include "Helpers/tuples.hpp"

#include "SireBase/packedarray2d.hpp"

using boost::python::register_tuple;

void register_SireQt_containers()
{
    register_list< QVector<quint32> >();
    register_list< QList<quint32> >();

    register_list< QVector<qint32> >();
    register_list< QList<qint32> >();
    
    register_list< QVector<quint64> >();
    register_list< QList<quint64> >();
    
    register_list< QVector<qint64> >();
    register_list< QList<qint64> >();

    register_list< QVector<float> >();
    register_list< QList<float> >();

    register_list< QVector<double> >();
    register_list< QList<double> >();

    register_list< QVector< QVector<double> > >();
    register_list< QVector< QVector< QVector<double> > > >();
    register_list< QVector< QVector< QVector< QVector<double> > > > >();

    register_list< QVector<bool> >();
    register_list< QVector< QVector<bool> > >();

    register_list< QVector<QString> >();
    register_list< QList<QString> >();
    register_list< QStringList >();

    register_list< QVector<QByteArray> >();
    register_list< QList<QByteArray> >();

    register_list< QList<QUuid> >();

    register_tuple< boost::tuple<double,double> >();
    register_tuple< boost::tuple<double,double,double> >();

    #if QT_VERSION >= 0x402000
    register_set< QSet<QString> >();
    
    register_dict< QMap<QString,QVariant> >();
    
    #else
    register_set< QSet<QString>, QString >();
    
    register_dict< QMap<QString,QVariant>, QString, QVariant >();
    
    #endif
}
