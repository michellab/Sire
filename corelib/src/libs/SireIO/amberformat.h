/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#ifndef SIREIO_AMBERFORMAT_H
#define SIREIO_AMBERFORMAT_H

#include "SireIO/amberprm.h"

SIRE_BEGIN_HEADER

namespace SireIO
{

namespace detail
{

// The partial charges in the top file are not in electrons
static const double AMBERCHARGECONV = 18.2223;

static const double AMBER14COUL = 1.0 / 1.2 ;
static const double AMBER14LJ = 0.50 ;

/** Internal class used by the Amber IO helpers to hold format information

    @author Christopher Woods
*/
class AmberFormat
{
public:
    AmberFormat();
    
    AmberFormat(SireIO::AmberPrm::FLAG_TYPE flag, int num,
                int field, int point=0);
    
    AmberFormat(const QString &line);
    
    ~AmberFormat();
    
    QString toString() const;
    
    QString toAmberString() const;
    
    int width() const;
    
    int numValues() const;
    
    int pointWidth() const;
    
    int numValues( const QString &line ) const;

    SireIO::AmberPrm::FLAG_TYPE flagType() const;

private:
    SireIO::AmberPrm::FLAG_TYPE flag_type;
    int num_values;
    int field_width;
    int point_width;
};

QStringList writeIntData(const QVector<qint64> &data, AmberFormat format,
                         QStringList *errors=0, bool include_header=true);

QVector<qint64> readIntData(const QVector<QString> &lines, AmberFormat format,
                            const QPair<qint64,qint64> &index,
                            double *total_score, QStringList *errors=0);

QStringList writeFloatData(const QVector<double> &data, AmberFormat format,
                           QStringList *errors=0, bool include_header=true,
                           char float_format='E');

QVector<double> readFloatData(const QVector<QString> &lines, AmberFormat format,
                              const QPair<qint64,qint64> &index,
                              double *total_score, QStringList *errors);

QStringList writeStringData(const QVector<QString> &data, AmberFormat format,
                            QStringList *errors=0, bool include_header=true);

QVector<QString> readStringData(const QVector<QString> &lines, AmberFormat format,
                                const QPair<qint64,qint64> &index,
                                double *total_score, QStringList *errors);

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Internal function used to collapse an array of tuples of arrays ointo
    a tuple of arrays */
template<class T, class U>
std::tuple< QVector<T>, QVector<U> > collapseTuples(
                   const QVector< std::tuple< QVector<T>,QVector<U> > > &arrays )
{
    int nvals = 0;
    
    for (const auto array : arrays)
    {
        nvals += std::get<0>(array).count();
    }
    
    if (nvals == 0)
    {
        return std::make_tuple( QVector<T>(), QVector<U>() );
    }
    
    QVector<T> tvals;
    QVector<U> uvals;
    
    tvals.reserve(nvals);
    uvals.reserve(nvals);
    
    for (const auto array : arrays)
    {
        tvals += std::get<0>(array);
        uvals += std::get<1>(array);
    }
    
    return std::make_tuple(tvals, uvals);
}

/** Internal function used to collapse an array of tuples of arrays ointo
    a tuple of arrays */
template<class T, class U, class V, class W>
std::tuple< QVector<T>, QVector<U>, QVector<V>, QVector<W> > collapseTuples(
                   const QVector< std::tuple< QVector<T>,QVector<U>,
                                              QVector<V>,QVector<W> > > &arrays )
{
    int nvals = 0;
    
    for (const auto array : arrays)
    {
        nvals += std::get<0>(array).count();
    }
    
    if (nvals == 0)
    {
        return std::make_tuple( QVector<T>(), QVector<U>(), QVector<V>(), QVector<W>() );
    }
    
    QVector<T> tvals;
    QVector<U> uvals;
    QVector<V> vvals;
    QVector<W> wvals;
    
    tvals.reserve(nvals);
    uvals.reserve(nvals);
    vvals.reserve(nvals);
    wvals.reserve(nvals);
    
    for (const auto array : arrays)
    {
        tvals += std::get<0>(array);
        uvals += std::get<1>(array);
        vvals += std::get<2>(array);
        wvals += std::get<3>(array);
    }
    
    return std::make_tuple(tvals, uvals, vvals, wvals);
}

/** Internal function used to collapse an array of arrays of type T into
    a single array of type T */
template<class T>
QVector<T> collapse(const QVector< QVector<T> > &arrays)
{
    int nvals = 0;
    
    for (const auto array : arrays)
    {
        nvals += array.count();
    }
    
    if (nvals == 0)
    {
        return QVector<T>();
    }
    
    QVector<T> values;
    values.reserve(nvals);
    
    for (const auto array : arrays)
    {
        values += array;
    }
    
    return values;
}

/** The width of each field (number of columns) */
SIRE_ALWAYS_INLINE int AmberFormat::width() const
{
    return field_width;
}

/** The maximum number of items per line */
SIRE_ALWAYS_INLINE int AmberFormat::numValues() const
{
    return num_values;
}

/** The number of values after the decimal point for float values */
SIRE_ALWAYS_INLINE int AmberFormat::pointWidth() const
{
    return point_width;
}

/** The fortran flag format type */
SIRE_ALWAYS_INLINE SireIO::AmberPrm::FLAG_TYPE AmberFormat::flagType() const
{
    return flag_type;
}

/** Return the number of values to read from the passed line */
SIRE_ALWAYS_INLINE int AmberFormat::numValues( const QString &line ) const
{
    return qMin( num_values, line.length() / field_width );
}

#endif

} // end of namespace detail

} // end of namespace detail

SIRE_END_HEADER

#endif
