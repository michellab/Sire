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

#include "amberformat.h"

#include <QRegularExpression>

#include "SireIO/amberprm.h"

#include "SireIO/errors.h"

using namespace SireIO;
using namespace SireIO::detail;
using namespace SireMM;

AmberFormat::AmberFormat() : flag_type(SireIO::AmberPrm::UNKNOWN),
                num_values(0), field_width(0), point_width(0)
{}

AmberFormat::AmberFormat(SireIO::AmberPrm::FLAG_TYPE flag, int num,
            int field, int point)
        : flag_type(flag), num_values(num), field_width(field), point_width(point)
{}

/** Construct a format based on the Amber string description of that format. This
    is in the form "%FORMAT(20a4)" or "%FORMAT(10I8)" or "%FORMAT(5E16.8)" */
AmberFormat::AmberFormat(const QString &line)
{
    QRegularExpression re("%FORMAT\\((\\d+)(\\w)(\\d+)\\.?(\\d+)?\\)");

    auto m = re.match(line);

    if (not m.hasMatch())
    {
        throw SireIO::parse_error( QObject::tr(
                "Could not extract the format from the line '%1'. "
                "Expected to read a line similar to '%FORMAT(5E16.8)'.")
                    .arg(line.simplified()), CODELOC );
    }

    num_values = m.captured(1).toInt();
    field_width = m.captured(3).toInt();

    QString typ = m.captured(2).toLower();

    if (typ == "a")
        flag_type = AmberPrm::STRING;
    else if (typ == "i")
        flag_type = AmberPrm::INTEGER;
    else if (typ == "e")
        flag_type = AmberPrm::FLOAT;
    else
        flag_type = AmberPrm::UNKNOWN;

    if (not m.captured(4).isNull())
    {
        point_width = m.captured(4).toInt();
    }
}

AmberFormat::~AmberFormat()
{}

QString AmberFormat::toString() const
{
    switch(flag_type)
    {
        case AmberPrm::STRING:
            return QString("AmberFormat( %1 x string[width = %2] )")
                        .arg(num_values).arg(field_width);
        case AmberPrm::INTEGER:
            return QString("AmberFormat( %1 x integer[width = %2] )")
                        .arg(num_values).arg(field_width);
        case AmberPrm::FLOAT:
            return QString("AmberFormat( %1 x float[width = %2, precision = %3] )")
                        .arg(num_values).arg(field_width).arg(point_width);
        default:
            return QString("AmberFormat( UNKNOWN )");
    }
}

/** Return this format as an Amber format string, i.e. "%FORMAT(5E16.8)" */
QString AmberFormat::toAmberString() const
{
    switch(flag_type)
    {
        case AmberPrm::STRING:
            return QString("%FORMAT(%1a%2)").arg(num_values).arg(field_width);
        case AmberPrm::INTEGER:
            return QString("%FORMAT(%1I%2)").arg(num_values).arg(field_width);
        case AmberPrm::FLOAT:
            return QString("%FORMAT(%1E%2.%3)").arg(num_values).arg(field_width)
                                               .arg(point_width);
        default:
            return QString("%FORMAT(%1U%2)").arg(num_values).arg(field_width);
    }
}

namespace SireIO
{
namespace detail
{

/** Internal function used to write out an array of integers using the passed AmberFormat.
    This will include the Amber format string as a header if 'include_header'
    is true, and will append any errors encountered during writing to 'errors'
    if this is non-zero. This returns the text lines as a QStringList */
QStringList writeIntData(const QVector<qint64> &data, AmberFormat format,
                         QStringList *errors, bool include_header)
{
    QStringList lines;

    if (include_header)
    {
        lines.append( format.toAmberString() );
    }

    int n = 0;

    QString line;

    for (int i=0; i<data.count(); ++i)
    {
        const qint64 value = data.constData()[i];

        QString sval = QString("%1").arg(value, format.width());

        if (sval.length() > format.width())
        {
            //we couldn't fit this number into the specified width!
            if (errors)
                errors->append( QObject::tr(
                    "Could not write the integer at index %1, value '%2' into "
                    "the specified format %3.")
                        .arg(i).arg(value).arg(format.toString()) );

            sval = QString("%1").arg(0, format.width());
        }

        line += sval;
        n += 1;

        if (n == format.numValues())
        {
            lines.append(line);
            line = QString();
            n = 0;
        }
    }

    if (n > 0)
    {
        lines.append(line);
    }

    if (data.count() == 0)
    {
        lines.append(" ");
    }

    return lines;
}

/** Function to read and return an array of integers from the passed lines,
    according to the passed AmberFormat. The passed index specifies the
    indicies of the first and last lines in the list of lines to read,
    and the accumlated score and set of errors are updated assuming these
    are nonzero */
QVector<qint64> readIntData(const QVector<QString> &lines, AmberFormat format,
                            const QPair<qint64,qint64> &index,
                            double *total_score, QStringList *errors)
{
    QVector<qint64> data;

    double score = 0;

    //read in all of the lines...
    const int strt = index.first;
    const int end = index.first + index.second;

    const QString *l = lines.constData();

    data.reserve( (end-strt) * format.numValues() );

    for (int i=strt; i<end; ++i)
    {
        const QString &line = l[i];

        int nvalues = format.numValues(line);

        if (errors and nvalues < format.numValues() and i != end-1)
        {
            //one of the data lines has too little data
            errors->append( QObject::tr("Too few data values on line %1: "
                "Expected %2 values but only saw %3.")
                    .arg(i+1).arg(format.numValues()).arg(nvalues) );
        }

        double read_score = 0;

        for (int j=0; j<nvalues; ++j)
        {
            auto ref = line.midRef( j*format.width(), format.width() );

            if (not ref.isNull())
            {
                bool ok = true;
                qint64 value = ref.toLong(&ok);

                if (ok)
                {
                    data.append(value);
                    read_score += 1.0;
                }
                else if (errors)
                    errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 (%3) into an integer!")
                            .arg(i+1).arg(j+1).arg(ref.toString()) );
            }
            else if (errors)
            {
                errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 as the string is null!")
                            .arg(i+1).arg(j+1) );
            }
        }

        if (nvalues > 0)
        {
            score += read_score / nvalues;
        }
    }

    if (total_score)
        *total_score += score;

    return data;
}

/** Internal function used to write out an array of doubles using the passed AmberFormat.
    This will include the Amber format string as a header if 'include_header'
    is true, and will append any errors encountered during writing to 'errors'
    if this is non-zero. This returns the text lines as a QStringList */
QStringList writeFloatData(const QVector<double> &data, AmberFormat format,
                           QStringList *errors, bool include_header, char float_format)
{
    QStringList lines;

    if (include_header)
    {
        lines.append( format.toAmberString() );
    }

    int n = 0;

    QString line;

    for (int i=0; i<data.count(); ++i)
    {
        const double value = data.constData()[i];

        QString sval = QString("%1").arg(value, format.width(), float_format, format.pointWidth());

        if (sval.length() > format.width())
        {
            //we couldn't fit this number into the specified width!
            if (errors)
                errors->append( QObject::tr(
                    "Could not write the float at index %1, value '%2' into "
                    "the specified format %3.")
                        .arg(i).arg(value).arg(format.toString()) );

            sval = QString("%1").arg(0.0, format.width(), float_format, format.pointWidth());
        }

        line += sval;
        n += 1;

        if (n == format.numValues())
        {
            lines.append(line);
            line = QString();
            n = 0;
        }
    }

    if (n > 0)
    {
        lines.append(line);
    }

    if (data.count() == 0)
    {
        lines.append(" ");
    }

    return lines;
}

/** Function to read and return an array of doubles from the passed lines,
    according to the passed AmberFormat. The passed index specifies the
    indicies of the first and last lines in the list of lines to read,
    and the accumlated score and set of errors are updated assuming these
    are nonzero */
QVector<double> readFloatData(const QVector<QString> &lines, AmberFormat format,
                              const QPair<qint64,qint64> &index,
                              double *total_score, QStringList *errors)
{
    QVector<double> data;

    double score = 0;

    //read in all of the lines...
    const int strt = index.first;
    const int end = index.first + index.second;

    data.reserve( (end-strt) * format.numValues() );

    const QString *l = lines.constData();

    for (int i=strt; i<end; ++i)
    {
        const QString &line = l[i];

        int nvalues = format.numValues(line);

        if (errors and nvalues < format.numValues() and i != end-1)
        {
            //one of the data lines has too little data
            errors->append( QObject::tr("Too few data values on line %1: "
                "Expected %2 values but only saw %3.")
                    .arg(i+1).arg(format.numValues()).arg(nvalues) );
        }

        double read_score = 0;

        for (int j=0; j<nvalues; ++j)
        {
            auto ref = line.midRef( j*format.width(), format.width() );

            if (not ref.isNull())
            {
                bool ok = true;
                double value = ref.toDouble(&ok);

                if (ok)
                {
                    data.append(value);
                    read_score += 1;
                }
                else if (errors)
                    errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 (%3) into an integer!")
                            .arg(i+1).arg(j+1).arg(ref.toString()) );
            }
            else if (errors)
            {
                errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 as the string is null!")
                            .arg(i+1).arg(j+1) );
            }
        }

        if (nvalues > 0)
        {
            score += read_score / nvalues;
        }
    }

    if (total_score)
        *total_score += score;

    return data;
}

/** Internal function used to write out an array of strings using the passed AmberFormat.
    This will include the Amber format string as a header if 'include_header'
    is true, and will append any errors encountered during writing to 'errors'
    if this is non-zero. This returns the text lines as a QStringList */
QStringList writeStringData(const QVector<QString> &data, AmberFormat format,
                            QStringList *errors, bool include_header)
{
    QStringList lines;

    if (include_header)
    {
        lines.append( format.toAmberString() );
    }

    int n = 0;

    QString line;

    for (int i=0; i<data.count(); ++i)
    {
        const QString value = data.constData()[i];

        if (value.length() < format.width())
        {
            line += QString("%1").arg(value, -format.width());
        }
        else if (value.length() == format.width())
        {
            line += value;
        }
        else
        {
            if (errors)
                errors->append( QObject::tr(
                    "Could not fully write the string data at index %1, value '%2', "
                    "as it can't fit into the format %3.")
                        .arg(i).arg(value).arg(format.toString()) );

            line += value.mid(0,format.width());
        }

        n += 1;

        if (n == format.numValues())
        {
            lines.append(line);
            line = QString();
            n = 0;
        }
    }

    if (n > 0)
    {
        lines.append(line);
    }

    if (data.count() == 0)
    {
        lines.append(" ");
    }

    return lines;
}

/** Function to read and return an array of strings from the passed lines,
    according to the passed AmberFormat. The passed index specifies the
    indicies of the first and last lines in the list of lines to read,
    and the accumlated score and set of errors are updated assuming these
    are nonzero */
QVector<QString> readStringData(const QVector<QString> &lines, AmberFormat format,
                                const QPair<qint64,qint64> &index,
                                double *total_score, QStringList *errors)
{
    QVector<QString> data;

    double score = 0;

    //read in all of the lines...
    const int strt = index.first;
    const int end = index.first + index.second;

    data.reserve( (end-strt) * format.numValues() );

    const QString *l = lines.constData();

    for (int i=strt; i<end; ++i)
    {
        const QString &line = l[i];

        int nvalues = format.numValues(line);

        if (errors and nvalues < format.numValues() and i != end-1)
        {
            //one of the data lines has too little data
            errors->append( QObject::tr("Too few data values on line %1: "
                "Expected %2 values but only saw %3.")
                    .arg(i+1).arg(format.numValues()).arg(nvalues) );
        }

        double read_score = 0;

        for (int j=0; j<nvalues; ++j)
        {
            auto ref = line.midRef( j*format.width(), format.width() );

            if (not ref.isNull())
            {
                data.append( ref.toString() );
                read_score += 1;
            }
            else if (errors)
            {
                errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 as the string is null!")
                            .arg(i+1).arg(j+1) );
            }
        }

        if (nvalues > 0)
        {
            score += read_score / nvalues;
        }
    }

    if (total_score)
        *total_score += score;

    return data;
}

} // end of namespace detail
} // end of namespace SireMM
