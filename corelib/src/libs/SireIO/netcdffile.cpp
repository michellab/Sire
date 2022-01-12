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

#include "netcdffile.h"

#include <QFileInfo>
#include <QDebug>

#include "SireBase/unittest.h"

#include "SireError/errors.h"

#ifdef SIRE_USE_NETCDF
    #include "netcdf.h"  // CONDITIONAL_INCLUDE
#endif

using namespace SireIO;

/////////////
///////////// Implemenetation of NetCDFDataInfo
/////////////


#ifdef SIRE_USE_NETCDF
    static QString nc_type_to_string(nc_type typ)
    {
        switch(typ)
        {
        case NC_BYTE:
            return QString("NC_BYTE");
        case NC_UBYTE:
            return QString("NC_UBYTE");
        case NC_CHAR:
            return QString("NC_CHAR");
        case NC_SHORT:
            return QString("NC_SHORT");
        case NC_USHORT:
            return QString("NC_USHORT");
        case NC_INT:
            return QString("NC_INT");
        case NC_UINT:
            return QString("NC_UINT");
        case NC_INT64:
            return QString("NC_INT64");
        case NC_UINT64:
            return QString("NC_UINT64");
        case NC_FLOAT:
            return QString("NC_FLOAT");
        case NC_DOUBLE:
            return QString("NC_DOUBLE");
        case NC_STRING:
            return QString("NC_STRING");
        default:
            return QString("NC_UNKNOWN");
        }
    }

    static int nc_type_to_size(nc_type typ)
    {
        switch(typ)
        {
            case NC_BYTE:
                return 1;
            case NC_UBYTE:
                return 1;
            case NC_CHAR:
                return 1;
            case NC_SHORT:
                return sizeof(short);
            case NC_USHORT:
                return sizeof(short);
            case NC_INT:
                return 4;
            case NC_UINT:
                return 4;
            case NC_INT64:
                return 8;
            case NC_UINT64:
                return 8;
            case NC_FLOAT:
                return 4;
            case NC_DOUBLE:
                return 8;
            case NC_STRING:
                return 1;
            default:
                return 0;
        }
    }

    static nc_type string_to_nc_type(const QString &typ)
    {
        if (typ == "NC_BYTE") return NC_BYTE;
        else if (typ == "NC_UBYTE") return NC_UBYTE;
        else if (typ == "NC_CHAR") return NC_CHAR;
        else if (typ == "NC_SHORT") return NC_SHORT;
        else if (typ == "NC_USHORT") return NC_USHORT;
        else if (typ == "NC_INT") return NC_INT;
        else if (typ == "NC_UINT") return NC_UINT;
        else if (typ == "NC_INT64") return NC_INT64;
        else if (typ == "NC_UINT64") return NC_UINT64;
        else if (typ == "NC_FLOAT") return NC_FLOAT;
        else if (typ == "NC_DOUBLE") return NC_DOUBLE;
        else if (typ == "NC_STRING") return NC_STRING;
        else
        {
            throw SireError::io_error( QObject::tr(
                "Unrecognised NetCDF type - %1").arg(typ), CODELOC );

            return 0;
        }
    }

    static nc_type qvariant_to_nc_type(const QVariant &typ)
    {
        const QString name = typ.typeName();

        if (name == "QString")
        {
            return NC_CHAR;
        }
        else if (name == "double")
        {
            return NC_DOUBLE;
        }
        else if (name == "float")
        {
            return NC_FLOAT;
        }
        else
        {
            throw SireError::io_error( QObject::tr(
                    "Unable to convert the QVariant type %1 to an NC_TYPE")
                        .arg(name), CODELOC );
        }

        return 0;
    }
#else
    static QString nc_type_to_string(int typ)
    {
        return QString("NC_UNKNOWN");
    }

    static int nc_type_to_size(int typ)
    {
        return 0;
    }

    static int string_to_nc_type(const QString &typ)
    {
        throw SireError::io_error( QObject::tr(
            "Unrecognised NetCDF type - %1").arg(typ), CODELOC );
    }

    static int qvariant_to_nc_type(const QVariant &typ)
    {
        throw SireError::io_error( QObject::tr(
                "Unable to convert the QVariant type %1 to an NC_TYPE")
                    .arg(typ.typeName()), CODELOC );
        return 0;
    }
#endif

static QVector<QVariant> extract_values(const QByteArray &memdata, int nc_type, int nvals)
{
    //ensure there is enough space in the data
    #ifdef SIRE_USE_NETCDF
    if (nvals > (memdata.count() / nc_type_to_size(nc_type)))
    {
        throw SireError::invalid_arg( QObject::tr(
                "You cannot read %1 values, as the array has only allocated space (%2) "
                "for %3 values of type %4.")
                    .arg(nvals).arg(memdata.count())
                    .arg(memdata.count() / nc_type_to_size(nc_type))
                    .arg(nc_type_to_string(nc_type)), CODELOC );
    }
    #endif

    QVector<QVariant> vals;

    if (nvals > 0)
    {
    #ifdef SIRE_USE_NETCDF
        vals = QVector<QVariant>(nvals);

        const char *data = memdata.constData();

        switch(nc_type)
        {
            case NC_BYTE:
            case NC_UBYTE:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( memdata[i] );
                }
                break;
            case NC_CHAR:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( memdata[i] );
                }
                break;
            case NC_SHORT:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( qint32( *(reinterpret_cast<const short*>(data) + i) ) );
                }
                break;
            case NC_USHORT:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( quint32( *(reinterpret_cast<const unsigned short*>(data) + i) ) );
                }
                break;
            case NC_INT:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( qint32( *(reinterpret_cast<const qint32*>(data) + i) ) );
                }
                break;
            case NC_UINT:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( quint32( *(reinterpret_cast<const quint32*>(data) + i) ) );
                }
                break;
            case NC_INT64:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( qint64( *(reinterpret_cast<const qint64*>(data) + i) ) );
                }
                break;
            case NC_UINT64:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( quint32( *(reinterpret_cast<const quint64*>(data) + i) ) );
                }
                break;
            case NC_FLOAT:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( float( *(reinterpret_cast<const float*>(data) + i) ) );
                }
                break;
            case NC_DOUBLE:
                for (int i=0; i<nvals; ++i)
                {
                    vals.append( double( *(reinterpret_cast<const double*>(data) + i) ) );
                }
                break;
            case NC_STRING:
            default:
                break;
        }

    #endif
    }

    return vals;
}

static QVariant extract_value(const QByteArray &memdata, int nc_type)
{
    //see if this is a single value or an array
    #ifdef SIRE_USE_NETCDF
    int nvals = memdata.count() / nc_type_to_size(nc_type);
    #else
    int nvals = 0;
    #endif

    if (nvals > 0)
    {
    #ifdef SIRE_USE_NETCDF

        const char *data = memdata.constData();

        switch(nc_type)
        {
            case NC_BYTE:
            case NC_UBYTE:
                return QVariant(memdata);
            case NC_CHAR:
                return QVariant( QString::fromUtf8(data) );
            case NC_SHORT:
                if (nvals > 1)
                {
                    QList<QVariant> vals;
                    for (int i=0; i<nvals; ++i)
                    {
                        vals.append( qint32( *(reinterpret_cast<const short*>(data) + i) ) );
                    }

                    return QVariant(vals);
                }
                else
                {
                    return QVariant( qint32( *(reinterpret_cast<const short*>(data)) ) );
                }
            case NC_USHORT:
                if (nvals > 1)
                {
                    QList<QVariant> vals;
                    for (int i=0; i<nvals; ++i)
                    {
                        vals.append( quint32( *(
                                reinterpret_cast<const unsigned short*>(data) + i) ) );
                    }

                    return QVariant(vals);
                }
                else
                {
                    return QVariant( quint32( *(reinterpret_cast<const unsigned short*>(data)) ) );
                }
            case NC_INT:
                if (nvals > 1)
                {
                    QList<QVariant> vals;
                    for (int i=0; i<nvals; ++i)
                    {
                        vals.append( quint32( *(reinterpret_cast<const qint32*>(data) + i) ) );
                    }

                    return QVariant(vals);
                }
                else
                {
                    return QVariant( *(reinterpret_cast<const qint32*>(data)) );
                }
            case NC_UINT:
                if (nvals > 1)
                {
                    QList<QVariant> vals;
                    for (int i=0; i<nvals; ++i)
                    {
                        vals.append( quint32( *(reinterpret_cast<const quint32*>(data) + i) ) );
                    }

                    return QVariant(vals);
                }
                else
                {
                    return QVariant( *(reinterpret_cast<const quint32*>(data)) );
                }
            case NC_INT64:
                if (nvals > 1)
                {
                    QList<QVariant> vals;
                    for (int i=0; i<nvals; ++i)
                    {
                        vals.append( qint64( *(reinterpret_cast<const qint64*>(data) + i) ) );
                    }

                    return QVariant(vals);
                }
                else
                {
                    return QVariant( *(reinterpret_cast<const qint64*>(data)) );
                }
            case NC_UINT64:
                if (nvals > 1)
                {
                    QList<QVariant> vals;
                    for (int i=0; i<nvals; ++i)
                    {
                        vals.append( quint64( *(reinterpret_cast<const quint64*>(data) + i) ) );
                    }

                    return QVariant(vals);
                }
                else
                {
                    return QVariant( *(reinterpret_cast<const quint64*>(data)) );
                }
            case NC_FLOAT:
                if (nvals > 1)
                {
                    QList<QVariant> vals;
                    for (int i=0; i<nvals; ++i)
                    {
                        vals.append( float( *(reinterpret_cast<const float*>(data) + i) ) );
                    }

                    return QVariant(vals);
                }
                else
                {
                    return QVariant( *(reinterpret_cast<const float*>(data)) );
                }
            case NC_DOUBLE:
                if (nvals > 1)
                {
                    QList<QVariant> vals;
                    for (int i=0; i<nvals; ++i)
                    {
                        vals.append( double( *(reinterpret_cast<const double*>(data) + i) ) );
                    }

                    return QVariant(vals);
                }
                else
                {
                    return QVariant( *(reinterpret_cast<const double*>(data)) );
                }
            case NC_STRING:
            default:
                return QVariant();
        }

    #endif
    }

    return QVariant();
}

/** Constructor - completely null data type */
NetCDFDataInfo::NetCDFDataInfo() : idnum(-1), xtyp(-1)
{}

/** Internal constructor used by NetCDFFile to construct from the passed data */
NetCDFDataInfo::NetCDFDataInfo(int idn, QString name, int tp,
                               QStringList dim_ns, QList<int> dim_sz,
                               QStringList att_ns, QList<int> att_ts,
                               QList<QVariant> att_vs)
{
    idnum = idn;
    nme = name;
    xtyp = tp;
    dim_names = dim_ns;
    dim_sizes = dim_sz;
    att_names = att_ns;
    att_types = att_ts;
    att_values = att_vs;

    if (idnum < 0 or xtyp < 0)
    {
        throw SireError::invalid_arg( QObject::tr(
                "You cannot construct a NetCDFDataInfo with a negative ID number (%1) "
                "or negative xtype (%2).")
                    .arg(idnum).arg(xtyp), CODELOC );
    }

    if (dim_names.count() != dim_sizes.count())
    {
        throw SireError::invalid_arg( QObject::tr(
                "The number of dimension names (%1) must equal the number of "
                "dimension sizes (%2)").arg(dim_names.count())
                                       .arg(dim_sizes.count()), CODELOC );
    }

    if (att_names.count() != att_types.count() or att_names.count() != att_values.count())
    {
        throw SireError::invalid_arg( QObject::tr(
                "The number of attribute names (%1), types (%2) and values (%3) must "
                "all be equal!").arg(att_names.count())
                                .arg(att_types.count())
                                .arg(att_values.count()), CODELOC );
    }
}

NetCDFDataInfo::NetCDFDataInfo(int idnum, QString name, const QString &xtyp,
                               QStringList dim_names, QList<int> dim_sizes,
                               QStringList att_names, QList<int> att_types,
                               QList<QVariant> att_values)
               : NetCDFDataInfo(idnum, name, string_to_nc_type(xtyp),
                                dim_names, dim_sizes, att_names, att_types,
                                att_values)
{}

/** Copy constructor */
NetCDFDataInfo::NetCDFDataInfo(const NetCDFDataInfo &other)
               : nme(other.nme), dim_names(other.dim_names), dim_sizes(other.dim_sizes),
                 att_names(other.att_names), att_types(other.att_types),
                 att_values(other.att_values),
                 idnum(other.idnum), xtyp(other.xtyp)
{}

/** Destructor */
NetCDFDataInfo::~NetCDFDataInfo()
{}

/** Return whether or not this data type info is null */
bool NetCDFDataInfo::isNull() const
{
    return idnum == -1;
}

/** Return a string representation of this data info */
QString NetCDFDataInfo::toString() const
{
    #ifdef SIRE_USE_NETCDF
    if (isNull())
        return QObject::tr("NetCDFDataInfo::null");
    else if (dim_names.isEmpty())
    {
        if (att_names.isEmpty())
        {
            return QObject::tr("NetCDFDataInfo( %1 = %2[%3]() )")
                    .arg(idnum).arg(nme).arg(this->type());
        }
        else
        {
            QStringList atts;
            for (int i=0; i<att_names.count(); ++i)
            {
                atts.append( QString("%1=[%2 - %3]").arg(att_names[i])
                                                    .arg(nc_type_to_string(att_types[i]))
                                                    .arg(att_values[i].toString()) );
            }

            return QObject::tr("NetCDFDataInfo( %1 = %2[%3](), attributes:{ %4 } )")
                    .arg(idnum).arg(nme).arg(this->type()).arg(atts.join(", "));
        }
    }
    else
    {
        QStringList dims;
        for (int i=0; i<dim_names.count(); ++i)
        {
            dims.append( QString("%1:%2").arg(dim_names[i]).arg(dim_sizes[i]) );
        }

        if (att_names.isEmpty())
        {
            return QObject::tr("NetCDFDataInfo( %1 = %2[%3](%4) )")
                    .arg(idnum).arg(nme).arg(this->type()).arg(dims.join(","));
        }
        else
        {
            QStringList atts;
            for (int i=0; i<att_names.count(); ++i)
            {
                atts.append( QString("%1=[%2 - %3]").arg(att_names[i])
                                                    .arg(nc_type_to_string(att_types[i]))
                                                    .arg(att_values[i].toString()) );
            }

            return QObject::tr("NetCDFDataInfo( %1 = %2[%3](%4), attributes:{ %5 } )")
                    .arg(idnum).arg(nme).arg(this->type()).arg(dims.join(","))
                    .arg(atts.join(", "));
        }
    }
    #else
        return QObject::tr("NetCDFDataInfo::null (netcdf is not supported)");
    #endif
}

/** Return the data type of this piece of data. This is a string
    version of the NC_TYPE, e.g. NC_FLOAT, NC_STRING etc. */
QString NetCDFDataInfo::type() const
{
    #ifdef SIRE_USE_NETCDF
    return nc_type_to_string(xtyp);
    #else
    return QString();
    #endif
}

/** Return the size in bytes of a variable of this type */
int NetCDFDataInfo::typeSize() const
{
    #ifdef SIRE_USE_NETCDF
    return nc_type_to_size(xtyp);
    #else
    return 0;
    #endif
}

/** Return the number of values that should be held by this data */
int NetCDFDataInfo::nValues() const
{
    int base = 1;

    for (const auto &sz : dim_sizes)
    {
        base *= sz;
    }

    return base;
}

/** Assert that the number of values that can be held by this data is 'n' */
void NetCDFDataInfo::assertNValuesEquals(int n) const
{
    if (n != this->nValues())
    {
        throw SireError::incompatible_error( QObject::tr(
                "Cannot store %1 items in this piece of NetCDF data, as it is dimensioned "
                "to store %2 values. %3")
                    .arg(n).arg(this->nValues()).arg(this->toString()), CODELOC );
    }
}

/** Return the total size of the data to be loaded, in bytes */
int NetCDFDataInfo::dataSize() const
{
    return typeSize() * nValues();
}

/** Return all of the names of the attributes */
QStringList NetCDFDataInfo::attributeNames() const
{
    return att_names;
}

/** Return the value of the attribute called 'name'. This returns QVariant::null
    if this attribute doesn't exist */
QVariant NetCDFDataInfo::attribute(const QString &name) const
{
    for (int i=0; i<att_names.count(); ++i)
    {
        if (name == att_names[i])
        {
            return att_values[i];
        }
    }

    return QVariant();
}

/** Return the type of the attribute called 'name'. This returns a string version
    of the NC_TYPE of the attribute, i.e. NC_DOUBLE or NC_CHAR. This returns
    an empty string if there is not attribute with this name */
QString NetCDFDataInfo::attributeType(const QString &name) const
{
    #ifdef SIRE_USE_NETCDF
    for (int i=0; i<att_names.count(); ++i)
    {
        if (name == att_names[i])
        {
            return nc_type_to_string(att_types[i]);
        }
    }
    #endif

    return QString();
}

/** Return a hash of all of the values of all attributes */
QHash<QString,QVariant> NetCDFDataInfo::attributes() const
{
    QHash<QString,QVariant> atts;

    if (not att_names.isEmpty())
    {
        atts.reserve(att_names.count());

        for (int i=0; i<att_names.count(); ++i)
        {
            atts.insert(att_names[i], att_values[i]);
        }
    }

    return atts;
}

/** Return a hash of all of the attribute types of the attributes */
QHash<QString,QString> NetCDFDataInfo::attributeTypes() const
{
    QHash<QString,QString> atts;

    #ifdef SIRE_USE_NETCDF
    if (not att_names.isEmpty())
    {
        atts.reserve(att_names.count());

        for (int i=0; i<att_names.count(); ++i)
        {
            atts.insert(att_names[i], nc_type_to_string(att_types[i]));
        }
    }
    #endif

    return atts;
}

/////////////
///////////// Implemenetation of NetCDFData
/////////////

/** Constructor */
NetCDFData::NetCDFData() : NetCDFDataInfo()
{}

/** Copy constructor */
NetCDFData::NetCDFData(const NetCDFData &other)
           : NetCDFDataInfo(other), memdata(other.memdata)
{}

/** Destructor */
NetCDFData::~NetCDFData()
{}

/** Internal constructor used by NetCDFFile */
NetCDFData::NetCDFData(const NetCDFDataInfo &info)
           : NetCDFDataInfo(info)
{}

/** Internal function used by NetCDFFile to set the data */
void NetCDFData::setData(const QByteArray &data)
{
    memdata = data;
}

/** Internal function to extract attribute information */
QStringList NetCDFData::get_attribute_names(const QHash<QString,QVariant> &attributes)
{
    if (attributes.isEmpty())
        return QStringList();

    QStringList names = attributes.keys();
    std::sort(names.begin(), names.end());
    return names;
}

/** Internal function to extract attribute information */
QList<int> NetCDFData::get_attribute_types(const QHash<QString,QVariant> &attributes)
{
    if (attributes.isEmpty())
        return QList<int>();

    QList<int> typs;

    for (const auto &name : get_attribute_names(attributes))
    {
        typs.append( qvariant_to_nc_type( attributes.value(name) ) );
    }

    return typs;
}

/** Internal function to extract attribute information */
QList<QVariant> NetCDFData::get_attribute_values(const QHash<QString,QVariant> &attributes)
{
    if (attributes.isEmpty())
        return QList<QVariant>();

    QList<QVariant> vals;

    for (const auto &name : get_attribute_names(attributes))
    {
        vals.append( attributes.value(name) );
    }

    return vals;
}

/** Return the data as an array of QVariants */
QVector<QVariant> NetCDFData::toArray() const
{
    return extract_values(memdata, xtyp, this->nValues());
}

/** Return the data cast as an array of floats */
QVector<float> NetCDFData::toFloatArray() const
{
    #ifdef SIRE_USE_NETCDF
    const int nvals = this->nValues();
    QVector<float> values( nvals );

    if (xtyp == NC_FLOAT)
    {
        const char *data = memdata.constData();

        for (int i=0; i<nvals; ++i)
        {
            values[i] = *(reinterpret_cast<const float*>(data) + i);
        }
    }
    else if (xtyp == NC_DOUBLE)
    {
        const char *data = memdata.constData();

        for (int i=0; i<nvals; ++i)
        {
            values[i] = *(reinterpret_cast<const double*>(data) + i);
        }
    }
    else
    {
        //need to go via the QVariant list
        const auto vars = this->toArray();

        for (int i=0; i<this->nValues(); ++i)
        {
            values[i] = vars[i].toFloat();
        }
    }

    return values;
    #else
    return QVector<float>();
    #endif
}

/** Return the data cast as an array of doubles */
QVector<double> NetCDFData::toDoubleArray() const
{
    #ifdef SIRE_USE_NETCDF
    const int nvals = this->nValues();
    QVector<double> values( nvals );

    if (xtyp == NC_FLOAT)
    {
        const char *data = memdata.constData();

        for (int i=0; i<nvals; ++i)
        {
            values[i] = *(reinterpret_cast<const float*>(data) + i);
        }
    }
    else if (xtyp == NC_DOUBLE)
    {
        const char *data = memdata.constData();

        for (int i=0; i<nvals; ++i)
        {
            values[i] = *(reinterpret_cast<const double*>(data) + i);
        }
    }
    else
    {
        //need to go via the QVariant list
        const auto vars = this->toArray();

        for (int i=0; i<this->nValues(); ++i)
        {
            values[i] = vars[i].toDouble();
        }
    }

    return values;
    #else
    return QVector<double>();
    #endif
}

/////////////
///////////// Implemenetation of NetCDFFile
/////////////

static void assert_no_netcdf_error(int errnum)
{
    #ifdef SIRE_USE_NETCDF
        QString err;

        switch(errnum)
        {
            case NC_NOERR:
                return;
            case NC_EHDFERR:
                err = QObject::tr("HDF5 error!");
                break;
            case NC_EDIMMETA:
                err = QObject::tr("NetCDF-4 dimension metadata error!");
                break;
            case NC_EBADID:
                err = QObject::tr("Not a netcdf id");
                break;
            case NC_ENFILE:
                err = QObject::tr("Too many netcdfs open");
                break;
            case NC_EEXIST:
                err = QObject::tr("netcdf file exists");
                break;
            case NC_EINVAL:
                err = QObject::tr("Invalid Argument");
                break;
            case NC_EPERM:
                err = QObject::tr("Write to read only");
                break;
            case NC_ENOTINDEFINE:
                err = QObject::tr("Operation not allowed in data mode");
                break;
            case NC_EINDEFINE:
                err = QObject::tr("Operation not allowed in define mode");
                break;
            case NC_EINVALCOORDS:
                err = QObject::tr("Index exceeds dimension bound");
                break;
            case NC_EMAXDIMS:
                err = QObject::tr("NC_MAX_DIMS exceeded");
                break;
            case NC_ENAMEINUSE:
                err = QObject::tr("String match to name in use");
                break;
            case NC_ENOTATT:
                err = QObject::tr("Attribute not found");
                break;
            case NC_EMAXATTS:
                err = QObject::tr("NC_MAX_ATTRS exceeded");
                break;
            case NC_EBADTYPE:
                err = QObject::tr("Not a netcdf data type");
                break;
            case NC_EBADDIM:
                err = QObject::tr("Invalid dimension id or name");
                break;
            case NC_EUNLIMPOS:
                err = QObject::tr("NC_UNLIMITED in the wrong index");
                break;
            case NC_EMAXVARS:
                err = QObject::tr("NC_MAX_VARS exceeded");
                break;
            case NC_ENOTVAR:
                err = QObject::tr("Variable not found");
                break;
            case NC_EGLOBAL:
                err = QObject::tr("Action prohibited on NC_GLOBAL varid");
                break;
            case NC_ENOTNC:
                err = QObject::tr("Not a netcdf file");
                break;
            case NC_ESTS:
                err = QObject::tr("In Fortran, string too short");
                break;
            case NC_EMAXNAME:
                err = QObject::tr("NC_MAX_NAME exceeded");
                break;
            case NC_EUNLIMIT:
                err = QObject::tr("NC_UNLIMITED size already in use");
                break;
            case NC_ENORECVARS:
                err = QObject::tr("nc_rec op when there are no record vars");
                break;
            case NC_ECHAR:
                err = QObject::tr("Attempt to convert between text & numbers");
                break;
            case NC_EEDGE:
                err = QObject::tr("Edge+start exceeds dimension bound");
                break;
            case NC_ESTRIDE:
                err = QObject::tr("Illegal stride");
                break;
            case NC_EBADNAME:
                err = QObject::tr("Attribute or variable name contains illegal characters");
                break;
            case NC_ERANGE:
                err = QObject::tr("Math result not representable");
                break;
            case NC_ENOMEM:
                err = QObject::tr("Memory allocation (malloc) failure");
                break;
            case NC_EVARSIZE:
                err = QObject::tr("One or more variable sizes violate format constraints");
                break;
            case NC_EDIMSIZE:
                err = QObject::tr("Invalid dimension size");
                break;
            case NC_ETRUNC:
                err = QObject::tr("File likely truncated or possibly corrupted");
                break;
            default:
                err = QObject::tr("NetCDF experienced an unknown error! %1").arg(errnum);
        }

        throw SireError::io_error( QObject::tr(
            "NetCDF experienced an error: %1 (%2)").arg(err).arg(errnum), CODELOC );
    #else
        throw SireError::io_error( QObject::tr(
            "NetCDF experienced an error as it is not compiled and supported with "
            "this version of Sire."), CODELOC );
    #endif
}

/** Constructor */
NetCDFFile::NetCDFFile() : hndl(-1)
{}

/** Function used to call and check the output of netcdf operations */
int NetCDFFile::call_netcdf_function(std::function<int()> func, int ignored_error) const
{
    QMutexLocker lkr( const_cast<QMutex*>(&mutex) );
    int err = func();

    if (err != ignored_error)
    {
        assert_no_netcdf_error(err);
    }

    return err;
}

/** Construct to open the file 'filename' in read-only mode */
NetCDFFile::NetCDFFile(const QString &filename) : fname(filename)
{
    #ifdef SIRE_USE_NETCDF
        QByteArray c_filename = filename.toUtf8();
        call_netcdf_function(
            [&](){ return nc_open(c_filename.constData(), NC_NOWRITE, &hndl); }
                            );
    #else
        throw SireError::unsupported( QObject::tr(
                "Software is missing NetCDF support, so cannot read the NetCDF file '%1'")
                        .arg(filename), CODELOC );
    #endif
}

/** Construct to create the file 'filename' in write-only mode. If 'overwrite_file'
    is true, then this will overwrite any existing file. If use_64bit_offset is
    true, then a 64bit offset format file is created. If use_netcdf4 is true,
    then a NetCDF4 file is created (otherwise, NetCDF3 is used) */
NetCDFFile::NetCDFFile(const QString &filename, bool overwrite_file,
                       bool use_64bit_offset, bool use_netcdf4) : fname(filename)
{
    #ifdef SIRE_USE_NETCDF
        QFileInfo file(filename);

        if (file.exists())
        {
            if (not overwrite_file)
            {
                throw SireError::io_error( QObject::tr(
                        "Cannot create the NetCDF file '%1' as it already exists, and "
                        "the software is not allowed to overwrite an existing file!")
                            .arg(filename), CODELOC );
            }

            if (file.isDir())
            {
                throw SireError::io_error( QObject::tr(
                        "Cannot create the NetCDF file '%1' as it already exists and "
                        "is a directory!")
                            .arg(filename), CODELOC );
            }

            if (not file.isWritable())
            {
                throw SireError::io_error( QObject::tr(
                        "Cannot create the NetCDF file '%1' as it already exists and "
                        "is read-only")
                            .arg(filename), CODELOC );
            }
        }

        int flags = NC_WRITE;

        if (not overwrite_file)
        {
            flags |= NC_NOCLOBBER;
        }

        if (use_64bit_offset)
        {
            flags |= NC_64BIT_OFFSET;
        }

        if (use_netcdf4)
        {
            flags |= NC_NETCDF4;
        }

        QByteArray c_filename = file.absoluteFilePath().toUtf8();
        call_netcdf_function(
            [&](){ return nc_create(c_filename.constData(), flags, &hndl); }
                            );

    #else
        throw SireError::unsupported( QObject::tr(
                "Software is missing NetCDF support, so cannot write the NetCDF file '%1'")
                    .arg(filename), CODELOC );
    #endif
}

/** Destructor - this will close the NetCDFFile */
NetCDFFile::~NetCDFFile()
{
    #ifdef SIRE_USE_NETCDF
        if (hndl != -1)
        {
            nc_close(hndl);
            hndl = -1;
        }
    #endif
}

/** Return the full set of names and data types for all of the variables
    in the file */
QHash<QString,NetCDFDataInfo> NetCDFFile::getVariablesInfo() const
{
    QHash<QString,NetCDFDataInfo> vars;

    #ifdef SIRE_USE_NETCDF
    if (hndl != -1)
    {
        char *tmp_name = new char[NC_MAX_NAME+1];
        int *dim_ids = new int[NC_MAX_VAR_DIMS];

        try
        {
            int i=0;

            while (true)
            {
                nc_type var_type;
                int ndims;
                int natts;

                if (call_netcdf_function([&](){
                    return nc_inq_var(hndl, i, tmp_name, &var_type, &ndims, dim_ids, &natts); },
                    NC_ENOTVAR) == NC_ENOTVAR)
                {
                    break;
                }

                QString var_name = QString::fromUtf8(tmp_name);

                QStringList dim_names;
                QList<int> dim_sizes;

                for (int j=0; j<ndims; ++j)
                {
                    size_t dim_len;
                    call_netcdf_function([&](){
                        return nc_inq_dim(hndl, dim_ids[j], tmp_name, &dim_len);});

                    dim_names.append( QString::fromUtf8(tmp_name) );
                    dim_sizes.append(dim_len);
                }

                QStringList att_names;
                QList<int> att_types;
                QList<QVariant> att_values;

                if (natts > 0)
                {
                    for (int j=0; j<natts; ++j)
                    {
                        //first read in the name of the attribute
                        call_netcdf_function( [&]()
                                { return nc_inq_attname(hndl, i, j, tmp_name); } );

                        QString attname = QString::fromUtf8(tmp_name);

                        //now read in metadata about the attribute
                        nc_type xtype;
                        size_t len;
                        call_netcdf_function( [&]()
                            { return nc_inq_att(hndl, i, tmp_name, &xtype, &len); } );

                        //now read in the value of the attribute
                        QByteArray memdata;

                        if (xtype == NC_CHAR)
                        {
                            memdata.fill('\0', len+1);
                        }
                        else
                        {
                            memdata.resize( nc_type_to_size(xtype) * len );
                        }

                        call_netcdf_function( [&]()
                            {return nc_get_att(hndl, i, tmp_name, memdata.data()); } );

                        QVariant val = extract_value(memdata, xtype);

                        att_names.append(attname);
                        att_types.append( int(xtype) );
                        att_values.append(val);
                    }
                }

                vars.insert( var_name, NetCDFDataInfo(i,var_name,var_type,
                                                      dim_names,dim_sizes,
                                                      att_names,att_types,att_values) );

                i += 1;
            }
        }
        catch(...)
        {
            delete[] tmp_name;
            delete[] dim_ids;
            throw;
        }

        delete[] tmp_name;
        delete[] dim_ids;
    }
    #endif

    return vars;
}

/** Write all of the passed data to the file */
void NetCDFFile::writeData(const QHash<QString,QString> &globals,
                           const QHash<QString,NetCDFData> &variable_data)
{
    if (hndl != -1)
    {
    #ifdef SIRE_USE_NETCDF

        //always write the data in alphabetical order, so that
        //we get the same file every time
        QStringList variables = variable_data.keys();
        std::sort(variables.begin(), variables.end());

        //map of dimension names to IDs
        QHash<QString,int> dimension_ids;

        //map of variable info names to IDs
        QHash<QString,int> var_ids;

        //first we have to set up the NetCDF file, so get all of the
        //dimensions, variables and attributes
        {
            QHash<QString,int> dimensions;

            for (const auto &variable : variables)
            {
                const auto vardata = variable_data[variable];

                const auto dims = vardata.dimensions();
                const auto dim_sizes = vardata.dimensionSizes();

                SireBase::assert_equal( dims.count(), dim_sizes.count(), CODELOC );

                for (int i=0; i<dims.count(); ++i)
                {
                    if (not dimensions.contains(dims[i]))
                    {
                        dimensions.insert(dims[i], dim_sizes[i]);
                    }
                    else
                    {
                        SireBase::assert_equal( dim_sizes[i], dimensions[dims[i]], CODELOC );
                    }
                }
            }

            //now write all of the dimensions to the file, saving the ID of each dimension
            QStringList dims = dimensions.keys();
            std::sort(dims.begin(), dims.end());

            dimension_ids.reserve(dims.count());

            for (const auto &dim : dims)
            {
                int idp;
                const QByteArray c_dim = dim.toUtf8();

                if (c_dim.length() > NC_MAX_NAME)
                {
                    throw SireError::io_error( QObject::tr(
                            "The length of the name of the dimension '%1' (%2) cannot "
                            "be greater than NC_MAX_NAME, which is %2.")
                                .arg(dim).arg(c_dim.length())
                                .arg(NC_MAX_NAME), CODELOC );
                }

                call_netcdf_function( [&](){ return nc_def_dim(hndl, c_dim.constData(),
                                                               dimensions[dim], &idp); } );

                dimension_ids.insert(dim, idp);
            }
        }

        //now go through and save info about all of the variables
        var_ids.reserve(variables.count());

        for (const auto &variable : variables)
        {
            const auto vardata = variable_data[variable];

            //get the IDs of all of the dimensions
            const auto dims = vardata.dimensions();

            QVarLengthArray<int,8> dim_ids;

            for (const auto &dim : dims)
            {
                dim_ids.append( dimension_ids[dim] );
            }

            int ndims = dims.count();

            int idp;
            const QByteArray c_var = variable.toUtf8();

            if (c_var.length() > NC_MAX_NAME)
            {
                throw SireError::io_error( QObject::tr(
                        "The length of the name of the variable '%1' (%2) cannot "
                        "be greater than NC_MAX_NAME, which is %2.")
                            .arg(variable).arg(c_var.length())
                            .arg(NC_MAX_NAME), CODELOC );
            }

            //get the type of the data
            nc_type xtyp = string_to_nc_type( vardata.type() );

            //now write the variable info to the netcdf file, saving the variable ID
            call_netcdf_function( [&](){ return nc_def_var(hndl, c_var.constData(),
                                                           xtyp, ndims, dim_ids.constData(),
                                                           &idp); } );

            var_ids.insert(variable, idp);
        }

        //now write all of the global attributes
        QStringList global_attributes = globals.keys();
        std::sort(global_attributes.begin(), global_attributes.end());

        for (const auto &global_attribute : global_attributes)
        {
            QByteArray c_name = global_attribute.toUtf8();
            QByteArray c_att = globals.value(global_attribute).toUtf8();

            //save the attribute to the file
            call_netcdf_function( [&]()
                { return nc_put_att_text(hndl, NC_GLOBAL, c_name.constData(),
                                         c_att.count(), c_att.constData()); } );
        }

        //now write all of the attributes
        for (const auto &variable : variables)
        {
            const auto vardata = variable_data[variable];

            const int id = var_ids.value(variable, -1);

            if (id < 0)
            {
                throw SireError::program_bug( QObject::tr(
                        "How can the variable '%1' have a negative ID? %2")
                            .arg(variable).arg(id), CODELOC );
            }

            //find all of the attributes of this variable
            for (const auto &attribute : vardata.attributeNames())
            {
                const auto att_value = vardata.attribute(attribute);
                const auto att_type = vardata.attributeType(attribute);

                const QByteArray c_attname = attribute.toUtf8();

                nc_type xtyp = string_to_nc_type(att_type);

                switch(xtyp)
                {
                    case NC_CHAR:
                    {
                        const QByteArray c_att = att_value.toString().toUtf8();
                        size_t len = c_att.count();

                        call_netcdf_function( [&](){
                            return nc_put_att(hndl, id, c_attname.constData(),
                                              xtyp, len, c_att.constData()); } );
                        break;
                    }
                    case NC_DOUBLE:
                    {
                        double val = att_value.toDouble();
                        size_t len = 1;

                        call_netcdf_function( [&](){
                            return nc_put_att(hndl, id, c_attname.constData(),
                                              xtyp, len, &val); } );
                        break;
                    }
                    case NC_FLOAT:
                    {
                        float val = att_value.toFloat();
                        size_t len = 1;

                        call_netcdf_function( [&](){
                            return nc_put_att(hndl, id, c_attname.constData(),
                                              xtyp, len, &val); } );
                        break;
                    }
                    default:
                        throw SireError::unsupported( QObject::tr(
                                "Not yet able to write a NetCDF attribute of type '%1' "
                                "for variable %2. (%3)")
                                    .arg(att_type).arg(variable).arg(xtyp), CODELOC );
                }
            }
        }

        //we have finished writing the metadata about the file
        call_netcdf_function( [&](){ return nc_enddef(hndl); } );

        //now that the metadata has been written, we can now write the actual data
        for (const auto &variable : variables)
        {
            const auto vardata = variable_data[variable];

            const int id = var_ids.value(variable, -1);

            call_netcdf_function( [&](){
                return nc_put_var( hndl, id, vardata.memdata.constData() ); } );
        }

        //finished writing the file :-)
        call_netcdf_function( [&](){ return nc_close(hndl); } );
        hndl = -1;

    #endif
    }
}

/** Return the names and sizes of all of the dimensions in the file */
QHash<QString,int> NetCDFFile::getDimensions() const
{
    QHash<QString,int> dims;

    #ifdef SIRE_USE_NETCDF
    if (hndl != -1)
    {
        int ndims;
        call_netcdf_function( [&](){ return nc_inq_ndims(hndl, &ndims); } );

        if (ndims <= 0)
            return dims;

        char *dim_name = new char[NC_MAX_NAME+1];

        try
        {
            for (int i=0; i<ndims; ++i)
            {
                size_t dim_len;
                call_netcdf_function( [&](){ return nc_inq_dim(hndl, i, dim_name, &dim_len); } );
                dims.insert( QString::fromUtf8(dim_name), dim_len );
            }
        }
        catch(...)
        {
            delete[] dim_name;
            throw;
        }

        delete[] dim_name;

    }
    #endif

    return dims;
}

/** Write the passed NetCDFData to the file 'filename'. This will overwrite
    the file if 'overwrite_file' is true, will write the file using a 64bit offset
    if 'use_64bit_offset' is true, and will write using NetCDF 4 if 'use_netcdf4'
    is true (otherwise it will write using NetCDF 3) */
QString NetCDFFile::write(const QString &filename,
                          const QHash<QString,QString> &globals,
                          const QHash<QString,NetCDFData> &data,
                          bool overwrite_file, bool use_64bit_offset,
                          bool use_netcdf4)
{
    #ifdef SIRE_USE_NETCDF
        QFileInfo file(filename);

        if (file.exists())
        {
            if (not overwrite_file)
            {
                throw SireError::io_error( QObject::tr(
                        "Cannot create the NetCDF file '%1' as it already exists, and "
                        "the software is not allowed to overwrite an existing file!")
                            .arg(filename), CODELOC );
            }

            if (file.isDir())
            {
                throw SireError::io_error( QObject::tr(
                        "Cannot create the NetCDF file '%1' as it already exists and "
                        "is a directory!")
                            .arg(filename), CODELOC );
            }

            if (not file.isWritable())
            {
                throw SireError::io_error( QObject::tr(
                        "Cannot create the NetCDF file '%1' as it already exists and "
                        "is read-only")
                            .arg(filename), CODELOC );
            }
        }

        QString absfile = file.absoluteFilePath();

        NetCDFFile netcdf(absfile, overwrite_file, use_64bit_offset, use_netcdf4);
        netcdf.writeData(globals, data);

        return absfile;
    #else
        throw SireError::unsupported( QObject::tr(
            "This version of Sire does not have support for writing NetCDF files!"),
                CODELOC );
    #endif
}

/** Return the value of the string attribute 'name'.*/
QString NetCDFFile::getStringAttribute(const QString &name) const
{
    if (hndl != -1)
    {
    #ifdef SIRE_USE_NETCDF

        QByteArray c_name = name.toUtf8();

        //get the size of the attribute
        size_t vsize;

        call_netcdf_function( [&]()
            { return nc_inq_attlen(hndl, NC_GLOBAL, c_name.constData(), &vsize); } );

        //get the attribute
        char *c_value = new char[vsize+1];

        try
        {
            call_netcdf_function( [&]()
                { return nc_get_att_text(hndl, NC_GLOBAL, c_name.constData(), c_value); } );

            c_value[vsize] = '\0';

            QString value = QString::fromUtf8(c_value);

            delete[] c_value;

            return value;
        }
        catch(...)
        {
            delete[] c_value;
            throw;
        }
    #endif
    }

    throw SireError::invalid_key( QObject::tr(
            "There is not string attribute called '%1' in the NetCDF file '%2'")
                .arg(name).arg(fname), CODELOC );

    return QString();
}

/** Read in and return the NetCDFData associated with the passed NetCDFDataInfo */
NetCDFData NetCDFFile::read(const NetCDFDataInfo &variable) const
{
    NetCDFData data(variable);

    int data_size = data.dataSize();

    if (hndl != -1 and data_size > 0)
    {
    #ifdef SIRE_USE_NETCDF
        QByteArray memdata;
        memdata.fill('\0', data_size);
        call_netcdf_function( [&](){ return nc_get_var(hndl, variable.ID(), memdata.data()); } );
        data.setData(memdata);
    #endif
    }

    return data;
}
