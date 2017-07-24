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

#ifndef SIREIO_NETCDFFILE_H
#define SIREIO_NETCDFFILE_H

#include "sireglobal.h"

#include <QString>
#include <QStringList>
#include <QVariant>
#include <QMutex>

#include <boost/noncopyable.hpp>

SIRE_BEGIN_HEADER

namespace SireIO
{

/** This class provides information about a data variable in 
    a NetCDF file
    
    @author Christopher Woods
*/
class SIREIO_EXPORT NetCDFDataInfo
{

friend class NetCDFFile;

public:
    NetCDFDataInfo();
    NetCDFDataInfo(const NetCDFDataInfo &other);
    
    ~NetCDFDataInfo();
    
    int ID() const;
    QString name() const;
    QString type() const;
    
    int typeSize() const;
    int dataSize() const;
    
    QStringList dimensions() const;
    QList<int> dimensionSizes() const;
    
    int nValues() const;
    
    int nAttributes() const;
    
    QStringList attributeNames() const;

    QVariant attribute(const QString &name) const;
    QString attributeType(const QString &name) const;

    QHash<QString,QVariant> attributes() const;
    QHash<QString,QString> attributeTypes() const;
    
    QString toString() const;

    bool isNull() const;
    
protected:
    NetCDFDataInfo(int idnum, QString name, int xtyp,
                   QStringList dim_names, QList<int> dim_sizes,
                   QStringList att_names, QList<int> att_types,
                   QList<QVariant> att_values);

    /** The name of the variable */
    QString nme;
    
    /** The names of each of the dimensions of the variable */
    QStringList dim_names;
    
    /** The size of each of the dimensions */
    QList<int> dim_sizes;
    
    /** The names of all attributes associated with the variable */
    QStringList att_names;
    /** The types of all of the attributes */
    QList<int> att_types;
    /** The values of all of the attributes */
    QList<QVariant> att_values;
    
    /** The ID number of the variable in the data file */
    int idnum;

    /** The type of the data in the data file */
    int xtyp;
};

/** This class holds the actual data read from a NetCDF file

    @author Christopher Woods
*/
class SIREIO_EXPORT NetCDFData : public NetCDFDataInfo
{

friend class NetCDFFile;

public:
    NetCDFData();
    NetCDFData(const NetCDFData &other);
    
    ~NetCDFData();

    QVector<QVariant> toArray() const;

    QVector<float> toFloatArray() const;
    QVector<double> toDoubleArray() const;
    
    QVector<qint32> toInt32Array() const;
    QVector<qint64> toInt64Array() const;

protected:
    NetCDFData(const NetCDFDataInfo &info);
    
    void setData(const QByteArray &data);
    
private:
    /** Raw memory containing the data */
    QByteArray memdata;
};

/** This class provides an internal interface to NetCDF files 

    @author Christopher Woods
*/
class SIREIO_EXPORT NetCDFFile : public boost::noncopyable
{
public:
    NetCDFFile();

    NetCDFFile(const QString &filename, bool overwrite_file,
               bool use_64bit_offset=true, bool use_netcdf4=false);

    NetCDFFile(const QString &filename);
    
    ~NetCDFFile();
    
    QString getStringAttribute(const QString &name) const;
    
    QHash<QString,int> getDimensions() const;
    
    QHash<QString,NetCDFDataInfo> getVariablesInfo() const;
    
    NetCDFData read(const NetCDFDataInfo &variable) const;

    void writeData(const QHash<QString,NetCDFData> &data);

private:
    int call_netcdf_function( std::function<int()> func,
                              int ignored_error = 0) const;

    /** The name of the file */
    QString fname;

    /** Handle to the NetCDF file */
    int hndl;
    
    /** Mutex to serialise all file IO operations */
    QMutex mutex;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the ID number of this piece of data */
inline int NetCDFDataInfo::ID() const
{
    return idnum;
}

/** Return the name of this piece of data */
inline QString NetCDFDataInfo::name() const
{
    return nme;
}

/** Return the names of the dimensions of this data */
inline QStringList NetCDFDataInfo::dimensions() const
{
    return dim_names;
}

/** Return the number of values for each of the dimensions of this data */
inline QList<int> NetCDFDataInfo::dimensionSizes() const
{
    return dim_sizes;
}

/** Return the number of attributes of this data in the file */
inline int NetCDFDataInfo::nAttributes() const
{
    return att_names.count();
}

#endif

}

SIRE_END_HEADER

#endif
