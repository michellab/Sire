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

#include "SireError/errors.h"

#ifdef SIRE_USE_NETCDF
    #include "netcdf.h"  // CONDITIONAL_INCLUDE
#endif

using namespace SireIO;

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
    #endif
}

/** Constructor */
NetCDFFile::NetCDFFile() : hndl(-1)
{}

/** Construct to open the file 'filename' */
NetCDFFile::NetCDFFile(const QString &filename) : fname(filename)
{
    #ifdef SIRE_USE_NETCDF
        QFileInfo file(filename);
    
        if (not file.exists())
        {
            throw SireError::io_error( QObject::tr(
                    "Cannot open '%1' as it does not appear to exist.")
                        .arg(filename), CODELOC );
        }
    
        int errnum;
        QByteArray c_filename = file.absoluteFilePath().toUtf8();
        errnum = nc_open(c_filename.constData(), NC_NOWRITE, &hndl);
        assert_no_netcdf_error(errnum);

    #else
        throw SireError::unsupported( QObject::tr(
                "Software is missing NetCDF support, so cannot read Amber Rst files!"),
                    CODELOC );
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

/** Return the value of the string attribute 'name'.*/
QString NetCDFFile::getStringAttribute(const QString &name) const
{
    if (hndl != -1)
    {
    #ifdef SIRE_USE_NETCDF
    
        QByteArray c_name = name.toUtf8();
    
        //get the size of the attribute
        size_t vsize;
        
        int err = nc_inq_attlen(hndl, NC_GLOBAL, c_name.constData(), &vsize);
        assert_no_netcdf_error(err);
        
        //get the attribute
        char *c_value = new char[vsize+1];
        
        err = nc_get_att_text(hndl, NC_GLOBAL, c_name.constData(), c_value);
        assert_no_netcdf_error(err);
        c_value[vsize] = '\0';
        
        QString value = QString::fromUtf8(c_value);
        
        delete[] c_value;
        
        return value;
    #endif
    }
    
    throw SireError::invalid_key( QObject::tr(
            "There is not string attribute called '%1' in the NetCDF file '%2'")
                .arg(name).arg(fname), CODELOC );
    
    return QString();
}
