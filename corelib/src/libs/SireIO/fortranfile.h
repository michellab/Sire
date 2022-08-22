/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREIO_FORTRANFILE_H
#define SIREIO_FORTRANFILE_H

#include "sireglobal.h"

#include <QByteArray>

SIRE_BEGIN_HEADER

namespace SireIO
{

/** This represents a single Fortran record, from which you
 *  can extract data
 */
class FortranRecord
{
public:
    FortranRecord();
    FortranRecord(const QByteArray &data, bool is_little_endian);
    FortranRecord(const FortranRecord &other);

    ~FortranRecord();

    FortranRecord& operator=(const FortranRecord &other);

    QString readChar(int n) const;

    QVector<double> readDouble(int n) const;
    QVector<float> readFloat(int n) const;
    QVector<qint32> readInt32(int n) const;
    QVector<qint64> readInt64(int n) const;

private:
    QByteArray data;
    bool is_little_endian;
};

/** This class is used to read and write fortran binary
 *  unformatted files (written as record files, not
 *  streaming files).
 *
 *  This automatically detects the endianness of the file
 */
class FortranFile
{
public:
    FortranFile();
    FortranFile(const QString &filename);
    FortranFile(const FortranFile &other);
    ~FortranFile();

    FortranFile& operator=(const FortranFile &other);

    int nRecords() const;

    FortranRecord operator[](int i) const;

private:
    bool try_read();

    QString abs_filename;

    QVector<qint64> record_pointers;
    QVector<qint64> record_sizes;

    int int_size;

    bool is_little_endian;
};

}

SIRE_END_HEADER

#endif
