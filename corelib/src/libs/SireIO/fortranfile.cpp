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

#include "fortranfile.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include <QFile>
#include <QFileInfo>
#include <QDataStream>
#include <QtEndian>
#include <QDebug>

using namespace SireIO;

//////
////// Implementatin of FortranFile
//////

FortranFile::FortranFile()
            : int_size(4), is_little_endian(true)
{}

bool FortranFile::try_read()
{
    record_pointers.clear();
    record_sizes.clear();

    QFile file(abs_filename);

    if (not file.open(QIODevice::ReadOnly))
    {
        throw SireError::io_error(QObject::tr(
            "Could not open file %1. Please check it exists and is readable.")
                .arg(abs_filename), CODELOC);
    }

    QDataStream ds(&file);

    QByteArray start_buffer(int_size, 0);
    QByteArray end_buffer(int_size, 0);

    qint64 read_count = 0;

    // each fortran record starts and ends with an integer that
    // gives the size in bytes of the record. We will now scan
    // through the file, and if we can consistently get all of the
    // sizes, then we much have the right int_size and endianness
    while (not ds.atEnd())
    {
        int read_size = ds.readRawData(start_buffer.data(), int_size);

        if (read_size != int_size)
        {
            return false;
        }

        qint64 start_size;

        if (int_size == 4)
        {
            if (is_little_endian)
                start_size = qFromLittleEndian<qint32>(start_buffer.data());
            else
                start_size = qFromBigEndian<qint32>(start_buffer.data());
        }
        else
        {
            if (is_little_endian)
                start_size = qFromLittleEndian<qint64>(start_buffer.data());
            else
                start_size = qFromBigEndian<qint64>(start_buffer.data());
        }

        read_count += int_size;

        read_size = ds.skipRawData(start_size);

        if (read_size != start_size)
        {
            //could not read this much!
            return false;
        }

        read_size = ds.readRawData(end_buffer.data(), int_size);

        if (read_size != int_size)
        {
            return false;
        }

        qint64 end_size;

        if (int_size == 4)
        {
            if (is_little_endian)
                end_size = qFromLittleEndian<qint32>(end_buffer.data());
            else
                end_size = qFromBigEndian<qint32>(end_buffer.data());
        }
        else
        {
            if (is_little_endian)
                end_size = qFromLittleEndian<qint64>(end_buffer.data());
            else
                end_size = qFromBigEndian<qint64>(end_buffer.data());
        }

        if (start_size != end_size)
        {
            //disagreement - cannot be a valid record
            return false;
        }

        record_pointers.append(read_count);
        record_sizes.append(start_size);

        read_count += start_size + int_size;
    }

    return true;
}

FortranFile::FortranFile(const QString &filename)
            : int_size(4), is_little_endian(true)
{
    abs_filename = QFileInfo(filename).absoluteFilePath();

    // try to read using 4 byte header and native endian
    int_size = 4;
    is_little_endian = true;

    if (try_read())
        return;

    // try to read using 8 byte header and native endian
    int_size = 8;
    is_little_endian = true;

    if (try_read())
        return;

    // try to read using 4 byte header and swapped endian
    int_size = 4;
    is_little_endian = false;

    if (try_read())
        return;

    // try to read using 8 byte header and swapped endian
    int_size = 8;
    is_little_endian = false;

    if (try_read())
        return;

    throw SireError::io_error(QObject::tr(
        "Could not read a consistent set of records from %1. "
        "It could not be read as a record-based unformatted "
        "Fortran binary file.").arg(abs_filename), CODELOC);
}

FortranFile::FortranFile(const FortranFile &other)
            : abs_filename(other.abs_filename),
              record_pointers(other.record_pointers),
              record_sizes(other.record_sizes),
              int_size(other.int_size),
              is_little_endian(other.is_little_endian)
{}

FortranFile::~FortranFile()
{}

FortranFile& FortranFile::operator=(const FortranFile &other)
{
    if (this != &other)
    {
        abs_filename = other.abs_filename;
        record_pointers = other.record_pointers;
        record_sizes = other.record_sizes;
        int_size = other.int_size;
        is_little_endian = other.is_little_endian;
    }

    return *this;
}

int FortranFile::nRecords() const
{
    return record_pointers.count();
}

FortranRecord FortranFile::operator[](int i) const
{
    i = SireID::Index(i).map(this->nRecords());

    QFile file(this->abs_filename);

    if (!file.open(QIODevice::ReadOnly))
        throw SireError::io_error( QObject::tr(
            "Problem opening file '%1'. Please make sure it is readable."
                    ).arg(abs_filename), CODELOC);

    QDataStream ds(&file);

    qint64 pointer = this->record_pointers[i];

    int skipped = ds.skipRawData(pointer);

    if (pointer != skipped)
    {
        throw SireError::io_error(QObject::tr(
            "Problem reading record %1. Needed to skip %2 bytes, but could "
            "only skip %3. Is the file corrupted?")
                .arg(abs_filename).arg(pointer).arg(skipped), CODELOC);
    }

    qint64 size = this->record_sizes[i];

    QByteArray array(size, 0);

    int read = ds.readRawData(array.data(), size);

    if (size != read)
    {
        throw SireError::io_error(QObject::tr(
            "Problem reading record %1. Needed to read %3 bytes, but could "
            "only read %3. Is the file corrupted?")
                .arg(abs_filename).arg(size).arg(read), CODELOC);
    }

    return FortranRecord(array, is_little_endian);
}

//////
////// Implementatin of FortranRecord
//////

FortranRecord::FortranRecord() : cursor(0)
{
    is_little_endian = true;
}

FortranRecord::FortranRecord(const QByteArray &d, bool le)
              : data(d), cursor(0), is_little_endian(le)
{}

FortranRecord::FortranRecord(const FortranRecord &other)
              : data(other.data),
                cursor(other.cursor),
                is_little_endian(other.is_little_endian)
{}

FortranRecord::~FortranRecord()
{}

FortranRecord& FortranRecord::operator=(const FortranRecord &other)
{
    if (this != &other)
    {
        data = other.data;
        cursor = other.cursor;
        is_little_endian = other.is_little_endian;
    }

    return *this;
}

int FortranRecord::size() const
{
    return data.count();
}

void FortranRecord::_assertValid(int n, int size) const
{
    if (cursor + (n*size) > data.count())
    {
        throw SireError::io_error(QObject::tr(
            "Cannot read %1 x %2 bytes of data as only %3 bytes remain!")
                .arg(n).arg(size).arg(data.count()-cursor), CODELOC);
    }
}

void FortranRecord::_assertPosValid(int pos, int size) const
{
    if (pos < 0)
    {
        throw SireError::io_error(QObject::tr(
            "Cannot read %1 bytes at position %2 as position is negative!")
                .arg(size).arg(pos), CODELOC);

    }
    else if (pos + size > data.count())
    {
        throw SireError::io_error(QObject::tr(
            "Cannot read %1 bytes at position %2 as only %3 bytes remain!")
                .arg(size).arg(pos).arg(data.count()-pos), CODELOC);
    }
}

char FortranRecord::readCharAt(int pos) const
{
    _assertPosValid(pos, 1);
    return data.constData()[pos];
}

double FortranRecord::readFloat64At(int pos) const
{
    _assertPosValid(pos, 8);

    double ret;

    #if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        if (is_little_endian)
        {
            // just copy the data
            memcpy(&ret, data.constData()+pos, 8);
        }
        else
        {
            //need to reverse the data
            throw SireError::incomplete_code(QObject::tr(
                "Haven't written code to deal with big endian doubles..."),
                    CODELOC);
        }
    #else
        if (is_little_endian)
        {
            //need to reverse the data
            throw SireError::incomplete_code(QObject::tr(
                "Haven't written code to deal with little endian doubles..."),
                    CODELOC);
        }
        else
        {
            // just copy the data
            memcpy(&ret, data.constData()+pos, 8);
        }
    #endif

    return ret;
}

float FortranRecord::readFloat32At(int pos) const
{
    _assertPosValid(pos, 4);

    float ret;

    #if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        if (is_little_endian)
        {
            // just copy the data
            memcpy(&ret, data.constData()+pos, 4);
        }
        else
        {
            //need to reverse the data
            throw SireError::incomplete_code(QObject::tr(
                "Haven't written code to deal with big endian doubles..."),
                    CODELOC);
        }
    #else
        if (is_little_endian)
        {
            //need to reverse the data
            throw SireError::incomplete_code(QObject::tr(
                "Haven't written code to deal with little endian doubles..."),
                    CODELOC);
        }
        else
        {
            // just copy the data
            memcpy(&ret, data.constData()+pos, 4);
        }
    #endif

    return ret;
}

qint32 FortranRecord::readInt32At(int pos) const
{
    _assertPosValid(pos, 4);

    qint32 ret;

    #if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        if (is_little_endian)
        {
            // just copy the data
            memcpy(&ret, data.constData()+pos, 4);
        }
        else
        {
            ret = qFromBigEndian<qint32>(data.constData()+pos);
        }
    #else
        if (is_little_endian)
        {
            ret = qFromLittleEndian<qint32>(data.constData()+pos);
        }
        else
        {
            // just copy the data
            memcpy(&ret, data.constData()+pos, 4);
        }
    #endif

    return ret;
}

qint64 FortranRecord::readInt64At(int pos) const
{
    _assertPosValid(pos, 8);

    qint64 ret;

    #if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        if (is_little_endian)
        {
            // just copy the data
            memcpy(&ret, data.constData()+pos, 8);
        }
        else
        {
            ret = qFromBigEndian<qint64>(data.constData()+pos);
        }
    #else
        if (is_little_endian)
        {
            ret = qFromLittleEndian<qint64>(data.constData()+pos);
        }
        else
        {
            // just copy the data
            memcpy(&ret, data.constData()+pos, 8);
        }
    #endif

    return ret;
}

QString FortranRecord::readChar(int n)
{
    if (n <= 0)
        return QString();

    _assertValid(n, 1);

    auto ret = QString::fromUtf8(data.constData()+cursor, n);

    cursor += n;

    return ret;
}

QVector<double> FortranRecord::readFloat64(int n)
{
    if (n <= 0)
        return QVector<double>();

    if (sizeof(double) != 8)
        throw SireError::incomplete_code(QObject::tr(
            "Haven't written code to deal with non-64bit doubles..."),
                CODELOC);

    _assertValid(n, 8);

    QVector<double> ret(n);

    #if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        if (is_little_endian)
        {
            // just copy the data
            memcpy(ret.data(), data.constData()+cursor, n*8);
            cursor += n*8;
        }
        else
        {
            //need to reverse the data
            throw SireError::incomplete_code(QObject::tr(
                "Haven't written code to deal with big endian doubles..."),
                    CODELOC);
        }
    #else
        if (is_little_endian)
        {
            //need to reverse the data
            throw SireError::incomplete_code(QObject::tr(
                "Haven't written code to deal with little endian doubles..."),
                    CODELOC);
        }
        else
        {
            // just copy the data
            memcpy(ret.data(), data.constData()+cursor, n*8);
            cursor += n*8;
        }
    #endif

    return ret;
}

QVector<float> FortranRecord::readFloat32(int n)
{
    if (n <= 0)
        return QVector<float>();

    if (sizeof(float) != 4)
        throw SireError::incomplete_code(QObject::tr(
            "Haven't written code to deal with non-32bit floats..."),
                CODELOC);

    _assertValid(n, 4);

    QVector<float> ret(n);

    #if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        if (is_little_endian)
        {
            // just copy the data
            memcpy(ret.data(), data.constData()+cursor, n*4);
            cursor += n*4;
        }
        else
        {
            //need to reverse the data
            throw SireError::incomplete_code(QObject::tr(
                "Haven't written code to deal with big endian floats..."),
                    CODELOC);
        }
    #else
        if (is_little_endian)
        {
            //need to reverse the data
            throw SireError::incomplete_code(QObject::tr(
                "Haven't written code to deal with little endian floats..."),
                    CODELOC);
        }
        else
        {
            // just copy the data
            memcpy(ret.data(), data.constData()+cursor, n*4);
            cursor += n*4;
        }
    #endif

    return ret;
}

QVector<qint32> FortranRecord::readInt32(int n)
{
    if (n <= 0)
        return QVector<qint32>();

    _assertValid(n, 4);

    QVector<qint32> ret(n);

    #if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        if (is_little_endian)
        {
            // just copy the data
            memcpy(ret.data(), data.constData()+cursor, n*4);
            cursor += n*4;
        }
        else
        {
            //need to reverse the data
            for (int i=0; i<n; ++i)
            {
                ret[i] = qFromBigEndian<qint32>(data.constData()+cursor);
                cursor += 4;
            }
        }
    #else
        if (is_little_endian)
        {
            //need to reverse the data
            for (int i=0; i<n; ++i)
            {
                ret[i] = qFromLittleEndian<qint32>(data.constData()+cursor);
                cursor += 4;
            }
        }
        else
        {
            // just copy the data
            memcpy(ret.data(), data.constData()+cursor, n*4);
            cursor += n*4;
        }
    #endif

    return ret;
}

QVector<qint64> FortranRecord::readInt64(int n)
{
    if (n <= 0)
        return QVector<qint64>();

    _assertValid(n, 8);

    QVector<qint64> ret(n);

    #if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        if (is_little_endian)
        {
            // just copy the data
            memcpy(ret.data(), data.constData()+cursor, n*8);
            cursor += n*8;
        }
        else
        {
            //need to reverse the data
            for (int i=0; i<n; ++i)
            {
                ret[i] = qFromBigEndian<qint64>(data.constData()+cursor);
                cursor += 8;
            }
        }
    #else
        if (is_little_endian)
        {
            //need to reverse the data
            for (int i=0; i<n; ++i)
            {
                ret[i] = qFromLittleEndian<qint64>(data.constData()+cursor);
                cursor += 8;
            }
        }
        else
        {
            // just copy the data
            memcpy(ret.data(), data.constData()+cursor, n*8);
            cursor += n*8;
        }
    #endif

    return ret;
}
