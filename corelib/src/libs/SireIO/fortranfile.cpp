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
    qDebug() << "\nTRY_READ" << int_size << is_little_endian;

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
        qDebug() << "NEW RECORD!";
        int read_size = ds.readRawData(start_buffer.data(), int_size);

        if (read_size != int_size)
        {
            qDebug() << CODELOC << read_size << int_size;
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
            qDebug() << CODELOC << read_size << start_size;
            return false;
        }

        read_size = ds.readRawData(end_buffer.data(), int_size);

        if (read_size != int_size)
        {
            qDebug() << CODELOC << read_size << int_size;
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
            qDebug() << CODELOC << start_size << end_size;
            return false;
        }

        record_pointers.append(read_count);
        record_sizes.append(start_size);

        read_count += int_size;
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

FortranRecord::FortranRecord()
{
    is_little_endian = true;
}

FortranRecord::FortranRecord(const QByteArray &d, bool le)
              : data(d), is_little_endian(le)
{}

FortranRecord::FortranRecord(const FortranRecord &other)
              : data(other.data), is_little_endian(other.is_little_endian)
{}

FortranRecord::~FortranRecord()
{}

FortranRecord& FortranRecord::operator=(const FortranRecord &other)
{
    if (this != &other)
    {
        data = other.data;
        is_little_endian = other.is_little_endian;
    }

    return *this;
}

QString FortranRecord::readChar(int n) const
{
    return QString();
}

QVector<double> FortranRecord::readDouble(int n) const
{
    return QVector<double>();
}

QVector<float> FortranRecord::readFloat(int n) const
{
    return QVector<float>();
}

QVector<qint32> FortranRecord::readInt32(int n) const
{
    return QVector<qint32>();
}

QVector<qint64> FortranRecord::readInt64(int n) const
{
    return QVector<qint64>();
}
