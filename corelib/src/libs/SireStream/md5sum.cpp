/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "md5sum.h"

#include <QString>
#include <QRegExp>
#include <QFile>
#include <QByteArray>
#include <QDataStream>

#include "ThirdParty/md5.h"      // CONDITIONAL_INCLUDE

using namespace SireStream;

QDataStream SIRESTREAM_EXPORT &operator<<(QDataStream &ds, const MD5Sum &md5sum)
{
    for (int i=0; i<16; ++i)
    {
        ds << quint8(md5sum.dgst[i]);
    }
    
    return ds;
}

QDataStream SIRESTREAM_EXPORT &operator>>(QDataStream &ds, MD5Sum &md5sum)
{
    for (int i=0; i<16; ++i)
    {
        quint8 c;
        ds >> c;
        
        md5sum.dgst[i] = c;
    }
    
    return ds;
}

MD5Sum::MD5Sum()
{
    //generate the null digest
    generate(0,0);
}

MD5Sum::MD5Sum(const QByteArray &buffer)
{
    generate(buffer.data(),buffer.size());
}

MD5Sum::MD5Sum(const QString &file)
{
    //get the string of bytes in 'file'
    QFile f(file);
    if (!f.open(QIODevice::ReadOnly))
    {
        //generate the null digest
        generate(0,0);
    }

    //get the string of bytes
    QByteArray bytes = f.readAll();
    f.close();

    generate(bytes.data(),bytes.size());
}

MD5Sum::MD5Sum(const char *buffer, unsigned int sz)
{
    generate(buffer,sz);
}

MD5Sum::MD5Sum(const MD5Sum &other)
{
    for (int i=0; i<16; i++)
        dgst[i] = other.dgst[i];
}

MD5Sum::~ MD5Sum()
{}

void MD5Sum::generate(const char* buffer, unsigned int sz)
{
    //use L. Peter Deutsch's free implementation of
    //the md5 algorithm
    md5_state_t state;
    //initialise the md5 engine
    md5_init(&state);
    //make it decode the bits of the image
    md5_append(&state,(unsigned char*)buffer,sz);

    //actually calculate the digest
    md5_finish(&state,dgst);
}

const MD5Sum& MD5Sum::operator=(const MD5Sum &other)
{
    if (*this == other)
        return *this;
    else
    {
        for (int i=0; i<16; i++)
            dgst[i] = other.dgst[i];

        return *this;
    }
}

QString MD5Sum::toString() const
{
    return QString().sprintf("%02x%02x%02x%02x-%02x%02x%02x%02x-%02x%02x%02x%02x-%02x%02x%02x%02x",
                             dgst[0],dgst[1],dgst[2],dgst[3],
                             dgst[4],dgst[5],dgst[6],dgst[7],
                             dgst[8],dgst[9],dgst[10],dgst[11],
                             dgst[12],dgst[13],dgst[14],dgst[15]);
}

const md5_byte_t* MD5Sum::digest() const
{
    return dgst;
} 
 
bool MD5Sum::operator==(const MD5Sum &other) const
{
    return (dgst[0] == other.dgst[0] and
            dgst[1] == other.dgst[1] and
            dgst[2] == other.dgst[2] and
            dgst[3] == other.dgst[3] and
            dgst[4] == other.dgst[4] and
            dgst[5] == other.dgst[5] and
            dgst[6] == other.dgst[6] and
            dgst[7] == other.dgst[7] and
            dgst[8] == other.dgst[8] and
            dgst[9] == other.dgst[9] and
            dgst[10] == other.dgst[10] and
            dgst[11] == other.dgst[11] and
            dgst[12] == other.dgst[12] and
            dgst[13] == other.dgst[13] and
            dgst[14] == other.dgst[14] and
            dgst[15] == other.dgst[15]);
}

bool MD5Sum::operator!=(const MD5Sum &other) const
{
    return (dgst[0] != other.dgst[0] or
            dgst[1] != other.dgst[1] or
            dgst[2] != other.dgst[2] or
            dgst[3] != other.dgst[3] or
            dgst[4] != other.dgst[4] or
            dgst[5] != other.dgst[5] or
            dgst[6] != other.dgst[6] or
            dgst[7] != other.dgst[7] or
            dgst[8] != other.dgst[8] or
            dgst[9] != other.dgst[9] or
            dgst[10] != other.dgst[10] or
            dgst[11] != other.dgst[11] or
            dgst[12] != other.dgst[12] or
            dgst[13] != other.dgst[13] or
            dgst[14] != other.dgst[14] or
            dgst[15] != other.dgst[15]);
}
