/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIRESTREAM_MD5SUM_H
#define SIRESTREAM_MD5SUM_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

class QFile;
class QString;
class QByteArray;

namespace SireStream
{
class MD5Sum;
}

SIRESTREAM_EXPORT QDataStream& operator<<(QDataStream&, const SireStream::MD5Sum&);
SIRESTREAM_EXPORT QDataStream& operator>>(QDataStream&, SireStream::MD5Sum&);

typedef unsigned char md5_byte_t;

namespace SireStream
{

/**
 * These functions are used to simplify the generation of md5sums of 
 * buffers and files
 *
 * @author Christopher Woods
 */

class SIRESTREAM_EXPORT MD5Sum
{

friend SIRESTREAM_EXPORT QDataStream& ::operator<<(QDataStream&, const MD5Sum&);
friend SIRESTREAM_EXPORT QDataStream& ::operator>>(QDataStream&, MD5Sum&);

public:
  /** Construct a null MD5Sum */
  MD5Sum();
  /** Construct an MD5Sum from a buffer */
  MD5Sum(const QByteArray &buffer);
  /**
   * Construct an MD5Sum from a const char* buffer.
   * @param buffer The buffer containing the data.
   * @param sz This size of the buffer. Don't lie otherwise it will crash!
   * @return A constructed, valid MD5Sum object. 
   */
  MD5Sum(const char *buffer, unsigned int sz);
  /** Construct the MD5Sum for a file */
  MD5Sum(const QString &file);
  /** Copy constructor */
  MD5Sum(const MD5Sum &sum);
  
  ~MD5Sum();
  
  /** Assignment operator */
  const MD5Sum& operator=(const MD5Sum &other);
  
  /** Equality operators */
  bool operator==(const MD5Sum &other) const;
  bool operator!=(const MD5Sum &other) const;
  
  /** Return a pointer to the digest array (16 element array) */
  const md5_byte_t* digest() const;
  
  /** Return a string representation of the digest */
  QString toString() const;

private:
  /**
   * Generate the md5sum digest.
   * @param buffer The buffer
   * @param sz The size of the buffer
   */
  void generate(const char *buffer, unsigned int sz);
  
  /** Pointer to the storage of the md5 digest */
  md5_byte_t dgst[16];
};

}

SIRE_EXPOSE_CLASS( SireStream::MD5Sum )

SIRE_END_HEADER

#endif
