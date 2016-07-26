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

#include "magic_error.h"

#include "datastream.h"

#include "version_error.h"

#include "SireError/errors.h"

#include <QStringList>
#include <QDataStream>

using namespace SireStream;
using namespace SireError;

// Define the streaming operators for SireError here
// This is to remove the circular dependency of SireError on SireStream

/** Implementation of SireError::exception */
static const RegisterMetaType<exception> r_exception( MAGIC_ONLY,
                                                      "SireError::exception" );

/** Serialise to a binary data stream */
QDataStream SIRESTREAM_EXPORT &operator<<(QDataStream &ds, const exception &e)
{
    writeHeader(ds, r_exception, 1)
       << e.err << e.plce << e.bt << e.pidstr;

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIRESTREAM_EXPORT &operator>>(QDataStream &ds, exception &e)
{
    VersionID v = readHeader(ds, r_exception);

    if (v == 1)
    {
        ds >> e.err >> e.plce >> e.bt >> e.pidstr;
    }
    else
        throw SireStream::version_error(v, "1", r_exception, CODELOC);

    return ds;
}

const char* magic_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<magic_error>() );
}

static const RegisterMetaType<SireError::program_bug> r_program_bug;
static const RegisterMetaType<SireError::unsupported> r_unsupported;
static const RegisterMetaType<SireError::id_error> r_id_error;
static const RegisterMetaType<SireError::invalid_key> r_invalid_key;
static const RegisterMetaType<SireError::invalid_index> r_invalid_index;
static const RegisterMetaType<SireError::invalid_cast> r_invalid_cast;
static const RegisterMetaType<SireError::noncopyable_error> r_noncopyable_error;
static const RegisterMetaType<SireError::nullptr_error> r_nullptr_error;
static const RegisterMetaType<SireError::lock_error> r_lock_error;
static const RegisterMetaType<SireError::assertation_failed> r_assertation_failed;
static const RegisterMetaType<SireError::incompatible_error> r_incompatible_error;
static const RegisterMetaType<SireError::file_error> r_file_error;
static const RegisterMetaType<SireError::io_error> r_io_error;
static const RegisterMetaType<SireError::process_error> r_process_error;
static const RegisterMetaType<SireError::invalid_arg> r_invalid_arg;
static const RegisterMetaType<SireError::invalid_state> r_invalid_state;
static const RegisterMetaType<SireError::invalid_operation> r_invalid_operation;
static const RegisterMetaType<SireError::unavailable_resource> r_unavailable_resource;
static const RegisterMetaType<SireError::incomplete_code> r_incomplete_code;
static const RegisterMetaType<SireError::std_exception> r_std_exception;
static const RegisterMetaType<SireError::unknown_exception> r_unknown_exception;
static const RegisterMetaType<SireError::unknown_type> r_unknown_type;
static const RegisterMetaType<SireError::dependency_error> r_dependency_error;
static const RegisterMetaType<SireError::version_error> r_version_error;

static RegisterMetaType<magic_error> r_magic_error;
