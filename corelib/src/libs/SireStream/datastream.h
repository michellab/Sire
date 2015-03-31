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

#ifndef SIRESTREAM_DATASTREAM_H
#define SIRESTREAM_DATASTREAM_H

// This file contains functions and definitions that
// will (hopefully!) make writing binary streaming functions
// more easy and more robust.
//
// @author Christopher Woods

#include "SireStream/version_error.h"

class QDataStream;

SIRE_BEGIN_HEADER

namespace SireStream
{

QDataStream& writeHeader(QDataStream &ds, const RegisterMetaTypeBase &r_type,
                         VersionID version);

QDataStream& writeHeader(QDataStream &ds, MagicID magicid, VersionID version);

VersionID readHeader(QDataStream &ds, const RegisterMetaTypeBase &r_type);

VersionID readHeader(QDataStream &ds, MagicID magicid, const char *type_name);

}

SIRE_END_HEADER

#endif
