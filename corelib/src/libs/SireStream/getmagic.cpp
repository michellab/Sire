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

#include "sireglobal.h"

#include "ThirdParty/md5.h"

#include <QByteArray>

using namespace Sire;

/** This function is used to obtain a reasonably unique
    MagicID number for each class type (with typename 'name') */
MagicID SIRE_EXPORT Sire::getMagic(const char *name)
{
    //use L. Peter Deutsch's free implementation of
    //the md5 algorithm
    md5_state_t state;
    //initialise the md5 engine
    md5_init(&state);
    //make it encode the name...
    md5_append(&state, (unsigned char*)name, qstrlen(name));
    //this variable will hold the returned digest
    md5_byte_t dgst[16];
    //actually calculate the digest
    md5_finish(&state,dgst);

    //only use the first 4 bytes of the digest to calculate the
    //MagicID
    MagicID magic = dgst[0] + 256*dgst[1] + 65536*dgst[2] + 16777216*dgst[3];

    return magic;
}
