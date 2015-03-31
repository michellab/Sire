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

#include "functionsignature.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

static const RegisterMetaType<FunctionSignature> r_sig(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const FunctionSignature &sig)
{
    writeHeader(ds, r_sig, 1) << sig.name() << sig.args();

    return ds;
}
/** Deserialise from a binary data stream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, FunctionSignature &sig)
{
    VersionID v = readHeader(ds, r_sig);

    if (v == 1)
    {
        QString name;
        QSet<SireCAS::SymbolID> args;

        ds >> name >> args;

        sig = SireCAS::FunctionSignature(name,args);
    }
    else
        throw version_error(v, "1", r_sig, CODELOC);

    return ds;
}

const char* FunctionSignature::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FunctionSignature>() );
}
