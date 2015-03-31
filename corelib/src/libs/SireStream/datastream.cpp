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

#include <QDataStream>

#include "datastream.h"
#include "magic_error.h"
#include "registeralternativename.h"
#include "version_error.h"

#include <QDebug>

namespace SireStream
{
    /** Write a header to the data stream that describes the type and version
        of the object that is about to be written */
    QDataStream SIRESTREAM_EXPORT &writeHeader(QDataStream &ds,
                                               const RegisterMetaTypeBase &r_type,
                                               VersionID version)
    {
        ds << r_type.magicID() << version;
        return ds;
    }

    /** Write a header to the data stream that contains the magic and version
        of the object that is about to be written */
    QDataStream SIRESTREAM_EXPORT &writeHeader(QDataStream &ds,
                                               MagicID magicid,
                                               VersionID version)
    {
        ds << magicid << version;
        return ds;
    }

    static QHash< QString,QSet<MagicID> > alt_ids;
    static QHash< QString,QSet<QString> > alt_names;

    namespace detail
    {
        /** Internal function used to get all of the alternative names for a class */
        QSet<QString> getAlternativeNames(QString type_name)
        {
            QHash< QString,QSet<QString> >::const_iterator it = alt_names.constFind(type_name);

            if (it == alt_names.constEnd())
                return QSet<QString>();
            else
                return it.value();
        }

        bool matchesMagic(const RegisterMetaTypeBase &r_type, MagicID id)
        {
            if (id == r_type.magicID())
                return true;
            
            QHash< QString,QSet<MagicID> >::const_iterator it
                        = alt_ids.constFind( QLatin1String(r_type.typeName()) );
            
            if (it != alt_ids.constEnd())
                return it->contains(id);
            else
                return false;
        }
    }

    /** Function used to register alternative names for a class. This allows a class
        to change names, but to still be loaded using its old name */
    void SIRESTREAM_EXPORT registerAlternativeName(const char *name, const char *alternative)
    {
        QLatin1String str(name);
        
        if (alt_ids.contains(str))
        {
            alt_ids[str].insert( getMagic(alternative) );
        }
        else
        {
            QSet<MagicID> alts;
            alts.insert( getMagic(alternative) );
            alt_ids.insert(str, alts);
        }
        
        QLatin1String altname(alternative);
        
        if (alt_names.contains(altname))
        {
            alt_names[altname].insert( QLatin1String(name) );
        }
        else
        {
            QSet<QString> names;
            names.insert( QLatin1String(name) );
            alt_names.insert(altname, names);
        }
    }

    /** Read the header of the binary object to check that the type is correct
        and to obtain the binary data version */
    VersionID SIRESTREAM_EXPORT readHeader(QDataStream &ds,
                                           const RegisterMetaTypeBase &r_type)
    {
        MagicID id;
        VersionID v;

        ds >> id >> v;

        if (not (id == r_type.magicID() or detail::matchesMagic(r_type, id)))
            throw SireStream::magic_error( id, r_type, CODELOC );

        return v;
    }

    /** Read the header of the binary object to check that the type is correct
        and to obtain the binary data version */
    VersionID SIRESTREAM_EXPORT readHeader(QDataStream &ds,
                                           MagicID magicid,
                                           const char *type_name)
    {
        MagicID id;
        VersionID v;

        ds >> id >> v;

        if (id != magicid)
            throw SireStream::magic_error( id, magicid, type_name, CODELOC );

        return v;
    }

} //end of namespace SireStream


