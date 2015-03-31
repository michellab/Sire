/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREBASE_TEMPDIR_H
#define SIREBASE_TEMPDIR_H

#include <QDir>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireBase
{

/** This class provides a temporary directory that is created in the 
    constructor, and deleted in the destructor.
    
    @author Christopher Woods
*/
class SIREBASE_EXPORT TempDir
{
public:
    TempDir();
    TempDir(const QString &temp_root);

    ~TempDir();
    
    static const char* typeName()
    {
        return "SireBase::TempDir";
    }
    
    const char* what() const
    {
        return TempDir::typeName();
    }
    
    QString toString() const;
    
    QString path() const;
    
    void doNotDelete();
    
private:
    TempDir(const TempDir&)
    {}
    
    TempDir& operator=(const TempDir&)
    {
        return *this;
    }

    void createDirectory(const QString &temp_root);

    /** The temporary directory */
    QDir tmpdir;
    
    /** Whether or not to delete the directory at exit */
    bool do_not_delete;
};

}

SIRE_EXPOSE_CLASS( SireBase::TempDir )

SIRE_END_HEADER

#endif
