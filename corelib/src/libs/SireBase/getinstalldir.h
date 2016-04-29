/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#ifndef SIREBASE_GETINSTALLDIR_H
#define SIREBASE_GETINSTALLDIR_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
    QString getInstallDir();
    void setInstallDir(QString dir);
    
    QString getBinDir();
    QString getLibDir();
    QString getBundledLibDir();
    QString getShareDir();
    
    QString getSireDir(const QString &path, bool assert_exists=true);
    
    QString getRepositoryURL();
    QString getRepositoryVersion();
    QString getRepositoryBranch();
    bool getRepositoryVersionIsClean();
    
    QString getReleaseVersion();
}

SIRE_EXPOSE_FUNCTION(SireBase::getInstallDir)
SIRE_EXPOSE_FUNCTION(SireBase::setInstallDir)
SIRE_EXPOSE_FUNCTION(SireBase::getBinDir)
SIRE_EXPOSE_FUNCTION(SireBase::getLibDir)
SIRE_EXPOSE_FUNCTION(SireBase::getBundledLibDir)
SIRE_EXPOSE_FUNCTION(SireBase::getShareDir)
SIRE_EXPOSE_FUNCTION(SireBase::getSireDir)
SIRE_EXPOSE_FUNCTION(SireBase::getRepositoryURL)
SIRE_EXPOSE_FUNCTION(SireBase::getRepositoryVersion)
SIRE_EXPOSE_FUNCTION(SireBase::getRepositoryBranch)
SIRE_EXPOSE_FUNCTION(SireBase::getRepositoryVersionIsClean)
SIRE_EXPOSE_FUNCTION(SireBase::getReleaseVersion)

SIRE_END_HEADER

#endif
