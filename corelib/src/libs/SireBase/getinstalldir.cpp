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

#include "sire_config.h"
#include "sire_version.h"
#include "getinstalldir.h"

#include <QDir>
#include <QFileInfo>

#include <QDebug>

#ifdef Q_OS_MAC
    extern "C" { int _NSGetExecutablePath(char* buf, uint32_t* bufsize); }
#endif

#include "SireError/errors.h"

namespace SireBase
{
    static QString install_dir;

    QString stripDir(QString last_part, QString path)
    {
        if (last_part.isEmpty() or path.isEmpty())
            return path;
    
        int j = path.lastIndexOf(last_part);
        
        if (j == path.length() - last_part.length())
        {
            path.remove(j, last_part.length());
            return path;
        }
        else
        {
            qDebug() << "CANNOT REMOVE" << last_part << "FROM" << path
                 << "AS IT IS NOT IN THE PATH?";

            return path;
        }
    }

    /** This function returns the URL of the source repository that contains
        the core Sire library */
    QString SIREBASE_EXPORT getRepositoryURL()
    {
        return QString(SIRE_REPOSITORY_URL);
    }
    
    /** This function returns the version number(s) of this copy of the corelib
        from the online source repository */
    QString SIREBASE_EXPORT getRepositoryVersion()
    {
        return QString(SIRE_REPOSITORY_VERSION);
    }

    /** This function is used to set the path to the installation directory.
        This can be used to override the path found from the running
        executable. This is useful when Sire is loaded as an external
        module or library from another executable */
    void SIREBASE_EXPORT setInstallDir(QString dir)
    {
        QDir d(dir);
        
        if (not d.exists())
            throw SireError::file_error( QObject::tr(
                "You cannot set the installation directory of Sire to a value "
                "that doesn't exist (%1).")
                    .arg(dir), CODELOC );
        
        install_dir = d.absolutePath();
    }

    /** This function locates the directory where Sire is installed. This
        returns the install prefix of the installation, e.g. the base
        from which the bin, lib etc. directories can be found */
    QString SIREBASE_EXPORT getInstallDir()
    {
        if (not install_dir.isEmpty())
            return install_dir;
    
        //first, find the full path to the running executable. We assume that
        //we are using a Sire executable
        
        //we follow the instructions from a stackoverflow answer by mark40
        //<http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe/1024937#1024937>
        //Mac OS X: _NSGetExecutablePath() (man 3 dyld)
        //Linux: readlink /proc/self/exe
        //Solaris: getexecname()
        //FreeBSD: sysctl CTL_KERN KERN_PROC KERN_PROC_PATHNAME -1
        //BSD with procfs: readlink /proc/curproc/file
        //Windows: GetModuleFileName() with hModule = NULL
        
        #ifdef Q_OS_MAC
            char pathbuf[PATH_MAX + 1];
            uint32_t bufsize = PATH_MAX+1;
            int ok = _NSGetExecutablePath((char*)pathbuf, &bufsize);
        
            if (ok != 0)
                throw SireError::program_bug( QObject::tr(
                        "For some reason, _NSGetExecutablePath has not worked! (%1)")
                            .arg(ok), CODELOC );
        
            QFileInfo f(pathbuf);

            //sometimes we may use a python executable in python.app/Contents/MacOS.
            //we will use a special case to remove this from the path
            QString path = f.canonicalPath();
            if (path.endsWith("python.app/Contents/MacOS"))
            {
                path.chop(25);
                setInstallDir(path);
            }
            else
            {
                setInstallDir( stripDir(SIRE_BIN_DIR,path) );
            }

            return install_dir;
        #else
        #ifdef Q_OS_LINUX
            QFileInfo f( "/proc/self/exe" );
        
            if (not f.exists())
                throw SireError::program_bug( QObject::tr(
                        "For some reason /proc/self/exe does not exist for your "
                        "version of Linux. Let the developers know and we will get it "
                        "working for you."), CODELOC );
        
            setInstallDir( stripDir(SIRE_BIN_DIR,f.canonicalPath()) );
            return install_dir;
        #else
            throw SireError::unavailable_code( QObject::tr(
                    "Ask the Sire developers to write the \"getInstallDir\" function "
                    "for your platform. Sorry that it has yet to be written."), CODELOC );
            return QString::null;
        #endif
        #endif
    }
    
    /** This function is used to get the full path to the file or directory 'path'
        that is contained in the Sire directory, e.g. "bin" would
        return SIRE_INSTALL_PREFIX/bin. An exception will be raised
        if this file or directory doesn't exist, and 'assert_exists'
        is true */
    QString SIREBASE_EXPORT getSireDir(const QString &path, bool assert_exists)
    {
        QDir dir(getInstallDir());
        
        if (assert_exists)
        {
            if (not dir.exists(path))
                throw SireError::file_error( QObject::tr(
                        "Cannot find the file or directory \"%2\" in the Sire install "
                        "directory \"%1\" (%3).")
                            .arg(dir.absolutePath(), path)
                            .arg(dir.absoluteFilePath(path)), CODELOC );
        }
        
        return dir.absoluteFilePath(path);
    }
    
    /** This returns the directory containing the Sire executables */
    QString SIREBASE_EXPORT getBinDir()
    {
        return getSireDir(SIRE_BIN_DIR);
    }
    
    /** This returns the directory containing the Sire libraries */
    QString SIREBASE_EXPORT getLibDir()
    {
        return getSireDir(SIRE_LIBS_DIR);
    }
    
    /** This returns the directory containing the Sire bundled libraries */
    QString SIREBASE_EXPORT getBundledLibDir()
    {
        return getSireDir(SIRE_BUNDLED_LIBS_DIR);
    }
    
    /** This returns the directory containing the Sire support files */
    QString SIREBASE_EXPORT getShareDir()
    {
        return getSireDir(SIRE_SHARE_DIR);
    }
    
    /** This returns the release version of Sire */
    QString SIREBASE_EXPORT getReleaseVersion()
    {
        return "2014.2";
    }
}
