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

#include <QUuid>
#include <QMutex>

#include "tempdir.h"

#include "tostring.h"

#include "SireError/errors.h"

#include <QDebug>

using namespace SireBase;

static QString getUserName()
{
    #ifdef Q_OS_UNIX
        return std::getenv("USER");
    #else
        return "USER";
    #endif
}

static QMutex tmpdir_mutex;

static QDir createDirectory(const QString &temp_root, int ntries)
{
    QString dirname = QString("%1/%2_sire_%3")
                            .arg(temp_root, getUserName(),
                                 QUuid::createUuid().toString());
                                 
    QMutexLocker lkr(&tmpdir_mutex);
    
    QDir tmpdir;
    
    if (tmpdir.exists(dirname))
    {
        //this directory already exists - try again
        lkr.unlock();
        return ::createDirectory(temp_root, ntries);
    }

    if (not tmpdir.mkdir(dirname))
    {
        if (ntries > 0)
        {
            return ::createDirectory(temp_root, ntries-1);
        }
        else
            throw SireError::file_error( QObject::tr(
                "Could not create the temporary directory \"%1\". Please "
                "ensure that there is enough disk space and that you have "
                "permission to write to the directory \"%2\".")
                    .arg(dirname, temp_root), CODELOC );
    }
    
    return QDir(dirname);
}

/** Create a new temporary directory in QDir::tempPath() */
void TempDir::createDirectory(const QString &temp_root)
{
    tmpdir = ::createDirectory(temp_root, 5);
}

/** This creates a new temporary directory in QDir::tempPath() */
TempDir::TempDir() : do_not_delete(false)
{
    this->createDirectory( QDir::tempPath() );
}

/** This creates a new temporary directory in 'temp_root' */
TempDir::TempDir(const QString &temp_root) : do_not_delete(false)
{
    this->createDirectory(temp_root);
}

static void removeDirectory(QDir &dir)
{
    //get the list of all files in this directory
    // ( cannot use NoDotAndDotDot as this doesn't work on my mac)
    QStringList files = dir.entryList();
    
    foreach (const QString &filename, files)
    {
        if (filename == "." or filename == "..")
            continue;
    
        if (not dir.exists(filename))
            continue;
       
        QFileInfo fileinfo( dir.absoluteFilePath(filename) );
    
        if ( fileinfo.isDir() )
        {
            //remove this directory
            QDir subdir( fileinfo.absoluteFilePath() );
            ::removeDirectory(subdir);
        }
        else
        {
            dir.remove(filename);
        }
    }
    
    //now remove this directory
    dir.rmdir(dir.absolutePath());
}

/** The destructor deletes the temporary directory and
    anything that it contains */
TempDir::~TempDir()
{
    if (not do_not_delete)
        ::removeDirectory(tmpdir);
}

/** This returns the complete path to the temporary directory */
QString TempDir::path() const
{
    return tmpdir.absolutePath();
}

/** Return a string representation of this TempDir */
QString TempDir::toString() const
{
    return QString("TempDir(\"%1\")").arg(this->path());
}

/** Tell the TempDir not to delete the directory when this object
    is deleted - this can be used when you are debugging to prevent
    the directory from disappearing! */
void TempDir::doNotDelete()
{
    do_not_delete = true;
}
