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

#include "findexe.h"

#include <QRegExp>
#include <QStringList>
#include <QDir>
#include <QProcess>

#include "SireError/errors.h"

QFileInfo SIREBASE_EXPORT SireBase::findExe(const QString &exe)
{
    //does it exist in the current directory?
    QFileInfo progfile( exe );

    if (progfile.isExecutable())
        return progfile;

    //try to find exe in the application PATH
    QStringList env_variables = QProcess::systemEnvironment();

    QRegExp rexp("^PATH=(.*)");
    rexp.setCaseSensitivity(Qt::CaseInsensitive);

    QStringList path;

    foreach (QString env_variable, env_variables)
    {
        if (rexp.indexIn(env_variable) != -1)
        {
            path = rexp.cap(1).split(":", QString::SkipEmptyParts);
            break;
        }
    }

    foreach (QString dir, path)
    {
        QFileInfo new_exe( QDir(dir).filePath(progfile.fileName()) );

        if (new_exe.isExecutable())
        {
            //we have found the executable!
            return new_exe;
        }
    }

    throw SireError::process_error( QObject::tr(
              "Could not find the executable \"%1\" in the current directory, nor "
              "in the system PATH.").arg(exe), CODELOC );

    return QFileInfo();
}
