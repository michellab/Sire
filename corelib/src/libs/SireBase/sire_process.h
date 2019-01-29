/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREBASE_PROCESS_H
#define SIREBASE_PROCESS_H

#include "sireglobal.h"

#include <QStringList>

#include <boost/shared_ptr.hpp>
#include <boost/throw_exception.hpp>

SIRE_BEGIN_HEADER

namespace SireBase
{

namespace detail
{
class ProcessData;
}

/** This class provides a means to run an external process
    (executable). This provides the equivalent of
    std::system, but with added error handling and
    signal handling (which ensures that any child processes
    are killed when Sire exits)

    @author Christopher Woods, Lester Hedges
*/
class SIREBASE_EXPORT Process
{
public:
    Process();
    Process(const Process &other);

    ~Process();

    Process& operator=(const Process &other);

    bool operator==(const Process &other) const;
    bool operator!=(const Process &other) const;

    static const char* typeName()
    {
        return "SireBase::Process";
    }

    const char* what() const
    {
        return Process::typeName();
    }

    void wait();
    bool wait(int ms);

    static Process run(const QString &command);
    static Process run(const QString &command, const QString &stdout_file,
        const QString &stderr_file);

    static Process run(const QString &command, const QString& arg);
    static Process run(const QString &command, const QString& arg,
        const QString &stdout_file, const QString &stderr_file);

    static Process run(const QString &command, const QStringList& arguments);
    static Process run(const QString &command, const QStringList& arguments,
        const QString &stdout_file, const QString &stderr_file);

    static void killAll();

    bool isRunning();
    bool hasFinished();

    bool isError();
    bool wasKilled();

    void kill();

private:
    void cleanUpJob(int status, int child_exit_status);

    void testWait();

    /** PIMPL pointer to the data for this process */
    boost::shared_ptr<detail::ProcessData> d;
};

}

SIRE_EXPOSE_CLASS( SireBase::Process )

SIRE_END_HEADER

#endif
