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

#include <QMutex>
#include <QList>
#include <QTime>

#include <boost/weak_ptr.hpp>

#include "sire_process.h"

#include "SireError/errors.h"

#include <unistd.h>  // CONDITIONAL_INCLUDE
#include <signal.h>  // CONDITIONAL_INCLUDE
#include <errno.h>   // CONDITIONAL_INCLUDE
#include <string.h>  // CONDITIONAL_INCLUDE

#include <sys/types.h>  // CONDITIONAL_INCLUDE
#include <sys/wait.h>   // CONDITIONAL_INCLUDE

#include <QDebug>

#ifndef Q_OS_UNIX
    #error This file can only be compiled on a UNIX-like machine!
#endif

using namespace SireBase;
using boost::shared_ptr;
using boost::weak_ptr;

namespace SireBase
{
namespace detail
{

/** Private implementation of Process */
class ProcessData
{
public:
    ProcessData() : pid(-1), is_running(false), is_error(false), was_killed(false)
    {}
    
    ~ProcessData()
    {}
    
    QMutex datamutex;

    /** The command being run */
    QString command;
    
    /** The arguments to the command */
    QStringList arguments;
    
    /** The process ID of the process */
    pid_t pid;
   
    /** Whether or not the process is running */
    bool is_running;
    
    /** Whether or not the job exited with an error */
    bool is_error;
    
    /** Whether or not the job was killed */
    bool was_killed;
};

} // end of namespace detail
} // end of namespace SireBase

using namespace SireBase::detail;

Q_GLOBAL_STATIC( QList< weak_ptr<ProcessData> >, processRegistry );
Q_GLOBAL_STATIC( QMutex, registryMutex );

/** Null constructor */
Process::Process()
{}

/** Copy constructor */
Process::Process(const Process &other) : d(other.d)
{}

/** Destructor */
Process::~Process()
{
    if (d.unique())
    {
        this->kill();
    }
}

/** Copy assignment operator */
Process& Process::operator=(const Process &other)
{
    if (d.get() != other.d.get())
        d = other.d;
    
    return *this;
}

/** Comparison operator */
bool Process::operator==(const Process &other) const
{
    return d.get() == other.d.get();
}

/** Comparison operator */
bool Process::operator!=(const Process &other) const
{
    return d.get() != other.d.get();
}

/** Internal function used to clean up a running job that has 
    just finished - only call this function if you are
    holding d->datamutex */
void Process::cleanUpJob(int status, int child_exit_status)
{
    if (d.get() == 0)
        return;
        
    if (WEXITSTATUS(child_exit_status) != 0)
    {
        //something went wrong with the job
        d->is_error = true;
    }

    if (WIFSIGNALED(child_exit_status))
    {
        if (WTERMSIG(child_exit_status) == SIGKILL or
            WTERMSIG(child_exit_status) == SIGHUP)
        {
            d->was_killed = true;
            d->is_error = true;
        }
        else
        {
            d->was_killed = false;
            d->is_error = true;
        }
    }
    else if (WIFSTOPPED(child_exit_status))
    {
        if (WSTOPSIG(child_exit_status) == SIGKILL or
            WSTOPSIG(child_exit_status) == SIGHUP)
        {
            d->was_killed = true;
            d->is_error = true;
        }
        else
        {
            d->was_killed = false;
            d->is_error = true;
        }
    }
    
    //make sure that all of the child processes have finished
    //by killing the child's process group
    killpg(d->pid, SIGKILL);
    
    d->is_running = false;
}

/** From the return value in 'child_exit_status' work out 
    and return whether or not the child process is still running */
static bool processRunning(int child_exit_status)
{
    if (WIFEXITED(child_exit_status) == 0)
    {
        if (WIFSIGNALED(child_exit_status) != 0 or
            WIFSTOPPED(child_exit_status) != 0)
        {
            //the job was killed or stopped by a signal
            return false;
        }
        else
            return true;
    }
    
    //the job has exited normally
    return false;
}

/** Wait until the process has finished */
void Process::wait()
{
    if (d.get() == 0)
        return;
    
    QMutexLocker lkr( &(d->datamutex) );
    
    if (not d->is_running)
        return;
    
    int child_exit_status;
    int status = waitpid(d->pid, &child_exit_status, 0);
    
    if (status == -1)
    {
        qDebug() << "waitpid exited with status -1!" << strerror(errno);
        return;
    }
    else if (status != d->pid)
    {
        qDebug() << "waitpid exited with the wrong PID (" << status
                 << "vs." << d->pid << ")" << strerror(errno);
        return;
    }
    
    if ( processRunning(child_exit_status) )
    {
        qDebug() << "child process has not finished!";
    }
    else
    {
        this->cleanUpJob(status, child_exit_status);
    }
}

#include <QThread>

class Sleeper : protected QThread
{
public:
    static void msleep(unsigned long msecs)
    {
        QThread::msleep(msecs);
    }
};

/** Wait until the process has finished, or until 'ms' milliseconds have passed.
    This returns whether or not the process has finished */
bool Process::wait(int ms)
{
    if (d.get() == 0)
        return true;
        
    QTime t;
    t.start();
    
    #if QT_VERSION >= 0x040300
    if (d->datamutex.tryLock(ms))
    #else
    if (d->datamutex.tryLock())
    #endif
    {
        if (not d->is_running)
        {
            d->datamutex.unlock();
            return true;
        }
    
        try
        {
            int child_exit_status;
            int status;
        
            while (t.elapsed() < ms)
            {
                status = waitpid(d->pid, &child_exit_status, WNOHANG);

                if (status == 0)
                {
                    Sleeper::msleep(25);
                    continue;
                }
                else if (status == -1)
                {
                    qDebug() << "waitpid exited with status -1!" << strerror(errno);
                    return true;
                }
                else if (status != d->pid)
                {
                    qDebug() << "waitpid exited with the wrong PID (" << status
                             << "vs." << d->pid << ")" << strerror(errno);
                    return true;
                }
                
                if ( processRunning(child_exit_status) )
                {
                    //the job is still running - wait then
                    //try again
                    Sleeper::msleep(25);
                }
                else
                {
                    //the job has finished - process the finished job
                    this->cleanUpJob(status, child_exit_status);
                    d->datamutex.unlock();
                    return true;
                }
            }
            
            d->datamutex.unlock();
            return false;
        }
        catch(...)
        {
            d->datamutex.unlock();
            return true;
        }
    }
    else
        return false;
}

/** Run the command 'command' with the arguments 'arguments', and
    return a Process object that can be used to query and control the 
    job */
Process Process::run(const QString &command,  const QStringList &arguments)
{
    //fork to run the command
    pid_t pid = fork();
    
    if (pid == -1)
    {
        throw SireError::unsupported( QObject::tr(
            "It is not possible to use fork on this platform. Running "
            "the external command \"%1\" is therefore not possible.")
                .arg(command), CODELOC );
    }
    else if (pid == 0)
    {
        //this is the child process!

        //move this child (and all of its children) into a new 
        //process ID group - the new group will have the same
        //process ID as this child
        setpgrp();
        
        QByteArray cmd = command.toUtf8();
        QList<QByteArray> args;
        char** char_args = new char*[ arguments.count() + 2 ];
        
        char_args[0] = cmd.data();
        
        for (int i=0; i<arguments.count(); ++i)
        {
            args.append( arguments[i].toUtf8() );
            char_args[i+1] = args[i].data();
        }
        
        char_args[arguments.count()+1] = 0;
       
        //now run the command
        int status = execvp( char_args[0], char_args );
        
        delete[] char_args;
        
        if (status != 0)
        {
            qDebug() << "Problem starting the command" << command 
                     << "with arguments" << arguments.join(" ")
                     << ". Status ==" << status << ", error =" 
                     << strerror(errno);
                     
            exit(-1);
        }
        
        exit(0);
    }
    else
    {
        //parent
        Process p;
        
        p.d.reset( new ProcessData() );
        p.d->pid = pid;
        
        p.d->is_running = true;

        p.d->command = command;
        p.d->arguments = arguments;

        //record this process in the list of running processes
        QMutexLocker lkr( registryMutex() );
        processRegistry()->append( weak_ptr<ProcessData>(p.d) );

        return p;
    }
    
    return Process();
}
           
/** Run the command 'command' and return a Process object that can be
    used to monitor the command */        
Process Process::run(const QString &command)
{
    return Process::run(command, QStringList());
}

/** Run the command 'command' with the solitary argument 'arg' */
Process Process::run(const QString &command, const QString &arg)
{
    QStringList args;
    args.append(arg);
    
    return Process::run(command, args);
}

/** Kill this process */
void Process::kill()
{
    if (d.get() == 0)
        return;
    
    //kill the job
    if (d->is_running)
    {
        qDebug() << "Killing job " << d->command.toUtf8().constData()
                 << d->arguments.join(" ").toUtf8().constData();
        killpg(d->pid, SIGKILL);
    }
    
    //now wait for it to finish
    if (not this->wait(1000))
    {
        //it still hasn't finished - print a warning
        qDebug() << "...job still not dead. You may want to check it "
                 << "is not still running.";
    }
}

/** Use this function to kill all of the jobs that are currently running */
void Process::killAll()
{
    QMutexLocker lkr( registryMutex() );
    
    QList< weak_ptr<ProcessData> > &process_list = *(processRegistry());

    for (QList< weak_ptr<ProcessData> >::iterator it = process_list.begin();
         it != process_list.end();
         ++it)
    {
        Process p;
        p.d = it->lock();
        
        p.kill();
    }
    
    process_list.clear();
}

/** Return whether or not the job is running */
bool Process::isRunning()
{
    if (d.get() == 0)
        return false;
        
    QMutexLocker lkr( &(d->datamutex) );
    return d->is_running;
}

/** Return whether or not this process has finished running */
bool Process::hasFinished()
{
    if (d.get() == 0)
        return true;
        
    QMutexLocker lkr( &(d->datamutex) );
    return not d->is_running;
}

/** Return whether or not the process exited in error */
bool Process::isError()
{
    if (d.get() == 0)
        return false;

    QMutexLocker lkr( &(d->datamutex) );
    return d->is_error;
}

/** Return whether or not the process was killed */
bool Process::wasKilled()
{
    if (d.get() == 0)
        return false;

    QMutexLocker lkr( &(d->datamutex) );
    return d->was_killed;
}
