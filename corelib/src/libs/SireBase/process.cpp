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

#include "sireglobal.h"

#include <QMutex>
#include <QList>
#include <QElapsedTimer>

#include <boost/weak_ptr.hpp>

#include "sire_process.h"

#include "SireError/errors.h"

#ifndef Q_OS_WIN
#include <fcntl.h>      // CONDITIONAL_INCLUDE
#include <unistd.h>     // CONDITIONAL_INCLUDE
#include <signal.h>     // CONDITIONAL_INCLUDE
#include <sys/wait.h>   // CONDITIONAL_INCLUDE
#else
#define WIN32_LEAN_AND_MEAN
#define _WIN32_WINNT 0x0500
#include <Windows.h>
#include <Objbase.h>
#include <QVector>

static QVector<wchar_t> toWCharVec(const QString &str)
{
    int s;
    QVector<wchar_t> str_w(str.size() + 1);
    s = str.toWCharArray(str_w.data());
    str_w[s] = L'\0';
    str_w.resize(s + 1);
    return str_w;
}
#endif
#include <errno.h>      // CONDITIONAL_INCLUDE
#include <string>     // CONDITIONAL_INCLUDE

#include <sys/stat.h>   // CONDITIONAL_INCLUDE
#include <sys/types.h>  // CONDITIONAL_INCLUDE

#include <QDebug>

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
    ProcessData() :
        #ifndef Q_OS_WIN
        pid(-1),
        #else
        proc_info({}),
        startup_info({}),
        job_handle(NULL),
        #endif
        is_running(false), is_error(false), was_killed(false)
    {}

    ~ProcessData()
    {}

    QMutex datamutex;

    /** The command being run */
    QString command;

    /** The arguments to the command */
    QStringList arguments;

    #ifndef Q_OS_WIN
    /** The process ID of the process */
    pid_t pid;

    /** Kill the whole process group */
    void kill_process_group()
    {
        killpg(pid, SIGKILL);
    }

    /** Close all process handles */
    void close_handles()
    {
        is_running = false;
        pid = 0;
    }
    #else
    enum exitcode
    {
        EXIT_CODE_KILLED = 30001
    };

    /** The PROCESS_INFORMATION of the process */
    PROCESS_INFORMATION proc_info;

    /** The STARTUPINFOW of the process */
    STARTUPINFOW startup_info;

    /** The HANDLE of the job to which the process is assigned */
    HANDLE job_handle;

    /** Kill the whole process group */
    void kill_process_group()
    {
        if (job_handle)
            TerminateJobObject(job_handle, EXIT_CODE_KILLED);
    }

    /** Close all process handles */
    void close_handles()
    {
        if (proc_info.hProcess)
            CloseHandle(proc_info.hProcess);
        proc_info.hProcess = NULL;
        if (proc_info.hThread)
            CloseHandle(proc_info.hThread);
        proc_info.hThread = NULL;
        if (job_handle)
            CloseHandle(job_handle);
        job_handle = NULL;
        if (startup_info.hStdOutput)
            CloseHandle(startup_info.hStdOutput);
        startup_info.hStdOutput = NULL;
        if (startup_info.hStdError)
            CloseHandle(startup_info.hStdError);
        startup_info.hStdError = NULL;
        is_running = false;
    }
    #endif

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
{
}

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
void Process::cleanUpJob(int, int child_exit_status)
{
    if (d.get() == 0)
        return;

    #ifndef Q_OS_WIN
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
    #else
    if (child_exit_status)
    {
        //something went wrong with the job
        d->is_error = true;
    }
    if (child_exit_status == ProcessData::EXIT_CODE_KILLED)
        d->was_killed = true;
    #endif

    // make sure that all of the child processes have finished
    // by killing the child's process group
    d->kill_process_group();
    d->close_handles();
}

/** From the return value in 'child_exit_status' work out
    and return whether or not the child process is still running */
static bool processRunning(int child_exit_status)
{
    #ifndef Q_OS_WIN
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
    #else
    if (child_exit_status == STILL_ACTIVE)
        return true;
    #endif

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

    #ifndef Q_OS_WIN
    int child_exit_status;
    int status = waitpid(d->pid, &child_exit_status, 0);

    if (status == -1)
    {
        return;
    }
    else if (status != d->pid)
    {
        qDebug() << "waitpid exited with the wrong PID (" << status
                 << "vs." << d->pid << ")" << strerror(errno);
        return;
    }
    #else
    DWORD status = WaitForSingleObject(d->proc_info.hProcess, INFINITE);
    if (status == WAIT_FAILED)
    {
        qDebug() << "WaitForSingleObject() exited with error code" << GetLastError();
        return;
    }
    DWORD child_exit_status;
    status = GetExitCodeProcess(d->proc_info.hProcess, &child_exit_status);
    if (!status)
    {
        qDebug() << "GetExitCodeProcess() exited with error code" << GetLastError();
        return;
    }
    #endif

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

/** Test the process to see if it's still running. If it has just finished, clean up
    and set appropriate flags. WARNING: You must hold datamutex while calling this
    function. */
void Process::testWait()
{
    if (not d->is_running)
        return;

    #ifndef Q_OS_WIN
    int child_exit_status;
    int status;

    status = waitpid(d->pid, &child_exit_status, WNOHANG);

    if (status == 0)
        return;

    if (status == -1)
    {
        qDebug() << "waitpid exited with status -1!" << strerror(errno);

        d->is_running = false;
        d->pid = 0;
        return;
    }
    else if (status != d->pid)
    {
        qDebug() << "waitpid exited with the wrong PID (" << status
                    << "vs." << d->pid << ")" << strerror(errno);

        d->is_running = false;
        d->pid = 0;
        return;
    }
    #else
    DWORD child_exit_status;
    BOOL status = GetExitCodeProcess(d->proc_info.hProcess, &child_exit_status);
    if (!status)
    {
        qDebug() << "GetExitCodeProcess() exited with error code" << GetLastError();

        d->is_running = false;
        return;
    }
    #endif

    if ( processRunning(child_exit_status) )
    {
        return;
    }
    else
    {
        //the job has finished - process the finished job
        this->cleanUpJob(status, child_exit_status);
        return;
    }
}

/** Wait until the process has finished, or until 'ms' milliseconds have passed.
    This returns whether or not the process has finished */
bool Process::wait(int ms)
{
    if (d.get() == 0)
        return true;

    QElapsedTimer t;
    t.start();

    #if QT_VERSION >= QT_VERSION_CHECK(4, 3, 0)
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
            while (t.elapsed() < ms)
            {
                testWait();

                if (d->is_running)
                {
                    Sleeper::msleep(25);
                }
                else
                {
                    d->datamutex.unlock();
                    return true;
                }
            }

            d->datamutex.unlock();
            return false;
        }
        catch(...)
        {
            bool is_running = d->is_running;
            d->datamutex.unlock();
            return is_running;
        }
    }
    else
        return false;
}

/** Run the command 'command' with the arguments 'arguments', and
    return a Process object that can be used to query and control the
    job */
Process Process::run(const QString &command, const QStringList &arguments)
{
    return Process::run(command, arguments, QString(), QString());
}

/** Run the command 'command' with the arguments 'arguments', and
    return a Process object that can be used to query and control the
    job. Stdout and stderr of the running process are redirected to
    the user specified files. */
Process Process::run(const QString &command, const QStringList &arguments,
    const QString &stdout_file, const QString &stderr_file)
{
    QString error_msg = QString("Problem running the command %1 "
                                "with arguments %2.").arg(command).arg(
                                arguments.join(" "));

    #ifndef WIN32
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

        // Redirect stdout.
        if (not stdout_file.isNull())
        {
            int out = open(stdout_file.toLatin1().data(),
                O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IRGRP | S_IWGRP | S_IWUSR);
            dup2(out, STDOUT_FILENO);
            close(out);
        }

        // Redirect stderr.
        if (not stderr_file.isNull())
        {
            int err = open(stderr_file.toLatin1().data(),
                O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IRGRP | S_IWGRP | S_IWUSR);
            dup2(err, STDERR_FILENO);
            close(err);
        }

        //now run the command
        int status = execvp( char_args[0], char_args );

        delete[] char_args;

        if (status != 0)
        {
            qDebug() << error_msg
                     << "Status ==" << status << ", error ="
                     << strerror(errno);

            exit(-1);
        }

        exit(0);
    }
    else
    #else
    QVector<wchar_t> command_w = toWCharVec(command);
    QVector<wchar_t> args_w = toWCharVec(QString("%1 %2").arg(command).arg(arguments.join(" ")));
    HANDLE h_nul_write = INVALID_HANDLE_VALUE;
    if (stdout_file.isNull() || stderr_file.isNull())
    {
        h_nul_write = CreateFileW(L"nul", GENERIC_WRITE, FILE_SHARE_WRITE,
            NULL, OPEN_EXISTING, 0, NULL);
        if (h_nul_write == INVALID_HANDLE_VALUE) {
            qDebug() << error_msg << "CreateFileW(\"nul\") for writing failed.";
            return Process();
        }
    }
    SECURITY_ATTRIBUTES out_sa_attr = {};
    out_sa_attr.nLength = static_cast<DWORD>(sizeof(SECURITY_ATTRIBUTES));
    out_sa_attr.bInheritHandle = TRUE;
    SECURITY_ATTRIBUTES err_sa_attr(out_sa_attr);
    PROCESS_INFORMATION proc_info = {};
    HANDLE job_handle = NULL;
    STARTUPINFOW startup_info = {};

    // Redirect stdin from NUL.
    startup_info.hStdInput = CreateFileW(L"nul", GENERIC_READ,
        FILE_SHARE_READ, NULL, OPEN_EXISTING, 0, NULL);
    if (startup_info.hStdInput == INVALID_HANDLE_VALUE) {
        qDebug() << error_msg << "CreateFileW(\"nul\") for reading failed.";
        return Process();
    }

    // Redirect stdout.
    if (!stdout_file.isNull())
    {
        QVector<wchar_t> stdout_file_w = toWCharVec(stdout_file);
        startup_info.hStdOutput = CreateFileW(stdout_file_w.data(),
            GENERIC_WRITE, FILE_SHARE_WRITE, &out_sa_attr,
            CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
        if (startup_info.hStdOutput == INVALID_HANDLE_VALUE) {
            qDebug() << error_msg << "CreateFileW(" << stdout_file << ") failed.";
            return Process();
        }
    }
    else
        startup_info.hStdOutput = h_nul_write;

    // Redirect stderr.
    if (stdout_file == stderr_file)
        startup_info.hStdError = startup_info.hStdOutput;
    else if (!stderr_file.isNull())
    {
        QVector<wchar_t> stderr_file_w = toWCharVec(stderr_file);
        startup_info.hStdError = CreateFileW(stderr_file_w.data(),
            GENERIC_WRITE, FILE_SHARE_WRITE, &err_sa_attr,
            CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
        if (startup_info.hStdError == INVALID_HANDLE_VALUE) {
            qDebug() << error_msg << "CreateFileW(" << stderr_file << ") failed.";
            return Process();
        }
    }
    else
        startup_info.hStdError = h_nul_write;

    startup_info.cb = sizeof(startup_info);
    startup_info.dwFlags = STARTF_USESTDHANDLES;

    // Generate a GUID to name the job object.
    GUID guid = {};
    if (CoCreateGuid(&guid) != S_OK)
    {
        qDebug() << error_msg << "CoCreateGuid() failed.";
        return Process();
    }
    static const int MAX_GUID_LEN = 40;
    QVector<wchar_t> guid_w(MAX_GUID_LEN, L'\0');
    QByteArray guid_a(MAX_GUID_LEN, '\0');
    size_t guid_len = StringFromGUID2(guid, guid_w.data(), guid_w.size());
    if (!guid_len)
    {
        qDebug() << error_msg << "StringFromGUID2() failed.";
        return Process();
    }
    guid_w.resize(static_cast<int>(guid_len));
    guid_len = WideCharToMultiByte(CP_ACP, 0, guid_w.data(), -1,
        guid_a.data(), guid_a.size(), NULL, NULL);
    if (!guid_len)
    {
        qDebug() << error_msg << "WideCharToMultiByte() failed.";
        return Process();
    }
    guid_a.resize(static_cast<int>(guid_len));

    // Create the job object with the generated GUID as name.
    job_handle = CreateJobObject(NULL, guid_a.data());
    if (!job_handle)
    {
        qDebug() << error_msg << "CreateJobObject() failed.";
        return Process();
    }

    // Configure the job object to kill its children when it dies.
    JOBOBJECT_EXTENDED_LIMIT_INFORMATION jeli = {};
    jeli.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
    if (!SetInformationJobObject(job_handle,
        JobObjectExtendedLimitInformation, &jeli, sizeof(jeli)))
    {
        qDebug() << error_msg << "SetInformationJobObject() failed.";
        return Process();
    }

    // Now create a process to run the command.
    if (!CreateProcessW(command_w.data(), args_w.data(), NULL, NULL, TRUE,
        CREATE_NEW_PROCESS_GROUP | CREATE_NO_WINDOW, NULL, NULL, &startup_info, &proc_info)) {
        qDebug() << error_msg << "CreateProcessW() failed with error code" << GetLastError();
        return Process();
    }

    // Assign the process to a new job object.
    if (!AssignProcessToJobObject(job_handle, proc_info.hProcess))
    {
        qDebug() << error_msg << "AssignProcessToJobObject() failed with error code" << GetLastError();
        return Process();
    }

    // Close the stdin handle as we are not using it.
    if (startup_info.hStdInput)
        CloseHandle(startup_info.hStdInput);
    startup_info.hStdInput = NULL;
    #endif
    {
        //parent
        Process p;

        p.d.reset( new ProcessData() );
        #ifndef Q_OS_WIN
        p.d->pid = pid;
        #else
        p.d->proc_info = proc_info;
        p.d->startup_info = startup_info;
        p.d->job_handle = job_handle;
        #endif

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
    return Process::run(command, QStringList(), QString(), QString());
}

/** Run the command 'command' and return a Process object that can be
    used to monitor the command. Stdout and stderr of the running
    process are redirected to the user specified files. */
Process Process::run(const QString &command, const QString &stdout_file, const QString &stderr_file)
{
    return Process::run(command, QStringList(), stdout_file, stderr_file);
}

/** Run the command 'command' with the solitary argument 'arg' */
Process Process::run(const QString &command, const QString &arg)
{
    QStringList args;
    args.append(arg);

    return Process::run(command, args, QString(), QString());
}

/** Run the command 'command' with the solitary argument 'arg'.
    Stdout and stderr of the running process are redirected to
    the user specified files. */
Process Process::run(const QString &command, const QString &arg,
    const QString &stdout_file, const QString &stderr_file)
{
    QStringList args;
    args.append(arg);

    return Process::run(command, args, stdout_file, stderr_file);
}

/** Kill this process */
void Process::kill()
{
    if (d.get() == 0)
        return;

    //kill the job
    if (d->is_running)
    {
        testWait();

        if (d->is_running)
        {
            qDebug() << "Killing job " << d->command.toUtf8().constData()
                     << d->arguments.join(" ").toUtf8().constData();
            d->kill_process_group();
        }
        else
            return;
    }

    //now wait for it to finish
    if (not this->wait(1000))
    {
        //it still hasn't finished - print a warning
        qDebug() << "...job still not dead. You may want to check it "
                 << "is not still running.";
    }

    QMutexLocker lkr(&d->datamutex);

    d->was_killed = true;
    d->close_handles();
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
    testWait();
    return d->is_running;
}

/** Return whether or not this process has finished running */
bool Process::hasFinished()
{
    if (d.get() == 0)
        return true;

    QMutexLocker lkr( &(d->datamutex) );
    testWait();
    return not d->is_running;
}

/** Return whether or not the process exited in error */
bool Process::isError()
{
    if (d.get() == 0)
        return false;

    QMutexLocker lkr( &(d->datamutex) );
    testWait();
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
