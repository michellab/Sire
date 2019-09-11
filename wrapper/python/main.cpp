
#include <Python.h>

#include <cstdio>

#include <QDir>

#include "tbb/task_scheduler_init.h"

#include "SireError/errors.h"
#include "SireError/printerror.h"

#include "SireBase/sire_process.h"
#include "SireBase/getinstalldir.h"
#include "SireBase/cpuid.h"

#include "SireCluster/cluster.h"

#include "sire_config.h"
#include "sire_python_config.h"

#include <boost/scoped_array.hpp>

using std::printf;

using namespace SireBase;
using namespace SireCluster;

#include <QDebug>

#ifdef Q_OS_UNIX

#include <signal.h>

//handle CTRL-C signal - this should kill the calculation
// - with thanks to
//  <http://www.gnu.org/software/libtool/manual/libc/Termination-in-Handler.html#Termination-in-Handler>

volatile sig_atomic_t fatal_error_in_progress = 0;

void fatal_error_signal (int sig)
{
    // Since this handler is established for more than one kind of signal,
    // it might still get invoked recursively by delivery of some other kind
    // of signal.  Use a static variable to keep track of that.
    if (fatal_error_in_progress)
        raise (sig);

    fatal_error_in_progress = 1;

    printf("You're killing me!!!\n");

    // Kill any child processes
    SireBase::Process::killAll();

    printf("\nI, and all of my children are now dead. Adieu...\n");

    // Now reraise the signal.  We reactivate the signal's
    // default handling, which is to terminate the process.
    // We could just call exit or abort,
    // but reraising the signal sets the return status
    // from the process correctly.
    signal (sig, SIG_DFL);
    raise (sig);
}

#endif // Q_OS_UNIX

int main(int argc, char **argv)
{
    int status = 0;

    try
    {
        //run through the command line arguments and filter out the ones we want
        //(we add 5 to the length as we may add arguments to the list)
        boost::scoped_array<wchar_t*> python_argv( new wchar_t*[argc+5] );
        int python_argc = 0;

        bool ignore_pythonpath = true;
        bool ignore_ipython = false;

        // by default, use all of the cores on the node
        CPUID cpuid;
        int ppn = cpuid.numCores();

        // See if the OMP_NUM_THREADS environment variable has been set.
        QString num_threads = qgetenv("OMP_NUM_THREADS");

        // If the variable has been set, try to convert to an int.
        if (not num_threads.isEmpty())
            ppn = num_threads.toInt();

        // See if the SIRE_NUM_THREADS environment variable has been set.
        // This takes precedence over OMP_NUM_THREADS.
        num_threads = qgetenv("SIRE_NUM_THREADS");

        // If the variable has been set, try to convert to an int.
        if (not num_threads.isEmpty())
            ppn = num_threads.toInt();

        if (ppn <= 0)
        {
            throw SireError::invalid_arg( QObject::tr(
                "Invalid OMP_NUM_THREADS or SIRE_NUM_THREADS environment "
                "variable! Must be a positive integer."));
        }

        QList< std::wstring > warg_strings;

        for (int i=0; i<argc; ++i)
        {
            QString arg = argv[i];
            //qDebug() << "ARG" << i << arg;

            // Command-line arg takes precedence over all environment variables.
            if (arg.startsWith("--ppn"))
            {
                QStringList parts = arg.split("=", QString::SkipEmptyParts);

                if (parts.count() > 1)
                {
                    bool ok;
                    int num = parts.last().toInt(&ok);
                    if (ok and num > 0)
                        ppn = num;
                }
            }
            else if (arg == "--include-pythonpath")
            {
                ignore_pythonpath = false;
            }
            else if (arg == "--ignore-ipython")
            {
                ignore_ipython = true;
            }
            else
            {
                warg_strings.append(arg.toStdWString());
                python_argv[python_argc] = const_cast<wchar_t*>(warg_strings.last().data());
                python_argc += 1;
            }
        }

        #ifdef Q_OS_UNIX
            signal(SIGINT, fatal_error_signal);
            signal(SIGTERM, fatal_error_signal);
        #endif // Q_OS_UNIX

        QDir site_packages( QString("%1/%2").arg( getInstallDir(), SIRE_PYTHON2_DIR ) );

        /*
        //if (not site_packages.exists())
        //    throw SireError::file_error( QObject::tr(
        //        "Cannot find the directory containing the Sire python modules (%1). "
        //        "Please check your installation of Sire in directory %2.")
        //            .arg(site_packages.absolutePath()).arg(getInstallDir()), CODELOC );
        */

        QString pythonpath;

        if (not ignore_pythonpath)
            pythonpath = qgetenv("PYTHONPATH");

        if (pythonpath.isEmpty())
            pythonpath = site_packages.canonicalPath();
        else
            pythonpath = QString("%1:%2").arg(site_packages.canonicalPath()).arg(pythonpath);

        QDir python_home( QString("%1/%2/..").arg( getInstallDir(), SIRE_BUNDLED_LIBS_DIR ) );

        /**
        if (not python_home.exists())
        throw SireError::file_error( QObject::tr(
            "Cannot find the directory containing the bundled files (%1). "
            "Please check your installations of Sire in directory %2.")
                .arg(python_home.absolutePath()).arg(getInstallDir()), CODELOC );

        qputenv("PYTHONPATH", pythonpath.toUtf8());
        qputenv("PYTHONHOME", python_home.canonicalPath().toUtf8());
        */

        //now look at the name of the executable. If there is a script with this
        //name in share/scripts then run that script
        QDir scripts_dir( QString("%1/scripts").arg(getShareDir()) );

        if (scripts_dir.exists())
        {
            QFileInfo my_script( scripts_dir, QString("%1.py").arg( QString(argv[0]).split("/").last() ) );

            if (my_script.exists())
            {
                //there is a matching script, so automatically run this script
                for (int i=python_argc; i>1; --i)
                {
                    python_argv[i] = python_argv[i-1];
                }

                warg_strings.append(my_script.absoluteFilePath().toStdWString());
                python_argv[1] = const_cast<wchar_t*>(warg_strings.last().data());
                python_argc += 1;

                //we must not now run ipython
                ignore_ipython = true;
            }
        }

        Cluster::start(ppn);
        printf("Starting %ls: number of threads equals %d\n", python_argv[0], ppn);

        // parallel implementation
        tbb::task_scheduler_init init(ppn);

        if (not ignore_ipython)
        {
            //if ipython is installed in sire.app/bundled/bin/ipython3 then automatically
            //run that as the first script. This will provide a nice environment for running
            //sire scripts
            QFileInfo ipython_file( python_home, "bin/ipython3" );

            if (ipython_file.exists())
            {
                for (int i=python_argc; i>1; --i)
                {
                    python_argv[i] = python_argv[i-1];
                }

                warg_strings.append(ipython_file.absoluteFilePath().toStdWString());
                python_argv[1] = const_cast<wchar_t*>(warg_strings.last().data());
                python_argc += 1;
            }
        }

        //name this process and thread
        SireError::setProcessString("master");
        SireError::setThreadString("main");

        // run the standard python interpreter
        status = Py_Main(python_argc, python_argv.get());
    }
    catch(const SireError::exception &e)
    {
        SireError::printError(e);
        status = -1;
    }
    catch(const std::exception &e)
    {
        SireError::printError( SireError::std_exception(e) );
        status = -1;
    }
    catch(...)
    {
        SireError::printError(SireError::unknown_exception(
                                 QObject::tr("An unknown error occurred!"), CODELOC ) );

        status = -1;
    }

    Cluster::shutdown();
    SireBase::Process::killAll();

    return status;
}
