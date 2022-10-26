
#ifdef SIRE_USE_MPI
    //mpich requires that mpi.h is included first
    #include <mpi.h>
#endif

#include <QFile>
#include <QFileInfo>
#include <QRegExp>
#include <QByteArray>
#include <QDateTime>

#include "SireError/errors.h"
#include "SireError/printerror.h"

#include "SireCluster/cluster.h"
#include "SireCluster/nodes.h"
#include "SireCluster/node.h"
#include "SireCluster/promise.h"

#include "SireBase/sire_process.h"

#include "SireSystem/system.h"
#include "SireMove/suprasystem.h"

#include "SireMove/moves.h"
#include "SireMove/supramoves.h"

#include "SireMove/simpacket.h"
#include "SireMove/suprasimpacket.h"

#include "SireStream/streamdata.hpp"

#include "sire_version.h"

using namespace SireCluster;
using namespace SireMove;
using namespace SireSystem;
using namespace SireStream;

#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
    // need this for older Qt
    namespace Qt
    {
        static auto endl = ::endl;
    }
#endif

void printOut(const QString &line)
{
    QTextStream stream(stdout);

    stream << line;
    Qt::endl(stream);
}

static QString repeated(const QString &s, int n)
{
    #if QT_VERSION >= QT_VERSION_CHECK(4, 5, 0)
        return s.repeated(n);
    #else
        QString r = s;

        for (int i=1; i<n; ++i)
        {
            r += s;
        }

        return r;
    #endif
}

void printBox(const QString &line, QTextStream &stream)
{
    QStringList lines = line.split("\n");

    int maxlength = 0;

    for (QStringList::iterator it = lines.begin(); it != lines.end(); ++it)
    {
        *it = it->simplified();

        if (it->length() > maxlength)
            maxlength = it->length();
    }

    const int max_maxlength = 80;

    if (maxlength > max_maxlength)
        maxlength = max_maxlength;

    QString hashline = ::repeated( "-", maxlength + 2 );

    Qt::endl(stream);
    stream << "*" << hashline << "*";
    Qt::endl(stream);

    foreach (const QString &l, lines)
    {
        if (l.length() > max_maxlength)
        {
            for (int j=0; j<l.length(); j+=max_maxlength)
            {
                stream << "| " << l.mid(j,max_maxlength).leftJustified(maxlength)
                       << " |";

                Qt::endl(stream);
            }
        }
        else
        {
            stream << "| " << l.leftJustified(maxlength) << " |";
            Qt::endl(stream);
        }
    }

    stream << "*" << hashline << "*";
    Qt::endl(stream);
    Qt::endl(stream);
}

void printBox(const QString &line)
{
    QTextStream stream(stdout);
    printBox(line, stream);
}

void printError(const QString &line)
{
    QTextStream stream(stderr);

    printBox( QObject::tr("Command line error!"), stream );

    stream << line;
    Qt::endl(stream);
    Qt::endl(stream);
}

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

    printBox( QObject::tr("You're killing me!!!") );

    // Kill any child processes
    SireBase::Process::killAll();

    // Now do the clean up actions:
    Cluster::shutdown();

    printBox( QObject::tr("I, and all of my children are now dead. Adieu...") );

    // Now reraise the signal.  We reactivate the signal's
    // default handling, which is to terminate the process.
    // We could just call exit or abort,
    // but reraising the signal sets the return status
    // from the process correctly.
    signal (sig, SIG_DFL);
    raise (sig);
}

#endif // Q_OS_UNIX

void throwIncompatibleError(const Property *p0, const Property *p1)
{
    throw SireError::incompatible_error( QObject::tr(
        "We can only perform a restart simulation "
        "with a System and Move(s) object, or with a SupraSystem and "
        "SupraMove(s) object. It cannot work with a %1 and %2.")
            .arg(p0->what()).arg(p1->what()), CODELOC );
}

/** Load up all of the passed restart file to see if it is valid */
void testLoad(const QString &filename)
{
    QFile f(filename);

    if (not f.open( QIODevice::ReadOnly) )
        throw SireError::file_error(f, CODELOC);

    QByteArray restart_data = f.readAll();

    if (restart_data.isEmpty())
        throw SireError::file_error( QObject::tr(
            "There was an error reading data from the file %1. Either "
            "the file is empty, or some read error has occured.")
                .arg(filename), CODELOC );

    //read the header
    {
        FileHeader header = SireStream::getDataHeader(restart_data);

        printBox( header.toString() );
    }

    //unpack the binary data
    QList< boost::tuple<boost::shared_ptr<void>,QString> > objects
                         = SireStream::load(restart_data);

    for (int i=0; i<objects.count(); ++i)
    {
        printOut( QObject::tr("Read object %1 of %2. It is a %3.")
                    .arg(i+1).arg(objects.count())
                    .arg(objects.at(i).get<1>()) );
    }
}

/** This reads a simulation restart file and creates a workpacket
    to run the next step in the simulation */
WorkPacket createWorkPacket(const QString &filename,
                            int nmoves, bool record_stats)
{
    //read the contents of the restart file into memory
    printOut( QObject::tr("Opening %1...").arg(filename) );
    QFile f(filename);

    if (not f.open( QIODevice::ReadOnly) )
        throw SireError::file_error(f, CODELOC);

    printOut( QObject::tr("Reading the file...") );
    QByteArray restart_data = f.readAll();

    if (restart_data.isEmpty())
        throw SireError::file_error( QObject::tr(
            "There was an error reading data from the file %1. Either "
            "the file is empty, or some read error has occured.")
                .arg(filename), CODELOC );

    //sanity check the header
    {
        printOut( QObject::tr("Reading the file header...") );
        FileHeader header = SireStream::getDataHeader(restart_data);

        printOut( QObject::tr("File Header ==\n%1").arg(header.toString()) );

        if (header.dataTypes().count() != 2)
        {
            throw SireError::incompatible_error( QObject::tr(
                "A RestartPacket can only process a restart file that contains "
                "two objects; the first should be the System or SupraSystem to "
                "be simulated, while the second should be the moves to be "
                "applied. The types available in the passed restart file "
                "are [ %1 ].").arg(header.dataTypes().join(", ")), CODELOC );
        }
    }

    printOut( QObject::tr("Unpacking the data...") );

    //unpack the binary data
    QList< boost::tuple<boost::shared_ptr<void>,QString> > objects
                         = SireStream::load(restart_data);

    printOut( QObject::tr("Number of loaded objects equals %1").arg(objects.count()) );

    if (objects.count() != 2)
        throw SireError::file_error( QObject::tr(
            "The restart file may be corrupted as despite the header claiming "
            "there were two objects, the number of objects is actually equal "
            "to %1.").arg(objects.count()), CODELOC );

    //the objects must both be derived from Property - if they are not
    //then this will cause a segfault
    Property *p0 = static_cast<Property*>(objects[0].get<0>().get());
    Property *p1 = static_cast<Property*>(objects[1].get<0>().get());

    //the first object should be derived from System or SupraSystem
    if (p0->isA<System>())
    {
        printOut( QObject::tr("Creating a SimPacket simulation...") );

        //the second object must be a 'Move' or 'Moves'
        if (p1->isA<Moves>())
        {
            return SimPacket(p0->asA<System>(), p1->asA<Moves>(),
                             nmoves, record_stats);
        }
        else if (p1->isA<Move>())
        {
            return SimPacket(p0->asA<System>(), SameMoves(p1->asA<Move>()),
                             nmoves, record_stats);
        }
        else
            ::throwIncompatibleError(p0, p1);
    }
    else if (p0->isA<SupraSystem>())
    {
        printOut( QObject::tr("Creating a SupraSimPacket simulation...") );

        //the second object must be a 'SupraMove' or 'SupraMoves'
        if (p1->isA<SupraMoves>())
        {
            return SupraSimPacket(p0->asA<SupraSystem>(), p1->asA<SupraMoves>(),
                                  nmoves, record_stats);
        }
        else if (p1->isA<SupraMove>())
        {
            return SupraSimPacket(p0->asA<SupraSystem>(),
                                  SameSupraMoves(p1->asA<SupraMove>()),
                                  nmoves, record_stats);
        }
        else
            ::throwIncompatibleError(p0, p1);
    }
    else
        ::throwIncompatibleError(p0, p1);

    return WorkPacket();
}

/** Return the name of the next restart file in the sequence. Restart
    files are name [filename]_[number].[extension] where [number] has to
    increment */
QString getNextName(QString old_name)
{
    QFileInfo filename(old_name);

    QString base_name = filename.baseName();

    QRegExp regexp("([\\w\\d]+)_(\\d+)");

    if (regexp.indexIn(base_name) != -1)
    {
        base_name = QString("%1_%2").arg(regexp.cap(1))
                                    .arg(regexp.cap(2).toInt() + 1);
    }
    else
    {
        base_name += QString("_1");
    }

    if (filename.completeSuffix().isEmpty())
    {
        return QString("%1/%2")
                    .arg(filename.path(), base_name);
    }
    else
    {
        return QString("%1/%2.%3")
                    .arg(filename.path(), base_name, filename.completeSuffix());
    }
}

/** Save the results of the simulation run in the passed WorkPacket
    to the file with name 'filename' */
void saveRestartFile(const WorkPacket &workpacket, const QString &filename)
{
    try
    {
        if (workpacket.isA<SimPacket>())
        {
            SimPacket simpacket = workpacket.asA<SimPacket>();

            SireStream::saveToFile( simpacket.system(),
                                    simpacket.moves().read(),
                                    filename );
        }
        else if (workpacket.isA<SupraSimPacket>())
        {
            SupraSimPacket simpacket = workpacket.asA<SupraSimPacket>();

            SireStream::saveToFile( simpacket.system(),
                                    simpacket.moves(),
                                    filename );
        }
        else
        {
            throw SireError::program_bug( QObject::tr(
                    "This code only works with SimPacket and SupraSimPacket "
                    "WorkPackets. It won't work with a %1.")
                        .arg(workpacket.what()), CODELOC );
        }
    }
    catch(const SireError::exception &e)
    {
        SireError::printError(e);
    }
}

struct RestartFile
{
    QString filename;
    int nmoves;
    bool record_stats;
};

void printLicense()
{
    Q_INIT_RESOURCE(sire);

    QFile f(":COPYING");

    f.open( QIODevice::ReadOnly );

    QString copying = f.readAll();

    printBox( QObject::tr("This program is licensed under the terms of the\n"
                          "GNU General Public License (GPL) Version 2 or greater.\n\n"
                          "%1" ).arg(copying) );
}

void printHelp()
{
    printOut( QObject::tr("Usage: sire\n"
                          "            (-h / --help)\n"
                          "            (-v / --version)\n"
                          "            (-T / --test)\n"
                          "            (-r number / --repeat=number)\n"
                          "            (-t hh:mm:ss / --time=hh:mm:ss)\n"
                          "            (-n number / --nmoves=number)\n"
                          "            (-s true|false / --statistics=true|false)\n"
                          "            restart_file\n"
                          "            (restart_file2 restart_file3 ...)\n\n"
                          "where:\n"
                          "  -h / --help       : Print this help text.\n"
                          "  -v / --version    : Print the version of this program.\n"
                          "  -T / --test       : Just test loading the restart files.\n"
                          "  -r / --repeat     : The number of times to repeat the "
                          "simulation.\n"
                          "  -t / --time       : The maximum run time for the simulation "
                          " (in hours, minutes and seconds)\n"
                          "  -n / --nmoves     : The number of moves to perform for each "
                          " subsequently listed restart file.\n"
                          "  -s / --statistics : Whether or not to record statistics "
                          " during simulations of the subsequently listed restart "
                          "files.\n\n"
                          "e.g. sire -n 5000 -t 2:00:00 -s true restart_file.s3\n\n"
                          "would run the simulation in \"restart_file.s3\" for a "
                          "maximum of two hours, performing 5000 moves, recording "
                          "statistics during the simulation.\n\n"
                          "sire -n 1000 restart_1.s3 -n 2000 restart_2.s3\n\n"
                          "would run 1000 moves of the system in \"restart_1.s3\" and "
                          "(possibly in parallel) 2000 moves of the system in "
                          "\"restart_2.s3\".\n") );
}

void printVersion()
{
    printOut( QObject::tr("Sire version: %1.%2.%3 (%4)")
                .arg(SIRE_VERSION_MAJOR).arg(SIRE_VERSION_MINOR)
                .arg(SIRE_VERSION_PATCH).arg(SIRE_REPOSITORY_VERSION) );
}

int readTime(const QString &time, bool &ok)
{
    QStringList parts = time.split(":");

    int times[3];

    while (parts.count() < 3)
    {
        parts.prepend("0");
    }

    for (int i=0; i<3; ++i)
    {
        times[i] = parts.at(i).toInt(&ok);

        if (times[i] < 0)
            ok = false;

        if (not ok)
            printError( QObject::tr("Could not interpret the time from \"%1\". "
                         "The correct format is hh:mm:ss, where hh is hours, "
                         "mm is minutes and ss is seconds, e.g. 2 hours, 14 minutes "
                         "and 52 seconds would be 02:14:52.").arg(time) );
    }

    return 3600 * times[0] + 60 * times[1] + times[2];
}

QString timeToString(int time)
{
    int hours = time / 3600;
    int mins = (time - hours*3600) / 60;
    int secs = time - hours*3600 - mins*60;

    QStringList parts;

    if (hours == 1)
        parts.append( QObject::tr("1 hour") );
    else if (hours > 0)
        parts.append( QObject::tr("%1 hours").arg(hours) );

    if (mins == 1)
        parts.append( QObject::tr("1 minute") );
    else if (mins > 0)
        parts.append( QObject::tr("%1 minutes").arg(mins) );

    if (secs == 1)
        parts.append( QObject::tr("1 second") );
    else if (secs > 0)
        parts.append( QObject::tr("%1 seconds").arg(secs) );

    if (parts.isEmpty())
        return QObject::tr("Unspecified time");

    else if (parts.count() == 1)
    {
        return parts.at(0);
    }
    else if (parts.count() == 2)
    {
        return QObject::tr("%1 and %2").arg(parts.at(0), parts.at(1));
    }
    else
    {
        return QObject::tr("%1, %2 and %3")
                    .arg(parts.at(0), parts.at(1), parts.at(2));
    }
}

bool readBool(const QString &flag, bool &ok)
{
    QString lower_flag = flag.toLower();

    if (lower_flag == QLatin1String("true") or
        lower_flag == QLatin1String("on") or
        lower_flag == QLatin1String("yes") or
        lower_flag == QLatin1String("1"))
    {
        return true;
    }
    else if (lower_flag == QLatin1String("false") or
             lower_flag == QLatin1String("off") or
             lower_flag == QLatin1String("no") or
             lower_flag == QLatin1String("0"))
    {
        return false;
    }
    else
    {
        printError( QObject::tr("Could not read a true|false value from \"%1\". "
                                "You should use true, on, yes or 1 to mean true, or "
                                "false, off, no or 0 to mean false (case insensitive).")
                                    .arg(flag) );

        ok = false;

        return false;
    }
}

QString boolToString(bool flag)
{
    if (flag)
        return QObject::tr("true");
    else
        return QObject::tr("false");
}

int readRepeat(const QString &repeat, bool &ok)
{
    int r = repeat.toInt(&ok);

    if (r <= 0)
        ok = false;

    if (not ok)
        printError( QObject::tr( "Could not read the number of repetitions "
                                 "from \"%1\". This should be an integer (whole number).")
                                    .arg(repeat) );

    return r;
}

int readNMoves(const QString &nmoves, bool &ok)
{
    int n = nmoves.toInt(&ok);

    if (n <= 0)
        ok = false;

    if (not ok)
        printError( QObject::tr( "Could not read the number of moves from \"%1\". "
                                 "This should be an integer (whole number).")
                                    .arg(nmoves) );

    return n;
}

QList<RestartFile> parseCommandLine(int argc, char **argv,
                                    int &repeat, int &time, bool &ok, bool &test)
{
    //  sire --repeat=10 --time=02:00:00 --nmoves=100 --statistics=true
    //                    restart_file --nmoves=150 restart_file2

    //  sire -r 10 -t 02:00:00 -n 100 -s true restart_file -n 150 restart_file2

    //  sire --help

    ok = true;

    repeat = 1;
    time = -1;

    int nmoves = 1;
    bool record_stats = true;
    test = false;

    QList<RestartFile> restart_files;

    if (argc == 1)
    {
        printHelp();
        return restart_files;
    }

    int i = 1;

    while (i < argc)
    {
        if (not ok)
            return QList<RestartFile>();

        const QString arg(argv[i]);

        if (arg.startsWith("--"))
        {
            if (arg == QLatin1String("--help"))
            {
                printHelp();
                return QList<RestartFile>();
            }
            else if (arg == QLatin1String("--license"))
            {
                printLicense();
                return QList<RestartFile>();
            }
            else if (arg == QLatin1String("--version"))
            {
                printVersion();
                return QList<RestartFile>();
            }
            else if (arg == QLatin1String("--test"))
            {
                test = true;
            }
            else if (arg.startsWith("--time"))
            {
                QStringList times = arg.split("=");

                if (times.count() != 2)
                {
                    ok = false;
                    printHelp();

                    printError( QObject::tr(
                        "Cannot interpret the time %1. The time should be set "
                        "using the argument \"--time=hh:mm:ss\" or \"-t hh:mm:ss\".")
                            .arg(arg) );
                }
                else
                {
                    time = readTime(times[1], ok);
                }
            }
            else if (arg.startsWith("--repeat"))
            {
                QStringList repeats = arg.split("=");

                if (repeats.count() != 2)
                {
                    ok = false;
                    printHelp();

                    printError( QObject::tr(
                        "Cannot interpret the number of repetitions from %1. "
                        "This should be set using the argument "
                        "\"--repeat=n\" or \"-r n\" where \"n\" is the number "
                        "of times you wish to repeat the simulation.").arg(arg) );
                }
                else
                {
                    repeat = readRepeat(repeats[1], ok);
                }
            }
            else if (arg.startsWith("--nmoves"))
            {
                QStringList moves = arg.split("=");

                if (moves.count() != 2)
                {
                    ok = false;
                    printHelp();

                    printError( QObject::tr(
                        "Cannot interpret the number of moves from %1. "
                        "This should be set using the argument "
                        "\"--nmoves=N\" or \"-n N\" where \"N\" is the number "
                        "of moves to perform.").arg(arg) );
                }
                else
                {
                    nmoves = readNMoves(moves[1], ok);
                }
            }
            else if (arg.startsWith("--statistics"))
            {
                QStringList stats = arg.split("=");

                if (stats.count() != 2)
                {
                    ok = false;
                    printHelp();

                    printError( QObject::tr(
                        "Cannot interpret whether or not to record statistics "
                        "from %1. This should be set using the argument "
                        "\"--statistics=true|false\" or \"-s true|false\".")
                            .arg(arg) );
                }
                else
                {
                    record_stats = readBool(stats[1], ok);
                }
            }
            else
            {
                ok = false;
                printHelp();

                printError( QObject::tr("Unrecognised argument \"%1\".").arg(arg) );
            }
        }
        else if (arg.startsWith("-"))
        {
            if (arg == QLatin1String("-h"))
            {
                printHelp();
                return QList<RestartFile>();
            }
            else if (arg == QLatin1String("-l"))
            {
                printLicense();
                return QList<RestartFile>();
            }
            else if (arg == QLatin1String("-v"))
            {
                printVersion();
                return QList<RestartFile>();
            }
            else if (arg == QLatin1String("-T"))
            {
                test = true;
            }
            else if (arg == QLatin1String("-t"))
            {
                ++i;

                if (i == argc)
                {
                    ok = false;
                    printHelp();

                    printError( QObject::tr(
                        "Cannot interpret the time. The time should be set using "
                        "the argument \"--time=hh:nm:ss\" or \"-t hh:mm:ss\".") );

                    return QList<RestartFile>();
                }

                time = readTime(argv[i], ok );

                if (not ok)
                    return QList<RestartFile>();
            }
            else if (arg == QLatin1String("-r"))
            {
                ++i;

                if (i == argc)
                {
                    ok = false;
                    printHelp();

                    printError( QObject::tr(
                        "Cannot interpret the number of repetitions. "
                        "This should be set using the argument "
                        "\"--repeat=n\" or \"-r n\" where \"n\" is the number "
                        "of times you wish to repeat the simulation.") );
                }
                else
                {
                    repeat = readRepeat(argv[i], ok);
                }
            }
            else if (arg == QLatin1String("-n"))
            {
                ++i;

                if (i == argc)
                {
                    ok = false;
                    printHelp();

                    printError( QObject::tr(
                        "Cannot interpret the number of moves from %1. "
                        "This should be set using the argument "
                        "\"--nmoves=N\" or \"-n N\" where \"N\" is the number "
                        "of moves to perform.") );
                }
                else
                {
                    nmoves = readNMoves(argv[i], ok);
                }
            }
            else if (arg == QLatin1String("-s"))
            {
                ++i;

                if (i == argc)
                {
                    ok = false;
                    printHelp();

                    printError( QObject::tr(
                        "Cannot interpret whether or not to record statistics. "
                        "This should be set using the argument "
                        "\"--statistics=true|false\" or \"-s true|false\".") );
                }
                else
                {
                    record_stats = readBool(argv[i], ok);
                }
            }
            else
            {
                ok = false;
                printHelp();

                printError( QObject::tr("Unrecognised argument \"1\".").arg(arg) );
            }
        }
        else
        {
            //this must be a restart file name - add it to the list
            //with the current number of moves and state of recording statistics
            RestartFile restart_file;

            restart_file.filename = arg;
            restart_file.nmoves = nmoves;
            restart_file.record_stats = record_stats;

            restart_files.append(restart_file);
        }

        ++i;
    }

    if (ok and restart_files.isEmpty())
        printHelp();

    return restart_files;
}

int main(int argc, char **argv)
{
    QString hostname, username;

    #ifdef Q_OS_UNIX
        signal(SIGINT, fatal_error_signal);
        signal(SIGTERM, fatal_error_signal);

        username = std::getenv("USER");

        char buffer[128];
        gethostname(buffer, 128);
        hostname = buffer;

    #else
        username = "unknown";
        hostname = username;
    #endif // Q_OS_UNIX

    int status = 0;
    const int ppn = 1;  //the number of processes per node

    try
    {
        #ifdef SIRE_USE_MPI
            //start MPI - ABSOLUTELY must use multi-threaded MPI
            int level;
            MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &level);
        #endif

        //are we the first node in the cluster?
        if (Cluster::getRank() == 0)
        {
            //we need to know when the program started so that we know
            //how much time left we have in our allocation
            const QDateTime program_start = QDateTime::currentDateTime();

            printBox( QObject::tr(
                "Sire Simulation Runner, Copyright (C) 2009 Christopher Woods. \n\n"
                "This program comes with ABSOLUTELY NO WARRANTY.\n\n"
                "This is free software and you are welcome to redistribute it\n"
                "under certain conditions; type \"sire -l\" or \"sire --license\"\n"
                "for warranty and licensing conditions.\n\n"
                "For more information and to contact the authors please visit\n\n"
                "http://siremol.org") );

            printBox( QObject::tr(
                    "%4@%5: Starting primary node (%1 of %2): nThreads()=%3")
                       .arg(Cluster::getRank()).arg(Cluster::getCount())
                       .arg(ppn)
                       .arg(username,hostname) );

            //name this process and thread
            SireError::setProcessString("primary");
            SireError::setThreadString("main");

            //start the cluster - on the primary we need one extra
            //thread for the Python interpreter
            #ifdef SIRE_USE_MPI
                MPI_Barrier( MPI_COMM_WORLD );
            #endif

            Cluster::start(ppn);

            #ifdef SIRE_USE_MPI
                MPI_Barrier( MPI_COMM_WORLD );
            #endif

            int time, repeat;
            bool ok, test_load;

            QList<RestartFile> restart_files = parseCommandLine(argc, argv,
                                                                repeat, time,
                                                                ok, test_load);

            int nrestarts = restart_files.count();

            if (not ok)
            {
                status = -1;
            }
            else if (test_load and nrestarts > 0)
            {
                //just test loading all of the restart files
                for (int i=0; i<nrestarts; ++i)
                {
                    printOut( QObject::tr("Loading restart file %1...")
                                .arg(restart_files.at(i).filename) );

                    try
                    {
                        testLoad(restart_files.at(i).filename);
                    }
                    catch(const SireError::exception &e)
                    {
                        printError(e);
                    }
                }
            }
            else if (nrestarts > 0)
            {
                printOut( QObject::tr("Number of repetitions equals %1.")
                            .arg(repeat) );

                if (time > 0)
                    printOut( QObject::tr("Total allowed runtime equals %1.")
                                  .arg(timeToString(time)) );

                printBox( QObject::tr("Restart files to run:") );

                for (int i=0; i<nrestarts; ++i)
                {
                    const RestartFile &r = restart_files.at(i);

                    printOut( QObject::tr(
                                "%1 : %2 ( nmoves = %3, recording statistics = %4 )")
                                .arg(i+1).arg(r.filename).arg(r.nmoves)
                                .arg(boolToString(r.record_stats)) );
                }

                printOut("\n");

                QVector<int> error_count(nrestarts, 0);

                for (int ir=0; ir<repeat; ++ir)
                {
                    const QDateTime time_start = QDateTime::currentDateTime();

                    if (repeat > 1)
                        printOut( QObject::tr("\nRunning repetition %1 of %2...")
                                    .arg(ir+1).arg(repeat) );

                    printOut( QObject::tr("Running %1 simulation(s)...").arg(nrestarts) );

                    Nodes nodes = Cluster::getNodes(nrestarts);
                    printOut( QObject::tr("Number of nodes in the cluster is %1")
                                    .arg(nodes.count()) );

                    ThisThread this_thread = nodes.borrowThisThread();

                    QList<Promise> promises;
                    QList<int> restart_idxs;
                    QStringList running_files;

                    int nskipped = 0;

                    printOut( QObject::tr("About to loop over simulations to submit") );

                    //submit all of the simulations
                    for (int i=0; i<nrestarts; ++i)
                    {
                        const RestartFile &r = restart_files.at(i);

                        if (error_count[i] > 5)
                        {
                            printOut( QObject::tr(
                                "Skipping %1 as it has failed five times in a row!")
                                    .arg(r.filename) );

                            ++nskipped;

                            continue;
                        }

                        try
                        {
                            //create a workpacket for this simulation
                            printOut( QObject::tr("Creating workpacket %1...")
                                                .arg(i+1) );
                            WorkPacket workpacket = createWorkPacket(r.filename,
                                                                     r.nmoves,
                                                                     r.record_stats);

                            printOut( QObject::tr("Running simulation %1 of %2: %3")
                                        .arg(i+1).arg(nrestarts).arg(r.filename) );

                            Node node = nodes.getNode();

                            printOut( QObject::tr("Starting the job...") );
                            promises.append( node.startJob(workpacket) );
                            running_files.append( r.filename );
                            restart_idxs.append( i );
                        }
                        catch (const SireError::exception &e)
                        {
                            printOut( QObject::tr("There was a problem when reading %1.")
                                        .arg(r.filename) );

                            SireError::printError(e);

                            error_count[i] += 1;
                        }
                    }

                    if (nskipped == nrestarts)
                    {
                        printBox( QObject::tr(
                                    "All restart files have finished in error!") );
                        break;
                    }

                    printOut( QObject::tr(
                                    "Waiting for all submitted jobs to finish...") );

                    //wait for them all to finish
                    for (int i=0; i < promises.count(); ++i)
                    {
                        promises[i].wait();
                    }

                    printOut( QObject::tr("All submitted jobs have finished") );

                    //save the output for all processes that completed successfully
                    for (int i=0; i < promises.count(); ++i)
                    {
                        if (not promises[i].isError())
                        {
                            printOut( QObject::tr(
                                        "The simulation connected with %1 has finished.")
                                            .arg(running_files[i]) );

                            //get the name of the new binary restart file
                            QString newfile = getNextName(running_files[i]);

                            //save the files to the new binary restart file
                            saveRestartFile(promises[i].result(), newfile);

                            //update the name of this restart file (for repetitions)
                            restart_files[ restart_idxs[i] ].filename = newfile;

                            error_count[ restart_idxs[i] ] = 0;
                        }
                    }

                    //did any simulations finish in error?
                    for (int i=0; i < promises.count(); ++i)
                    {
                        if (promises[i].isError())
                        {
                            printOut( QObject::tr("There was a problem when running %1.")
                                        .arg(running_files[i]) );

                            error_count[ restart_idxs[i] ] += 1;

                            try
                            {
                                promises[i].throwError();
                            }
                            catch(const SireError::exception &e)
                            {
                                SireError::printError(e);
                            }
                        }
                    }

                    const QDateTime time_end = QDateTime::currentDateTime();

                    int work_time = time_start.secsTo(time_end);
                    int elapsed_time = program_start.secsTo(time_end);

                    printOut( QObject::tr("Finished repetition %1 of %2. This "
                          "took %3.")
                            .arg(ir+1).arg(repeat)
                            .arg(timeToString(work_time)) );

                    if (ir+1 != repeat and time > 0)
                    {
                        printOut( QObject::tr("Remaining time: %1.")
                                    .arg(timeToString(time - elapsed_time)) );

                        if (elapsed_time + work_time + 600 >= time)
                        {
                            printOut( QObject::tr("There is not sufficient remaning "
                                "time to complete another repetition. Exiting...") );

                            break;
                        }
                    }
                } // end of loop over repetitions
            }
        }
        else
        {
            //this is one of the compute nodes...
            printBox( QObject::tr(
                        "%4@%5: Starting one of the compute nodes (%1 of %2): nThreads()=%3")
                            .arg(Cluster::getRank()).arg(Cluster::getCount())
                            .arg(ppn)
                            .arg(username,hostname) );

            //name this process
            SireError::setProcessString( QString("compute%1").arg(Cluster::getRank()) );
            SireError::setThreadString( "main" );

            //exec the Cluster - this starts the cluster and then
            //blocks while it is running
            printOut( QObject::tr("compute%1 waiting to start...")
                                .arg(Cluster::getRank()) );

            #ifdef SIRE_USE_MPI
                MPI_Barrier( MPI_COMM_WORLD );
            #endif

            printOut( QObject::tr("compute%1 starting...").arg(Cluster::getRank()) );

            Cluster::start(ppn);

            printOut( QObject::tr("compute%1 waiting to wait...")
                                .arg(Cluster::getRank()) );

            #ifdef SIRE_USE_MPI
                MPI_Barrier( MPI_COMM_WORLD );
            #endif

            printOut( QObject::tr("compute%1 waiting...").arg(Cluster::getRank()) );

            Cluster::wait();
            status = 0;

            printOut( QObject::tr("compute%1 finished").arg(Cluster::getRank()) );
        }
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

    //shutdown the cluster
    #ifdef SIRE_USE_MPI
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (rank == 0)
        {
            printOut( QObject::tr("Shutting down the cluster...") );
            Cluster::shutdown();
        }

        //wait for all of the MPI jobs to finish
        MPI_Barrier( MPI_COMM_WORLD );

        if (rank == 0)
        {
            printOut( QObject::tr("The entire cluster has now shutdown.") );
        }

        MPI_Finalize();
    #else
        printOut( QObject::tr("Shutting down the cluster...") );
        Cluster::shutdown();
        printOut( QObject::tr("The entire cluster has now shutdown.") );
    #endif

    return status;
}
