
#include "sire_process.h"
#include "SireError/errors.h"

using namespace SireBase;

Process::Process()
{
    throw SireError::incomplete_code( QObject::tr(
        "SireBase::Process has not yet been ported to Windows..."), CODELOC );
}

Process::Process(const Process &other)
{}

Process::~Process()
{}

Process& Process::operator=(const Process &other)
{
    return *this;
}

bool Process::operator==(const Process &other) const
{
    return true;
}

bool Process::operator!=(const Process &other) const
{
    return false;
}

void Process::wait()
{}

bool Process::wait(int ms)
{
    return true;
}

Process Process::run(const QString &command)
{
    return Process();
}

Process Process::run(const QString &command, const QString &arg)
{
    return Process();
}

Process Process::run(const QString &command, const QStringList &arguments) 
{
    return Process();
}

Process Process::run(const QString &command, const QString &stdout_file,
                     const QString &stderr_file)
{
    return Process();
}

Process Process::run(const QString &command, const QString& arg,
                     const QString &stdout_file, const QString &stderr_file)
{
    return Process();
}

Process Process::run(const QString &command, const QStringList& arguments,
                     const QString &stdout_file, const QString &stderr_file)
{
    return Process();
}

void Process::killAll()
{
    return;
}

bool Process::isRunning() 
{
    return false;
}

bool Process::hasFinished()
{
    return true;
}

bool Process::isError()
{
    return true;
}

bool Process::wasKilled()
{
    return true;
}

void Process::kill()
{
}
