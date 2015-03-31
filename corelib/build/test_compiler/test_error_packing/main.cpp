
#include "SireError/errors.h"

#include <QByteArray>
#include <QDataStream>

#include <QDebug>

using namespace SireError;

int main(int argc, const char **argv)
{
    try
    {
        unknown_exception u = unknown_exception( "This is an unknown error!", CODELOC );

        QByteArray data = u.pack();

        qDebug() << "Data has been packed into " << data.count() << "bytes";   

        try
        {
            qDebug() << "Unpacking and throwing this exception...";
            exception::unpackAndThrow(data);
        }
        catch(const SireError::unknown_exception &e)
        {
            qDebug() << "THIS IS CORRECT AND EXPECTED BEHAVIOUR";
            qDebug() << "Caught unknown_error!!!" << e.toString();
            qDebug() << "THIS IS CORRECT AND EXPECTED BEHAVIOUR";
            return 0;
        }
        catch(const SireError::exception &e)
        {
            qDebug() << "Something went wrong...!";
            qDebug() << e.toString();
        }
        catch(...)
        {
            qDebug() << "SOMETHING WENT WRONG! CANNOT DETECT EXCEPTION TYPE!";
        }
    
    }
    catch(const SireError::exception &e)
    {
        qDebug() << "Something went really wrong!";
        qDebug() << e.toString();
    }
    catch(...)
    {
        qDebug() << "SOMETHING WENT WRONG - WE COULDN'T DETECT";
        qDebug() << "THE TYPE OF EXCEPTION!!!";
    }

    return -1;
}

