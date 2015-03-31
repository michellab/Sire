
#include "SireStream/streamdata.hpp"

#include "SireMol/molecule.h"
#include "SireMol/mover.hpp"

#include "SireBase/property.h"

#include "SireMol/molecule.h"

#include "SireError/errors.h"

#include <QByteArray>
#include <QDataStream>

#include <QDebug>

using boost::tuple;
using boost::shared_ptr;

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

int main(int argc, const char **argv)
{
    try
    {
        VariantProperty v( 42 );

        qDebug() << "Saving a number...";

        QByteArray data = save(v);

        qDebug() << "Loading a number...";

        VariantProperty new_v = loadType<VariantProperty>(data);      

        qDebug() << new_v.value<int>();
        qDebug() << new_v.toInt();

        qDebug() << "Done!";
        
        Molecule mol;
        
        qDebug() << "Saving a molecule...";
        
        data = save(mol);
        
        qDebug() << "Loading a molecule...";
        
        Molecule new_mol = loadType<Molecule>(data);
        
        qDebug() << "Done!";
    }
    catch(const SireError::exception &e)
    {
        qDebug() << "ERROR!!!";
        qDebug() << e.toString();
    }

    return 0;
}

