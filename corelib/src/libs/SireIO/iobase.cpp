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

#include "iobase.h"

#include "SireMol/molecule.h"
#include "SireMol/molidx.h"
#include "SireMol/mover.hpp"
#include "SireMol/cuttingfunction.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include <QFile>
#include <QDebug>

using namespace SireIO;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

////////////
//////////// Implementation of IOParameterBase
////////////

IOParametersBase::IOParametersBase()
{}

IOParametersBase::~IOParametersBase()
{}

PropertyName IOParametersBase::coords_property("coordinates");
PropertyName IOParametersBase::elements_property("element");

PropertyName IOParametersBase::cutting_function( "cutting-function",
                                                 CuttingFunction::null() );

////////////
//////////// Implementation of IOBase
////////////

/** Constructor */
IOBase::IOBase() : Property()
{}

/** Copy constructor */
IOBase::IOBase(const IOBase &other) : Property(other)
{}

/** Destructor */
IOBase::~IOBase()
{}

/** Read all of the molecules contained in the data 'data', using 
    the (optional) passed properties in 'map', and returning a 
    molecule group containing the molecules in the same order as 
    they appear in the file */
MoleculeGroup IOBase::read(const QByteArray &data, const PropertyMap &map) const
{
    if (data.isEmpty())
        return MoleculeGroup();

    return this->readMols(data, map);
}

/** Read all of the molecules contained in the file 'filename', using
    the (optional) passed properties in 'map', and returning a molecule
    group containing the molecules in the same order as they appear
    in the file, and with the molecule group name being 'filename' */
MoleculeGroup IOBase::read(const QString &filename, const PropertyMap &map) const
{
    //open the file for reading
    QFile fle(filename); 
    
    if (not fle.exists())
        throw SireError::file_error( QObject::tr(
            "Cannot find the file \"%1\".").arg(filename), CODELOC );
    
    if (not fle.open(QIODevice::ReadOnly|QIODevice::Text))
    {
        throw SireError::file_error(fle, CODELOC);
    }

    //get all of the data out of the file
    QByteArray contents = fle.readAll();
    fle.close();

    if (contents.isEmpty())
        return MoleculeGroup();

    MoleculeGroup molecules = this->readMols(contents, map);
    
    molecules.setName(filename);

    return molecules;
}

/** Simple overload of IOBase::read(QString) */
MoleculeGroup IOBase::read(const char *filename, const PropertyMap &map) const
{
    return this->read(QLatin1String(filename), map);
}

/** Read all of the molecules contained in the IO device 'dev', using
    the (optional) passed properties in 'map', and returning a molecule
    group containing the molecules in the same order as they appear
    on the device, and with the molecule group name being taken
    from the device */
MoleculeGroup IOBase::read(QIODevice &dev, const PropertyMap &map) const
{
    if (not dev.isReadable())
    {
        throw SireError::io_error(QObject::tr("Cannot read molecules from an unreadble "
                                              "device!"), CODELOC );
    }

    //get all of the data from the device
    QByteArray contents = dev.readAll();

    if (contents.isEmpty())
        return MoleculeGroup();

    MoleculeGroup molecules = this->readMols(contents, map);

    //maybe one day try to name the group...

    return molecules;
}

/** Read a single molecule from the passed data - this returns only
    the first molecule from the data */
Molecule IOBase::readMolecule(const QByteArray &data, const PropertyMap &map) const
{
    MoleculeGroup molecules = this->read(data, map);
    return molecules.at(MolIdx(0)).molecule();
}

/** Read a single molecule from the passed file - this returns only
    the first molecule from the file */
Molecule IOBase::readMolecule(const QString &filename, const PropertyMap &map) const
{
    MoleculeGroup molecules = this->read(filename, map);
    return molecules.at(MolIdx(0)).molecule();
}

/** Simple overload designed to prevent confusion with QByteArray function */
Molecule IOBase::readMolecule(const char *filename, const PropertyMap &map) const
{
    return this->readMolecule( QLatin1String(filename), map );
}

/** Read a single molecule from the passed IO device - this returns
    only the first molecule from the device */
Molecule IOBase::readMolecule(QIODevice &dev, const PropertyMap &map) const
{
    MoleculeGroup molecules = this->read(dev, map);
    return molecules.at(MolIdx(0)).molecule();
}

/** Write the molecules in the passed group to memory and return
    that data.
    
    This writes the molecules in the same order as they appear in the
    passed group. 
*/
QByteArray IOBase::write(const MoleculeGroup &molgroup, const PropertyMap &map) const
{
    return this->writeMols(molgroup, map);
}

/** Write the molecules in the passed group to the file called 'filename'.
    This writes the molecules in the same order as they appear in the
    passed group. */
void IOBase::write(const MoleculeGroup &molgroup, const QString &filename,
                   const PropertyMap &map) const
{
    //write the group to a binary blob
    QByteArray data = this->writeMols(molgroup, map);

    //open a file into which to write the molecules
    QFile fle(filename);
    
    if (not fle.open(QIODevice::WriteOnly|QIODevice::Text))
    {
        throw SireError::file_error(fle);
    }
    
    //dump the blob into the file
    if (fle.write(data) == -1)
    {
        throw SireError::file_error(fle);
    }
}

/** Write the molecules to memory, which is returned */
QByteArray IOBase::write(const Molecules &molecules, const PropertyMap &map) const
{
    return this->writeMols(molecules, map);
}

/** Write the molecules in the passed group to the file called 'filename'. */
void IOBase::write(const Molecules &molecules, const QString &filename,
                   const PropertyMap &map) const
{
    //write the group to a binary blob
    QByteArray data = this->writeMols(molecules, map);

    //open a file into which to write the molecules
    QFile fle(filename);
    
    if (not fle.open(QIODevice::WriteOnly|QIODevice::Text))
    {
        throw SireError::file_error(fle);
    }
    
    //dump the blob into the file
    if (fle.write(data) == -1)
    {
        throw SireError::file_error(fle);
    }
}

/** Write the passed molecule to memory, which is returned */
QByteArray IOBase::write(const MoleculeView &molecule, const PropertyMap &map) const
{
    return this->write( Molecules(molecule), map );
}

/** Write the passed molecule to the file called 'filename' */
void IOBase::write(const MoleculeView &molecule, const QString &filename,
                   const PropertyMap &map) const
{
    this->write( Molecules(molecule), filename, map );
}

/** Write the molecules in the passed group to the IO device 'dev'.
    This writes the molecules in the same order as they appear in the
    passed group. */
void IOBase::write(const MoleculeGroup &molgroup, QIODevice &dev,
                   const PropertyMap &map) const
{
    //write the group to a binary blob
    QByteArray data = this->writeMols(molgroup, map);

    //dump the blob to the device
    if (dev.write(data) == -1)
    {
        throw SireError::file_error(dev.errorString(), CODELOC);
    }
}

/** Write the molecules in the passed group to the IO device 'dev'. */
void IOBase::write(const Molecules &molecules, QIODevice &dev,
                   const PropertyMap &map) const
{
    //write the group to a binary blob
    QByteArray data = this->writeMols(molecules, map);

    //dump the blob into the file
    if (dev.write(data) == -1)
    {
        throw SireError::file_error(dev.errorString(), CODELOC);
    }
}

/** Write the passed molecule to the IO device 'dev' */
void IOBase::write(const MoleculeView &molecule, QIODevice &dev,
                   const PropertyMap &map) const
{
    this->write( Molecules(molecule), dev, map );
}

Q_GLOBAL_STATIC( NullIO, globalNullIOBase );

/** Return the global null IOBase object (a PDB writer) */
NullIO IOBase::null()
{
    return *(globalNullIOBase());
}

////////////
//////////// Implementation of NullIO
////////////

static const RegisterMetaType<NullIO> r_nullio;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const NullIO &nullio)
{
    writeHeader(ds, r_nullio, 1);
    
    ds << static_cast<const IOBase&>(nullio);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, NullIO &nullio)
{
    VersionID v = readHeader(ds, r_nullio);
    
    if (v == 1)
    {
        ds >> static_cast<IOBase&>(nullio);
    }
    else
        throw version_error(v, "1", r_nullio, CODELOC);
        
    return ds;
}

/** Constructor */
NullIO::NullIO() : ConcreteProperty<NullIO,IOBase>()
{}

/** Copy constructor */
NullIO::NullIO(const NullIO &other) : ConcreteProperty<NullIO,IOBase>(other)
{}

/** Destructor */
NullIO::~NullIO()
{}

/** Copy assignment operator */
NullIO& NullIO::operator=(const NullIO &other)
{
    IOBase::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullIO::operator==(const NullIO&) const
{
    return true;
}

/** Comparison operator */
bool NullIO::operator!=(const NullIO&) const
{
    return false;
}

/** Nothing will be read */
MoleculeGroup NullIO::readMols(const QByteArray&, const PropertyMap&) const
{
    return MoleculeGroup();
}

/** Nothing will be written */
QByteArray NullIO::writeMols(const MoleculeGroup&, const PropertyMap&) const
{
    return QByteArray();
}

/** Nothing will be written */
QByteArray NullIO::writeMols(const Molecules&, const PropertyMap&) const
{
    return QByteArray();
}

const char* NullIO::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullIO>() );
}
