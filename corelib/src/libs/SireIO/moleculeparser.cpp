/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#include "moleculeparser.h"

#include "SireError/errors.h"

#include "SireBase/parallel.h"

#include "SireSystem/system.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QFile>
#include <QTextStream>
#include <QElapsedTimer>
#include <QDebug>

using namespace SireIO;
using namespace SireSystem;
using namespace SireBase;
using namespace SireStream;

//////////////
////////////// Implementation of MoleculeParser
//////////////

static const RegisterMetaType<MoleculeParser> r_parser( MAGIC_ONLY, MoleculeParser::typeName() );

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const MoleculeParser &parser)
{
    writeHeader(ds, r_parser, 1);
    
    SharedDataStream sds(ds);
    sds << parser.lnes << static_cast<const Property&>(parser);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, MoleculeParser &parser)
{
    VersionID v = readHeader(ds, r_parser);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> parser.lnes >> static_cast<Property&>(parser);
    }
    else
        throw version_error(v,"1", r_parser, CODELOC);
    
    return ds;
}

/** Constructor */
MoleculeParser::MoleculeParser() : Property()
{}

/** Construct the parser, parsing in all of the lines in the file
    with passed filename */
MoleculeParser::MoleculeParser(const QString &filename)
               : Property()
{
    QFile file(filename);
    
    if (not file.open(QIODevice::ReadOnly | QIODevice::Unbuffered))
    {
        throw SireError::file_error(file, CODELOC);
    }
    
    QTextStream ts(&file);
    
    QStringList l;
    
    while (not ts.atEnd())
    {
        l.append( ts.readLine() );
    }
    
    file.close();
    
    lnes = l.toVector();
}

/** Copy constructor */
MoleculeParser::MoleculeParser(const MoleculeParser &other)
               : Property(other), lnes(other.lnes)
{}

/** Destructor */
MoleculeParser::~MoleculeParser()
{}

const char* MoleculeParser::typeName()
{
    return "SireIO::MoleculeParser";
}

MoleculeParser& MoleculeParser::operator=(const MoleculeParser &other)
{
    lnes = other.lnes;
    Property::operator=(other);
    return *this;
}

bool MoleculeParser::operator==(const MoleculeParser &other) const
{
    return lnes == other.lnes and Property::operator==(other);
}

bool MoleculeParser::operator!=(const MoleculeParser &other) const
{
    return not MoleculeParser::operator==(other);
}

/** Return whether or not this is a lead parser. The lead parser is responsible
    for starting the process of turning the parsed file into the System. There
    must be one and one-only lead parser in a set of parsers creating a System */
bool MoleculeParser::isLead() const
{
    return false;
}

/** Write the parsed data back to the file called 'filename'. This will
    overwrite the file if it exists already, so be careful! */
void MoleculeParser::write(const QString &filename) const
{
    if (not this->isTextFile())
        throw SireError::program_bug( QObject::tr(
            "Dear programmer - please override the MoleculeParser::save function "
            "to work with your binary file format. Text writing is not supported "
            "for the parser %1.").arg(this->what()), CODELOC );

    QElapsedTimer t;
    t.start();

    QFile f(filename);
    
    if (not f.open( QIODevice::WriteOnly | QIODevice::Text ))
    {
        throw SireError::file_error(f, CODELOC);
    }
    
    QTextStream ts(&f);
    
    for (const QString &line : lnes)
    {
        ts << line << '\n';
    }
    
    f.close();
    
    qint64 ns = t.nsecsElapsed();
    
    qDebug() << "File write took" << (0.000001*ns) << "ms";
}

/** Parse the passed file, returning the resulting Parser. This employs a lot
    of magic to automatically work out the format of the file and whether or
    not this is parseable by Sire... This raises an exception if the file
    cannot be recognised, or if there is an error in parsing. */
MoleculeParserPtr MoleculeParser::parse(const QString &filename,
                                        const PropertyMap &map)
{
    return MoleculeParserPtr();
}

/** Parse the passed set of files, returning the resulting Parsers */
QList<MoleculeParserPtr> MoleculeParser::parse(const QStringList &filenames,
                                               const PropertyMap &map)
{
    QList<MoleculeParserPtr> result;

    if (filenames.count() == 1)
    {
        result.append( MoleculeParser::parse(filenames[0]) );
    }
    else
    {
        QVector<MoleculeParserPtr> parsers(filenames.count());
    
        //parse the files in parallel - we use a grain size of 1
        //as each file can be pretty big, and there won't be many of them
        tbb::parallel_for( tbb::blocked_range<int>(0,filenames.count(),1),
                           [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                parsers[i] = MoleculeParser::parse(filenames[i], map);
            }
        }, tbb::simple_partitioner());
        
        result = parsers.toList();
    }
    
    return result;
}

/** Read the passed file called 'filename', returning the System contained therein */
System MoleculeParser::read(const QString &filename, const PropertyMap &map)
{
    MoleculeParserPtr parser = MoleculeParser::parse(filename, map);
    return parser.read().toSystem(map);
}

/** Read the two passed files, returning the System contained therein. The two
    files must refer to the same System, i.e. they could be a parameter + coordinate file */
System MoleculeParser::read(const QString &file1, const QString &file2,
                            const PropertyMap &map)
{
    MoleculeParserPtr parser1, parser2;
    
    tbb::parallel_invoke( [&](){ parser1 = MoleculeParser::parse(file1,map); },
                          [&](){ parser2 = MoleculeParser::parse(file2,map); } );
    
    return parser1.read().toSystem(parser2.read(),map);
}

/** Read the files with passed filenames, returning the System contained therein.
    Note that all of the files must be connected to the same system 
    (i.e. it could be the Amber Parm and Rst file) */
System MoleculeParser::read(const QStringList &filenames, const PropertyMap &map)
{
    QList<MoleculeParserPtr> parsers = MoleculeParser::parse(filenames, map);
    
    if (parsers.isEmpty())
        return System();
    
    MoleculeParserPtr parser = parsers.takeFirst();
    
    return parser.read().toSystem(parsers, map);
}

/** Return the System that is constructed from the data in this parser */
System MoleculeParser::toSystem(const PropertyMap &map) const
{
    return this->startSystem(map);
}

/** Return the System that is constructed from the data in the two
    passed parsers (i.e. representing a topology and a coordinate file) */
System MoleculeParser::toSystem(const MoleculeParser &other, const PropertyMap &map) const
{
    if (this->isLead())
    {
        if (other.isLead())
            throw SireError::io_error( QObject::tr(
                    "Cannot construct a System from two lead parsers: %1 and %2")
                        .arg(this->toString()).arg(other.toString()), CODELOC );
    
        System system = this->startSystem(map);
        other.addToSystem(system, map);
        return system;
    }
    else
    {
        if (not other.isLead())
            throw SireError::io_error( QObject::tr(
                    "Cannot construct a System when you have no lead parsers: %1 and %2")
                        .arg(this->toString()).arg(other.toString()), CODELOC );
        
        System system = other.startSystem(map);
        this->addToSystem(system,map);
        return system;
    }
}

/** Return the System that is constructed from the information in the passed
    parsers. This will parse the information in order, meaning that data contained
    in earlier parsers may be overwritten by data from later parsers */
System MoleculeParser::toSystem(const QList<MoleculeParserPtr> &others,
                                const PropertyMap &map) const
{
    if (this->isLead())
    {
        //make sure that there is not more than one lead parser
        for (auto other : others)
        {
            if (other.read().isLead())
                throw SireError::io_error( QObject::tr(
                    "Cannot construct a System when you have more than one lead parser: %1 and %2")
                        .arg(this->toString()).arg(other.read().toString()), CODELOC );
        }
        
        System system = this->startSystem(map);
        
        for (auto other : others)
        {
            other.read().addToSystem(system, map);
        }
        
        return system;
    }
    else
    {
        //find the first lead parser in the set
        int lead_idx = -1;
        
        for (int i=0; i<others.count(); ++i)
        {
            if (others[i].read().isLead())
            {
                lead_idx = i;
                break;
            }
        }
        
        if (lead_idx >= 0)
        {
            QList<MoleculeParserPtr> followers(others);
        
            MoleculeParserPtr lead = followers.takeAt(lead_idx);
            followers.push_front(*this);
            
            return lead.read().toSystem(followers, map);
        }
        else
        {
            QStringList p;
            
            p.append(this->toString());
            
            for (auto other : others)
            {
                p.append( other.read().toString() );
            }
        
            throw SireError::io_error( QObject::tr(
                "Cannot construct a System when there are no lead parsers! [ %1 ]")
                    .arg(p.join(", ")), CODELOC );
            
            return System();
        }
    }
}

/** Start creating a new System using the information contained in this parser,
    using the (optional) property map to name the properties */
System MoleculeParser::startSystem(const PropertyMap &map) const
{
    throw SireError::io_error( QObject::tr(
            "There is not enough information in this parser (%1) to start "
            "the creation of a new System. You need to use a more detailed input file.")
                .arg(this->toString()), CODELOC );
}

/** Continue adding data to the passed System using the information contained in 
    this parser, using the (optional) property map to name the properties */
void MoleculeParser::addToSystem(System &system, const PropertyMap &map) const
{
    throw SireError::io_error( QObject::tr(
            "This parser (%1) cannot be used to add additional information to a "
            "System. It can only be used to create a new System from scratch.")
                .arg(this->toString()), CODELOC );
}

Q_GLOBAL_STATIC( NullParser, nullParser )

const NullParser& MoleculeParser::null()
{
    return *(nullParser());
}

//////////////
////////////// Implementation of NullParser
//////////////

static const RegisterMetaType<NullParser> r_null;

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const NullParser &parser)
{
    writeHeader(ds, r_null, 1);
    ds << static_cast<const MoleculeParser&>(parser);
    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, NullParser &parser)
{
    VersionID v = readHeader(ds, r_null);
    
    if (v == 1)
    {
        ds >> static_cast<MoleculeParser&>(parser);
    }
    else
        throw version_error(v, "1", r_parser, CODELOC);
    
    return ds;
}

NullParser::NullParser() : ConcreteProperty<NullParser,MoleculeParser>()
{}

NullParser::NullParser(const NullParser &other)
           : ConcreteProperty<NullParser,MoleculeParser>(other)
{}

NullParser::~NullParser()
{}

NullParser& NullParser::operator=(const NullParser &other)
{
    MoleculeParser::operator=(other);
    return *this;
}

bool NullParser::operator==(const NullParser &other) const
{
    return MoleculeParser::operator==(other);
}

bool NullParser::operator!=(const NullParser &other) const
{
    return MoleculeParser::operator!=(other);
}

const char* NullParser::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullParser>() );
}

System NullParser::toSystem(const PropertyMap&) const
{
    return System();
}

System NullParser::toSystem(const MoleculeParser &other,
                            const PropertyMap&) const
{
    if (not other.isA<NullParser>())
        throw SireError::incompatible_error( QObject::tr(
                "Null parsers cannot be combined with other parsers (%1)")
                    .arg(other.toString()), CODELOC );
    
    return System();
}

System NullParser::toSystem(const QList<MoleculeParserPtr> &others,
                            const PropertyMap&) const
{
    for (auto other : others)
    {
        if (not other.isNull())
        {
            throw SireError::incompatible_error( QObject::tr(
                    "Null parsers cannot be combined with other parsers (%1)")
                        .arg(other.read().toString()), CODELOC );
        }
    }
    
    return System();
}
