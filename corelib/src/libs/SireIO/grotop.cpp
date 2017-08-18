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

#include "SireIO/grotop.h"

#include "SireSystem/system.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/booleanproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QRegularExpression>
#include <QFileInfo>

using namespace SireIO;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireStream;

const RegisterParser<GroTop> register_grotop;
static const RegisterMetaType<GroTop> r_grotop;

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const GroTop &grotop)
{
    writeHeader(ds, r_grotop, 1);
    
    SharedDataStream sds(ds);
    
    sds << grotop.include_path << grotop.included_files
        << static_cast<const MoleculeParser&>(grotop);
    
    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, GroTop &grotop)
{
    VersionID v = readHeader(ds, r_grotop);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
    
        sds >> grotop.include_path >> grotop.included_files
            >> static_cast<MoleculeParser&>(grotop);
    }
    else
        throw version_error(v, "1", r_grotop, CODELOC);

    return ds;
}

//first thing is to parse in the gromacs files. These use #include, #define, #if etc.
//so we need to pull all of them together into a single set of lines

/** Constructor */
GroTop::GroTop() : ConcreteProperty<GroTop,MoleculeParser>()
{}

/** This function gets the gromacs include path from the passed property map,
    as well as the current system environment */
void GroTop::getIncludePath(const PropertyMap &map)
{
    QStringList path;

    //now, see if the path is given in "GROMACS_PATH" in map
    try
    {
        const auto p = map["GROMACS_PATH"];
        
        if (p.hasValue())
        {
            path += p.value().asA<StringProperty>().toString().split(":", QString::SkipEmptyParts);
        }
        else if (p.source() != "GROMACS_PATH")
        {
            path += p.source().split(":", QString::SkipEmptyParts);
        }
    }
    catch(...)
    {}
    
    //now, see if the path is given in the "GROMACS_PATH" environment variable
    QString val = QString::fromLocal8Bit( qgetenv("GROMACS_PATH") );
    
    if (not val.isEmpty())
    {
        path += val.split(":", QString::SkipEmptyParts);
    }
    
    include_path = path;
}

/** Construct to read in the data from the file called 'filename'. The 
    passed property map can be used to pass extra parameters to control
    the parsing */
GroTop::GroTop(const QString &filename, const PropertyMap &map)
       : ConcreteProperty<GroTop,MoleculeParser>(filename,map)
{
    this->getIncludePath(map);

    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.
    
    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    //parse the data in the parse function
    this->parseLines(map);
    
    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
GroTop::GroTop(const QStringList &lines, const PropertyMap &map)
       : ConcreteProperty<GroTop,MoleculeParser>(lines,map)
{
    this->getIncludePath(map);

    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.
    
    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    //parse the data in the parse function
    this->parseLines(map);
    
    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
GroTop::GroTop(const SireSystem::System &system, const PropertyMap &map)
       : ConcreteProperty<GroTop,MoleculeParser>(map)
{
    this->getIncludePath(map);

    //look through the system and extract out the information
    //that is needed to go into the file, based on the properties
    //found in 'map'. Do this to create the set of text lines that
    //will make up the file
    
    QStringList lines;  // = code used to generate the lines
    
    //now that you have the lines, reparse them back into a GroTop object,
    //so that the information is consistent, and you have validated that
    //the lines you have written are able to be correctly read in by
    //this parser. This will also implicitly call 'assertSane()'
    GroTop parsed( lines, map );
    
    this->operator=(parsed);
}

/** Copy constructor */
GroTop::GroTop(const GroTop &other)
       : ConcreteProperty<GroTop,MoleculeParser>(other),
         include_path(other.include_path), included_files(other.included_files)
{}

/** Destructor */
GroTop::~GroTop()
{}

/** Copy assignment operator */
GroTop& GroTop::operator=(const GroTop &other)
{
    if (this != &other)
    {
        include_path = other.include_path;
        included_files = other.included_files;
        MoleculeParser::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool GroTop::operator==(const GroTop &other) const
{
    return include_path == other.include_path and
           included_files == other.included_files and
           MoleculeParser::operator==(other);
}

/** Comparison operator */
bool GroTop::operator!=(const GroTop &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* GroTop::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GroTop>() );
}

/** Return the C++ name for this class */
const char* GroTop::what() const
{
    return GroTop::typeName();
}

/** Return the list of names of directories in which to search for 
    include files. The directories are either absolute, or relative
    to the current directory. If "absolute_paths" is true then
    the full absolute paths for directories that exist on this 
    machine will be returned */
QStringList GroTop::includePath(bool absolute_paths) const
{
    if (absolute_paths)
    {
        QStringList abspaths;
        
        for (const auto path : include_path)
        {
            QFileInfo file(path);
            
            if (file.exists())
                abspaths.append( file.absoluteFilePath() );
        }
        
        return abspaths;
    }
    else
        return include_path;
}

/** Return the list of names of files that were included when reading or
    writing this file. The files are relative. If "absolute_paths"
    is true then the full absolute paths for the files will be
    used */
QStringList GroTop::includedFiles(bool absolute_paths) const
{
    if (absolute_paths)
        return included_files.values();
    else
        return included_files.keys();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr GroTop::construct(const QString &filename,
                                  const PropertyMap &map) const
{
    return GroTop(filename,map);
}

/** Return the parser that has been constructed by reading in the passed  
    text lines using the passed properties */
MoleculeParserPtr GroTop::construct(const QStringList &lines,
                                  const PropertyMap &map) const
{
    return GroTop(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr GroTop::construct(const SireSystem::System &system,
                                    const PropertyMap &map) const
{
    return GroTop(system,map);
}

/** Return a string representation of this parser */
QString GroTop::toString() const
{
    return QObject::tr("GroTop( includePath() = [%1], includedFiles() = [%2] )")
                .arg(includePath().join(", "))
                .arg(includedFiles().join(", "));
}

/** Return the format name that is used to identify this file format within Sire */
QString GroTop::formatName() const
{
    return "GroTop";
}

/** Return a description of the file format */
QString GroTop::formatDescription() const
{
    return QObject::tr("Gromacs Topology format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList GroTop::formatSuffix() const
{
    static const QStringList suffixes = { "top" };
    return suffixes;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void GroTop::assertSane() const
{
    //check state, raise SireError::program_bug if we are in an invalid state
}

/** Return whether or not the gromacs preprocessor would change these lines */
static bool gromacs_preprocess_would_change(const QVector<QString> &lines,
                                            bool use_parallel,
                                            const QHash<QString,QString> &defines)
{
    //create the regexps that are needed to find all of the
    //data that may be #define'd
    QVector<QRegularExpression> regexps;
    
    if (not defines.isEmpty())
    {
        regexps.reserve(defines.count());

        for (const auto key : defines.keys())
        {
            regexps.append( QRegularExpression( QString("\\s+%1\\s*").arg(key) ) );
        }
    }

    //function that says whether or not an individual line would change
    auto lineWillChange = [&](const QString &line)
    {
        if (line.indexOf(QLatin1String(";")) != -1 or
            line.indexOf(QLatin1String("#include")) != -1 or
            line.indexOf(QLatin1String("#ifdef")) != -1 or
            line.indexOf(QLatin1String("#else")) != -1 or
            line.indexOf(QLatin1String("#endif")) != -1 or
            line.indexOf(QLatin1String("#define")) != -1 )
        {
            return true;
        }
        else
        {
            for (int i=0; i<regexps.count(); ++i)
            {
                if (line.contains( regexps.constData()[i] ))
                    return true;
            }
            
            return false;
        }
    };
    
    const auto lines_data = lines.constData();
    
    if (use_parallel)
    {
        QMutex mutex;
    
        bool must_change = false;
    
        tbb::parallel_for( tbb::blocked_range<int>(0, lines.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            if (not must_change)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    if (lineWillChange(lines_data[i]))
                    {
                        QMutexLocker lkr(&mutex);
                        must_change = true;
                        break;
                    }
                }
            }
        });
        
        return must_change;
    }
    else
    {
        for (int i=0; i<lines.count(); ++i)
        {
            if (lineWillChange(lines_data[i]))
                return true;
        }
    }

   return false;
}

/** Return the full path to the file 'filename' searching through the
    Gromacs file path. This throws an exception if the file is not found */
QString GroTop::findIncludeFile(QString filename, QString current_dir)
{
    QString dirfile = filename;
    
    if (current_dir != ".")
    {
        dirfile = QString("%1/%2").arg(current_dir).arg(filename);
    }

    //have we seen this file before?
    if (included_files.contains(dirfile))
    {
        return included_files.value(dirfile);
    }

    //new file, so first see if this filename is absolute
    QFileInfo file(dirfile);
    
    //does this exist from the current directory?
    if (not (file.exists() and file.isReadable()))
    {
        //otherwise search the GROMACS_PATH
        for (const auto path : include_path)
        {
            file = QFileInfo( QString("%1/%2").arg(path).arg(dirfile) );
            
            if (file.exists() and file.isReadable())
            {
                file = QFileInfo(file.absoluteFilePath());
                break;
            }
        }
    }
    
    if (not (file.exists() and file.isReadable()))
    {
        //nothing was found!
        throw SireError::io_error( QObject::tr(
                "Cannot find the file '%1' using GROMACS_PATH = [ %2 ], current directory '%3'. "
                "Please make "
                "sure the file exists and is readable within your GROMACS_PATH from the "
                "current directory '%3' (e.g. "
                "set the GROMACS_PATH environment variable to include the directory "
                "that contains '%1', or copy this file into one of the existing "
                "directories [ %2 ])")
                    .arg(dirfile).arg(include_path.join(", ")).arg(current_dir), CODELOC );
    }

    //the file has been found - get the absolute file path
    QString fullpath = file.absoluteFilePath();
    
    //the file has been found. Cache the full path
    included_files.insert(dirfile,fullpath);
    
    return fullpath;
}

/** This function will use the Gromacs search path to find and load the
    passed include file. This will load the file and return the 
    un-preprocessed text. The file, together with its QFileInfo, will
    be saved in the 'included_files' hash */
QVector<QString> GroTop::loadInclude(QString filename, QString current_dir)
{
    //try to find the file
    QString absfile = findIncludeFile(filename, current_dir);
    
    //now load the file
    return MoleculeParser::readTextFile(absfile);
}

/** This function scans through a set of gromacs file lines and expands all
    macros, removes all comments and includes all #included files */
QVector<QString> GroTop::preprocess(const QVector<QString> &lines,
                                    QHash<QString,QString> defines,
                                    const QString &current_directory)
{
    //first, scan through to see if anything needs changing
    if (not gromacs_preprocess_would_change(lines, usesParallel(), defines))
    {
        //nothing to do
        return lines;
    }
    
    //Ok, we have to change the lines...
    QVector<QString> new_lines;
    new_lines.reserve(lines.count());

    //regexps used to parse the files...
    QRegularExpression include_regexp("#include\\s+['\"]([\\w\\d/\\.]+)['\"]\\s*");
    
    //loop through all of the lines...
    QVectorIterator<QString> lines_it(lines);
    
    QList<bool> ifparse;
    
    while (lines_it.hasNext())
    {
        QString line = lines_it.next();
    
        //remove any comments
        if (line.indexOf(QLatin1String(";")) != -1)
        {
            line = line.mid(0, line.indexOf(QLatin1String(";"))).simplified();
            
            //this is just an empty line, so ignore it
            if (line.isEmpty())
            {
                continue;
            }
        }
        else if (line.startsWith("*"))
        {
            //the whole line is a comment
            continue;
        }
        else
        {
            //simplify the line to remove weirdness
            line = line.simplified();
        }
        
        //now look to see if there is an #ifdef
        if (line.startsWith("#ifdef"))
        {
            //we have an ifdef - has it been defined?
            auto symbol = line.split(" ", QString::SkipEmptyParts).last();

            //push the current parse state (whether we parse if or else)
            ifparse.append( defines.value(symbol,"0") != "0" );
            continue;
        }

        //now look to see if there is an #ifndef
        if (line.startsWith("#ifndef"))
        {
            //we have an ifndef - has it been defined?
            auto symbol = line.split(" ", QString::SkipEmptyParts).last();
            
            //push the current parse state (whether we parse if or else)
            ifparse.append( defines.value(symbol,"0") == "0" );
            continue;
        }
        
        if (line == "#else")
        {
            //switch the last ifdef state
            if (ifparse.isEmpty())
                throw SireIO::parse_error( QObject::tr(
                    "Unmatched '#else' in the GROMACS file!"), CODELOC );

            ifparse.last() = not ifparse.last();
            continue;
        }
        
        if (line == "#endif")
        {
            //pop off the last 'ifdef' state
            if (ifparse.isEmpty())
                throw SireIO::parse_error( QObject::tr(
                    "Unmatched '#endif' in the GROMACS file!"), CODELOC );
        
            ifparse.removeLast();
            continue;
        }
        
        if (not ifparse.isEmpty())
        {
            //are we allowed to read this?
            if (not ifparse.last())
            {
                //no, this is blocked out
                continue;
            }
        }
        
        //now look for any #define lines
        if (line.startsWith("#define"))
        {
            auto words = line.split(" ", QString::SkipEmptyParts);
            
            if (words.count() == 1)
                throw SireIO::parse_error( QObject::tr(
                    "Malformed #define line in Gromacs file? %1").arg(line), CODELOC );
            
            if (words.count() == 2)
            {
                defines.insert(words[1], "1");
            }
            else
            {
                defines.insert(words[1], words[2]);
            }

            continue;
        }
        
        //now try to substitute any 'defines' in the line with their defined values
        for (auto it = defines.constBegin(); it != defines.constEnd(); ++it)
        {
            if (line.indexOf(it.key()) != -1)
            {
                auto words = line.split(" ", QString::SkipEmptyParts);
                
                for (int i=0; i<words.count(); ++i)
                {
                    if (words[i] == it.key())
                    {
                        words[i] = it.value();
                    }
                }
                
                line = words.join(" ");
            }
        }
        
        //now look for #include lines
        if (line.startsWith("#include"))
        {
            //now insert the contents of any included files
            auto m = include_regexp.match(line);
        
            if (not m.hasMatch())
            {
                throw SireIO::parse_error( QObject::tr(
                    "Malformed #include line in Gromacs file? %1").arg(line), CODELOC );
            }
            
            //we have to include a file
            auto filename = m.captured(1);
            
            //now find the absolute path to the file...
            auto absfile = findIncludeFile(filename, current_directory);

            //now load the file
            auto included_lines = MoleculeParser::readTextFile(absfile);
            
            //now get the absolute path to the included file
            auto parts = absfile.split("/");
            parts.removeLast();
            
            //fully preprocess these lines using the current set of defines
            included_lines = preprocess(included_lines, defines,
                                        QFileInfo(filename).path());
            
            //add these included lines to the set
            new_lines.reserve( new_lines.count() + included_lines.count() );
            new_lines += included_lines;
            continue;
        }
        
        //finally, make sure that we have not missed any '#' directives...
        if (line.startsWith("#"))
        {
            throw SireIO::parse_error( QObject::tr(
                "Unrecognised directive on Gromacs file line '%1'").arg(line), CODELOC );
        }
        
        //skip empty lines
        if (not line.isEmpty())
        {
            //otherwise this is a normal line, so append this to the set of new_lines
            new_lines.append(line);
        }
    }
    
    if (not ifparse.isEmpty())
    {
        throw SireIO::parse_error( QObject::tr(
            "Unmatched #ifdef or #ifndef in Gromacs file!"), CODELOC );
    }
    
    return new_lines;
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void GroTop::parseLines(const PropertyMap &map)
{
    //first, go through an expand any macros and include the contents of any
    //included files
    {
        QVector<QString> preprocessed_lines = preprocess(lines());
        
        if (preprocessed_lines != lines())
        {
            //the lines have changed
            setLines(preprocessed_lines);
        }
    }

    for (const auto line : lines())
    {
        qDebug() << line;
    }

    //now we know that there are no macros to expand, no other files to
    //include, and everything should be ok... ;-)

    this->setScore(0);
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map' */
System GroTop::startSystem(const PropertyMap &map) const
{
    //you should now take the data that you have parsed into the intermediate
    //format and use it to start a new System and create all molecules,
    //which are added to this System
    
    System system;
    
    
    return system;
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in 
    the passed System that are missing coordinate data. */
void GroTop::addToSystem(System &system, const PropertyMap &map) const
{
    //you should loop through each molecule in the system and work out
    //which ones are described in the file, and then add data from the file
    //to thise molecules.
}
