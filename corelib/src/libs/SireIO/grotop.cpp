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

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QRegularExpression>
#include <QFileInfo>
#include <QDir>

using namespace SireIO;
using namespace SireUnits;
using namespace SireMol;
using namespace SireMM;
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
        << grotop.expanded_lines
        << grotop.atom_types
        << grotop.nb_func_type
        << grotop.combining_rule << grotop.fudge_lj
        << grotop.fudge_qq << grotop.generate_pairs
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
            >> grotop.expanded_lines
            >> grotop.atom_types
            >> grotop.nb_func_type
            >> grotop.combining_rule >> grotop.fudge_lj
            >> grotop.fudge_qq >> grotop.generate_pairs
            >> static_cast<MoleculeParser&>(grotop);
    }
    else
        throw version_error(v, "1", r_grotop, CODELOC);

    return ds;
}

//first thing is to parse in the gromacs files. These use #include, #define, #if etc.
//so we need to pull all of them together into a single set of lines

/** Constructor */
GroTop::GroTop()
       : ConcreteProperty<GroTop,MoleculeParser>(),
         nb_func_type(0), combining_rule(0), fudge_lj(0), fudge_qq(0),
         generate_pairs(false)
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
    
    //now go through each path and convert it into an absolute path based on the
    //current directory
    for (const auto p : path)
    {
        include_path.append( QFileInfo(p).canonicalFilePath() );
    }
}

/** Construct to read in the data from the file called 'filename'. The 
    passed property map can be used to pass extra parameters to control
    the parsing */
GroTop::GroTop(const QString &filename, const PropertyMap &map)
       : ConcreteProperty<GroTop,MoleculeParser>(filename,map),
         nb_func_type(0), combining_rule(0), fudge_lj(0), fudge_qq(0),
         generate_pairs(false)
{
    this->getIncludePath(map);

    //parse the data in the parse function, passing in the absolute path
    //to the directory that contains this file
    this->parseLines( QFileInfo(filename).absolutePath(), map);
    
    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
GroTop::GroTop(const QStringList &lines, const PropertyMap &map)
       : ConcreteProperty<GroTop,MoleculeParser>(lines,map),
         nb_func_type(0), combining_rule(0), fudge_lj(0), fudge_qq(0),
         generate_pairs(false)
{
    this->getIncludePath(map);

    //parse the data in the parse function, assuming the file has
    //come from the current directory
    this->parseLines(QDir::current().absolutePath(), map);
    
    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
GroTop::GroTop(const SireSystem::System &system, const PropertyMap &map)
       : ConcreteProperty<GroTop,MoleculeParser>(map),
         nb_func_type(0), combining_rule(0), fudge_lj(0), fudge_qq(0),
         generate_pairs(false)
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
         include_path(other.include_path), included_files(other.included_files),
         expanded_lines(other.expanded_lines),
         atom_types(other.atom_types),
         nb_func_type(other.nb_func_type), combining_rule(other.combining_rule),
         fudge_lj(other.fudge_lj), fudge_qq(other.fudge_qq),
         generate_pairs(other.generate_pairs)
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
        expanded_lines = other.expanded_lines;
        atom_types = other.atom_types;
        nb_func_type = other.nb_func_type;
        combining_rule = other.combining_rule;
        fudge_lj = other.fudge_lj;
        fudge_qq = other.fudge_qq;
        generate_pairs = other.generate_pairs;
        MoleculeParser::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool GroTop::operator==(const GroTop &other) const
{
    return include_path == other.include_path and
           included_files == other.included_files and
           expanded_lines == other.expanded_lines and
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
    //first, go through the list of included files
    QStringList files;
    
    for (auto it = included_files.constBegin(); it != included_files.constEnd(); ++it)
    {
        files += it.value();
    }

    if (absolute_paths)
    {
        //these are already absolute filenames
        return files;
    }
    else
    {
        //subtract any paths that relate to the current directory or GROMACS_PATH
        QString curpath = QDir::current().absolutePath();
        
        for (auto it = files.begin(); it != files.end(); ++it)
        {
            if ( it->startsWith(curpath) )
            {
                *it = it->mid(curpath.length()+1);
            }
            else
            {
                for (const auto path : include_path)
                {
                    if (it->startsWith(path))
                    {
                        *it = it->mid(path.length()+1);
                    }
                }
            }
        }
        
        return files;
    }
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

/** Return the atom type data for the passed atom type. This returns
    null data if it is not present */
GromacsAtomType GroTop::atomType(const QString &atm) const
{
    return atom_types.value(atm, GromacsAtomType());
}

/** Return the ID string for the bond atom types 'atm0' 'atm1'. This
    creates the string 'atm0;atm1' or 'atm1;atm0' depending on which
    of the atoms is lower. The ';' character is used as a separator
    as it cannot be in the atom names, as it is used as a comment 
    character in the Gromacs Top file */
static QString get_bond_id(const QString &atm0, const QString &atm1)
{
    if (atm0 < atm1)
    {
        return QString("%1;%2").arg(atm0,atm1);
    }
    else
    {
        return QString("%1;%2").arg(atm1,atm0);
    }
}

/** Return the bond potential data for the passed pair of atoms */
GromacsBond GroTop::bond(const QString &atm0, const QString &atm1) const
{
    return bond_potentials.value( get_bond_id(atm0,atm1), GromacsBond() );
}

/** Return the atom types loaded from this file */
QHash<QString,GromacsAtomType> GroTop::atomTypes() const
{
    return atom_types;
}

/** Return the bond potentials loaded from this file */
QHash<QString,SireMM::GromacsBond> GroTop::bonds() const
{
    return bond_potentials;
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
            
            if (line.trimmed().endsWith("\\"))
            {
                //this is a continuation line
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
    //new file, so first see if this filename is absolute
    QFileInfo file(filename);

    //is the filename absolute?
    if (file.isAbsolute())
    {
        if (not (file.exists() and file.isReadable()))
        {
            throw SireError::io_error( QObject::tr(
                "Cannot find the file '%1'. Please make sure that this file exists "
                "and is readable").arg(filename), CODELOC );
        }
        
        return filename;
    }
    
    //does this exist from the current directory?
    file = QFileInfo( QString("%1/%2").arg(current_dir).arg(filename) );
    
    if (file.exists() and file.isReadable())
        return file.absoluteFilePath();
    
    //otherwise search the GROMACS_PATH
    for (const auto path : include_path)
    {
        file = QFileInfo( QString("%1/%2").arg(path).arg(filename) );
            
        if (file.exists() and file.isReadable())
        {
            return file.absoluteFilePath();
        }
    }

    //nothing was found!
    throw SireError::io_error( QObject::tr(
            "Cannot find the file '%1' using GROMACS_PATH = [ %2 ], current directory '%3'. "
            "Please make "
            "sure the file exists and is readable within your GROMACS_PATH from the "
            "current directory '%3' (e.g. "
            "set the GROMACS_PATH environment variable to include the directory "
            "that contains '%1', or copy this file into one of the existing "
            "directories [ %2 ])")
                .arg(filename).arg(include_path.join(", ")).arg(current_dir), CODELOC );

    return QString();
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
                                    const QString &current_directory,
                                    const QString &parent_file)
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
        
        //now look to see if the line should be joined to the next line
        while (line.endsWith("\\"))
        {
            if (not lines_it.hasNext())
            {
                throw SireIO::parse_error( QObject::tr(
                    "Continuation line on the last line of the Gromacs file! '%1'")
                        .arg(line), CODELOC );
            }
            
            line += lines_it.next();
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
                                        parts.join("/"), absfile);
            
            //add these included lines to the set
            new_lines.reserve( new_lines.count() + included_lines.count() );
            new_lines += included_lines;

            //finally, record that this file depends on the included file
            included_files[parent_file].append(absfile);

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

/** Return the non-bonded function type for the molecules in this file */
int GroTop::nonBondedFunctionType() const
{
    return nb_func_type;
}

/** Return the combining rules to use for the molecules in this file */
int GroTop::combiningRules() const
{
    return combining_rule;
}

/** Return the Lennard Jones fudge factor for the molecules in this file */
double GroTop::fudgeLJ() const
{
    return fudge_lj;
}

/** Return the electrostatic fudge factor for the molecules in this file */
double GroTop::fudgeQQ() const
{
    return fudge_qq;
}

/** Return whether or not the non-bonded pairs should be automatically generated
    for the molecules in this file */
bool GroTop::generateNonBondedPairs() const
{
    return generate_pairs;
}

/** Return the expanded set of lines (after preprocessing) */
const QVector<QString>& GroTop::expandedLines() const
{
    return expanded_lines;
}

/** Internal function to return a LJParameter from the passed W and V values
    for the passed Gromacs combining rule */
static LJParameter toLJParameter(double v, double w, int rule)
{
    if (rule == 2 or rule == 3)
    {
        // v = sigma in nm, and w = epsilon in kJ mol-1
        return LJParameter::fromSigmaAndEpsilon( v * nanometer, w * kJ_per_mol );
    }
    else
    {
        // v = 4 epsilon sigma^6 in kJ mol-1 nm^6, w = 4 epsilon sigma^12 in kJ mol-1 nm^12
        // so sigma = (w/v)^1/6 and epsilon = v^2 / 4w
        return LJParameter::fromSigmaAndEpsilon( std::pow(w/v, 1.0/6.0) * nanometer,
                                                 (v*v / (4.0*w)) * kJ_per_mol );
    }
}

/** Internal function to convert a LJParameter to V and W based on the passed gromacs 
    combining rule */
static std::tuple<double,double> fromLJParameter(const LJParameter &lj, int rule)
{
    const double sigma = lj.sigma().to(nanometer);
    const double epsilon = lj.epsilon().to(kJ_per_mol);

    if (rule == 2 or rule == 3)
    {
        return std::make_tuple(sigma, epsilon);
    }
    else
    {
        double sig6 = SireMaths::pow(sigma, 6);
        double v = 4.0 * epsilon * sig6;
        double w = v * sig6;
        
        return std::make_tuple(v, w);
    }
}

/** Internal function, called by ::interpret() that processes all of the data
    from all of the directives, returning a set of warnings */
QStringList GroTop::processDirectives(const QMap<int,QString> &taglocs,
                                      const QHash<QString,int> &ntags)
{
    //internal function that returns the lines associated with the
    //specified directive
    auto getLines = [&](const QString &directive, int n) -> QStringList
    {
        if (n >= ntags.value(directive,0))
        {
            return QStringList();
        }
        
        int start = 0;
        int end = expandedLines().count();
        
        //find the tag
        for (auto it = taglocs.constBegin(); it != taglocs.constEnd(); ++it)
        {
            if (it.value() == directive)
            {
                if (n == 0)
                {
                    start = it.key()+1;
                    
                    ++it;
                    
                    if (it != taglocs.constEnd())
                    {
                        end = it.key();
                    }
                    
                    break;
                }
            }
        }
        
        QStringList lines;
        
        for (int i=start; i<end; ++i)
        {
            lines.append( expandedLines().constData()[i] );
        }
        
        return lines;
    };

    //return all of the lines associated with all copies of the passed directive
    auto getAllLines = [&](const QString &directive) -> QStringList
    {
        QStringList lines;
        
        for (int i=0; i<ntags.value(directive,0); ++i)
        {
            lines += getLines(directive, i);
        }
        
        return lines;
    };

    //interpret a bool from the passed string
    auto gromacs_toBool = [&](const QString &word, bool *ok)
    {
        QString w = word.toLower();

        if (ok) *ok = true;
        
        if (w == "yes" or w == "y" or w == "true" or w == "1")
        {
            return true;
        }
        else if (w == "no" or w == "n" or w == "false" or w == "0")
        {
            return false;
        }
        else
        {
            if (ok) *ok = false;
            return false;
        }
    };

    //internal function to process the defaults lines
    auto processDefaults = [&]()
    {
        QStringList warnings;
    
        //there should only be one defaults line
        const auto lines = getLines("defaults", 0);
        
        if (lines.isEmpty())
            throw SireIO::parse_error( QObject::tr(
                "The required data for the '[defaults]' directive in Gromacs is "
                "not supplied. This is not a valid Gromacs topology file!"), CODELOC );
        
        auto words = lines[0].split(" ", QString::SkipEmptyParts);
        
        //there should be five words; non-bonded function type, combinination rule,
        //                            generate pairs, fudge LJ and fudge QQ
        if (words.count() < 5)
        {
            throw SireIO::parse_error( QObject::tr(
                "There is insufficient data for the '[defaults]' line '%1'. This is "
                "not a valid Gromacs topology file!").arg(lines[0]), CODELOC );
        }

        bool ok;
        int nbtyp = words[0].toInt(&ok);
        
        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The first value for the '[defaults]' line '%1' is not an integer. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );
        
        int combrule = words[1].toInt(&ok);
        
        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The second value for the '[defaults]' line '%1' is not an integer. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );
        
        bool gen_pairs = gromacs_toBool(words[2], &ok);
        
        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The third value for the '[defaults]' line '%1' is not a yes/no. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );
        
        double lj = words[3].toDouble(&ok);
        
        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The fourth value for the '[defaults]' line '%1' is not a double. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );
        
        double qq = words[4].toDouble(&ok);
        
        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The fifth value for the '[defaults]' line '%1' is not a double. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );

        //validate and then save these values
        if (nbtyp <= 0 or nbtyp > 2)
        {
            warnings.append( QObject::tr("A non-supported non-bonded function type (%1) "
              "is requested.").arg(nbtyp) );
        }
        
        if (combrule <= 0 or combrule > 3)
        {
            warnings.append( QObject::tr("A non-supported combinig rule (%1) is requested!")
                .arg(combrule) );
        }
        
        if (lj < 0 or lj > 1)
        {
            warnings.append( QObject::tr("An invalid value of fudge_lj (%1) is requested!")
                .arg(lj) );
            
            if (lj < 0) lj = 0;
            else if (lj > 1) lj = 1;
        }
        
        if (qq < 0 or qq > 1)
        {
            warnings.append( QObject::tr("An invalid value of fudge_qq (%1) is requested1")
                .arg(qq) );
            
            if (qq < 0) qq = 0;
            else if (qq > 1) qq = 1;
        }
        
        nb_func_type = nbtyp;
        combining_rule = combrule;
        fudge_lj = lj;
        fudge_qq = qq;
        generate_pairs = gen_pairs;
        
        return warnings;
    };

    //internal function to process the atomtypes lines
    auto processAtomTypes = [&]()
    {
        QStringList warnings;
        
        //get all 'atomtypes' lines
        const auto lines = getAllLines("atomtypes");
    
        //the database of all atom types
        QHash<QString,GromacsAtomType> typs;
    
        //now parse each atom
        for (const auto line : lines)
        {
            const auto words = line.split(" ", QString::SkipEmptyParts);
            
            //should either have 2 words (atom type, mass) or
            //have 6 words; atom type, mass, charge, type, V, W or
            //have 7 words; atom type, atom number, mass, charge, type, V, W
            if (words.count() < 2)
            {
                warnings.append( QObject::tr( "There is not enough data for the "
                  "atomtype data '%1'. Skipping this line." ).arg(line) );
                continue;
            }
            
            GromacsAtomType typ;

            if (words.count() < 6)
            {
                //only getting the atom type and mass
                bool ok_mass;
                double mass = words[1].toDouble( &ok_mass );
                
                if (not ok_mass)
                {
                    warnings.append( QObject::tr( "Could not interpret the atom type data "
                       "from line '%1'. Skipping this line.").arg(line) );
                    continue;
                }
                
                typ = GromacsAtomType(words[0], mass*g_per_mol);
            }
            else if (words.count() < 7)
            {
                bool ok_mass, ok_charge, ok_ptyp, ok_v, ok_w;
                double mass = words[1].toDouble(&ok_mass);
                double chg = words[2].toDouble(&ok_charge);
                auto ptyp = GromacsAtomType::toParticleType(words[3], &ok_ptyp);
                double v = words[4].toDouble(&ok_v);
                double w = words[5].toDouble(&ok_w);
                
                if (not (ok_mass and ok_charge and ok_ptyp and ok_v and ok_w))
                {
                    warnings.append( QObject::tr( "Could not interpret the atom type data "
                      "from line '%1'. Skipping this line.").arg(line) );
                    continue;
                }
                
                typ = GromacsAtomType(words[0], mass*g_per_mol, chg*mod_electron,
                                      ptyp, ::toLJParameter(v,w,combining_rule));
            }
            else if (words.count() < 8)
            {
                bool ok_mass, ok_elem, ok_charge, ok_ptyp, ok_v, ok_w;
                int nprotons = words[1].toInt(&ok_elem);
                double mass = words[2].toDouble(&ok_mass);
                double chg = words[3].toDouble(&ok_charge);
                auto ptyp = GromacsAtomType::toParticleType(words[4], &ok_ptyp);
                double v = words[5].toDouble(&ok_v);
                double w = words[6].toDouble(&ok_w);
                
                if (not (ok_elem and ok_mass and ok_charge and ok_ptyp and ok_v and ok_w))
                {
                    warnings.append( QObject::tr( "Could not interpret the atom type data "
                      "from line '%1'. Skipping this line.").arg(line) );
                    continue;
                }
                
                typ = GromacsAtomType(words[0], mass*g_per_mol, chg*mod_electron,
                                      ptyp, ::toLJParameter(v,w,combining_rule),
                                      Element(nprotons));
            }
            
            if (typs.contains(typ.atomType()))
            {
                //only replace if the new type is fully specified
                if (typ.hasMassOnly())
                    continue;
            
                warnings.append( QObject::tr( "The data for atom type '%1' exists already! "
                 "This will now be replaced with new data from line '%2'")
                    .arg(typ.atomType()).arg(line) );
            }
            
            typs.insert( typ.atomType(), typ );
        }
    
        //save the database of types
        atom_types = typs;
    
        return warnings;
    };

    //internal function to process the bondtypes lines
    auto processBondTypes = [&]()
    {
        QStringList warnings;
        
        //get all 'bondtypes' lines
        const auto lines = getAllLines("bondtypes");

        //save into a database of bonds
        QHash<QString,GromacsBond> bnds;

        for (const auto line : lines)
        {
            //each line should contain the atom types of the two atoms, then
            //the function type, then the parameters for the function
            const auto words = line.split(" ", QString::SkipEmptyParts);
            
            if (words.count() < 3)
            {
                warnings.append( QObject::tr("There is not enough data on the "
                  "line '%1' to extract a Gromacs bond parameter. Skipping line.")
                    .arg(line) );
                continue;
            }
            
            const auto atm0 = words[0];
            const auto atm1 = words[1];
            
            bool ok;
            int func_type = words[2].toInt(&ok);
            
            if (not ok)
            {
                warnings.append( QObject::tr("Unable to determine the function type "
                  "for the bond on line '%1'. Skipping line.")
                    .arg(line) );
                continue;
            }
            
            //now read in all of the remaining values as numbers...
            QList<double> params;
            
            for (int i=3; i<words.count(); ++i)
            {
                double param = words[i].toDouble(&ok);
                
                if (ok) params.append(param);
            }
            
            GromacsBond bond;
            
            try
            {
                bond = GromacsBond(func_type, params);
            }
            catch(const SireError::exception &e)
            {
                warnings.append( QObject::tr("Unable to extract the correct information "
                  "to form a bond from line '%1'. Error is '%2'")
                    .arg(line).arg(e.error()) );
                continue;
            }
            
            QString key = get_bond_id(atm0,atm1);
            
            if (bnds.contains(key))
            {
                warnings.append( QObject::tr("Duplicate bond entry for %1 : %2 on line "
                   "'%3'. Replacing the original entry.")
                        .arg(atm0).arg(atm1).arg(line) );
            }
            
            bnds.insert(key, bond);
        }

        bond_potentials = bnds;

        return warnings;
    };

    //internal function to process the pairtypes lines
    auto processPairTypes = [&]()
    {
        return QStringList();
    };

    //internal function to process the angletypes lines
    auto processAngleTypes = [&]()
    {
        return QStringList();
    };

    //internal function to process the dihedraltypes lines
    auto processDihedralTypes = [&]()
    {
        return QStringList();
    };

    //internal function to process the constrainttypes lines
    auto processConstraintTypes = [&]()
    {
        return QStringList();
    };

    //internal function to process the nonbond_params lines
    auto processNonBondParams = [&]()
    {
        return QStringList();
    };
    
    //process the defaults data first, as this affects the rest of the parsing
    auto warnings = processDefaults();

    //now we can process the other tags
    const QVector< std::function<QStringList()> > funcs =
                 { processAtomTypes, processBondTypes, processPairTypes,
                   processAngleTypes, processDihedralTypes, processConstraintTypes,
                   processNonBondParams
                 };

    if (usesParallel())
    {
        QMutex mutex;
        
        tbb::parallel_for( tbb::blocked_range<int>(0, funcs.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            QStringList local_warnings;
            
            for (int i=r.begin(); i<r.end(); ++i)
            {
                local_warnings += funcs[i]();
            }
            
            if (not local_warnings.isEmpty())
            {
                QMutexLocker lkr(&mutex);
                warnings += local_warnings;
            }
        });
    }
    else
    {
        for (int i=0; i<funcs.count(); ++i)
        {
            warnings += funcs[i]();
        }
    }
    
    return warnings;
}

/** Interpret the fully expanded set of lines to extract all of the necessary data */
void GroTop::interpret()
{
    //first, go through and find the line numbers of all tags
    const QRegularExpression re("\\[\\s*([\\w\\d]+)\\s*\\]");

    //map giving the type and line number of each directive tag
    QMap<int,QString> taglocs;

    const int nlines = expandedLines().count();
    const auto lines = expandedLines().constData();

    //run through this file to find all of the directives
    if (usesParallel())
    {
        QMutex mutex;
    
        tbb::parallel_for( tbb::blocked_range<int>(0,nlines),
                           [&](const tbb::blocked_range<int> &r)
        {
            QMap<int,QString> mylocs;
        
            for (int i=r.begin(); i<r.end(); ++i)
            {
                auto m = re.match(lines[i]);
                
                if (m.hasMatch())
                {
                    auto tag = m.captured(1);
                    mylocs.insert(i,tag);
                }
            }
            
            if (not mylocs.isEmpty())
            {
                QMutexLocker lkr(&mutex);
                
                for (auto it=mylocs.constBegin(); it!=mylocs.constEnd(); ++it)
                {
                    taglocs.insert(it.key(), it.value());
                }
            }
        });
    }
    else
    {
        for (int i=0; i<nlines; ++i)
        {
            auto m = re.match(lines[i]);
            
            if (m.hasMatch())
            {
                auto tag = m.captured(1);
                taglocs.insert(i,tag);
            }
        }
    }

    //now, validate that this looks like a gromacs top file. Rules are taken
    //from page 138 of the Gromacs 5.1 PDF reference manual
    
    //first, count up the number of each tag
    QHash<QString,int> ntags;
    
    for (auto it = taglocs.constBegin(); it != taglocs.constEnd(); ++it)
    {
        if (not ntags.contains(it.value()))
        {
            ntags.insert(it.value(), 1);
        }
        else
        {
            ntags[it.value()] += 1;
        }
    }

    //there should be only one 'defaults' tag
    if (ntags.value("defaults", 0) != 1)
    {
        throw SireIO::parse_error( QObject::tr(
            "This is not a valid GROMACS topology file. Such files contain one, and one "
            "only 'defaults' directive. The number of such directives in this file is %1.")
                .arg(ntags.value("defaults",0)), CODELOC );
    }

    //now process all of the directives
    auto warnings = this->processDirectives(taglocs, ntags);

    if (not warnings.isEmpty())
    {
        qDebug() << warnings.join("\n");
    }

    this->setScore(100);
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void GroTop::parseLines(const QString &path, const PropertyMap &map)
{
    //first, see if there are any GROMACS defines in the passed map
    //and then preprocess the lines to create the fully expanded file to parse
    {
        QHash<QString,QString> defines;

        try
        {
            const auto p = map["GROMACS_DEFINE"];
            
            QStringList d;
            
            if (p.hasValue())
            {
                d = p.value().asA<StringProperty>().toString().split(":", QString::SkipEmptyParts);
            }
            else if (p.source() != "GROMACS_DEFINE")
            {
                d = p.source().split(":", QString::SkipEmptyParts);
            }
            
            for (const auto define : d)
            {
                auto words = define.split("=");
                
                if (words.count() == 1)
                {
                    defines.insert( words[0].simplified(), "1" );
                }
                else
                {
                    defines.insert( words[0].simplified(), words[1].simplified() );
                }
            }
        }
        catch(...)
        {}

        //now go through an expand any macros and include the contents of any
        //included files
        expanded_lines = preprocess(lines(), defines, path, ".");
    }

    //now we know that there are no macros to expand, no other files to
    //include, and everything should be ok... ;-)
    this->interpret();
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
