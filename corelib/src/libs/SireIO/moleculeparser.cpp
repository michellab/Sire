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
#include "supplementary.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireBase/parallel.h"
#include "SireBase/booleanproperty.h"
#include "SireBase/stringproperty.h"

#include "SireFF/ffdetail.h"
#include "SireMM/mmdetail.h"

#include "SireMol/molecule.h"
#include "SireMol/trajectory.h"

#include "SireSystem/system.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QFile>
#include <QTextStream>
#include <QElapsedTimer>
#include <QFileInfo>
#include <QMutex>
#include <QDebug>

using namespace SireIO;
using namespace SireSystem;
using namespace SireFF;
using namespace SireBase;
using namespace SireStream;

//////////////
////////////// Implementation of ParserFactory and ParserFactoryHelper
//////////////

namespace SireIO
{
namespace detail
{
    /** Null constructor */
    ParserFactoryHelper::ParserFactoryHelper()
    {}

    /** Copy constructor */
    ParserFactoryHelper::ParserFactoryHelper(const ParserFactoryHelper &other)
                        : parser(other.parser)
    {}

    /** Destructor */
    ParserFactoryHelper::~ParserFactoryHelper()
    {}


    ParserFactoryHelper& ParserFactoryHelper::operator=(const ParserFactoryHelper &other)
    {
        parser = other.parser;
        return *this;
    }

    bool ParserFactoryHelper::operator<(const ParserFactoryHelper &other) const
    {
        if (isValid())
        {
            if (other.isValid())
            {
                return parser->formatName() < other.parser->formatName();
            }
            else
            {
                return true;
            }
        }
        else
        {
            return not other.isValid();
        }
    }

    bool ParserFactoryHelper::operator==(const ParserFactoryHelper &other) const
    {
        if (isValid())
        {
            if (other.isValid())
            {
                return parser->formatName() == other.parser->formatName();
            }
            else
                return false;
        }
        else
        {
            return not other.isValid();
        }
    }

    bool ParserFactoryHelper::operator>(const ParserFactoryHelper &other) const
    {
        return not (operator==(other) or operator<(other));
    }

    bool ParserFactoryHelper::operator!=(const ParserFactoryHelper &other) const
    {
        return not operator==(other);
    }

    bool ParserFactoryHelper::operator<=(const ParserFactoryHelper &other) const
    {
        return operator==(other) or operator<(other);
    }

    bool ParserFactoryHelper::operator>=(const ParserFactoryHelper &other) const
    {
        return not operator<(other);
    }

    /** Return whether or not this helper is valid */
    bool ParserFactoryHelper::isValid() const
    {
        return parser.get() != 0;
    }

    /** Return whether or not this is the supplementary parser */
    bool ParserFactoryHelper::isSupplementary() const
    {
        return parser.get() != 0 and parser->isA<Supplementary>();
    }

    /** Return the unique ID name of the parser in the program */
    QString ParserFactoryHelper::formatName() const
    {
        if (isValid())
        {
            return parser->formatName();
        }
        else
        {
            return QString();
        }
    }

    /** Return the description of the parser */
    QString ParserFactoryHelper::formatDescription() const
    {
        if (isValid())
        {
            return parser->formatDescription();
        }
        else
        {
            return QString();
        }
    }

    QString ParserFactoryHelper::toString() const
    {
        return QString("Parser( %1 : %2 )").arg(formatName()).arg(formatDescription());
    }

    /** Return all of the suffixes recognised by this parser, in their order
        of preference */
    QStringList ParserFactoryHelper::suffixes() const
    {
        if (isValid())
        {
            return parser->formatSuffix();
        }
        else
        {
            return QStringList();
        }
    }

    /** Return the preferred suffix for the parser */
    QString ParserFactoryHelper::preferredSuffix() const
    {
        const auto s = this->suffixes();

        if (not s.isEmpty())
        {
            return s[0];
        }
        else
        {
            return QString();
        }
    }

    /** Use this factory helper to construct a new parser that parses
        the file called 'filename' */
    MoleculeParserPtr ParserFactoryHelper::construct(const QString &filename,
                                                     const PropertyMap &map) const
    {
        if (isValid())
        {
            return parser->construct(filename, map);
        }
        else
            return MoleculeParserPtr();
    }

    /** Use this factory helper to construct a new parser that parses
        the data in the passed lines of text */
    MoleculeParserPtr ParserFactoryHelper::construct(const QStringList &lines,
                                                     const PropertyMap &map) const
    {
        if (isValid())
        {
            return parser->construct(lines, map);
        }
        else
            return MoleculeParserPtr();
    }

    /** Use this factory helper to construct a new parser from the information
        contained in the passed system */
    MoleculeParserPtr ParserFactoryHelper::construct(const SireSystem::System &system,
                                                     const PropertyMap &map) const
    {
        if (isValid())
        {
            return parser->construct(system, map);
        }
        else
            return MoleculeParserPtr();
    }

    /** The parser factory */
    class ParserFactory
    {
    public:
        ParserFactory()
        {}

        ~ParserFactory()
        {}

        void registerParser(const ParserFactoryHelper &helper)
        {
            if (not helper.isValid())
            {
                return;
            }

            QMutexLocker lkr(&mutex);

            helpers_by_id.insert( helper.formatName(), helper );

            for (const auto &suffix : helper.suffixes())
            {
                helpers_by_suffix.insert(suffix, helper);
            }
        }

        QList<ParserFactoryHelper> getFactories(const QStringList &parser_names)
        {
            QMutexLocker lkr(&mutex);
            QList<ParserFactoryHelper> helpers;
            QStringList missing;

            for (const auto &name : parser_names)
            {
                auto helper = helpers_by_id.value(name);

                if (not helper.isValid())
                {
                    // search for the helped in the secondary suffixes...
                    for (auto h : helpers_by_id.values())
                    {
                        for (auto suffix : h.suffixes())
                        {
                            if (name.toLower() == suffix.toLower())
                            {
                                helper = h;
                                break;
                            }
                        }

                        if (helper.isValid())
                            break;
                    }
                }

                helpers.append(helper);

                if (not helpers.last().isValid())
                {
                    missing.append(name);
                }
            }

            if (not missing.isEmpty())
            {
                lkr.unlock();
                throw SireError::io_error( QObject::tr(
                        "Cannot find parsers that support the following formats: %1.\n"
                        "Supported parsers are:\n%2")
                            .arg(missing.join(", "))
                            .arg(this->supportedFormats()), CODELOC );
            }

            return helpers;
        }

        QList<ParserFactoryHelper> factoriesForSuffix(const QString &suffix,
                                                      bool disable_supplementary)
        {
            QMutexLocker lkr(&mutex);
            auto helpers = helpers_by_suffix.values(suffix);
            std::sort(helpers.begin(), helpers.end());

            if (disable_supplementary)
            {
                QMutableListIterator<ParserFactoryHelper> it(helpers);

                while (it.hasNext())
                {
                    const auto &value = it.next();

                    if (value.isSupplementary())
                    {
                        it.remove();
                    }
                }
            }

            return helpers;
        }

        QList<ParserFactoryHelper> factoriesExcludingSuffix(const QString &suffix,
                                                            bool disable_supplementary)
        {
            QMutexLocker lkr(&mutex);

            if (suffix.isEmpty())
            {
                auto helpers = helpers_by_id.values();
                std::sort(helpers.begin(), helpers.end());

                if (disable_supplementary)
                {
                    QMutableListIterator<ParserFactoryHelper> it(helpers);

                    while (it.hasNext())
                    {
                        const auto &value = it.next();

                        if (value.isSupplementary())
                        {
                            it.remove();
                        }
                    }
                }

                return helpers;
            }

            QList<ParserFactoryHelper> helpers;

            for (const auto &helper : helpers_by_id)
            {
                if (not helper.suffixes().contains(suffix))
                {
                    helpers.append(helper);
                }
            }

            std::sort(helpers.begin(), helpers.end());

            if (disable_supplementary)
            {
                QMutableListIterator<ParserFactoryHelper> it(helpers);

                while (it.hasNext())
                {
                    const auto &value = it.next();

                    if (value.isSupplementary())
                    {
                        it.remove();
                    }
                }
            }

            return helpers;
        }

        ParserFactoryHelper factory(const QString &name)
        {
            QMutexLocker lkr(&mutex);
            return helpers_by_id.value(name);
        }

        QString supportedFormats()
        {
            QMutexLocker lkr(&mutex);

            auto keys = helpers_by_id.keys();
            std::sort(keys.begin(), keys.end());

            QStringList lines;

            for (const auto &key : keys)
            {
                const auto parser = helpers_by_id.value(key);

                if (not parser.isSupplementary())
                {
                    lines.append( QObject::tr("## Parser %1 ##").arg(key) );
                    lines.append( QObject::tr("Supports files: %1")
                                        .arg(parser.suffixes().join(", ")) );
                    lines.append( parser.formatDescription() );
                    lines += QString("#").repeated(13 + key.length()) + "\n";
                }
            }

            return lines.join("\n");
        }

    private:
        /** Mutex to serialise access to the factory */
        QMutex mutex;

        /** All of the factory helpers arranged by the suffix of
            file that they support */
        QMultiHash<QString,ParserFactoryHelper> helpers_by_suffix;

        /** All of the factor helpers arranged by their unique ID */
        QHash<QString,ParserFactoryHelper> helpers_by_id;
    };

} // end of namespace detail
} // end of namespace SireIO

Q_GLOBAL_STATIC( SireIO::detail::ParserFactory, getParserFactory );

/** This registers a ParserFactoryHelper with the ParserFactory for the
    specified parser */
SireIO::detail::ParserFactoryHelper::ParserFactoryHelper(MoleculeParser *p)
{
    parser.reset( p );
    getParserFactory()->registerParser(*this);
}

//////////////
////////////// Implementation of MoleculeParser
//////////////

static const RegisterMetaType<MoleculeParser> r_parser( MAGIC_ONLY, MoleculeParser::typeName() );

QDataStream &operator<<(QDataStream &ds, const MoleculeParser &parser)
{
    writeHeader(ds, r_parser, 2);

    SharedDataStream sds(ds);
    sds << parser.fname << parser.lnes << parser.scr << parser.run_parallel
        << static_cast<const Property&>(parser);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, MoleculeParser &parser)
{
    VersionID v = readHeader(ds, r_parser);

    if (v == 2)
    {
        SharedDataStream sds(ds);
        sds >> parser.fname >> parser.lnes >> parser.scr >> parser.run_parallel
            >> static_cast<Property&>(parser);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> parser.lnes >> parser.scr >> parser.run_parallel
            >> static_cast<Property&>(parser);

        parser.fname = QString();
    }
    else
        throw version_error(v,"1", r_parser, CODELOC);

    return ds;
}

/** Constructor */
MoleculeParser::MoleculeParser(const PropertyMap &map) : Property(), scr(0), run_parallel(true)
{
    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }
}

/** Internal function that provides a file cache */
class FileContentsCache
{
public:
    FileContentsCache()
    {}

    ~FileContentsCache()
    {}

    QVector<QString> read(QString filename)
    {
        QMutexLocker lkr(&cache_mutex);
        auto it = cache.constFind(filename);

        if (it != cache.constEnd())
        {
            return it.value();
        }
        else
        {
            return QVector<QString>();
        }
    }

    void save(QString filename, QVector<QString> filecontents)
    {
        QMutexLocker lkr(&cache_mutex);
        cache[filename] = filecontents;
    }

    void clear()
    {
        QMutexLocker lkr(&cache_mutex);
        cache.clear();
    }

private:
    QMutex cache_mutex;
    QHash< QString, QVector<QString> > cache;
};

Q_GLOBAL_STATIC( FileContentsCache, getFileCache );

/** Internal function that can be used by the parsers to read the contents
    of a text file into memory. This uses a cache to ensure that every file
    is read only once */
QVector<QString> MoleculeParser::readTextFile(QString filename)
{
    filename = QFileInfo(filename).absoluteFilePath();

    QVector<QString> lines = getFileCache()->read(filename);

    if (not lines.isEmpty())
        return lines;

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

    lines = l.toVector();

    if (not lines.isEmpty())
        getFileCache()->save(filename, lines);

    return lines;
}

/** Construct the parser, parsing in all of the lines in the file
    with passed filename */
MoleculeParser::MoleculeParser(const QString &filename,
                               const PropertyMap &map)
               : Property(), scr(0), run_parallel(true)
{
    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    fname = QFileInfo(filename).absoluteFilePath();
    lnes = readTextFile(fname);
}

/** Construct the parser, parsing in all of the passed text lines */
MoleculeParser::MoleculeParser(const QStringList &lines,
                               const PropertyMap &map)
               : Property(), scr(0), run_parallel(true)
{
    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    if (not lines.isEmpty())
    {
        lnes = lines.toVector();
    }
}

/** Copy constructor */
MoleculeParser::MoleculeParser(const MoleculeParser &other)
               : Property(other), fname(other.fname),
                 lnes(other.lnes), scr(other.scr),
                 run_parallel(other.run_parallel)
{}

/** Destructor */
MoleculeParser::~MoleculeParser()
{}

/** Remove any comment lines (those that start with 'comment_flag')
 *  from the file. This should make parsing easier
*/
void MoleculeParser::removeCommentLines(const QString &comment_flag)
{
    bool has_comments = false;

    for (int i=0; i<lnes.count(); ++i)
    {
        const auto &line = lnes.constData()[i];

        if (line.startsWith(comment_flag))
        {
            has_comments = true;
            break;
        }
    }

    if (has_comments)
    {
        QMutableVectorIterator<QString> it(lnes);

        while(it.hasNext())
        {
            const auto &line = it.next();
            if (line.startsWith(comment_flag))
            {
                it.remove();
            }
        }
    }
}

/** Function used by derived classes to set the lines */
void MoleculeParser::setLines( const QVector<QString> &lines )
{
    lnes = lines;
}

/** Functions used by derived classes to set the filename */
void MoleculeParser::setFilename(const QString &filename)
{
    fname = QFileInfo(filename).absoluteFilePath();
}

const char* MoleculeParser::typeName()
{
    return "SireIO::MoleculeParser";
}

MoleculeParser& MoleculeParser::operator=(const MoleculeParser &other)
{
    if (this != &other)
    {
        fname = other.fname;
        lnes = other.lnes;
        scr = other.scr;
        run_parallel = other.run_parallel;
        Property::operator=(other);
    }

    return *this;
}

bool MoleculeParser::operator==(const MoleculeParser &other) const
{
    return fname == other.fname and lnes == other.lnes and scr == other.scr and
           run_parallel == other.run_parallel and Property::operator==(other);
}

bool MoleculeParser::operator!=(const MoleculeParser &other) const
{
    return not MoleculeParser::operator==(other);
}

/** Return the name of the file that was parsed */
QString MoleculeParser::filename() const
{
    return fname;
}

/** Return whether or not this parser is broken */
bool MoleculeParser::isBroken() const
{
    return false;
}

/** Return the error report, if this parser is broken. If it isn't,
 *  then an empty string is returned. */
QString MoleculeParser::errorReport() const
{
    return QString();
}

/** Return the number of trajectory frames contained in this parser.
 *  Trajectory frames contain coordinates and/or velocities and/or
 *  forces data. It is possible for a parser to have zero frames,
 *  e.g. if it only contains topology information.
*/
int MoleculeParser::nFrames() const
{
    return 0;
}

/** Return the ith trajectory frame from this parser. Note that
 *  some parsers may have to re-read the file, so this may fail
 *  if the filename changes since the last time this parser
 *  was used
 */
SireMol::Frame MoleculeParser::getFrame(int i) const
{
    // this will raise an exception as we will only
    // be calling this function if the parser has zero frames
    i = SireID::Index(i).map(0);

    return SireMol::Frame();
}

/** Enable code to parse files in parallel */
void MoleculeParser::enableParallel()
{
    run_parallel = true;
}

/** Disable code to parse files in parallel - parsing will happen in serial */
void MoleculeParser::disableParallel()
{
    run_parallel = false;
}

/** Set whether or not to parse files in parallel or serial */
void MoleculeParser::setUseParallel(bool on)
{
    run_parallel = on;
}

/** Return whether or not this is a lead parser. The lead parser is responsible
    for starting the process of turning the parsed file into the System. There
    must be one and one-only lead parser in a set of parsers creating a System */
bool MoleculeParser::isLead() const
{
    return false;
}

/** Return whether or not this parser can follow a lead parser and add data
    to an existing molecular system. By default, all parsers can follow. */
bool MoleculeParser::canFollow() const
{
    return true;
}

/** Extract and return a FFDetail forcefield that is compatible with all of the
    molecules in this system, using the passed property map to find the property.
    Note that this will raise an incompatible_error exception if there is no
    forcefield that adequately covers all of the molecules */
PropertyPtr MoleculeParser::getForceField(const System &system, const PropertyMap &map) const
{
    const auto ffprop = map["forcefield"];

    //make sure that the user is not telling us a specific forcefield to use
    if (ffprop.hasValue())
    {
        if (not ffprop.value().isA<FFDetail>())
            throw SireError::incompatible_error( QObject::tr(
                "Cannot convert the passed mapped value '%1' to an object of type FFDetail. "
                "If you want to specify the focefield it must be an object derived from FFDetail.")
                    .arg(ffprop.value().toString()), CODELOC );

        return ffprop.value();
    }

    const auto molnums = system.molNums().toVector();

    if (molnums.isEmpty())
        return SireMM::MMDetail();

    PropertyPtr ffield;
    QStringList errors;

    for (int i=0; i<molnums.count(); ++i)
    {
        const auto mol = system[ molnums[i] ].molecule();
        PropertyPtr molff;

        if (mol.hasProperty(ffprop))
        {
            const auto &p = mol.property(ffprop);

            if (p.isA<FFDetail>())
                molff = p;
        }

        if (molff.isNull())
        {
            errors.append( QObject::tr("Molecule '%1' does not have a valid 'forcefield' "
               "property. Please make sure that it has a property called '%2' that is "
               "derived from FFDetail")
                    .arg(mol.toString()).arg(ffprop.source()) );
        }
        else if (ffield.isNull())
        {
            //this is the first valid forcefield
            ffield = molff;
        }
        else if (not ffield.read().asA<FFDetail>()
                           .isCompatibleWith(molff.read().asA<FFDetail>()))
        {
            //incompatible forcefields!
            errors.append( QObject::tr("The forcefield for molecule '%1' is not compatible "
              "with that for other molecules.\n%2\nversus\n%3.")
                    .arg(mol.toString()).arg(molff.read().toString())
                    .arg(ffield.read().toString()) );
        }
    }

    if (not errors.isEmpty())
        throw SireError::incompatible_error( QObject::tr("There were some problems "
          "extracting a valid forcefield object from all of the molecules.\n\n%1")
            .arg(errors.join("\n")), CODELOC );

    //we should have a valid FFDetail now...
    return ffield.read().asA<FFDetail>();
}

/** Write the parsed data back to the file called 'filename'. This will
    overwrite the file if it exists already, so be careful! */
void MoleculeParser::writeToFile(const QString &filename) const
{
    if (lnes.isEmpty())
        return;

    if (not this->isTextFile())
        throw SireError::program_bug( QObject::tr(
            "Dear programmer - please override the MoleculeParser::writeToFile function "
            "to work with your binary file format. Text writing is not supported "
            "for the parser %1.").arg(this->what()), CODELOC );

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
}

/** Internal function that actually tries to parse the supplied file with name
    'filename'. This will try to find a parser based on suffix, but if that fails,
    it will try all parsers. It will choose the parser that doesn't raise an error
    that scores highest. */
MoleculeParserPtr MoleculeParser::_pvt_parse(const QString &filename,
                                             const PropertyMap &map)
{
    bool disable_supplementary = false;

    if (map.specified("DISABLE_SUPPLEMENTARY"))
    {
        disable_supplementary = map["DISABLE_SUPPLEMENTARY"]
                                    .value().asA<BooleanProperty>().value();
    }

    QFileInfo info(filename);

    if (not (info.isFile() and info.isReadable()))
    {
        throw SireError::file_error( QObject::tr(
                "There is no file readable called '%1'.")
                    .arg(filename), CODELOC );
    }

    //try to find the right parser based on the suffix
    QString suffix = info.suffix();

    QStringList errors;
    QStringList suffix_errors;
    QMap<float,MoleculeParserPtr> parsers;

    if (not suffix.isEmpty())
    {
        for (auto factory : getParserFactory()->factoriesForSuffix(suffix, disable_supplementary))
        {
            try
            {
                const auto parser = factory.construct(info.absoluteFilePath(), map);

                if (parser.read().score() <= 0)
                {
                    suffix_errors.append( QObject::tr("*-- Failed to parse '%1' with parser '%2'.\n"
                       "The file is not recognised as being of the required format.")
                                    .arg(filename).arg(factory.formatName()) );
                }
                else
                {
                    parsers.insert(parser.read().score(), parser);
                }
            }
            catch(const SireError::exception &e)
            {
                suffix_errors.append( QObject::tr("*-- Failed to parse '%1' with parser '%2'.\n%3")
                                .arg(filename).arg(factory.formatName())
                                .arg(e.error()) );
            }
        }

        if (not parsers.isEmpty())
        {
            //return the parser with the highest score
            return parsers.last();
        }
    }

    //none of the tested parsers worked, so let's now try all of the parsers
    for (auto factory : getParserFactory()->factoriesExcludingSuffix(suffix, disable_supplementary))
    {
        try
        {
            const auto parser = factory.construct(filename, map);

            if (parser.read().score() <= 0)
            {
                errors.append( QObject::tr("Failed to parse '%1' with parser '%2' "
                       "as this file is not recognised as being of the required format.")
                                    .arg(filename).arg(factory.formatName()) );
            }
            else
            {
                parsers.insert(parser.read().score(), parser);
            }
        }
        catch(const SireError::exception &e)
        {
            errors.append( QObject::tr("Failed to parse '%1' with parser '%2'\n%3")
                            .arg(filename).arg(factory.formatName())
                            .arg(e.error()) );
        }
    }

    if (not parsers.isEmpty())
    {
        return parsers.last();
    }
    else
    {
        if (not suffix.isEmpty())
        {
            return MoleculeParserPtr(BrokenParser(filename, suffix,
                                                  suffix_errors));
        }
        else
        {
            return MoleculeParserPtr(BrokenParser(filename,
                                                  suffix_errors + errors));
        }

        return MoleculeParserPtr();
    }
}

/** Parse the passed file, returning the resulting Parser. This employs a lot
    of magic to automatically work out the format of the file and whether or
    not this is parseable by Sire... This raises an exception if the file
    cannot be recognised, or if there is an error in parsing. */
MoleculeParserPtr MoleculeParser::parse(const QString &filename,
                                        const PropertyMap &map)
{
    MoleculeParserPtr parser = MoleculeParser::_pvt_parse(filename, map);
    getFileCache()->clear();
    return parser;
}

/** Parse the passed set of files, returning the resulting Parsers */
QList<MoleculeParserPtr> MoleculeParser::parse(const QStringList &filenames,
                                               const PropertyMap &map)
{
    QList<MoleculeParserPtr> result;

    if (filenames.count() == 1)
    {
        result.append( MoleculeParser::_pvt_parse(filenames[0], map) );
    }
    else
    {
        QVector<MoleculeParserPtr> parsers(filenames.count());

        bool run_parallel = true;

        if (map["parallel"].hasValue())
        {
            run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
        }

        if (run_parallel)
        {
            //parse the files in parallel - we use a grain size of 1
            //as each file can be pretty big, and there won't be many of them
            tbb::parallel_for( tbb::blocked_range<int>(0,filenames.count(),1),
                               [&](tbb::blocked_range<int> r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    parsers[i] = MoleculeParser::_pvt_parse(filenames[i], map);
                }
            }, tbb::simple_partitioner());
        }
        else
        {
            for (int i=0; i<filenames.count(); ++i)
            {
                parsers[i] = MoleculeParser::_pvt_parse(filenames[i], map);
            }
        }

        result = parsers.toList();
    }

    getFileCache()->clear();

    return result;
}

/** Read the passed file called 'filename', returning the System contained therein */
System MoleculeParser::read(const QString &filename, const PropertyMap &map)
{
    MoleculeParserPtr parser = MoleculeParser::parse(filename, map);
    auto system = parser.read().toSystem(map);

    if (system.name().isEmpty())
    {
        system.setName(QFileInfo(filename).baseName());
    }

    return system;
}

/** Read the two passed files, returning the System contained therein. The two
    files must refer to the same System, i.e. they could be a parameter + coordinate file */
System MoleculeParser::read(const QString &file1, const QString &file2,
                            const PropertyMap &map)
{
    MoleculeParserPtr parser1, parser2;

    bool run_parallel = true;

    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    if (run_parallel)
    {
        tbb::parallel_invoke( [&](){ parser1 = MoleculeParser::parse(file1,map); },
                              [&](){ parser2 = MoleculeParser::parse(file2,map); } );
    }
    else
    {
        parser1 = MoleculeParser::parse(file1,map);
        parser2 = MoleculeParser::parse(file2,map);
    }

    auto system = parser1.read().toSystem(parser2.read(),map);

    if (system.name().isEmpty())
    {
        auto p1 = QFileInfo(file1).baseName();
        auto p2 = QFileInfo(file2).baseName();

        if (p1 == p2)
        {
            system.setName(p1);
        }
        else
        {
            system.setName(QString("%1:%2").arg(p1).arg(p2));
        }
    }

    return system;
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

    auto system = parser.read().toSystem(parsers, map);

    if (system.name().isEmpty())
    {
        QSet<QString> parts;
        for (const auto &filename : filenames)
        {
            parts.insert(QFileInfo(filename).baseName());
        }

        //QStringList names(parts.constBegin(), parts.constEnd());
        QStringList names = parts.values();

        system.setName(names.join(":"));
    }

    return system;
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QString &filename, const PropertyMap &map)
{
    return MoleculeParser::read(filename, map);
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QString &file1, const QString &file2, const PropertyMap &map)
{
    return MoleculeParser::read(file1, file2, map);
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QStringList &filenames, const PropertyMap &map)
{
    return MoleculeParser::read(filenames, map);
}

/** Parse the passed file, returning the resulting Parser. This employs a lot
    of magic to automatically work out the format of the file and whether or
    not this is parseable by Sire... This raises an exception if the file
    cannot be recognised, or if there is an error in parsing. */
MoleculeParserPtr MoleculeParser::parse(const QString &filename)
{
    return parse(filename, PropertyMap());
}

/** Parse the passed set of files, returning the resulting Parsers */
QList<MoleculeParserPtr> MoleculeParser::parse(const QStringList &filenames)
{
    return parse(filenames, PropertyMap());
}

/** Read the passed file called 'filename', returning the System contained therein */
System MoleculeParser::read(const QString &filename)
{
    return read(filename, PropertyMap());
}

/** Read the two passed files, returning the System contained therein. The two
    files must refer to the same System, i.e. they could be a parameter + coordinate file */
System MoleculeParser::read(const QString &file1, const QString &file2)
{
    return read(file1, file2, PropertyMap());
}

/** Read the files with passed filenames, returning the System contained therein.
    Note that all of the files must be connected to the same system
    (i.e. it could be the Amber Parm and Rst file) */
System MoleculeParser::read(const QStringList &filenames)
{
    return read(filenames, PropertyMap());
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QString &filename)
{
    return load(filename, PropertyMap());
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QString &file1, const QString &file2)
{
    return load(file1, file2, PropertyMap());
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QStringList &filenames)
{
    return load(filenames, PropertyMap());
}

/** Return the suffix (or suffixes) given to files that support this format.
    The first suffix is the preferred on to use */
QStringList MoleculeParser::formatSuffix() const
{
    //just use a lower-case version of the format name
    return QStringList( this->formatName().toLower() );
}

/** This returns a human readable set of lines describing the formats supported
    by MoleculeParser. Each line is formatted as "extension : description" where
    extension is the unique extension of the file used by MoleculeParser, and
    description is a description of the file format */
QString MoleculeParser::supportedFormats()
{
    return getParserFactory()->supportedFormats();
}

QStringList pvt_write(const System &system,
                      const QStringList &filenames,
                      const QStringList &fileformats,
                      const PropertyMap &map)
{
    if (filenames.count() != fileformats.count())
    {
        throw SireError::program_bug( QObject::tr(
            "Disagreement of the number of files... %1 vs %2")
                .arg(filenames.count()).arg(fileformats.count()), CODELOC );
    }

    QVector<QFileInfo> fileinfos(filenames.count());

    QStringList errors;

    for (int i=0; i<filenames.count(); ++i)
    {
        fileinfos[i] = QFileInfo(filenames[i]);

        if (fileinfos[i].exists())
        {
            if (fileinfos[i].isDir())
            {
                errors.append( QObject::tr("The file %1 is actually a directory, and not writable!")
                                    .arg(fileinfos[i].absoluteFilePath()) );
            }
            else if (not fileinfos[i].isWritable())
            {
                errors.append( QObject::tr("The file %1 exists and is not writable!")
                                    .arg(fileinfos[i].absoluteFilePath()) );
            }
        }
    }

    if (not errors.isEmpty())
    {
        throw SireError::io_error( QObject::tr(
                "Cannot write the files as the following errors occurred:\n%1")
                    .arg(errors.join("\n\n")), CODELOC );
    }

    //now get all of the parsers
    const auto factories = getParserFactory()->getFactories(fileformats);

    QVector<QString> written_files(filenames.count());

    //should we write the files in parallel?
    bool run_parallel = true;

    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    if (run_parallel)
    {
        tbb::spin_mutex error_mutex;

        tbb::parallel_for( tbb::blocked_range<int>(0,filenames.count(),1),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                const auto filename = fileinfos[i].absoluteFilePath();

                try
                {
                    written_files[i] = filename;
                    factories[i].construct(system, map).read().writeToFile(filename);
                }
                catch(const SireError::exception &e)
                {
                    tbb::spin_mutex::scoped_lock locker(error_mutex);
                    errors.append( QObject::tr("Failed to write the file '%1' using the parser "
                       "for fileformat '%2'. Errors reported by the parser are:\n%3")
                            .arg(filename).arg(factories[i].formatName())
                            .arg(e.error()) );
                }
            }
        }, tbb::simple_partitioner());
    }
    else
    {
        for (int i=0; i<filenames.count(); ++i)
        {
            const auto filename = fileinfos[i].absoluteFilePath();

            try
            {
                written_files[i] = filename;
                factories[i].construct(system, map).read().writeToFile(filename);
            }
            catch(const SireError::exception &e)
            {
                errors.append( QObject::tr("Failed to write the file '%1' using the parser "
                   "for fileformat '%2'. Errors reported by the parser are:\n%3")
                        .arg(filename).arg(factories[i].formatName())
                        .arg(e.error()) );
            }
        }
    }

    if (not errors.isEmpty())
    {
        throw SireError::io_error( QObject::tr(
                "Cannot write the (perhaps some of the ) files "
                "as the following errors occurred:\n%1")
                    .arg(errors.join("\n\n")), CODELOC );
    }

    return written_files.toList();
}

/** Save the passed System to the file called 'filename'. First, the 'fileformat'
    property is looked at in 'map'. This is used to set the format(s) of
    the files that are written (comma-separated list).

    If this does not exist, then the extension of the
    file is used to work out which format to use. If no extension is given,
    then the System will be queried to find out its preferred format (normally
    the format it was loaded with), via its 'fileformat' property
    (again, comma separated list).

    If their preferred format results in multiple files, then
    multiple files will be written. This returns the full pathnames to
    all of the files that are written
*/
QStringList MoleculeParser::write(const System &system, const QString &filename,
                                  const PropertyMap &map)
{
    if (filename.isEmpty() or QFileInfo(filename).baseName().isEmpty())
    {
        throw SireError::io_error( QObject::tr(
                "You must supply a valid filename. This '%1' is not sufficient.")
                    .arg(filename), CODELOC );
    }

    //build a list of filenames with their associated fileformats
    QStringList filenames;
    QStringList fileformats;

    const auto format_property = map["fileformat"];

    if (format_property.hasValue())
    {
        try
        {
            fileformats = format_property.value().asA<StringProperty>().toString().split(",");
        }
        catch(...)
        {
            fileformats.append( format_property.value().asA<MoleculeParser>()
                                               .formatName() );
        }

        QString basename = QFileInfo(filename).completeBaseName();

        for (const auto &format : fileformats)
        {
            filenames.append( QString("%1.%2").arg(basename,format.toLower()) );
        }
    }
    else
    {
        QString extension = QFileInfo(filename).completeSuffix();

        if (extension.isEmpty())
        {
            //we need to find the format from the system
            try
            {
                fileformats = system.property(format_property).asA<StringProperty>()
                                                              .toString().split(",");
            }
            catch(...)
            {
                throw SireError::io_error( QObject::tr(
                        "Cannot work out the fileformat to use to write the System to "
                        "file '%1'. You need to either supply the format using the "
                        "'fileformat' property in the passed map, add this to the System "
                        "as its 'fileformat' property, or pass a filename with an extension "
                        "whose fileformat can be determined. Supported fileformats are;\n%2")
                            .arg(filename).arg(MoleculeParser::supportedFormats()),
                                CODELOC );
            }

            for (const auto &format : fileformats)
            {
                filenames.append( QString("%1.%2").arg(filename).arg(format.toLower()) );
            }
        }
        else
        {
            filenames.append(filename);
            fileformats.append(extension.toUpper());
        }
    }

    //now we have a list of filenames and associated formats, actually
    //write the files
    return ::pvt_write(system, filenames, fileformats, map);
}

/** Extension of MoleculeParser::write which allows many filenames.
    The same rules to locate the fileformats are now used, except that now only
    the number of files written must match the number of filenames */
QStringList MoleculeParser::write(const System &system,
                                  const QStringList &files,
                                  const PropertyMap &map)
{
    if (files.isEmpty())
    {
        throw SireError::io_error( QObject::tr(
                "You must supply a valid filename. An empty list is not sufficient!"),
                    CODELOC );
    }
    else if (files.count() == 1)
    {
        return MoleculeParser::write(system, files[0], map);
    }

    //build a list of filenames with their associated fileformats
    QStringList filenames;
    QStringList fileformats;

    const auto format_property = map["fileformat"];

    if (format_property.hasValue())
    {
        try
        {
            fileformats = format_property.value().asA<StringProperty>().toString().split(",");
        }
        catch(...)
        {
            fileformats.append( format_property.value().asA<MoleculeParser>()
                                               .formatName() );
        }

        if (files.count() != fileformats.count())
        {
            throw SireError::io_error( QObject::tr(
                    "You must match up the number of filenames to fileformats when "
                    "specifying both the filenames [%1] and fileformats [%2].")
                        .arg(filenames.join(",")).arg(fileformats.join(",")), CODELOC );
        }

        for (int i=0; i<fileformats.count(); ++i)
        {
            const QString filename = files[i];

            if (filename.isEmpty() or QFileInfo(filename).baseName().isEmpty())
            {
                throw SireError::io_error( QObject::tr(
                        "You must supply a valid filename. This '%1' is not sufficient.")
                            .arg(filename), CODELOC );
            }

            QString basename = QFileInfo(filename).completeBaseName();
            filenames.append( QString("%1.%2").arg(basename,fileformats[i].toLower()) );
        }
    }
    else
    {
        //we may need to find the format from the system
        try
        {
            fileformats = system.property(format_property).asA<StringProperty>()
                                                          .toString().split(",");
        }
        catch(...)
        {}

        for (int i=0; i<files.count(); ++i)
        {
            const auto filename = files[i];

            QString extension = QFileInfo(filename).completeSuffix();

            if (extension.isEmpty())
            {
                if (i >= fileformats.count())
                {
                    throw SireError::io_error( QObject::tr(
                            "Cannot work out the fileformat to use to write the System to "
                            "file '%1'. You need to either supply the format using the "
                            "'fileformat' property in the passed map, add this to the System "
                            "as its 'fileformat' property, or pass a filename with an extension "
                            "whose fileformat can be determined. Supported fileformats are;\n%2")
                                .arg(filename).arg(MoleculeParser::supportedFormats()),
                                    CODELOC );
                }
                else
                {
                    filenames.append( QString("%1.%2").arg(filename)
                                                      .arg(fileformats[i].toLower()) );
                }
            }
            else if (i >= fileformats.count())
            {
                filenames.append(filename);
                fileformats.append(extension.toUpper());
            }
            else
            {
                filenames.append(filename);
                fileformats[i] = extension.toUpper();
            }
        }
    }

    //now we have a list of filenames and associated formats, actually
    //write the files
    return ::pvt_write(system, filenames, fileformats, map);
}

/** Extension of MoleculeParser::write which allows you to specify two filenames.
    The same rules to locate the fileformats are now used, except now only two
    files are permitted to be written */
QStringList MoleculeParser::write(const System &system, const QString &file1,
                                  const QString &file2, const PropertyMap &map)
{
    QStringList filenames;
    filenames.append(file1);
    filenames.append(file2);
    return MoleculeParser::write(system, filenames, map);
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system, const QString &filename,
                                 const PropertyMap &map)
{
    return MoleculeParser::write(system, filename, map);
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system,
                                 const QString &file1, const QString &file2,
                                 const PropertyMap &map)
{
    return MoleculeParser::write(system, file1, file2, map);
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system,
                                 const QStringList &filenames,
                                 const PropertyMap &map)
{
    return MoleculeParser::write(system, filenames, map);
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
    // Construct a list of parsers.
    QList<MoleculeParserPtr> parsers({*this, MoleculeParserPtr(other)});

    // A list of supplementary parsers.
    QList<MoleculeParserPtr> supplementary;

    // Sort the parsers: lead, then follower.
    sortParsers(parsers, supplementary, map);

    if (parsers.count() == 0)
        return System();

    if (parsers.first().read().isBroken())
    {
        //all of the parsers must be broken
        QTextStream cout(stdout, QIODevice::WriteOnly);

        cout << QObject::tr("Unable to read the file. Errors are below.\n\n");

        QStringList filenames;

        for (const auto &parser : parsers)
        {
            cout << "\n\n" << parser.read().errorReport();
            filenames.append(parser.read().filename());
        }

        throw SireIO::parse_error(QObject::tr(
            "Unable to load the file: %1")
                .arg(filenames.join(", ")), CODELOC);
    }

    // Instantiate an empty system.
    System system;

    // No supplementary data.
    if (supplementary.count() == 0)
    {
        // Construct the system: leader, then follower.
        system = parsers[0].read().startSystem(map);
        parsers[1].read().addToSystem(system, map);
    }
    else
    {
        system = parsers[0].read()
            .startSystem(supplementary[0].read().lines(), map);

        parsers[1].read().addToSystem(system, map);
    }

    return system;
}

/** Return the System that is constructed from the information in the passed
    parsers. This will parse the information in order, meaning that data contained
    in earlier parsers may be overwritten by data from later parsers */
System MoleculeParser::toSystem(const QList<MoleculeParserPtr> &others,
                                const PropertyMap &map) const
{
    // Make a copy of the list of parsers.
    auto parsers = others;

    // Add this parser to the list.
    parsers.append(*this);

    // A list of supplementary parers.
    QList<MoleculeParserPtr> supplementary;

    // Sort the parsers: leader, then followers.
    sortParsers(parsers, supplementary, map);

    if (parsers.count() == 0)
        return System();

    if (parsers.first().read().isBroken())
    {
        //all of the parsers must be broken
        QTextStream cout(stdout, QIODevice::WriteOnly);

        cout << QObject::tr("Unable to read some of the files. Errors are below.\n");

        QStringList filenames;

        for (const auto &parser : parsers)
        {
            cout << parser.read().errorReport() << "\n\n";
            filenames.append(parser.read().filename());
        }

        throw SireIO::parse_error(QObject::tr(
            "Unable to load some of the files: %1")
                .arg(filenames.join(", ")), CODELOC);
    }

    // Instantiate an empty system.
    System system;

    // No supplementary data.
    if (supplementary.count() == 0)
    {
        // Construct the initial system from the leader.
        system = parsers[0].read().startSystem(map);

        // Add to the system, using properties parsed by the followers.
        for (int i=1; i<parsers.count(); ++i)
            parsers[i].read().addToSystem(system, map);
    }
    else
    {
        // The list of supplementary record lines.
        QVector<QString> supp_lines;

        // Combine all of the supplementary data into a single
        // lists of record lines.
        for (const auto &supp : supplementary)
            supp_lines += supp.read().lines();

        // Start a system using the lead parser, passing the
        // supplementary data records.
        system = parsers[0].read().startSystem(supp_lines, map);

        // Add to the system, using properties parsed by the followers.
        for (int i=1; i<parsers.count(); ++i)
            parsers[i].read().addToSystem(system, map);
    }

    return system;
}

/** Save the passed System to the file called 'filename'. First, the 'fileformat'
    property is looked at in 'map'. This is used to set the format(s) of
    the files that are written (comma-separated list).

    If this does not exist, then the extension of the
    file is used to work out which format to use. If no extension is given,
    then the System will be queried to find out its preferred format (normally
    the format it was loaded with), via its 'fileformat' property
    (again, comma separated list).

    If their preferred format results in multiple files, then
    multiple files will be written. This returns the full pathnames to
    all of the files that are written
*/
QStringList MoleculeParser::write(const System &system, const QString &filename)
{
    return write(system, filename, PropertyMap());
}

/** Extension of MoleculeParser::write which allows many filenames.
    The same rules to locate the fileformats are now used, except that now only
    the number of files written must match the number of filenames */
QStringList MoleculeParser::write(const System &system,
                                  const QStringList &files)
{
    return write(system, files, PropertyMap());
}

/** Extension of MoleculeParser::write which allows you to specify two filenames.
    The same rules to locate the fileformats are now used, except now only two
    files are permitted to be written */
QStringList MoleculeParser::write(const System &system, const QString &file1,
                                  const QString &file2)
{
    return write(system, file1, file2, PropertyMap());
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system, const QString &filename)
{
    return save(system, filename, PropertyMap());
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system,
                                 const QString &file1, const QString &file2)
{
    return save(system, file1, file2, PropertyMap());
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system,
                                 const QStringList &filenames)
{
    return save(system, filenames, PropertyMap());
}

/** Return the System that is constructed from the data in this parser */
System MoleculeParser::toSystem() const
{
    return toSystem(PropertyMap());
}

/** Return the System that is constructed from the data in the two
    passed parsers (i.e. representing a topology and a coordinate file) */
System MoleculeParser::toSystem(const MoleculeParser &other) const
{
    return toSystem(other, PropertyMap());
}

/** Return the System that is constructed from the information in the passed
    parsers. This will parse the information in order, meaning that data contained
    in earlier parsers may be overwritten by data from later parsers */
System MoleculeParser::toSystem(const QList<MoleculeParserPtr> &others) const
{
    return toSystem(others, PropertyMap());
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

/** Start creating a new System using the information contained in this parser
    and the supplementary records contained in 'lines', using the (optional)
    property map to name the properties */
System MoleculeParser::startSystem(const QVector<QString> &lines, const PropertyMap &map) const
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

/** Sort the parsers: Lead first, then followers. */
void MoleculeParser::sortParsers(QList<MoleculeParserPtr> &parsers,
                                 QList<MoleculeParserPtr> &supplementary,
                                 const PropertyMap &map) const
{
    /* Parsers can be leaders or followers. Leaders are capable of
       constructing an entire molecular system on their own, whereas
       followers cannot. However, certain lead parsers (PDB2, MOl2)
       are also able to follow.

       Here we sort the parsers into order, ready to construct a
       system. If there is more than one lead, then we need to work
       out which are capable of following.
     */

    // The leaders. We should end up with one.
    QList<MoleculeParserPtr> leaders;

    // The follower parsers.
    QList<MoleculeParserPtr> followers;

    // The broken parsers.
    QList<MoleculeParserPtr> broken;

    // First pass: work out leaders and followers.
    for (auto parser : parsers)
    {
        // This is a lead parser.
        if (parser.read().isBroken())
        {
            broken.append(parser);
        }
        else if (parser.read().isLead())
        {
            leaders.append(parser);
        }
        else
        {
            // This parser can follow.
            if (parser.read().canFollow())
            {
                followers.append(parser);
            }
            // This parser must contain data that
            // is supplementary to the lead.
            else if (parser->isA<Supplementary>())
            {
                supplementary.append(parser);
            }
            else
            {
                throw SireIO::parse_error(QObject::tr(
                    "A parser cannot lead or follow? %1"
                        ).arg(parser.read().toString()), CODELOC);
            }
        }
    }

    if (broken.count() > 0)
    {
        //everything is broken!
        parsers = broken;
        supplementary.clear();
        return;
    }

    // No leaders.
    if (leaders.count() == 0)
    {
        if (supplementary.count() > 0)
        {
            // likely a Supplementary is hiding a broken parser - reparse it...
            for (const auto &parser : supplementary)
            {
                PropertyMap m2(map);
                m2.set("DISABLE_SUPPLEMENTARY", BooleanProperty(true));

                auto p = _pvt_parse(parser.read().filename(), m2);

                if (p.read().isBroken())
                {
                    broken.append(p);
                }
            }

            if (broken.count() > 0)
            {
                parsers = broken;
                supplementary.clear();
                return;
            }
        }

        throw SireIO::parse_error( QObject::tr(
            "Unable to load any molecules from the files as none "
            "contain the necessary molecular information to create a "
            "system. Only coordinate or trajectory information "
            "has been loaded. Structure or topology information, "
            "e.g. as would be found in a topology file, is missing."), CODELOC );
    }

    // If there are multiple leaders, then check whether any can follow.
    // If so, move them to the followers until only a single lead remains.
    if (leaders.count() > 1)
    {
        // Whether the parser is the new leader.
        QVector<bool> is_leader(leaders.count());
        is_leader.fill(true);

        // The number of leaders.
        int num_lead = leaders.count();

        for (int i=0; i<leaders.count(); ++i)
        {
            // Make sure we have one leader.
            if (num_lead > 1)
            {
                // Add the leader to followers and flag that it's been removed.
                if (leaders[i].read().canFollow())
                {
                    followers.append(leaders[i]);
                    is_leader[i] = false;
                    num_lead--;
                }
            }
        }

        // Can only have one leader.
        if (num_lead > 1)
        {
            throw SireError::incompatible_error(QObject::tr(
                "Cannot construct a System from multiple lead parsers if "
                "none can follow!"), CODELOC);
        }

        // Find the new lead parser..
        for (int i=0; i<leaders.count(); ++i)
        {
            // This parser is the new leader.
            if (is_leader[i])
            {
                // Copy the ptr to the new lead.
                MoleculeParserPtr leader = leaders[i];

                // Clear the leaders list.
                leaders.clear();

                // Append the new leader.
                leaders.append(leader);

                break;
            }
        }
    }

    // Make the sorted list of parsers: leader, then followers.
    parsers = leaders + followers;
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

QDataStream &operator<<(QDataStream &ds, const NullParser &parser)
{
    writeHeader(ds, r_null, 1);
    ds << static_cast<const MoleculeParser&>(parser);
    return ds;
}

QDataStream &operator>>(QDataStream &ds, NullParser &parser)
{
    VersionID v = readHeader(ds, r_null);

    if (v == 1)
    {
        ds >> static_cast<MoleculeParser&>(parser);
    }
    else
        throw version_error(v, "1", r_null, CODELOC);

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

QString NullParser::formatName() const
{
    return "NULL";
}

/** Return a description of the file format */
QString NullParser::formatDescription() const
{
    return QObject::tr("Null parser that should not be used for any real parsing.");
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

/** Return this parser constructed from the passed filename */
MoleculeParserPtr NullParser::construct(const QString &filename,
                                        const PropertyMap &map) const
{
    throw SireError::program_bug( QObject::tr(
            "The NullParser should not be used for an real file IO!"), CODELOC );

    return MoleculeParserPtr();
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr NullParser::construct(const QStringList &lines,
                                        const PropertyMap &map) const
{
    throw SireError::program_bug( QObject::tr(
            "The NullParser should not be used for an real file IO!"), CODELOC );

    return MoleculeParserPtr();
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr NullParser::construct(const SireSystem::System &system,
                                        const PropertyMap &map) const
{
    throw SireError::program_bug( QObject::tr(
            "The NullParser should not be used for an real file IO!"), CODELOC );

    return MoleculeParserPtr();
}

//////////////
////////////// Implementation of BrokenParser
//////////////

static const RegisterMetaType<BrokenParser> r_broken;

QDataStream &operator<<(QDataStream &ds, const BrokenParser &parser)
{
    writeHeader(ds, r_broken, 1);

    SharedDataStream sds(ds);

    sds << parser.error_report << parser.suffix
        << static_cast<const MoleculeParser&>(parser);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, BrokenParser &parser)
{
    VersionID v = readHeader(ds, r_broken);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> parser.error_report >> parser.suffix
            >> static_cast<MoleculeParser&>(parser);
    }
    else
        throw version_error(v, "1", r_broken, CODELOC);

    return ds;
}

BrokenParser::BrokenParser() : ConcreteProperty<BrokenParser,MoleculeParser>()
{}

BrokenParser::BrokenParser(const QString &filename, const PropertyMap &map)
             : ConcreteProperty<BrokenParser,MoleculeParser>(filename, map)
{}

BrokenParser::BrokenParser(const QString &filename, const QString &s,
                           const QStringList &errors)
             : ConcreteProperty<BrokenParser,MoleculeParser>(filename, PropertyMap())
{
    suffix = s;
    error_report = errors;
}

BrokenParser::BrokenParser(const QString &filename, const QStringList &errors)
             : ConcreteProperty<BrokenParser,MoleculeParser>(filename, PropertyMap())
{
    error_report = errors;
}

BrokenParser::BrokenParser(const QStringList &lines, const PropertyMap &map)
             : ConcreteProperty<BrokenParser,MoleculeParser>(lines, map)
{}

BrokenParser::BrokenParser(const System&, const PropertyMap &map)
             : ConcreteProperty<BrokenParser,MoleculeParser>(map)
{}

BrokenParser::BrokenParser(const BrokenParser &other)
             : ConcreteProperty<BrokenParser,MoleculeParser>(other),
               error_report(other.error_report), suffix(other.suffix)
{}

BrokenParser::~BrokenParser()
{}

BrokenParser& BrokenParser::operator=(const BrokenParser &other)
{
    if (this != &other)
    {
        error_report = other.error_report;
        suffix = other.suffix;
        MoleculeParser::operator=(other);
    }

    return *this;
}

bool BrokenParser::operator==(const BrokenParser &other) const
{
    return MoleculeParser::operator==(other);
}

bool BrokenParser::operator!=(const BrokenParser &other) const
{
    return MoleculeParser::operator!=(other);
}

const char* BrokenParser::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BrokenParser>() );
}

QString BrokenParser::formatName() const
{
    return "BROKEN";
}

/** Return a description of the file format */
QString BrokenParser::formatDescription() const
{
    return QObject::tr("Broken parser used to report an unparseable file.");
}

bool BrokenParser::isBroken() const
{
    return true;
}

QString BrokenParser::errorReport() const
{
    if (suffix.isEmpty())
    {
        return QObject::tr("== %1 ==\n\nThis file was not recognised by any of the file parsers!\n\n"
                           "%2\n")
                    .arg(this->filename())
                    .arg(error_report.join("\n\n"));
    }
    else
    {
        return QObject::tr("== %1 ==\n\nThis file could not be parsed by any of the file parsers! "
                           "It was recognised as a file of type %2, but all parsers failed "
                           "to parse this file. The errors from the parsers associated "
                           "with the suffix %2 are printed below:\n\n"
                           "%3\n")
                    .arg(this->filename()).arg(suffix)
                    .arg(error_report.join("\n\n"));
    }
}

System BrokenParser::toSystem(const PropertyMap&) const
{
    return System();
}

System BrokenParser::toSystem(const MoleculeParser &other,
                              const PropertyMap&) const
{
    if (not other.isA<BrokenParser>())
        throw SireError::incompatible_error( QObject::tr(
                "Broken parsers cannot be combined with other parsers (%1)")
                    .arg(other.toString()), CODELOC );

    return System();
}

System BrokenParser::toSystem(const QList<MoleculeParserPtr> &others,
                              const PropertyMap&) const
{
    for (const auto &other : others)
    {
        if (not other->isA<BrokenParser>())
        {
            throw SireError::incompatible_error( QObject::tr(
                    "Broken parsers cannot be combined with other parsers (%1)")
                        .arg(other.read().toString()), CODELOC );
        }
    }

    return System();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr BrokenParser::construct(const QString &filename,
                                          const PropertyMap &map) const
{
    return MoleculeParserPtr(BrokenParser(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr BrokenParser::construct(const QStringList &lines,
                                          const PropertyMap &map) const
{
    return MoleculeParserPtr(BrokenParser(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr BrokenParser::construct(const SireSystem::System &system,
                                          const PropertyMap &map) const
{
    return MoleculeParserPtr(BrokenParser(system, map));
}
