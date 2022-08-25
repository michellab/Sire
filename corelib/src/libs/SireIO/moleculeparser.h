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

#ifndef SIREIO_MOLECULEPARSER_H
#define SIREIO_MOLECULEPARSER_H

#include "SireBase/propertymap.h"

#include "SireMol/residuecutting.h"

#include <functional>
#include <memory>

SIRE_BEGIN_HEADER

namespace SireIO
{
class MoleculeParser;
class NullParser;
class BrokenParser;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::MoleculeParser&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::MoleculeParser&);

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::NullParser&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::NullParser&);

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::BrokenParser&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::BrokenParser&);

namespace SireSystem
{
class System;
}

namespace SireIO
{

using SireBase::PropertyMap;

typedef SireBase::PropPtr<MoleculeParser> MoleculeParserPtr;

namespace detail
{
    /** Base class of the parser registers */
    class SIREIO_EXPORT ParserFactoryHelper
    {
    public:
        ParserFactoryHelper();
        ParserFactoryHelper(MoleculeParser *parser);
        ParserFactoryHelper(const ParserFactoryHelper &other);

        ~ParserFactoryHelper();

        bool isValid() const;

        QString formatName() const;
        QString formatDescription() const;
        QString preferredSuffix() const;

        QStringList suffixes() const;

        QString toString() const;

        MoleculeParserPtr construct(const QString &filename,
                                    const PropertyMap &map) const;

        MoleculeParserPtr construct(const QStringList &lines,
                                    const PropertyMap &map) const;

        MoleculeParserPtr construct(const SireSystem::System &system,
                                    const PropertyMap &map) const;

        ParserFactoryHelper& operator=(const ParserFactoryHelper &other);

        bool operator<(const ParserFactoryHelper &other) const;
        bool operator>(const ParserFactoryHelper &other) const;

        bool operator==(const ParserFactoryHelper &other) const;
        bool operator!=(const ParserFactoryHelper &other) const;

        bool operator<=(const ParserFactoryHelper &other) const;
        bool operator>=(const ParserFactoryHelper &other) const;

    private:
        std::shared_ptr<MoleculeParser> parser;
    };
}

/** Class used as part of the parser registration system */
template<class T>
class SIREIO_EXPORT RegisterParser
{
public:
    RegisterParser()
    {
        detail::ParserFactoryHelper helper( new T() );
    }

    ~RegisterParser()
    {}
};

/** The base class of all molecule parsers */
class SIREIO_EXPORT MoleculeParser : public SireBase::Property
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const MoleculeParser&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, MoleculeParser&);

public:
    MoleculeParser(const PropertyMap &map = PropertyMap());
    MoleculeParser(const QString &filename, const PropertyMap &map);
    MoleculeParser(const QStringList &lines, const PropertyMap &map);

    MoleculeParser(const MoleculeParser &other);

    virtual ~MoleculeParser();

    static const char* typeName();

    static const NullParser& null();

    virtual MoleculeParser* clone() const=0;

    const QVector<QString>& lines() const;

    virtual MoleculeParserPtr construct(const QString &filename,
                                        const PropertyMap &map) const=0;

    virtual MoleculeParserPtr construct(const QStringList &lines,
                                        const PropertyMap &map) const=0;

    virtual MoleculeParserPtr construct(const SireSystem::System &system,
                                        const PropertyMap &map) const=0;

    static QString supportedFormats();

    static MoleculeParserPtr parse(const QString &filename,
                                   const PropertyMap &map);

    static QList<MoleculeParserPtr> parse(const QStringList &filenames,
                                          const PropertyMap &map);

    static SireSystem::System read(const QString &filename,
                                   const PropertyMap &map);
    static SireSystem::System read(const QString &file1, const QString &file2,
                                   const PropertyMap &map);
    static SireSystem::System read(const QStringList &filenames,
                                   const PropertyMap &map);

    static SireSystem::System load(const QString &filename,
                                   const PropertyMap &map);
    static SireSystem::System load(const QString &file1, const QString &file2,
                                   const PropertyMap &map);
    static SireSystem::System load(const QStringList &filenames,
                                   const PropertyMap &map);

    static QStringList write(const SireSystem::System &system,
                             const QString &filename,
                             const PropertyMap &map);

    static QStringList write(const SireSystem::System &system,
                             const QString &file1, const QString &file2,
                             const PropertyMap &map);

    static QStringList write(const SireSystem::System &system,
                             const QStringList &filenames,
                             const PropertyMap &map);

    static QStringList save(const SireSystem::System &system,
                            const QString &filename,
                            const PropertyMap &map);

    static QStringList save(const SireSystem::System &system,
                            const QString &file1, const QString &file2,
                            const PropertyMap &map);

    static QStringList save(const SireSystem::System &system,
                            const QStringList &filenames,
                            const PropertyMap &map);

    static MoleculeParserPtr parse(const QString &filename);

    static QList<MoleculeParserPtr> parse(const QStringList &filenames);

    static SireSystem::System read(const QString &filename);
    static SireSystem::System read(const QString &file1, const QString &file2);
    static SireSystem::System read(const QStringList &filenames);

    static SireSystem::System load(const QString &filename);
    static SireSystem::System load(const QString &file1, const QString &file2);
    static SireSystem::System load(const QStringList &filenames);

    static QStringList write(const SireSystem::System &system,
                             const QString &filename);

    static QStringList write(const SireSystem::System &system,
                             const QString &file1, const QString &file2);

    static QStringList write(const SireSystem::System &system,
                             const QStringList &filenames);

    static QStringList save(const SireSystem::System &system,
                            const QString &filename);

    static QStringList save(const SireSystem::System &system,
                            const QString &file1, const QString &file2);

    static QStringList save(const SireSystem::System &system,
                            const QStringList &filenames);

    virtual bool isLead() const;
    virtual bool canFollow() const;

    /** Return the unique name of this format. This should be a short
        string in capital letters that will identify the format within
        the program, e.g. "PRM7", "RST7", "PDB3" etc. */
    virtual QString formatName() const=0;
    virtual QString formatDescription() const=0;

    virtual QStringList formatSuffix() const;

    QString filename() const;

    double score() const;

    void enableParallel();
    void disableParallel();
    void setUseParallel(bool on);
    bool usesParallel() const;

    SireSystem::System toSystem(const PropertyMap &map) const;

    SireSystem::System toSystem(const MoleculeParser &other,
                                const PropertyMap &map) const;

    SireSystem::System toSystem(const QList<MoleculeParserPtr> &others,
                                const PropertyMap &map) const;

    SireSystem::System toSystem() const;

    SireSystem::System toSystem(const MoleculeParser &other) const;

    SireSystem::System toSystem(const QList<MoleculeParserPtr> &others) const;

    virtual void writeToFile(const QString &filename) const;

    virtual bool isTextFile() const;
    virtual bool isBinaryFile() const;

    virtual bool isEmpty() const;

protected:
    MoleculeParser& operator=(const MoleculeParser &other);

    bool operator==(const MoleculeParser &other) const;
    bool operator!=(const MoleculeParser &other) const;

    void setScore(double score);

    void removeCommentLines(const QString &comment_flag);

    virtual SireSystem::System startSystem(const PropertyMap &map) const;
    virtual SireSystem::System startSystem(const QVector<QString> &lines,
                                           const PropertyMap &map) const;
    virtual void addToSystem(SireSystem::System &system,
                             const PropertyMap &map) const;

    void setLines(const QVector<QString> &lines);

    static QVector<QString> readTextFile(QString filename);

    virtual SireBase::PropertyPtr getForceField(const SireSystem::System &system,
                                                const PropertyMap &map) const;

private:
    static MoleculeParserPtr _pvt_parse(const QString &filename, const PropertyMap &map);

    void sortParsers(QList<MoleculeParserPtr>& parsers,
                     QList<MoleculeParserPtr>& supplementary) const;

    /** The name of the file (if one was read) */
    QString fname;

    /** All of the lines in the file */
    QVector<QString> lnes;

    /** The score associated with the parser. The higher the score,
        the better the file was parsed */
    double scr;

    /** Whether or not to run in parallel */
    bool run_parallel;
};

/** This is a null parser */
class SIREIO_EXPORT NullParser : public SireBase::ConcreteProperty<NullParser,MoleculeParser>
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const NullParser&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, NullParser&);

public:
    NullParser();
    NullParser(const NullParser&);
    ~NullParser();

    NullParser& operator=(const NullParser&);

    bool operator==(const NullParser&) const;
    bool operator!=(const NullParser&) const;

    static const char* typeName();

    QString formatName() const;
    QString formatDescription() const;

    MoleculeParserPtr construct(const QString &filename,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const QStringList &lines,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const SireSystem::System &system,
                                const PropertyMap &map) const;

    SireSystem::System toSystem(const PropertyMap &map = PropertyMap()) const;
    SireSystem::System toSystem(const MoleculeParser &other,
                                const PropertyMap &map = PropertyMap()) const;
    SireSystem::System toSystem(const QList<MoleculeParserPtr> &others,
                                const PropertyMap &map = PropertyMap()) const;
};

/** This is a broken parser, returned when the file cannot be parsed */
class SIREIO_EXPORT BrokenParser : public SireBase::ConcreteProperty<BrokenParser,MoleculeParser>
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const BrokenParser&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, BrokenParser&);

public:
    BrokenParser();
    BrokenParser(const BrokenParser &other);
    ~BrokenParser();

    BrokenParser& operator=(const BrokenParser&);

    bool operator==(const BrokenParser&) const;
    bool operator!=(const BrokenParser&) const;

    static const char* typeName();

    QStringList errorReport() const;

    QString formatName() const;
    QString formatDescription() const;

    MoleculeParserPtr construct(const QString &filename,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const QStringList &lines,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const SireSystem::System &system,
                                const PropertyMap &map) const;

    SireSystem::System toSystem(const PropertyMap &map = PropertyMap()) const;
    SireSystem::System toSystem(const MoleculeParser &other,
                                const PropertyMap &map = PropertyMap()) const;
    SireSystem::System toSystem(const QList<MoleculeParserPtr> &others,
                                const PropertyMap &map = PropertyMap()) const;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Set the score from parsing */
SIRE_ALWAYS_INLINE void MoleculeParser::setScore(double score)
{
    scr = score;
}

/** Return the score from parsing. The higher the score, the better the
    parsing. You can use this to decide which parser thought it was better
    at parsing the file. The score should be a value of between
    nlines * [0 (failed to parse line) -> 1 (parsed line perfectly)],
    where the maximum score of nlines means that every line of the
    file was parsed perfectly. */
SIRE_ALWAYS_INLINE double MoleculeParser::score() const
{
    return scr;
}

/** Return the lines of the file. Note that this
    only returns something for text-based files */
SIRE_ALWAYS_INLINE const QVector<QString>& MoleculeParser::lines() const
{
    return lnes;
}

/** Return whether or not this parser works with text files */
SIRE_ALWAYS_INLINE bool MoleculeParser::isTextFile() const
{
    return true;
}

/** Return whether or not this parser works with binary files */
SIRE_ALWAYS_INLINE bool MoleculeParser::isBinaryFile() const
{
    return not this->isTextFile();
}

/** Return whether or not this parser is empty */
SIRE_ALWAYS_INLINE bool MoleculeParser::isEmpty() const
{
    return this->isTextFile() and lnes.isEmpty();
}

/** Return whether or not this parser runs in parallel */
SIRE_ALWAYS_INLINE bool MoleculeParser::usesParallel() const
{
    return run_parallel;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireIO::NullParser )
Q_DECLARE_METATYPE( SireIO::BrokenParser )

SIRE_EXPOSE_CLASS( SireIO::MoleculeParser )
SIRE_EXPOSE_CLASS( SireIO::BrokenParser )
SIRE_EXPOSE_CLASS( SireIO::NullParser )

SIRE_EXPOSE_PROPERTY( SireIO::MoleculeParserPtr, SireIO::MoleculeParser )


#endif
