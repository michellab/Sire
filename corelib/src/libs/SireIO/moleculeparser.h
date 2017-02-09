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

SIRE_BEGIN_HEADER

namespace SireIO
{
class MoleculeParser;
class NullParser;
}

QDataStream& operator<<(QDataStream&, const SireIO::MoleculeParser&);
QDataStream& operator>>(QDataStream&, SireIO::MoleculeParser&);

QDataStream& operator<<(QDataStream&, const SireIO::NullParser&);
QDataStream& operator>>(QDataStream&, SireIO::NullParser&);

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
    typedef std::function<MoleculeParserPtr (QString filename, const PropertyMap &map)>
            ParserFactory;

    /** Internal class used to help register parsers against extensions */
    class SIREIO_EXPORT RegisterParser
    {
    public:
        RegisterParser(const QString &extension, ParserFactory factory);
        ~RegisterParser();
    };
}

/** The base class of all molecule parsers */
class SIREIO_EXPORT MoleculeParser : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeParser&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeParser&);

public:
    MoleculeParser();
    MoleculeParser(const QString &filename);
    MoleculeParser(const MoleculeParser &other);

    virtual ~MoleculeParser();
    
    static const char* typeName();

    static const NullParser& null();

    virtual MoleculeParser* clone() const=0;

    const QVector<QString>& lines() const;

    static MoleculeParserPtr parse(const QString &filename,
                                   const PropertyMap &map = PropertyMap());
    
    static QList<MoleculeParserPtr> parse(const QStringList &filenames,
                                          const PropertyMap &map = PropertyMap());

    static SireSystem::System read(const QString &filename,
                                   const PropertyMap &map = PropertyMap());
    static SireSystem::System read(const QString &file1, const QString &file2,
                                   const PropertyMap &map = PropertyMap());
    static SireSystem::System read(const QStringList &filenames,
                                   const PropertyMap &map = PropertyMap());

    virtual bool isLead() const;

    double score() const;

    SireSystem::System toSystem(const PropertyMap &map = PropertyMap()) const;
    
    SireSystem::System toSystem(const MoleculeParser &other,
                                const PropertyMap &map = PropertyMap()) const;
    
    SireSystem::System toSystem(const QList<MoleculeParserPtr> &others,
                                const PropertyMap &map = PropertyMap()) const;
    
    virtual void write(const QString &filename) const;
    
    virtual bool isTextFile() const;
    virtual bool isBinaryFile() const;
    
    virtual bool isEmpty() const;
    
protected:
    MoleculeParser& operator=(const MoleculeParser &other);
    
    bool operator==(const MoleculeParser &other) const;
    bool operator!=(const MoleculeParser &other) const;

    void setScore(double score);

    virtual SireSystem::System startSystem(const PropertyMap &map) const;
    virtual void addToSystem(SireSystem::System &system,
                             const PropertyMap &map) const;

private:
    static MoleculeParserPtr _pvt_parse(const QString &filename, const PropertyMap &map);

    friend class detail::RegisterParser;
    static void registerParser(const QString &extension, detail::ParserFactory factory);

    /** All of the lines in the file */
    QVector<QString> lnes;
    
    /** The score associated with the parser. The higher the score,
        the better the file was parsed */
    double scr;
};

/** This is a null parser, returned when the file cannot be parsed */
class SIREIO_EXPORT NullParser : public SireBase::ConcreteProperty<NullParser,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const NullParser&);
friend QDataStream& ::operator>>(QDataStream&, NullParser&);

public:
    NullParser();
    NullParser(const NullParser&);
    ~NullParser();
    
    NullParser& operator=(const NullParser&);
    
    bool operator==(const NullParser&) const;
    bool operator!=(const NullParser&) const;
    
    static const char* typeName();

    SireSystem::System toSystem(const PropertyMap &map = PropertyMap()) const;
    SireSystem::System toSystem(const MoleculeParser &other,
                                const PropertyMap &map = PropertyMap()) const;
    SireSystem::System toSystem(const QList<MoleculeParserPtr> &others,
                                const PropertyMap &map = PropertyMap()) const;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Set the score from parsing */
inline void MoleculeParser::setScore(double score)
{
    scr = score;
}

/** Return the score from parsing. The higher the score, the better the
    parsing. You can use this to decide which parser thought it was better
    at parsing the file. The score should be a value of between
    nlines * [0 (failed to parse line) -> 1 (parsed line perfectly)], 
    where the maximum score of nlines means that every line of the 
    file was parsed perfectly. */
inline double MoleculeParser::score() const
{
    return scr;
}

/** Return the lines of the file. Note that this
    only returns something for text-based files */
inline const QVector<QString>& MoleculeParser::lines() const
{
    return lnes;
}

/** Return whether or not this parser works with text files */
inline bool MoleculeParser::isTextFile() const
{
    return true;
}

/** Return whether or not this parser works with binary files */
inline bool MoleculeParser::isBinaryFile() const
{
    return not this->isTextFile();
}

/** Return whether or not this parser is empty */
inline bool MoleculeParser::isEmpty() const
{
    return this->isTextFile() and lnes.isEmpty();
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

#define SIRE_REGISTER_PARSER( X, Y ) static SireIO::detail::RegisterParser \
_SireIO_detail_RegisterParser_##X( #X, [](QString filename, const PropertyMap &map) \
                  { return SireIO::MoleculeParserPtr(Y(filename,map));} );

Q_DECLARE_METATYPE( SireIO::NullParser )

SIRE_EXPOSE_CLASS( SireIO::MoleculeParser )
SIRE_EXPOSE_CLASS( SireIO::NullParser )

SIRE_EXPOSE_PROPERTY( SireIO::MoleculeParserPtr, SireIO::MoleculeParser )


#endif
