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

#ifndef SIREIO_GROTOP_H
#define SIREIO_GROTOP_H

#include "moleculeparser.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class GroTop;
}

QDataStream& operator<<(QDataStream&, const SireIO::GroTop&);
QDataStream& operator>>(QDataStream&, SireIO::GroTop&);

namespace SireIO
{

/** This class holds a parser for reading and writing Gromacs "top" topology files.

    @author Christopher Woods
*/
class SIREIO_EXPORT GroTop : public SireBase::ConcreteProperty<GroTop,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const GroTop&);
friend QDataStream& ::operator>>(QDataStream&, GroTop&);

public:
    GroTop();
    GroTop(const QString &filename,
           const PropertyMap &map = PropertyMap());

    GroTop(const QStringList &lines,
           const PropertyMap &map = PropertyMap());
    GroTop(const SireSystem::System &system,
           const PropertyMap &map = PropertyMap());
    
    GroTop(const GroTop &other);
    
    ~GroTop();
    
    GroTop& operator=(const GroTop &other);
    
    bool operator==(const GroTop &other) const;
    bool operator!=(const GroTop &other) const;
    
    static const char* typeName();
    
    const char* what() const;

    QStringList includePath(bool absolute_paths=false) const;
    QStringList includedFiles(bool absolute_paths=false) const;

    MoleculeParserPtr construct(const QString &filename,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const QStringList &lines,
                                const PropertyMap &map) const;
    
    MoleculeParserPtr construct(const SireSystem::System &system,
                                const PropertyMap &map) const;

    QString toString() const;

    QString formatName() const;
    QString formatDescription() const;
    QStringList formatSuffix() const;

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

    void getIncludePath(const PropertyMap &map);

    QString findIncludeFile(QString filename, QString current_directory);

    QVector<QString> loadInclude(QString filename, QString current_directory);

    QVector<QString> preprocess(const QVector<QString> &lines,
                                QHash<QString,QString> defines = QHash<QString,QString>(),
                                const QString &current_directory = ".");

private:
    void assertSane() const;
    void parseLines(const PropertyMap &map);
    
    /** This is the full search path of all directories that should
        be searched for Gromacs include files */
    QStringList include_path;
    
    /** This is the set of files that had to be included as part of parsing
        this file, together with the full path in which they were found */
    QHash<QString,QString> included_files;
};

}

Q_DECLARE_METATYPE( SireIO::GroTop )

SIRE_EXPOSE_CLASS( SireIO::GroTop )

SIRE_END_HEADER

#endif
