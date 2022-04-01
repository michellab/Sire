/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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


#ifndef SIREIO_PDB2_H
#define SIREIO_PDB2_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"
#include "SireMol/atomvelocities.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class SDF;
}

namespace SireMol
{
class Atom;
class MolEditor;
class MoleculeInfoData;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::SDF&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::SDF&);

namespace SireIO
{

/** This class holds a parser for reading and writing
    Structure Data File (SDF) molecular file formats

    @author Christopher Woods
*/
class SIREIO_EXPORT SDF : public SireBase::ConcreteProperty<SDF,MoleculeParser>
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const SDF&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, SDF&);

public:
    SDF();
    SDF(const QString &filename,
        const PropertyMap &map = PropertyMap());

    SDF(const QStringList &lines,
        const PropertyMap &map = PropertyMap());
    SDF(const SireSystem::System &system,
        const PropertyMap &map = PropertyMap());

    SDF(const SDF &other);

    ~SDF();

    SDF& operator=(const SDF &other);

    bool operator==(const SDF &other) const;
    bool operator!=(const SDF &other) const;

    static const char* typeName();

    const char* what() const;

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

    bool isLead() const;

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const PropertyMap &map);

    /** Any warnings that were raised when reading the file. */
    QStringList parse_warnings;
};

}

Q_DECLARE_METATYPE( SireIO::SDF )

SIRE_EXPOSE_CLASS( SireIO::SDF )

SIRE_END_HEADER

#endif
