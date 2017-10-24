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

#ifndef SIREIO_CHARMMPSF_H
#define SIREIO_CHARMMPSF_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
}

namespace SireMol
{
class Atom;
class MolEditor;
class MoleculeInfoData;
class MoleculeView;
class Residue;
}

namespace SireIO
{

/** This class holds a parser for reading and writing Charmm PSF files.

    @author Lester Hedges
*/
class SIREIO_EXPORT CharmmPSF : public SireBase::ConcreteProperty<CharmmPSF,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const CharmmPSF&);
friend QDataStream& ::operator>>(QDataStream&, CharmmPSF&);

public:
    CharmmPSF();
    CharmmPSF(const QString &filename,
         const PropertyMap &map = PropertyMap());

    CharmmPSF(const QStringList &lines,
         const PropertyMap &map = PropertyMap());
    CharmmPSF(const SireSystem::System &system,
         const PropertyMap &map = PropertyMap());

    CharmmPSF(const CharmmPSF &other);

    ~CharmmPSF();

    CharmmPSF& operator=(const CharmmPSF &other);

    bool operator==(const CharmmPSF &other) const;
    bool operator!=(const CharmmPSF &other) const;

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

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const PropertyMap &map);

    SireMol::MolStructureEditor getMolStructure(int imol,
        const SireBase::PropertyName &cutting) const;

    SireMol::MolEditor getMolecule(int imol,
        const PropertyMap &map = PropertyMap()) const;

    /** The name of the parsed file (if from a file). */
    QString filename;

    /** Any warnings that were raised when reading the file. */
    QStringList parse_warnings;
};

}

SIRE_END_HEADER

#endif
