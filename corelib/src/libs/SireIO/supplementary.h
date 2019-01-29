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

#ifndef SIREIO_SUPPLEMENTARY_H
#define SIREIO_SUPPLEMENTARY_H

#include "moleculeparser.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class Supplementary;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::Supplementary&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::Supplementary&);

namespace SireIO
{

/** A dunmmy parser class to hold data that is supplementary to a lead parser,
    e.g. data from a parameter file.

    @author Lester Hedges
*/
class SIREIO_EXPORT Supplementary : public SireBase::ConcreteProperty<Supplementary,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const Supplementary&);
friend QDataStream& ::operator>>(QDataStream&, Supplementary&);

public:
    Supplementary();
    Supplementary(const QString &filename,
         const PropertyMap &map = PropertyMap());

    Supplementary(const QStringList &lines,
         const PropertyMap &map = PropertyMap());

    Supplementary(const Supplementary &other);

    ~Supplementary();

    Supplementary& operator=(const Supplementary &other);

    bool operator==(const Supplementary &other) const;
    bool operator!=(const Supplementary &other) const;

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

    bool canFollow() const;

private:
    /** The name of the parsed file (if from a file). */
    QString filename;
};

}

Q_DECLARE_METATYPE( SireIO::Supplementary )

SIRE_EXPOSE_CLASS( SireIO::Supplementary )

SIRE_END_HEADER

#endif
