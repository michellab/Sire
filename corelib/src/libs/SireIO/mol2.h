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

#ifndef SIREIO_MOL2_H
#define SIREIO_MOL2_H

#include "moleculeparser.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class Mol2;
}

QDataStream& operator<<(QDataStream&, const SireIO::Mol2&);
QDataStream& operator>>(QDataStream&, SireIO::Mol2&);

namespace SireIO
{

/** This class holds a parser for reading and writing Tripos Mol2 files.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2 : public SireBase::ConcreteProperty<Mol2,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const Mol2&);
friend QDataStream& ::operator>>(QDataStream&, Mol2&);

public:
    Mol2();
    Mol2(const QString &filename,
         const PropertyMap &map = PropertyMap());

    Mol2(const QStringList &lines,
         const PropertyMap &map = PropertyMap());
    Mol2(const SireSystem::System &system,
         const PropertyMap &map = PropertyMap());
    
    Mol2(const Mol2 &other);
    
    ~Mol2();
    
    Mol2& operator=(const Mol2 &other);
    
    bool operator==(const Mol2 &other) const;
    bool operator!=(const Mol2 &other) const;
    
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
};

}

Q_DECLARE_METATYPE( SireIO::Mol2 )

SIRE_EXPOSE_CLASS( SireIO::Mol2 )

SIRE_END_HEADER

#endif
