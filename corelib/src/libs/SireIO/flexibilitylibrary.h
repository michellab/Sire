/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREIO_FLEXIBILITYLIBRARY_H
#define SIREIO_FLEXIBILITYLIBRARY_H

#include "iobase.h"

#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class FlexibilityLibrary;
class FlexibilityTemplate;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::FlexibilityLibrary&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::FlexibilityLibrary&);

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::FlexibilityTemplate&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::FlexibilityTemplate&);

namespace SireMol
{
class MoleculeView;
}

namespace SireMove
{
class Flexibility;
}

namespace SireIO
{  

using SireMol::BondID;
using SireMol::AngleID;
using SireMol::DihedralID;
using SireMol::MoleculeView;

using SireMove::Flexibility;

using SireUnits::Dimension::Length;
using SireUnits::Dimension::Angle;

/** Internal class used to store the data describing a single flexibility template
    
    @author Julien Michel
*/
class SIREIO_EXPORT FlexibilityTemplate
{

friend QDataStream& ::operator<<(QDataStream&, const FlexibilityTemplate&);
friend QDataStream& ::operator>>(QDataStream&, FlexibilityTemplate&);

public:
    FlexibilityTemplate();
    FlexibilityTemplate(const QString &name);
    
    FlexibilityTemplate(const FlexibilityTemplate &other);
    
    ~FlexibilityTemplate();

    static const char* typeName();

    FlexibilityTemplate& operator=(const FlexibilityTemplate &other);

    bool operator==(const FlexibilityTemplate &other) const;
    bool operator!=(const FlexibilityTemplate &other) const;

    const QString getName();
    
    void setRotation(const Angle &angle);
    void setTranslation(const Length &distance);
    void setMaximumBondVar(int maxvar);
    void setMaximumAngleVar(int maxvar);
    void setMaximumDihedralVar(int maxvar);

    void setBondDelta(const BondID &bondid, const Length &delta);
    void setAngleDelta(const AngleID &angleid, const Angle &delta);
    void setDihedralDelta(const DihedralID &dihedralid, const Angle &delta);
    
    Length getTranslation() const;
    Angle getRotation() const;
    int getMaximumBondVar() const;
    int getMaximumAngleVar() const;
    int getMaximumDihedralVar() const;
    
    Length getDelta(const BondID &bondid) const;
    Angle getDelta(const AngleID &angleid) const;
    Angle getDelta(const DihedralID &dihedralid) const;
    
    const QHash<BondID,Length>& getBondDeltas() const;
    const QHash<AngleID,Angle>& getAngleDeltas() const;
    const QHash<DihedralID,Angle>& getDihedralDeltas() const;

private:
    QString name;
    Length translation;
    Angle rotation;
    qint32 maxbondvar;
    qint32 maxanglevar;
    qint32 maxdihedralvar;
    
    QHash<BondID,Length> bonds;
    QHash<AngleID,Angle> angles;
    QHash<DihedralID,Angle> dihedrals;
};

/** This class is used to read templates describing how a 
    molecule will be moved by an InternalMove and to create 
    the property "flexibility" for a molecule whose template is available

    @author Julien Michel
*/

class SIREIO_EXPORT FlexibilityLibrary 
        : public SireBase::ConcreteProperty<FlexibilityLibrary,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const SireIO::FlexibilityLibrary&);
friend QDataStream& ::operator>>(QDataStream&, SireIO::FlexibilityLibrary&);

public:
    FlexibilityLibrary();
    FlexibilityLibrary(const QString &file);
    
    FlexibilityLibrary(const FlexibilityLibrary &other);
    
    ~FlexibilityLibrary();
    
    static const char* typeName();
    
    FlexibilityLibrary& operator=(const FlexibilityLibrary &other);

    bool operator==(const FlexibilityLibrary &other) const;
    bool operator!=(const FlexibilityLibrary &other) const;

    void loadTemplates(const QString &file);

    FlexibilityLibrary& operator+=(const FlexibilityLibrary &other);

    FlexibilityLibrary operator+(const FlexibilityLibrary &other) const;
    
    void add(const FlexibilityLibrary &other);
    
    const FlexibilityTemplate& getTemplate(const QString &key);
    
    void setTemplate(const QString &key, const FlexibilityTemplate &tmplate);
    
    Flexibility getFlexibility(const MoleculeView &molecule) const;

private:
    /** The flexibility templates, indexed by molecule name*/
    QHash<QString,FlexibilityTemplate> templates; 
};

} // end of namespace SireIO

Q_DECLARE_METATYPE( SireIO::FlexibilityLibrary )
Q_DECLARE_METATYPE( SireIO::FlexibilityTemplate )

SIRE_EXPOSE_CLASS( SireIO::FlexibilityLibrary )
SIRE_EXPOSE_CLASS( SireIO::FlexibilityTemplate )

SIRE_END_HEADER

#endif
