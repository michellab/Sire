/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SQUIRE_AM1BCC_H
#define SQUIRE_AM1BCC_H

#include "qmchargecalculator.h"
#include "mopac.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class AM1BCC;
}

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::AM1BCC&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::AM1BCC&);

namespace Squire
{

/** This charge calculator calculates AM1-BCC partial 
    atomic charges
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT AM1BCC 
        : public SireBase::ConcreteProperty<AM1BCC,QMChargeCalculator>
{

friend SQUIRE_EXPORT QDataStream& ::operator<<(QDataStream&, const AM1BCC&);
friend SQUIRE_EXPORT QDataStream& ::operator>>(QDataStream&, AM1BCC&);

public:
    AM1BCC();
    AM1BCC(const AM1BCC &other);
    
    ~AM1BCC();
    
    static const char* typeName();
    
    AM1BCC& operator=(const AM1BCC &other);
    
    bool operator==(const AM1BCC &other) const;
    bool operator!=(const AM1BCC &other) const;

    void setEnvironment(const QString &variable, const QString &value);

    const QHash<QString,QString>& environment() const;

    QString environment(const QString &variable) const;
    
    void setTotalCharge(int charge);
    int totalCharge() const;

    void setScaleFactor(double sclfac);

    AtomCharges operator()(const PartialMolecule &molecule,
                           const PropertyMap &map = PropertyMap()) const;

    bool mayChangeCharges(const PartialMolecule &oldmol,
                          const PartialMolecule &newmol,
                          const PropertyMap &map = PropertyMap()) const;

private:
    QString getAmberHome() const;

    AtomCharges convertAM1MullikenToAM1BCC(const AtomCharges &am1mulliken,
                                           const Molecule &molecule,
                                           const PropertyMap &map,
                                           const QString &amberhome) const;
    
    /** Currently, this uses Mopac to calculate the 
        AM1 Mulliken charges */
    Mopac mopac;
};

}

Q_DECLARE_METATYPE( Squire::AM1BCC )

SIRE_EXPOSE_CLASS( Squire::AM1BCC )

SIRE_END_HEADER

#endif
