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

#ifndef SQUIRE_QMCHARGECALCULATOR_H
#define SQUIRE_QMCHARGECALCULATOR_H

#include <QHash>

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireMol/atomcharges.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class QMChargeCalculator;
class NullQMChargeCalculator;
}

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::QMChargeCalculator&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::QMChargeCalculator&);

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::NullQMChargeCalculator&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::NullQMChargeCalculator&);

namespace SireMol
{
class PartialMolecule;
class Molecules;
class MolNum;
}

namespace Squire
{

using SireMol::PartialMolecule;
using SireMol::Molecules;
using SireMol::AtomCharges;
using SireMol::MolNum;

using SireBase::PropertyMap;

/** This is the base class of all functions which are used
    to calculate atomic partial charges of molecules from
    an underlying QM Hamiltonian
        
    @author Christopher Woods
*/
class SQUIRE_EXPORT QMChargeCalculator : public SireBase::Property
{

friend SQUIRE_EXPORT QDataStream& ::operator<<(QDataStream&, const QMChargeCalculator&);
friend SQUIRE_EXPORT QDataStream& ::operator>>(QDataStream&, QMChargeCalculator&);

public:
    QMChargeCalculator();
    
    QMChargeCalculator(const QMChargeCalculator &other);
    
    virtual ~QMChargeCalculator();
    
    static const char* typeName();
    
    virtual QMChargeCalculator* clone() const=0;

    virtual void setScaleFactor(double sclfactor);
    
    double scaleFactor() const;
    
    virtual AtomCharges operator()(const PartialMolecule &molecule,
                                   const PropertyMap &map = PropertyMap()) const=0;
    
    AtomCharges calculate(const PartialMolecule &molecule,
                          const PropertyMap &map = PropertyMap()) const;

    virtual bool mayChangeCharges(const PartialMolecule &oldmol,
                                  const PartialMolecule &newmol,
                                  const PropertyMap &map = PropertyMap()) const=0;

    static const NullQMChargeCalculator& null();
    
protected:
    QMChargeCalculator& operator=(const QMChargeCalculator &other);
    
    bool operator==(const QMChargeCalculator &other) const;
    bool operator!=(const QMChargeCalculator &other) const;

private:
    /** The scale factor that multiples the calculated charges */
    double sclfac;
};

/** This is a null charge calculator - this returns zero 
    charges for every molecule! */
class SQUIRE_EXPORT NullQMChargeCalculator
            : public SireBase::ConcreteProperty<NullQMChargeCalculator,QMChargeCalculator>
{

friend SQUIRE_EXPORT QDataStream& ::operator<<(QDataStream&, const NullQMChargeCalculator&);
friend SQUIRE_EXPORT QDataStream& ::operator>>(QDataStream&, NullQMChargeCalculator&);

public:
    NullQMChargeCalculator();
    NullQMChargeCalculator(const NullQMChargeCalculator &other);
    
    ~NullQMChargeCalculator();
    
    static const char* typeName();
    
    NullQMChargeCalculator& operator=(const NullQMChargeCalculator &other);
    
    bool operator==(const NullQMChargeCalculator &other) const;
    bool operator!=(const NullQMChargeCalculator &other) const;
        
    AtomCharges operator()(const PartialMolecule &molecule,
                           const PropertyMap &map = PropertyMap()) const;

    bool mayChangeCharges(const PartialMolecule&, const PartialMolecule&,
                          const PropertyMap&) const;
};

typedef SireBase::PropPtr<QMChargeCalculator> QMChargeCalculatorPtr;

}

Q_DECLARE_METATYPE( Squire::NullQMChargeCalculator )

SIRE_EXPOSE_CLASS( Squire::QMChargeCalculator )
SIRE_EXPOSE_CLASS( Squire::NullQMChargeCalculator )

SIRE_EXPOSE_PROPERTY( Squire::QMChargeCalculatorPtr, Squire::QMChargeCalculator )

SIRE_END_HEADER

#endif
