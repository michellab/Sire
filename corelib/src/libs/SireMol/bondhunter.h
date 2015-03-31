/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREMOL_BONDHUNTER_H
#define SIREMOL_BONDHUNTER_H

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class BondHunter;

class NullBondHunter;
class CovalentBondHunter;
class ChemicalBondHunter;
}

QDataStream& operator<<(QDataStream&, const SireMol::BondHunter&);
QDataStream& operator>>(QDataStream&, SireMol::BondHunter&);

QDataStream& operator<<(QDataStream&, const SireMol::NullBondHunter&);
QDataStream& operator>>(QDataStream&, SireMol::NullBondHunter&);

QDataStream& operator<<(QDataStream&, const SireMol::CovalentBondHunter&);
QDataStream& operator>>(QDataStream&, SireMol::CovalentBondHunter&);

QDataStream& operator<<(QDataStream&, const SireMol::ChemicalBondHunter&);
QDataStream& operator>>(QDataStream&, SireMol::ChemicalBondHunter&);

namespace SireMol
{

class Connectivity;
class MoleculeView;

using SireBase::PropertyMap;

/** Base class of all functions used to hunt for bonds in a molecule

    @author Christopher Woods
*/
class SIREMOL_EXPORT BondHunter : public SireBase::Property
{
public:
    BondHunter();
    BondHunter(const BondHunter &other);
    
    virtual ~BondHunter();
    
    static const char* typeName()
    {
        return "SireMol::BondHunter";
    }

    virtual BondHunter* clone() const=0;
    
    /** Return the connectivity of the molecule viewed in 'molview' 
        using this function to hunt for all of the bonded atoms.
        This only searches for the bonds between atoms that are
        part of this view */
    virtual Connectivity operator()(const MoleculeView &molview,
                                    const PropertyMap &map = PropertyMap()) const=0;

    static const NullBondHunter& null();
};

class SIREMOL_EXPORT CovalentBondHunterParameters
{
public:
    CovalentBondHunterParameters()
    {}
    
    ~CovalentBondHunterParameters()
    {}
    
    static const SireBase::PropertyName& coordinates()
    {
        return coords_param;
    }
    
    static const SireBase::PropertyName& element()
    {
        return elements_param;
    }

private:
    static SireBase::PropertyName coords_param;
    static SireBase::PropertyName elements_param;
};

/** This is a bond hunter that finds bonded atoms by comparing
    the distance between the atoms to the sum of their covalent radii.
    If the distance is less than the sum of the covalent radii then the atoms
    are bonded.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT CovalentBondHunter
           : public SireBase::ConcreteProperty<CovalentBondHunter,BondHunter>
{

friend QDataStream& ::operator<<(QDataStream&, const CovalentBondHunter&);
friend QDataStream& ::operator>>(QDataStream&, CovalentBondHunter&);

public:
    CovalentBondHunter(double tolerance=1.1);
    CovalentBondHunter(const CovalentBondHunter &other);
    
    ~CovalentBondHunter();
    
    static const char* typeName();

    CovalentBondHunter& operator=(const CovalentBondHunter &other);
    
    bool operator==(const CovalentBondHunter &other) const;
    bool operator!=(const CovalentBondHunter &other) const;
    
    static const CovalentBondHunterParameters& parameters()
    {
        return params;
    }
    
    Connectivity operator()(const MoleculeView &molview,
                            const PropertyMap &map = PropertyMap()) const;

private:
    static CovalentBondHunterParameters params;

    /** The tolerance added to the sum of vdw radii when hunting for bonds */
    double tol;
};

/** This is a null bond hunter. This finds no bonds in a molecule and is used
    when you want to create an empty connectivity object for a molecule
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT NullBondHunter
           : public SireBase::ConcreteProperty<NullBondHunter,BondHunter>
{
public:
    NullBondHunter();
    NullBondHunter(const NullBondHunter &other);
    
    ~NullBondHunter();
    
    static const char* typeName();

    const char* what() const;
    
    Connectivity operator()(const MoleculeView &molview,
                            const PropertyMap &map = PropertyMap()) const;
};

/** This is a bond hunter that hunts for bonds using the distance between
    atoms (and comparing that distance against the sum of atomic covalent
    radii), but then it runs over each atom and ensures that the atom does
    not contain too many bonds. If it does, then only the n closest bonds
    are retained.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ChemicalBondHunter 
           : public SireBase::ConcreteProperty<ChemicalBondHunter,CovalentBondHunter>
{
public:
    ChemicalBondHunter();
    ChemicalBondHunter(double tolerance);
    ChemicalBondHunter(const ChemicalBondHunter &other);
    
    ~ChemicalBondHunter();
    
    static const char* typeName();
    
    Connectivity operator()(const MoleculeView &molview,
                            const PropertyMap &map = PropertyMap()) const;
};

typedef SireBase::PropPtr<BondHunter> BondHunterPtr;

}

Q_DECLARE_METATYPE( SireMol::NullBondHunter )
Q_DECLARE_METATYPE( SireMol::CovalentBondHunter )
Q_DECLARE_METATYPE( SireMol::ChemicalBondHunter )

SIRE_EXPOSE_CLASS( SireMol::BondHunter )
SIRE_EXPOSE_CLASS( SireMol::NullBondHunter )
SIRE_EXPOSE_CLASS( SireMol::CovalentBondHunter )
SIRE_EXPOSE_CLASS( SireMol::ChemicalBondHunter )
SIRE_EXPOSE_CLASS( SireMol::CovalentBondHunterParameters )

SIRE_EXPOSE_PROPERTY( SireMol::BondHunterPtr, SireMol::BondHunter )

SIRE_END_HEADER

#endif
