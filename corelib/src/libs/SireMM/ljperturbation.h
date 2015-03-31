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

#ifndef SIREMOL_LJPERTURBATION_H
#define SIREMOL_LJPERTURBATION_H

#include "SireMol/perturbation.h"

#include "ljparameter.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class LJPerturbation;
}

QDataStream& operator<<(QDataStream&, const SireMM::LJPerturbation&);
QDataStream& operator>>(QDataStream&, SireMM::LJPerturbation&);

namespace SireMM
{

using SireBase::PropertyMap;

/** This is a perturbation that maps LJ parameters for a molecule
    from an initial to a final state
    
    @author Christopher Woods
*/
class SIREMM_EXPORT LJPerturbation 
         : public SireBase::ConcreteProperty<LJPerturbation,SireMol::Perturbation>
{

friend QDataStream& ::operator<<(QDataStream&, const LJPerturbation&);
friend QDataStream& ::operator>>(QDataStream&, LJPerturbation&);

public:
    enum MapType { MAP_SIGMA_AND_EPSILON = 1,
                   MAP_RMIN_AND_EPSILON = 2,
                   MAP_A_AND_B = 3 };

    LJPerturbation(const PropertyMap &map = PropertyMap());
    LJPerturbation(MapType maptype, const PropertyMap &map = PropertyMap());
    
    LJPerturbation(const SireCAS::Expression &mapping_function,
                   const PropertyMap &map = PropertyMap());
    
    LJPerturbation(const SireCAS::Expression &mapping_function,
                   MapType maptype,
                   const PropertyMap &map = PropertyMap());
    
    LJPerturbation(const SireCAS::Expression &sigma_mapping_function,
                   const SireCAS::Expression &epsilon_mapping_function,
                   const PropertyMap &map = PropertyMap());
    
    LJPerturbation(const SireCAS::Expression &sigma_mapping_function,
                   const SireCAS::Expression &epsilon_mapping_function,
                   MapType maptype,
                   const PropertyMap &map = PropertyMap());
    
    LJPerturbation(const LJPerturbation &other);
    
    ~LJPerturbation();
    
    static const char* typeName();
    
    LJPerturbation& operator=(const LJPerturbation &other);
    
    bool operator==(const LJPerturbation &other) const;
    bool operator!=(const LJPerturbation &other) const;

    QString toString() const;

    SireMol::PerturbationPtr recreate(const SireCAS::Expression &mapping_function) const;
    SireMol::PerturbationPtr recreate(const SireCAS::Expression &mapping_function,
                                      const PropertyMap &map) const;

    SireMol::PerturbationPtr substitute(const SireCAS::Identities &identities) const;

    const SireCAS::Expression& mappingFunction() const;

    const SireCAS::Expression& rMinMappingFunction() const;
    const SireCAS::Expression& sigmaMappingFunction() const;
    const SireCAS::Expression& epsilonMappingFunction() const;
    
    const SireCAS::Expression& A_MappingFunction() const;
    const SireCAS::Expression& B_MappingFunction() const;
    
    bool mapSigmaEpsilon() const;
    bool mapRMinEpsilon() const;
    bool mapAB() const;

    QSet<QString> requiredProperties() const;
    
    bool wouldChange(const SireMol::Molecule &molecule, 
                     const SireCAS::Values &values) const;

protected:
    void perturbMolecule(SireMol::MolEditor &molecule, 
                         const SireCAS::Values &values) const;

private:
    /** Mapping function for sigma or A value */
    SireCAS::Expression sigma_mapfunc;
    
    /** How the LJ parameter is mapped */
    MapType maptype;
};

}

Q_DECLARE_METATYPE( SireMM::LJPerturbation )

SIRE_EXPOSE_CLASS( SireMM::LJPerturbation )

SIRE_END_HEADER

#endif
