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

#ifndef SIREMM_PDBPARAMS_H
#define SIREMM_PDBPARAMS_H

#include "SireBase/propertymap.h"

#include "SireMol/molecule.h"
#include "SireMol/molviewproperty.h"

namespace SireMM
{
class PDBParams;
}

QDataStream& operator<<(QDataStream &ds, const SireMM::PDBParams &params);
QDataStream& operator>>(QDataStream &ds, SireMM::PDBParams &params);

namespace SireMM
{

/** This class holds all of the extra parameter data that has been
    read into the molecule from a PDB file
    
    @author Lester Hedges
*/
class SIREMM_EXPORT PDBParams
    : public SireBase::ConcreteProperty<PDBParams,SireMol::MoleculeProperty>
{

friend QDataStream& ::operator<<(QDataStream&, const PDBParams&);
friend QDataStream& ::operator>>(QDataStream&, PDBParams&);

public:
    PDBParams();
    PDBParams(const SireMol::MoleculeView &molecule,
              const SireBase::PropertyMap &map = SireBase::PropertyMap());
                
    PDBParams(const SireMol::MoleculeInfo &molinfo);
    PDBParams(const SireMol::MoleculeInfoData &molinfo);
    
    PDBParams(const PDBParams &other);

    ~PDBParams();

    static const char* typeName();

    PDBParams& operator=(const PDBParams &other);
    
    bool operator==(const PDBParams &other) const;
    bool operator!=(const PDBParams &other) const;

    SireMol::MoleculeInfo info() const;

    bool isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const;

    QString toString() const;
    
    void setPropertyMap(const SireBase::PropertyMap &map);
    const SireBase::PropertyMap& propertyMap() const;
  
    void updateFrom(const SireMol::MoleculeView &molview);

    QStringList validate() const;

private:
    void _pvt_createFrom(const SireMol::MoleculeData &moldata);
    void _pvt_updateFrom(const SireMol::MoleculeData &moldata);

    /** The molecule that this flexibility operates on */
    SireMol::MoleculeInfo molinfo;

    /** The property map used (if any) to identify the properties that
        hold the amber parameters */
    SireBase::PropertyMap propmap;
};

}

Q_DECLARE_METATYPE( SireMM::PDBParams )

SIRE_EXPOSE_CLASS( SireMM::PDBParams )

#endif
