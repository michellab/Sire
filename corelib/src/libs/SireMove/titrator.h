/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#ifndef SIREMOVE_TITRATOR_H
#define SIREMOVE_TITRATOR_H

#include "SireBase/property.h"
#include "SireBase/majorminorversion.h"

#include "SireMol/mgname.h"
#include "SireMol/mgnum.h"
#include "SireMol/molecule.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class Titrator;
}

QDataStream& operator<<(QDataStream&, const SireMove::Titrator&);
QDataStream& operator>>(QDataStream&, SireMove::Titrator&);

namespace SireSystem{ class System; }
namespace SireMol{ class MoleculeGroup; }

namespace SireMove
{

using SireSystem::System;
using SireMol::Molecule;
using SireMol::MoleculeGroup;
using SireBase::PropertyMap;

/** This property is used by the TitrationMove to maintain a list
    of which molecules can be titrated, and the list of titration
    states of each molecule.
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT Titrator
        : public SireBase::ConcreteProperty<Titrator,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const Titrator&);
friend QDataStream& ::operator>>(QDataStream&, Titrator&);

public:
    Titrator();
    Titrator(const MoleculeGroup &group);
    
    Titrator(const Titrator &other);
    
    ~Titrator();
    
    static const char* typeName();
    const char* what() const;
    
    Titrator& operator=(const Titrator &other);
    
    bool operator==(const Titrator &other) const;
    bool operator!=(const Titrator &other) const;

    void setMoleculeGroup(const MoleculeGroup &group);

    void setPositiveTemplate(const Molecule &positive_ion,
                             const QStringList &properties,
                             const PropertyMap &map = PropertyMap());
    
    void setNegativeTemplate(const Molecule &negative_ion,
                             const QStringList &properties,
                             const PropertyMap &map = PropertyMap());
    
    void setNeutralTemplate(const Molecule &neutral_mol,
                            const QStringList &properties,
                            const PropertyMap &map = PropertyMap());

    void setPositiveTemplate(const Molecule &positive_ion,
                             const PropertyMap &map = PropertyMap());
    void setNegativeTemplate(const Molecule &negative_ion,
                             const PropertyMap &map = PropertyMap());
    void setNeutralTemplate(const Molecule &neutral_mol,
                            const PropertyMap &map = PropertyMap());

    int nIons() const;
    int nNeutrals() const;

    int nPositiveIons() const;
    int nNegativeIons() const;

    int getIonIndex(int ion_index) const;
    int getPositiveIonIndex(int ion_index) const;
    int getNegativeIonIndex(int ion_index) const;
    int getNeutralIndex(int neutral_index) const;
    
    int getCharge(int i) const;
    void swapCharge(int i, int j);
    
    void randomiseCharge(int ncharges);
    void randomiseCharge(int npositive, int nnegative);
    
    double applyTo(System &system);

private:
    void initialiseFromGroup(const MoleculeGroup &group);
    void clearState();

    /** The name of the molecule group containing all of the 
        molecules that can swap */
    SireMol::MGName mgname;

    /** The number of the last molecule group to which this was applied */
    SireMol::MGNum mgnum;
    
    /** The version of the molecule group the last time this was applied */
    SireBase::Version mgversion;

    /** The current charge of each molecule in the group */
    QVector<qint32> chgs;
    
    /** The desired charge of each molecule in the group */
    QVector<qint32> desired_chgs;
    
    /** The template for the neutral molecule */
    Molecule neutral_template;
    
    /** The template for a negative ion */
    Molecule negative_template;
    
    /** The template for a positive ion */
    Molecule positive_template;
    
    /** The properties copied for a neutral molecule. If empty,
        everything is copied except for map["coordinates"] */
    QStringList neutral_properties;
    
    /** The properties copied for a negative ion. Again, if empty,
        everything is copied except for map["coordinates"] */
    QStringList negative_properties;
    
    /** The properties copied for a positive ion. Again, if empty,
        everything is copied except for map["coordinates"] */
    QStringList positive_properties;
    
    /** The property map used to extract properties from the neutral
        molecule */
    PropertyMap neutral_map;
    
    /** The property map used to extract properties from the negative ion */
    PropertyMap negative_map;
    
    /** The property map used to extract properties from the positive ion */
    PropertyMap positive_map;
    
    /** The property map used to put properties into the actual molecules */
    PropertyMap propmap;
    
    /** The charge of the positive ion template */
    qint32 pos_charge;
    
    /** The charge of the negative ion template */
    qint32 neg_charge;
};

} // end of namespace SireMove

Q_DECLARE_METATYPE( SireMove::Titrator )

SIRE_EXPOSE_CLASS( SireMove::Titrator )

SIRE_END_HEADER

#endif