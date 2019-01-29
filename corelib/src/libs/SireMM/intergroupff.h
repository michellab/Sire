/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIREMM_INTERGROUPFF_H
#define SIREMM_INTERGROUPFF_H

#include "cljgrid.h"
#include "cljfunction.h"
#include "cljgroup.h"
#include "multicljcomponent.h"

#include "SireBase/shareddatapointer.hpp"

#include "SireFF/g2ff.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class InterGroupFF;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::InterGroupFF&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::InterGroupFF&);

namespace SireMM
{

using SireBase::Property;
using SireBase::Properties;

namespace detail
{
    class InterGroupFFData;
}

/** This is a forcefield that calculates the intermolecular coulomb
    and Lennard Jones (LJ) energy between all molecules in group 0
    and all molecules in group 1.
 
    It also calculates the interactions between all molecules in group 0
    with any fixed atoms added to this forcefield
    
    @author Christopher Woods
*/
class SIREMM_EXPORT InterGroupFF : public SireBase::ConcreteProperty<InterGroupFF,SireFF::G2FF>
{

friend QDataStream& ::operator<<(QDataStream&, const InterGroupFF&);
friend QDataStream& ::operator>>(QDataStream&, InterGroupFF&);

public:
    InterGroupFF();
    InterGroupFF(const QString &name);
    
    InterGroupFF(const InterGroupFF &other);
    
    ~InterGroupFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    InterGroupFF& operator=(const InterGroupFF &other);
    
    bool operator==(const InterGroupFF &other) const;
    bool operator!=(const InterGroupFF &other) const;

    InterGroupFF* clone() const;

    const MultiCLJComponent& components() const;

    void setCLJFunction(const CLJFunction &cljfunc);
    const CLJFunction& cljFunction() const;

    void setCLJFunction(QString key, const CLJFunction &cljfunc);

    void removeCLJFunctionAt(QString key);
    void removeAllCLJFunctions();

    const CLJFunction& cljFunction(QString key) const;

    int nCLJFunctions() const;
    QStringList cljFunctionKeys() const;
    
    QHash<QString,CLJFunctionPtr> cljFunctions() const;

    void addFixedAtoms(const MoleculeView &molecule, const PropertyMap &map = PropertyMap());
    void addFixedAtoms(const Molecules &molecules, const PropertyMap &map = PropertyMap());
    void addFixedAtoms(const CLJAtoms &atoms);
    
    void setFixedAtoms(const MoleculeView &molecule, const PropertyMap &map = PropertyMap());
    void setFixedAtoms(const Molecules &molecules, const PropertyMap &map = PropertyMap());
    void setFixedAtoms(const CLJAtoms &atoms);

    void setFixedOnly(bool on);
    bool fixedOnly() const;

    void enableGrid();
    void disableGrid();
    void setUseGrid(bool on);
    bool usesGrid() const;
    
    void setGridBuffer(Length buffer);
    Length gridBuffer() const;
    
    void setGridSpacing(Length spacing);
    Length gridSpacing() const;

    GridInfo grid() const;

    void enableParallelCalculation();
    void disableParallelCalculation();
    void setUseParallelCalculation(bool on);
    bool usesParallelCalculation() const;

    void enableReproducibleCalculation();
    void disableReproducibleCalculation();
    void setUseReproducibleCalculation(bool on);
    bool usesReproducibleCalculation() const;

    bool setProperty(const QString &name, const Property &property);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;

    void mustNowRecalculateFromScratch();    

    void accept();
    bool needsAccepting() const;

private:
    void mustNowReallyRecalculateFromScratch();

    void recalculateEnergy();
    void rebuildProps();
    
    void regridAtoms();
    
    void _pvt_added(quint32 group_id,
                    const SireMol::PartialMolecule &mol,
                    const SireBase::PropertyMap &map);

    void _pvt_removed(quint32 group_id,
                      const SireMol::PartialMolecule &mol);

    void _pvt_changed(quint32 group_id, const SireMol::Molecule &molecule, bool auto_update);
    void _pvt_changed(quint32 group_id,
                      const QList<SireMol::Molecule> &molecules, bool auto_update);
    
    void _pvt_removedAll(quint32 group_id);
        
    bool _pvt_wouldChangeProperties(quint32 group_id, SireMol::MolNum molnum,
                                    const SireBase::PropertyMap &map) const;

    void _pvt_updateName();

    /** The CLJGroups containing all of the extracted molecules */
    CLJGroup cljgroup[2];

    /** Implicitly shared pointer to the (mostly) const data for this forcefield */
    SireBase::SharedDataPointer<detail::InterGroupFFData> d;
    
    /** Whether or not we need to 'accept' this move */
    bool needs_accepting;
};

}

Q_DECLARE_METATYPE( SireMM::InterGroupFF )

SIRE_EXPOSE_CLASS( SireMM::InterGroupFF )

SIRE_END_HEADER

#endif
