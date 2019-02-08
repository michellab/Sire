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

#ifndef SIREMM_INTRAGROUPFF_H
#define SIREMM_INTRAGROUPFF_H

#include "cljgroup.h"
#include "cljfunction.h"
#include "multicljcomponent.h"

#include "SireBase/chunkedhash.hpp"
#include "SireBase/shareddatapointer.hpp"

#include "SireFF/g2ff.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class IntraGroupFF;
namespace detail{ class IntraGroupFFMolData; }
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::IntraGroupFF&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::IntraGroupFF&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::detail::IntraGroupFFMolData&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::detail::IntraGroupFFMolData&);

namespace SireMM
{

namespace detail
{
class IntraGroupFFData;
}

/** This forcefield is used to calculate the intramolecular 
    coulomb and LJ energy of the contained molecules. Note
    that this is the coulomb and LJ energy of the non-bonded
    atoms *only*, i.e. it does not contain the scaled
    1-4 coulomb and LJ energies. These should be calculated
    separately, e.g. via additional terms added to InternalFF
    
    @author Christopher Woods
*/
class SIREMM_EXPORT IntraGroupFF : public SireBase::ConcreteProperty<IntraGroupFF,SireFF::G2FF>
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const IntraGroupFF&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, IntraGroupFF&);

public:
    IntraGroupFF();
    IntraGroupFF(const QString &name);
    
    IntraGroupFF(const IntraGroupFF &other);
    
    ~IntraGroupFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    IntraGroupFF& operator=(const IntraGroupFF &other);
    
    bool operator==(const IntraGroupFF &other) const;
    bool operator!=(const IntraGroupFF &other) const;

    IntraGroupFF* clone() const;

    const MultiCLJComponent& components() const;

    void setCLJFunction(const CLJIntraFunction &cljfunc);
    const CLJIntraFunction& cljFunction() const;

    void setCLJFunction(QString key, const CLJIntraFunction &cljfunc);

    void removeCLJFunctionAt(QString key);
    void removeAllCLJFunctions();

    const CLJIntraFunction& cljFunction(QString key) const;

    int nCLJFunctions() const;
    QStringList cljFunctionKeys() const;
    
    QHash<QString,CLJFunctionPtr> cljFunctions() const;

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
    
    MultiCLJEnergy calcEnergy(detail::IntraGroupFFMolData &moldata) const;
    
    void _pvt_added(quint32 group_idx,
                    const SireMol::PartialMolecule &mol,
                    const SireBase::PropertyMap &map);

    void _pvt_removed(quint32 group_idx,
                      const SireMol::PartialMolecule &mol);

    void _pvt_changed(quint32 group_idx,
                      const SireMol::Molecule &molecule, bool auto_update);
    void _pvt_changed(quint32 group_idx,
                      const QList<SireMol::Molecule> &molecules, bool auto_update);
    
    void _pvt_removedAll(quint32 group_idx);
        
    bool _pvt_wouldChangeProperties(quint32 group_idx,
                                    SireMol::MolNum molnum,
                                    const SireBase::PropertyMap &map) const;

    void _pvt_updateName();

    typedef SireBase::ChunkedHash< MolNum,
                        SireBase::SharedDataPointer<detail::IntraGroupFFMolData> > MolData;

    /** The CLJGroups for each of the molecules added to this forcefield */
    MolData moldata;

    /** Implicitly shared pointer to the (mostly) const data for this forcefield */
    SireBase::SharedDataPointer<detail::IntraGroupFFData> d;

    /** Whether or not we need to 'accept' this move */
    bool needs_accepting;
};

}

Q_DECLARE_METATYPE( SireMM::IntraGroupFF )

SIRE_EXPOSE_CLASS( SireMM::IntraGroupFF )

SIRE_END_HEADER

#endif

