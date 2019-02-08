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

#ifndef SIREMM_INTRAFF_H
#define SIREMM_INTRAFF_H

#include "cljgroup.h"
#include "cljfunction.h"
#include "multicljcomponent.h"

#include "SireBase/shareddatapointer.hpp"
#include "SireBase/chunkedhash.hpp"

#include "SireFF/g1ff.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class IntraFF;
namespace detail{ class IntraFFMolData; }
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::IntraFF&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::IntraFF&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::detail::IntraFFMolData&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::detail::IntraFFMolData&);

namespace SireMM
{

namespace detail
{
class IntraFFData;
}

/** This forcefield is used to calculate the intramolecular 
    coulomb and LJ energy of the contained molecules. Note
    that this is the coulomb and LJ energy of the non-bonded
    atoms *only*, i.e. it does not contain the scaled
    1-4 coulomb and LJ energies. These should be calculated
    separately, e.g. via additional terms added to InternalFF
    
    @author Christopher Woods
*/
class SIREMM_EXPORT IntraFF : public SireBase::ConcreteProperty<IntraFF,SireFF::G1FF>
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const IntraFF&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, IntraFF&);

public:
    IntraFF();
    IntraFF(const QString &name);
    
    IntraFF(const IntraFF &other);
    
    ~IntraFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    IntraFF& operator=(const IntraFF &other);
    
    bool operator==(const IntraFF &other) const;
    bool operator!=(const IntraFF &other) const;

    IntraFF* clone() const;

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
    
    MultiCLJEnergy calcEnergy(detail::IntraFFMolData &moldata) const;
    
    void _pvt_added(const SireMol::PartialMolecule &mol,
                    const SireBase::PropertyMap &map);

    void _pvt_removed(const SireMol::PartialMolecule &mol);

    void _pvt_changed(const SireMol::Molecule &molecule, bool auto_update);
    void _pvt_changed(const QList<SireMol::Molecule> &molecules, bool auto_update);
    
    void _pvt_removedAll();
        
    bool _pvt_wouldChangeProperties(SireMol::MolNum molnum,
                                    const SireBase::PropertyMap &map) const;

    void _pvt_updateName();

    typedef SireBase::ChunkedHash< MolNum,
                                   SireBase::SharedDataPointer<detail::IntraFFMolData> > MolData;

    /** The CLJGroups for each of the molecules added to this forcefield */
    MolData moldata;

    /** Implicitly shared pointer to the (mostly) const data for this forcefield */
    SireBase::SharedDataPointer<detail::IntraFFData> d;

    /** Whether or not we need to 'accept' this move */
    bool needs_accepting;
};

}

Q_DECLARE_METATYPE( SireMM::IntraFF )

SIRE_EXPOSE_CLASS( SireMM::IntraFF )

SIRE_END_HEADER

#endif
