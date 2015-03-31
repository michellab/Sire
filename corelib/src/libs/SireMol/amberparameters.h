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

#ifndef SIREMOL_AMBERPARAMETERS_H
#define SIREMOL_AMBERPARAMETERS_H

#include "SireBase/propertymap.h"
#include "SireBase/shareddatapointer.hpp"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"
#include "SireMol/atomidx.h"
#include "SireMol/mover.hpp"
#include "SireMol/molviewproperty.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class AmberParameters;
}

QDataStream& operator<<(QDataStream&, const SireMol::AmberParameters&);
QDataStream& operator>>(QDataStream&, SireMol::AmberParameters&);

namespace SireMol
{
 class BondID;
 class AngleID;
 class DihedralID;
 class Molecule;
 class MoleculeData;
 class MoleculeInfoData;
}

namespace SireMol
{
using SireBase::PropertyMap;

using SireMol::BondID;
using SireMol::AngleID;
using SireMol::DihedralID;
using SireMol::Molecule;
using SireMol::MoleculeData;
using SireMol::MoleculeInfoData;

/** This class stores AMBER bonded force field parameters for a collection of bonds, angles and dihedrals


    @author Julien Michel
 */
class SIREMOL_EXPORT AmberParameters
    : public SireBase::ConcreteProperty<AmberParameters,SireMol::MoleculeProperty>
{

friend QDataStream& ::operator<<(QDataStream&, const SireMol::AmberParameters&);
friend QDataStream& ::operator>>(QDataStream&, SireMol::AmberParameters&);

 public:
    AmberParameters();
    AmberParameters(const MoleculeData &molecule);
    AmberParameters(const AmberParameters &other);

    ~AmberParameters();

    static const char* typeName();

    AmberParameters& operator=(const AmberParameters &other);
    
    bool operator==(const AmberParameters &other) const;
    bool operator!=(const AmberParameters &other) const;

    const SireMol::MoleculeInfoData& info() const;

    bool isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const;

    void add(const BondID &bond, const double &k, const double &ro);
    void remove(const BondID &bond);
    QList<double> getParams(const BondID &bond);
    QList<BondID> getAllBonds();

    void add(const AngleID &angle, const double &k, const double &theta0);
    void remove(const AngleID &angle);
    QList<double> getParams(const AngleID &angle);
    QList<AngleID> getAllAngles();

    void add(const DihedralID &dihedral, const double &v, const double &periodicity, const double &phase);
    void remove(const DihedralID &dihedral);
    QList<double> getParams(const DihedralID &dihedral);
    QList<DihedralID> getAllDihedrals();

    void add(const ImproperID &improper, const double &v, const double &periodicity, const double &phase);
    void remove(const ImproperID &improper);
    QList<double> getParams(const ImproperID &improper);
    QList<ImproperID> getAllImpropers();

    void add14Pair(const BondID &pair, const double &cscl, const double &ljscl);
    void remove14Pair(const BondID &pair);
    QList<double> get14PairParams(const BondID &pair);
    QList<BondID> getAll14Pairs();
  
 private:
    /** The molecule that this flexibility operates on */
    SireBase::SharedDataPointer<SireMol::MoleculeInfoData> molinfo;

    /**A Hash of force constants and equilibrium bond lengths for bonds **/
    QHash< BondID, QList<double> > bonds;

    /**A Hash of force constants and equilibrium bond angles for angles **/
    QHash< AngleID, QList<double> > angles;

    /**A Hash of torsional parameters for dihedrals **/
    QHash< DihedralID, QList<double> > dihedrals;

    /**A Hash of torsional parameters for impropers **/
    QHash< ImproperID, QList<double> > impropers;    

    /**A Hash of coulombic and lennard jones scale factors for 1,4 pairs**/
    QHash< BondID, QList<double> > nb14pairs;
};

}

Q_DECLARE_METATYPE( SireMol::AmberParameters )

SIRE_EXPOSE_CLASS( SireMol::AmberParameters )

SIRE_END_HEADER

#endif
