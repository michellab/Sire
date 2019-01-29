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

#ifndef SIREIO_PERTURBATIONSLIBRARY_H
#define SIREIO_PERTURBATIONSLIBRARY_H

#include "iobase.h"

#include "SireCAS/expression.h"

#include "SireMM/ljparameter.h"
#include "SireMol/molecule.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"
#include "SireMol/improperid.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class PerturbationsLibrary;
class PerturbationsTemplate;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::PerturbationsLibrary&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::PerturbationsLibrary&);

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::PerturbationsTemplate&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::PerturbationsTemplate&);

namespace SireMM
{
  class LJParameter;
}

namespace SireIO
{  

  using SireCAS::Expression;
  using SireMol::BondID;
  using SireMol::AngleID;
  using SireMol::DihedralID;
  using SireMol::ImproperID;
  using SireMol::Molecule;

  using SireMM::LJParameter;
  using SireUnits::Dimension::Charge;

/** Internal class used to store the data describing a single perturbations template
    
    @author Julien Michel
*/
class SIREIO_EXPORT PerturbationsTemplate
{

friend QDataStream& ::operator<<(QDataStream&, const PerturbationsTemplate&);
friend QDataStream& ::operator>>(QDataStream&, PerturbationsTemplate&);

public:
    PerturbationsTemplate();
    PerturbationsTemplate(const QString &name);
    
    PerturbationsTemplate(const PerturbationsTemplate &other);
    
    ~PerturbationsTemplate();

    static const char* typeName();

    PerturbationsTemplate& operator=(const PerturbationsTemplate &other);

    bool operator==(const PerturbationsTemplate &other) const;
    bool operator!=(const PerturbationsTemplate &other) const;

    const QString getName();
    
    void setInitCharge(const QString &atomname, const Charge &atomcharge);
    void setFinalCharge(const QString &atomname, const Charge &atomcharge);

    void setInitLJ(const QString &atomname, const LJParameter &atomlj);
    void setFinalLJ(const QString &atomname, const LJParameter &atomlj);

    void setInitType(const QString &atomname, const QString &atype);
    void setFinalType(const QString &atomname, const QString &atype);

    void setInitBondK(const BondID &bond, const double &k);
    void setInitBondR(const BondID &bond, const double &r);
    void setFinalBondK(const BondID &bond, const double &k);
    void setFinalBondR(const BondID &bond, const double &r);

    void setInitAngleK(const AngleID &angle, const double &k);
    void setInitAngleT(const AngleID &angle, const double &r);
    void setFinalAngleK(const AngleID &angle, const double &k);
    void setFinalAngleT(const AngleID &angle, const double &r);

    //void setInitDihPotential(const DihedralID &dihedral, const Expression &pot);
    //void setFinalDihPotential(const DihedralID &dihedral, const Expression &pot);

    void setInitDihParams(const DihedralID &dihedral, const QList<double> &params);
    void setFinalDihParams(const DihedralID &dihedral, const QList<double> &params);

    //void setInitDihedralK0(const DihedralID &dihedral, const double &k0);
    //void setInitDihedralN(const DihedralID &dihedral, const int &n);
    //void setInitDihedralPhase(const DihedralID &dihedral, const double &phase);
    //void setFinalDihedralK0(const DihedralID &dihedral, const double &k0);
    //void setFinalDihedralN(const DihedralID &dihedral, const int &n);
    //void setFinalDihedralPhase(const DihedralID &dihedral, const double &phase);

    //void setInitImpPotential(const ImproperID &improper, const Expression &pot);
    //void setFinalImpPotential(const ImproperID &improper, const Expression &pot);

    void setInitImpParams(const ImproperID &improper, const QList<double> &params);
    void setFinalImpParams(const ImproperID &improper, const QList<double> &params);

    //void setInitImproperK0(const ImproperID &improper, const double &k0);
    //void setInitImproperN(const ImproperID &improper, const int &n);
    //void setInitImproperPhase(const ImproperID &improper, const double &phase);
    //void setFinalImproperK0(const ImproperID &improper, const double &k0);
    //void setFinalImproperN(const ImproperID &improper, const int &n);
    //void setFinalImproperPhase(const ImproperID &improper, const double &phase);

    Charge getInitCharge(const QString &atomname) const;
    Charge getFinalCharge(const QString &atomname) const;
    LJParameter getInitLJ(const QString &atomname) const;
    LJParameter getFinalLJ(const QString &atomname) const;
    QString getInitType(const QString &atomname) const;
    QString getFinalType(const QString &atomname) const;

    QList<BondID> getBonds() const;
    double getInitBondK(const BondID &bond) const;
    double getInitBondR(const BondID &bond) const;
    double getFinalBondK(const BondID &bond) const;
    double getFinalBondR(const BondID &bond) const;

    QList<AngleID> getAngles() const;
    double getInitAngleK(const AngleID &angle) const;
    double getInitAngleT(const AngleID &angle) const;
    double getFinalAngleK(const AngleID &angle) const;
    double getFinalAngleT(const AngleID &angle) const;

    QList<DihedralID> getDihedrals() const;
    //Expression getInitDihPotential(const DihedralID &dihedral) const;
    //Expression getFinalDihPotential(const DihedralID &dihedral) const;

    QList<double> getInitDihParams(const DihedralID &dihedral) const;
    QList<double> getFinalDihParams(const DihedralID &dihedral) const;
    
    //double getInitDihedralK0(const DihedralID &dihedral) const;
    //int getInitDihedralN(const DihedralID &dihedral) const;
    //double getInitDihedralPhase(const DihedralID &dihedral) const;
    //double getFinalDihedralK0(const DihedralID &dihedral) const;
    //int getFinalDihedralN(const DihedralID &dihedral) const;
    //double getFinalDihedralPhase(const DihedralID &dihedral) const;

    QList<ImproperID> getImpropers() const;
    //Expression getInitImpPotential(const ImproperID &improper) const;
    //Expression getFinalImpPotential(const ImproperID &improper) const;

    QList<double> getInitImpParams(const ImproperID &improper) const;
    QList<double> getFinalImpParams(const ImproperID &improper) const;
    
    //double getInitImproperK0(const ImproperID &improper) const;
    //int getInitImproperN(const ImproperID &improper) const;
    //double getInitImproperPhase(const ImproperID &improper) const;
    //double getFinalImproperK0(const ImproperID &improper) const;
    //int getFinalImproperN(const ImproperID &improper) const;
    //double getFinalImproperPhase(const ImproperID &improper) const;

private:
    QString name;
    // The atom charges 
    QHash<QString,Charge> initcharges;
    QHash<QString,Charge> finalcharges;
    // The atom LJs
    QHash<QString,LJParameter> initLJs;
    QHash<QString,LJParameter> finalLJs;
    // The atom types
    QHash<QString,QString> initatypes;
    QHash<QString,QString> finalatypes;
    // The bond parameters
    QHash<BondID,double> initbondsk;
    QHash<BondID,double> initbondsr;    
    QHash<BondID,double> finalbondsk;
    QHash<BondID,double> finalbondsr;
    // The angle parameters
    QHash<AngleID,double> initanglesk;
    QHash<AngleID,double> initanglest;    
    QHash<AngleID,double> finalanglesk;
    QHash<AngleID,double> finalanglest;
    // The dihedral parameters
    //QHash<DihedralID,Expression> initdihpotential;
    //QHash<DihedralID,Expression> finaldihpotential;

    QHash<DihedralID,QList<double> > initdihparams;
    QHash<DihedralID,QList<double> > finaldihparams;
    //QHash<DihedralID,double> initdihedralsk0;
    //QHash<DihedralID,double> initdihedralsn;    
    //QHash<DihedralID,double> initdihedralsphase;
    //QHash<DihedralID,double> finaldihedralsk0;
    //QHash<DihedralID,double> finaldihedralsn;    
    //QHash<DihedralID,double> finaldihedralsphase;
    // The improper parameters
    //QHash<ImproperID,Expression> initimppotential;
    //QHash<ImproperID,Expression> finalimppotential;

    QHash<ImproperID,QList<double> > initimpparams;
    QHash<ImproperID,QList<double> > finalimpparams;
    //QHash<ImproperID,double> initimpropersk0;
    //QHash<ImproperID,double> initimpropersn;    
    //QHash<ImproperID,double> initimpropersphase;
    //QHash<ImproperID,double> finalimpropersk0;
    //QHash<ImproperID,double> finalimpropersn;    
    //QHash<ImproperID,double> finalimpropersphase;

};

/** This class is used to read templates describing how a 
    molecule can be perturbed

    @author Julien Michel
*/

class SIREIO_EXPORT PerturbationsLibrary 
        : public SireBase::ConcreteProperty<PerturbationsLibrary,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const SireIO::PerturbationsLibrary&);
friend QDataStream& ::operator>>(QDataStream&, SireIO::PerturbationsLibrary&);

public:
    PerturbationsLibrary();
    PerturbationsLibrary(const QString &file);
    
    PerturbationsLibrary(const PerturbationsLibrary &other);
    
    ~PerturbationsLibrary();
    
    static const char* typeName();
    
    PerturbationsLibrary& operator=(const PerturbationsLibrary &other);

    bool operator==(const PerturbationsLibrary &other) const;
    bool operator!=(const PerturbationsLibrary &other) const;

    void loadTemplates(const QString &file);

    PerturbationsLibrary& operator+=(const PerturbationsLibrary &other);

    PerturbationsLibrary operator+(const PerturbationsLibrary &other) const;
    
    void add(const PerturbationsLibrary &other);
    
    const PerturbationsTemplate& getTemplate(const QString &key);
    
    void setTemplate(const QString &key, const PerturbationsTemplate &tmplate);
    
    Molecule applyTemplate(const Molecule &molecule) const;

private:
    /** The perturbations templates, indexed by molecule name*/
    QHash<QString,PerturbationsTemplate> templates; 
};

} // end of namespace SireIO

Q_DECLARE_METATYPE( SireIO::PerturbationsLibrary )
Q_DECLARE_METATYPE( SireIO::PerturbationsTemplate )

SIRE_EXPOSE_CLASS( SireIO::PerturbationsLibrary )
SIRE_EXPOSE_CLASS( SireIO::PerturbationsTemplate )

SIRE_END_HEADER

#endif
