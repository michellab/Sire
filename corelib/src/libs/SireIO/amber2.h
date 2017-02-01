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

#ifndef SIREIO_AMBER2_H
#define SIREIO_AMBER2_H

#include "SireBase/propertymap.h"
#include "SireBase/shareddatapointer.hpp"

#include "SireMol/residuecutting.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class MoleculeParser;
class AmberParm;
class AmberRst;
class Amber2;
}

QDataStream& operator<<(QDataStream&, const SireIO::MoleculeParser&);
QDataStream& operator>>(QDataStream&, SireIO::MoleculeParser&);

QDataStream& operator<<(QDataStream&, const SireIO::Amber2&);
QDataStream& operator>>(QDataStream&, SireIO::Amber2&);

QDataStream& operator<<(QDataStream&, const SireIO::AmberParm&);
QDataStream& operator>>(QDataStream&, SireIO::AmberParm&);

QDataStream& operator<<(QDataStream&, const SireIO::AmberRst&);
QDataStream& operator>>(QDataStream&, SireIO::AmberRst&);

namespace SireSystem
{
class System;
}

namespace SireIO
{

/** The base class of all molecule parsers */
class SIREIO_EXPORT MoleculeParser : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeParser&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeParser&);

public:
    MoleculeParser();
    MoleculeParser(const MoleculeParser &other);
    ~MoleculeParser();
    
    static const char* typeName();
    
protected:
    MoleculeParser& operator=(const MoleculeParser &other);
    
    bool operator==(const MoleculeParser &other) const;
    bool operator!=(const MoleculeParser &other) const;
};

/** This class represents an Amber-format restart/coordinate file (ascii),
    currently supporting these files from Amber7 to Amber16.
    
    The format of this file is described here;
    
    http://ambermd.org/formats.html
    
    (specifically the "AMBER coordinate/restart" file specification
 
    @author Christopher Woods
*/
class SIREIO_EXPORT AmberRst : public SireBase::ConcreteProperty<AmberRst,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const AmberRst&);
friend QDataStream& ::operator>>(QDataStream&, AmberRst&);

public:
    AmberRst();
    AmberRst(const AmberRst &other);
    
    ~AmberRst();
    
    AmberRst& operator=(const AmberRst &other);
    
    bool operator==(const AmberRst &other) const;
    bool operator!=(const AmberRst &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
};

/** This class represents an Amber-format parameter file, currently
    supporting top files produced from Amber7 until Amber16
    
    The format of this file is described here;
    
    http://ambermd.org/formats.html
    
    (specifically the "PARM" parameter/topology file specification)
 
    @author Christopher Woods
*/
class SIREIO_EXPORT AmberParm : public SireBase::ConcreteProperty<AmberParm,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const AmberParm&);
friend QDataStream& ::operator>>(QDataStream&, AmberParm&);

public:
    enum FLAG_TYPE { UNKNOWN = 0,
                     INTEGER = 1,
                     FLOAT = 2,
                     STRING = 3 };

    AmberParm();
    AmberParm(const QString &filename);
    AmberParm(const SireSystem::System &system,
               const SireBase::PropertyMap &map = SireBase::PropertyMap());
    
    AmberParm(const AmberParm &other);
    
    ~AmberParm();
    
    AmberParm& operator=(const AmberParm &other);
    
    bool operator==(const AmberParm &other) const;
    bool operator!=(const AmberParm &other) const;
  
    static const char* typeName();
  
    const char* what() const;
    
    QString toString() const;
    
    static AmberParm read(const QString &filename);
    static AmberParm write(const SireSystem::System &system,
                            const SireBase::PropertyMap &map = SireBase::PropertyMap());
    
    SireSystem::System toSystem() const;

    SireMol::Molecule getMolecule(int i) const;
    
    QStringList lines() const;
    QStringList lines(const QString &flag) const;
    
    QStringList flags() const;

    FLAG_TYPE flagType(const QString &flag) const;
    
    QList<qint64> intData(const QString &flag) const;
    QList<double> floatData(const QString &flag) const;
    QStringList stringData(const QString &flag) const;

    QString title() const;

    int nAtoms() const;
    
    int nTypes() const;
    
    int nBonds() const;
    int nBondsWithHydrogen() const;
    int nBondsNoHydrogen() const;
    
    int nAngles() const;
    int nAnglesWithHydrogen() const;
    int nAnglesNoHydrogen() const;
    
    int nDihedrals() const;
    int nDihedralsWithHydrogen() const;
    int nDihedralsNoHydrogen() const;
    
    int nExcluded() const;
    int nResidues() const;
    
    int nMolecules() const;
    
    int nAtoms(int molidx) const;
    
    void assertSane() const;

private:
    void rebuildAfterReload();

    SireMol::Molecule getMolecule(int start_idx, int natoms) const;

    QList< QPair<int,int> > moleculeIndicies() const;

    /** Function to process all flags */
    void processAllFlags();

    /** All of the lines contained in the parm file */
    QStringList lnes;
    
    /** A map showing the line number of all flags. This holds
        the start index and number of lines for each flag */
    QHash< QString,QPair<qint64,qint64> > flag_to_line;
    
    /** The raw int data for the integer flags */
    QHash< QString, QList<qint64> > int_data;
    
    /** The raw float data for the float flags */
    QHash< QString, QList<double> > float_data;
    
    /** The raw string data for the string flags */
    QHash< QString, QStringList > string_data;
    
    /** A copy of the POINTER data to prevent over-lookup */
    QList<qint64> pointers;
};

/** This class is used to read and write AMBER molecule files.
    The class will aim to support the range of file formats
    used by Amber. This class is based on the original Amber
    class written by Julien Michel
    
    @author Christopher Woods
*/
class SIREIO_EXPORT Amber2 : public SireBase::ConcreteProperty<Amber2,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const SireIO::Amber2&);
friend QDataStream& ::operator>>(QDataStream&, SireIO::Amber2&);
  
public:
    Amber2();
    Amber2(const Amber2 &other);
    ~Amber2();
  
    Amber2& operator=(const Amber2 &other);
    
    bool operator==(const Amber2 &other) const;
    bool operator!=(const Amber2 &other) const;
  
    static const char* typeName();
  
    const char* what() const;
    
    QString toString() const;
    
    void set14Factors(double coul_14, double lj_14);
    
    double coulomb14Factor() const;
    double lj14Factor() const;
  
    SireSystem::System readRstParm(const QString &rstfile,
                                   const QString &prmfile,
                                   const SireMol::CuttingFunction &cutting_function
                                            = SireMol::ResidueCutting(),
                                   const SireBase::PropertyMap &map
                                            = SireBase::PropertyMap()) const;
  
    SireSystem::System readRstParm(const QString &rstfile,
                                   const QString &prmfile,
                                   const SireBase::PropertyMap &map,
                                   const SireMol::CuttingFunction &cutting_function
                                           = SireMol::ResidueCutting()) const;

    void writeRstParm(const SireSystem::System &system,
                      const QString &rstfile, const QString &prmfile,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

private:
    double coul_14scl;
    double lj_14scl;
};

}

Q_DECLARE_METATYPE( SireIO::Amber2 )
Q_DECLARE_METATYPE( SireIO::AmberParm )
Q_DECLARE_METATYPE( SireIO::AmberRst )

SIRE_EXPOSE_CLASS( SireIO::MoleculeParser )
SIRE_EXPOSE_CLASS( SireIO::Amber2 )
SIRE_EXPOSE_CLASS( SireIO::AmberParm )
SIRE_EXPOSE_CLASS( SireIO::AmberRst )

SIRE_END_HEADER

#endif

