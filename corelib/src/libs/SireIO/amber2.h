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
class AmberParm7;
class Amber2;
}

QDataStream& operator<<(QDataStream&, const SireIO::Amber2&);
QDataStream& operator>>(QDataStream&, SireIO::Amber2&);

QDataStream& operator<<(QDataStream&, const SireIO::AmberParm7&);
QDataStream& operator>>(QDataStream&, SireIO::AmberParm7&);

namespace SireSystem
{
class System;
}

namespace SireIO
{

/** This class represents an Amber7-format parameter file, currently
    supporting top files produced from Amber7 until Amber16
    
    @author Christopher Woods
*/
class SIREIO_EXPORT AmberParm7
{

friend QDataStream& ::operator<<(QDataStream&, const AmberParm7&);
friend QDataStream& ::operator>>(QDataStream&, AmberParm7&);

public:
    enum FLAG_TYPE { UNKNOWN = 0,
                     INTEGER = 1,
                     FLOAT = 2,
                     STRING = 3 };

    AmberParm7();
    AmberParm7(const QString &filename);
    AmberParm7(const SireSystem::System &system,
               const SireBase::PropertyMap &map = SireBase::PropertyMap());
    
    AmberParm7(const AmberParm7 &other);
    
    ~AmberParm7();
    
    AmberParm7& operator=(const AmberParm7 &other);
    
    bool operator==(const AmberParm7 &other) const;
    bool operator!=(const AmberParm7 &other) const;
  
    static const char* typeName();
  
    const char* what() const;
    
    static AmberParm7 read(const QString &filename);
    static AmberParm7 write(const SireSystem::System &system,
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
class SIREIO_EXPORT Amber2
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
  
    void set14Factors(double coul_14, double lj_14);
    
    double coulomb14Factor() const;
    double lj14Factor() const;
  
    SireSystem::System readRst7Parm7(const QString &rstfile,
                                     const QString &prmfile,
                                     const SireMol::CuttingFunction &cutting_function
                                             = SireMol::ResidueCutting(),
                                     const SireBase::PropertyMap &map
                                             = SireBase::PropertyMap()) const;
  
    SireSystem::System readRst7Parm7(const QString &rstfile,
                                     const QString &prmfile,
                                     const SireBase::PropertyMap &map,
                                     const SireMol::CuttingFunction &cutting_function
                                             = SireMol::ResidueCutting()) const;

    void writeRst7Parm7(const SireSystem::System &system,
                        const QString &rstfile, const QString &prmfile,
                        const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

private:
    double coul_14scl;
    double lj_14scl;
};

}

Q_DECLARE_METATYPE( SireIO::Amber2 )
Q_DECLARE_METATYPE( SireIO::AmberParm7 )

SIRE_EXPOSE_CLASS( SireIO::Amber2 )
SIRE_EXPOSE_CLASS( SireIO::AmberParm7 )

SIRE_END_HEADER

#endif

