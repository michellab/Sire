/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIREMM_CLJNBPAIRS_H
#define SIREMM_CLJNBPAIRS_H

#include "atompairs.hpp"

#include "SireMol/connectivity.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class CoulombScaleFactor;
class LJScaleFactor;
class CLJScaleFactor;

class CoulombNBPairs;
class LJNBPairs;
class CLJNBPairs;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CoulombNBPairs&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CoulombNBPairs&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::LJNBPairs&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::LJNBPairs&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJNBPairs&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJNBPairs&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CoulombScaleFactor&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CoulombScaleFactor&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::LJScaleFactor&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::LJScaleFactor&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJScaleFactor&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJScaleFactor&);

namespace SireMM
{

/** This is the interatomic scale factor for the coulomb
    parameters for the intramolecular energy */
class SIREMM_EXPORT CoulombScaleFactor
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const CoulombScaleFactor&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, CoulombScaleFactor&);

public:
    CoulombScaleFactor(double scl=0);
    
    CoulombScaleFactor(const CoulombScaleFactor &other);
    
    ~CoulombScaleFactor();

    static const char* typeName();

    const char* what() const
    {
        return CoulombScaleFactor::typeName();
    }

    CoulombScaleFactor& operator=(const CoulombScaleFactor &other);

    bool operator==(const CoulombScaleFactor &other) const;
    bool operator!=(const CoulombScaleFactor &other) const;

    double coulomb() const;

private:
    /** The coulomb scale factor */
    double cscl;
};

/** This is the interatomic scale factor for the LJ
    parameters for the intramolecular energy */
class SIREMM_EXPORT LJScaleFactor
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const LJScaleFactor&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, LJScaleFactor&);

public:
    LJScaleFactor(double scl=0);
    
    LJScaleFactor(const LJScaleFactor &other);
    
    ~LJScaleFactor();

    static const char* typeName();

    const char* what() const
    {
        return LJScaleFactor::typeName();
    }

    LJScaleFactor& operator=(const LJScaleFactor &other);

    bool operator==(const LJScaleFactor &other) const;
    bool operator!=(const LJScaleFactor &other) const;

    double lj() const;

private:
    /** The LJ scale factor */
    double ljscl;
};

/** This is the interatomic scale factor for the coulomb and
    LJ parameters for the intramolecular energy. */
class SIREMM_EXPORT CLJScaleFactor : public CoulombScaleFactor,
                                     public LJScaleFactor
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const CLJScaleFactor&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, CLJScaleFactor&);

public:
    CLJScaleFactor(double scl=0);

    CLJScaleFactor(double scale_coul, double scale_lj);
    
    CLJScaleFactor(const CLJScaleFactor &other);
    
    ~CLJScaleFactor();

    static const char* typeName();

    const char* what() const
    {
        return CLJScaleFactor::typeName();
    }

    QString toString() const;

    CLJScaleFactor& operator=(const CLJScaleFactor &other);

    bool operator==(const CLJScaleFactor &other) const;
    bool operator!=(const CLJScaleFactor &other) const;
};

/** This class holds all of the non-bonded scale factors that are used
    to scale the intramolecular atom-atom coulomb
    interactions between pairs of atoms, e.g. for most MM forcefields,
    the scale factors for 1-1, 1-2 and 1-3 pairs are zero, the
    1-4 pairs are scaled by a coulomb factor (e.g. 0.5 for OPLS)
    and the 1-5 and above pairs are not scaled (i.e. the factors equal 1)

    @author Christopher Woods
*/
class SIREMM_EXPORT CoulombNBPairs 
        : public SireBase::ConcreteProperty< CoulombNBPairs, 
                                             AtomPairs<CoulombScaleFactor> >
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const CoulombNBPairs&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, CoulombNBPairs&);

public:
    typedef AtomPairs<CoulombScaleFactor>::CGPairs CGPairs;

    CoulombNBPairs();

    CoulombNBPairs(const MoleculeInfoData &molinfo,
                   const CoulombScaleFactor &default_scale = CoulombScaleFactor(1));

    CoulombNBPairs(const MoleculeView &molview,
                   const CoulombScaleFactor &default_scale = CoulombScaleFactor(1));

    CoulombNBPairs(const CLJNBPairs &cljpairs);

    CoulombNBPairs(const CoulombNBPairs &other);

    ~CoulombNBPairs();

    CoulombNBPairs& operator=(const CoulombNBPairs &other);

    CoulombNBPairs& operator=(const CLJNBPairs &cljpairs);

    static const char* typeName();
    
    bool operator==(const CoulombNBPairs &other) const;
    bool operator!=(const CoulombNBPairs &other) const;
};

/** This class holds all of the non-bonded scale factors that are used
    to scale the intramolecular atom-atom Lennard-Jones
    interactions between pairs of atoms, e.g. for most MM forcefields,
    the scale factors for 1-1, 1-2 and 1-3 pairs are zero, the
    1-4 pairs are scaled by a LJ factor (e.g. 0.5 for OPLS)
    and the 1-5 and above pairs are not scaled (i.e. the factors equal 1)

    @author Christopher Woods
*/
class SIREMM_EXPORT LJNBPairs 
        : public SireBase::ConcreteProperty< LJNBPairs, 
                                             AtomPairs<LJScaleFactor> >
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const LJNBPairs&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, LJNBPairs&);

public:
    typedef AtomPairs<LJScaleFactor>::CGPairs CGPairs;

    LJNBPairs();

    LJNBPairs(const MoleculeView &molview,
              const LJScaleFactor &default_scale = LJScaleFactor(1));

    LJNBPairs(const MoleculeInfoData &molinfo,
              const LJScaleFactor &default_scale = LJScaleFactor(1));

    LJNBPairs(const CLJNBPairs &cljpairs);

    LJNBPairs(const LJNBPairs &other);

    ~LJNBPairs();

    LJNBPairs& operator=(const LJNBPairs &other);
    LJNBPairs& operator=(const CLJNBPairs &cljpairs);

    static const char* typeName();
    
    bool operator==(const LJNBPairs &other) const;
    bool operator!=(const LJNBPairs &other) const;
};

/** This class holds all of the non-bonded scale factors that are used
    to scale the intramolecular atom-atom coulomb and Lennard-Jones
    interactions between pairs of atoms, e.g. for most MM forcefields,
    the scale factors for 1-1, 1-2 and 1-3 pairs are zero, the
    1-4 pairs are scaled by a coulomb and LJ factor (e.g. 0.5 for OPLS)
    and the 1-5 and above pairs are not scaled (i.e. the coulomb and
    LJ factors equal 1)

    @author Christopher Woods
*/
class SIREMM_EXPORT CLJNBPairs 
        : public SireBase::ConcreteProperty< CLJNBPairs, 
                                             AtomPairs<CLJScaleFactor> >
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const CLJNBPairs&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, CLJNBPairs&);

public:
    typedef AtomPairs<CLJScaleFactor>::CGPairs CGPairs;

    CLJNBPairs();

    CLJNBPairs(const MoleculeView &molview,
               const CLJScaleFactor &default_scale = CLJScaleFactor(1,1));

    CLJNBPairs(const MoleculeInfoData &molinfo,
               const CLJScaleFactor &default_scale = CLJScaleFactor(1,1));

    CLJNBPairs(const SireMol::MoleculeInfo &molinfo,
               const CLJScaleFactor &default_scale = CLJScaleFactor(1,1));

    CLJNBPairs(const SireMol::Connectivity &connectivity,
               const CLJScaleFactor &scale14);

    CLJNBPairs(const CLJNBPairs &other);

    ~CLJNBPairs();

    CLJNBPairs& operator=(const CLJNBPairs &other);

    static const char* typeName();
    
    bool operator==(const CLJNBPairs &other) const;
    bool operator!=(const CLJNBPairs &other) const;
    
    QString toString() const;
    
    int nExcludedAtoms() const;
    QVector<AtomIdx> excludedAtoms() const;
    
    int nExcludedAtoms(const AtomID &atomid) const;
    QVector<AtomIdx> excludedAtoms(const AtomID &atomid) const;
    
};

}

Q_DECLARE_METATYPE(SireMM::CoulombScaleFactor)
Q_DECLARE_METATYPE(SireMM::LJScaleFactor)
Q_DECLARE_METATYPE(SireMM::CLJScaleFactor)

Q_DECLARE_METATYPE(SireMM::CLJNBPairs)
Q_DECLARE_METATYPE(SireMM::CoulombNBPairs)
Q_DECLARE_METATYPE(SireMM::LJNBPairs)

Q_DECLARE_TYPEINFO( SireMM::CoulombScaleFactor, Q_MOVABLE_TYPE );
Q_DECLARE_TYPEINFO( SireMM::LJScaleFactor, Q_MOVABLE_TYPE );
Q_DECLARE_TYPEINFO( SireMM::CLJScaleFactor, Q_MOVABLE_TYPE );

SIRE_EXPOSE_CLASS( SireMM::CoulombScaleFactor )
SIRE_EXPOSE_CLASS( SireMM::LJScaleFactor )
SIRE_EXPOSE_CLASS( SireMM::CLJScaleFactor )

SIRE_EXPOSE_CLASS( SireMM::CoulombNBPairs )
SIRE_EXPOSE_ALIAS( SireMM::AtomPairs<SireMM::CoulombScaleFactor>,
                   SireMM::AtomPairs_CoulombScaleFactor_ )

SIRE_EXPOSE_CLASS( SireMM::LJNBPairs )
SIRE_EXPOSE_ALIAS( SireMM::AtomPairs<SireMM::LJScaleFactor>,
                   SireMM::AtomPairs_LJScaleFactor_ )

SIRE_EXPOSE_CLASS( SireMM::CLJNBPairs )
SIRE_EXPOSE_ALIAS( SireMM::AtomPairs<SireMM::CLJScaleFactor>,
                   SireMM::AtomPairs_CLJScaleFactor_ )

SIRE_END_HEADER

#endif
