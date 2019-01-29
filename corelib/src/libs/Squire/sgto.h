/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SQUIRE_SGTO_H
#define SQUIRE_SGTO_H

#include "gto.h"

#include "SireBase/trigarray2d.hpp"

#include "SireMaths/trigmatrix.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class S_GTO;
class SS_GTO;
}

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::S_GTO&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::S_GTO&);

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::SS_GTO&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::SS_GTO&);

namespace Squire
{

using SireMaths::TrigMatrix;
using SireBase::TrigArray2D;

class PointCharge;
class PointDipole;

/** This is a single S-type Gaussian Type Orbital shell */
class SQUIRE_EXPORT S_GTO : public SireBase::ConcreteProperty<S_GTO,GTO>
{
public:
    S_GTO();
    S_GTO(double alpha, double scale=1);
    
    S_GTO(const S_GTO &other);
    
    ~S_GTO();
    
    static const char* typeName();
    
    S_GTO& operator=(const S_GTO &other);
    
    bool operator==(const S_GTO &other) const;
    bool operator!=(const S_GTO &other) const;
    
    QString toString() const;
    
    int angularMomentum() const;
    int nOrbitals() const;
};

/** This is the combined SS shell pair composed from two S-type
    GTO shells */
class SQUIRE_EXPORT SS_GTO : public SireBase::ConcreteProperty<SS_GTO,GTOPair>
{

friend QDataStream& ::operator<<(QDataStream&, const SS_GTO&);
friend QDataStream& ::operator>>(QDataStream&, SS_GTO&);

public:
    SS_GTO();
    SS_GTO(const Vector &A, const S_GTO &a,
           const Vector &B, const S_GTO &b);

    SS_GTO(const SS_GTO &other);
    
    ~SS_GTO();
    
    static const char* typeName();
    
    SS_GTO& operator=(const SS_GTO &other);
    
    bool operator==(const SS_GTO &other) const;
    bool operator!=(const SS_GTO &other) const;

    int angularMomentum0() const;
    int angularMomentum1() const;
    
    int nOrbitals0() const;
    int nOrbitals1() const;
};    

/** This class is used to calculate integrals involving just SS pairs.
    It is used to aid in the calculation of integrals.

	@author Christopher Woods
*/
class SQUIRE_EXPORT SS_GTOs
{
public:
	SS_GTOs();
    
    SS_GTOs(const QVector<S_GTO> &s_gtos, 
        	const QVector<Vector> &centers);
    
    SS_GTOs(const SS_GTOs &other);
    
    ~SS_GTOs();
    
    SS_GTOs& operator=(const SS_GTOs &other);

	TrigMatrix overlap_integral() const;
    TrigMatrix kinetic_integral() const;
    
    TrigMatrix potential_integral(const QVector<PointCharge> &C) const;
    TrigMatrix potential_integral(const QVector<PointCharge> &C, int m) const;
    
    TrigMatrix coulomb_integral(const SS_GTOs &other) const;
    TrigMatrix exchange_integral(const SS_GTOs &other) const;
    
private:
	/** All of the orbital pairs */
    TrigArray2D<SS_GTO> orbs;
};

//////////
////////// Integrals involving only S-orbitals
//////////

SQUIRE_EXPORT double kinetic_integral(const SS_GTO &P);
SQUIRE_EXPORT double overlap_integral(const SS_GTO &P);

SQUIRE_EXPORT double potential_integral(const QVector<PointCharge> &C, const SS_GTO &P);
SQUIRE_EXPORT double potential_integral(const QVector<PointDipole> &C, const SS_GTO &P);

SQUIRE_EXPORT double potential_integral(const QVector<PointCharge> &C, const SS_GTO &P, int m);
SQUIRE_EXPORT double potential_integral(const QVector<PointDipole> &C, const SS_GTO &P, int m);

SQUIRE_EXPORT double potential_integral(const PointCharge &Q, const SS_GTO &P);
SQUIRE_EXPORT double potential_integral(const PointDipole &Q, const SS_GTO &P);

SQUIRE_EXPORT double potential_integral(const PointCharge &C, const SS_GTO &P, int m);
SQUIRE_EXPORT double potential_integral(const PointDipole &C, const SS_GTO &P, int m);

SQUIRE_EXPORT double electron_integral(const SS_GTO &P, const SS_GTO &Q);
SQUIRE_EXPORT double electron_integral(const SS_GTO &P, const SS_GTO &Q, int m);

}

Q_DECLARE_METATYPE( Squire::S_GTO )
Q_DECLARE_METATYPE( Squire::SS_GTO )

SIRE_EXPOSE_CLASS( Squire::S_GTO )
SIRE_EXPOSE_CLASS( Squire::SS_GTO )
SIRE_EXPOSE_CLASS( Squire::SS_GTOs )

SIRE_EXPOSE_FUNCTION( Squire::kinetic_integral )
SIRE_EXPOSE_FUNCTION( Squire::potential_integral )
SIRE_EXPOSE_FUNCTION( Squire::overlap_integral )
SIRE_EXPOSE_FUNCTION( Squire::electron_integral )

SIRE_END_HEADER

#endif
