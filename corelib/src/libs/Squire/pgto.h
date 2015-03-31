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

#ifndef SQUIRE_PGTO_H
#define SQUIRE_PGTO_H

#include "gto.h"
#include "sgto.h"

#include "SireMaths/vector.h"
#include "SireMaths/matrix.h"
#include "SireMaths/nmatrix.h"

#include "SireBase/trigarray2d.hpp"
#include "SireBase/array2d.hpp"

SIRE_BEGIN_HEADER

namespace Squire
{
class P_GTO;

class PP_GTO;
class PS_GTO;
}

QDataStream& operator<<(QDataStream&, const Squire::P_GTO&);
QDataStream& operator>>(QDataStream&, Squire::P_GTO&);

QDataStream& operator<<(QDataStream&, const Squire::PS_GTO&);
QDataStream& operator>>(QDataStream&, Squire::PS_GTO&);

QDataStream& operator<<(QDataStream&, const Squire::PP_GTO&);
QDataStream& operator>>(QDataStream&, Squire::PP_GTO&);

namespace Squire
{

class PointCharge;
class PointDipole;

using SireMaths::Vector;
using SireMaths::Matrix;
using SireMaths::TrigMatrix;
using SireMaths::NMatrix;

using SireBase::TrigArray2D;
using SireBase::Array2D;

/** This is a single P-type shell of Gaussian Type Orbitals */
class SQUIRE_EXPORT P_GTO : public SireBase::ConcreteProperty<P_GTO,GTO>
{

friend QDataStream& ::operator<<(QDataStream&, const P_GTO&);
friend QDataStream& ::operator>>(QDataStream&, P_GTO&);

public:
    P_GTO();
    P_GTO(double alpha, double scale=1);
    
    P_GTO(const P_GTO &other);
    
    ~P_GTO();
    
    static const char* typeName();
    
    P_GTO& operator=(const P_GTO &other);
    
    bool operator==(const P_GTO &other) const;
    bool operator!=(const P_GTO &other) const;
    
    QString toString() const;
    
    int angularMomentum() const;
    int nOrbitals() const;
};

/** This is a combined S-P GTO shell pair */
class SQUIRE_EXPORT PS_GTO : public SireBase::ConcreteProperty<PS_GTO,GTOPair>
{

friend QDataStream& ::operator<<(QDataStream&, const PS_GTO&);
friend QDataStream& ::operator>>(QDataStream&, PS_GTO&);

public:
    PS_GTO();
    PS_GTO(const Vector &A, const S_GTO &a,
           const Vector &B, const P_GTO &b);
    PS_GTO(const Vector &A, const P_GTO &a,
           const Vector &B, const S_GTO &b);
           
    PS_GTO(const PS_GTO &other);
    
    ~PS_GTO();
    
    static const char* typeName();
    
    PS_GTO& operator=(const PS_GTO &other);
    
    bool operator==(const PS_GTO &other) const;
    bool operator!=(const PS_GTO &other) const;
    
    const Vector& P_minus_A() const;
    const Vector& P_minus_B() const;
    
    const Vector& Q_minus_C() const;
    const Vector& Q_minus_D() const;

    double scale() const;

    int angularMomentum0() const;
    int angularMomentum1() const;
    
    int nOrbitals0() const;
    int nOrbitals1() const;

private:
    /** The vector from the center of the P-orbital to
        the 'center of mass' of the shell-pair */
    Vector p_minus_a;
    
    /** The extra factor for the scaling constant */
    double norm_scl;
};

/** This is a combined P-P GTO shell pair */
class SQUIRE_EXPORT PP_GTO : public SireBase::ConcreteProperty<PP_GTO,GTOPair>
{

friend QDataStream& ::operator<<(QDataStream&, const PP_GTO&);
friend QDataStream& ::operator>>(QDataStream&, PP_GTO&);

public:
    PP_GTO();
    PP_GTO(const Vector &A, const P_GTO &a,
           const Vector &B, const P_GTO &b);
           
    PP_GTO(const PP_GTO &other);
    
    ~PP_GTO();
    
    static const char* typeName();
    
    PP_GTO& operator=(const PP_GTO &other);
    
    bool operator==(const PP_GTO &other) const;
    bool operator!=(const PP_GTO &other) const;

    const Vector& P_minus_A() const;
    const Vector& P_minus_B() const;

    const Vector& Q_minus_C() const;
    const Vector& Q_minus_D() const;

    double scale() const;

    int angularMomentum0() const;
    int angularMomentum1() const;
    
    int nOrbitals0() const;
    int nOrbitals1() const;

private:
    /** The vector from the center of the first P orbital to the 
        center of mass of the gaussian */
    Vector p_minus_a;
    
    /** The vector from the center of the second P orbital to the
        center of mass of the gaussian */
    Vector p_minus_b;
    
    /** The extra factor for the scaling constant */
    double norm_scl;
};

/** This class is used to calculate integrals involving PS orbital
    pairs
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT PS_GTOs
{
public:
    PS_GTOs();

    PS_GTOs(const QVector<S_GTO> &s_gtos,
            const QVector<Vector> &s_centers,
            const QVector<P_GTO> &p_gtos,
            const QVector<Vector> &p_centers);

    PS_GTOs(const QVector<P_GTO> &p_gtos,
            const QVector<Vector> &p_centers,
            const QVector<S_GTO> &s_gtos,
            const QVector<Vector> &s_centers);
            
    PS_GTOs(const PS_GTOs &other);
    
    ~PS_GTOs();
    
    PS_GTOs& operator=(const PS_GTOs &other);
    
    bool operator==(const PS_GTOs &other) const;
    bool operator!=(const PS_GTOs &other) const;
    
    NMatrix overlap_integral() const;
    NMatrix kinetic_integral() const;
    
    NMatrix potential_integral(const QVector<PointCharge> &C) const;
    NMatrix potential_integral(const QVector<PointCharge> &C, int m) const;

private:
    void _pvt_create(const QVector<S_GTO> &s_gtos,
                     const QVector<Vector> &s_centers,
                     const QVector<P_GTO> &p_gtos,
                     const QVector<Vector> &p_centers);

    /** All of the orbital pairs */
    Array2D<PS_GTO> orbs;
};

/** This class is used to calculate integrals involving PP orbital pairs */
class SQUIRE_EXPORT PP_GTOs
{
public:
    PP_GTOs();

    PP_GTOs(const QVector<P_GTO> &p_gtos,
            const QVector<Vector> &p_centers);
            
    PP_GTOs(const PP_GTOs &other);
    
    ~PP_GTOs();
    
    PP_GTOs& operator=(const PP_GTOs &other);
    
    bool operator==(const PP_GTOs &other) const;
    bool operator!=(const PP_GTOs &other) const;
    
    TrigMatrix overlap_integral() const;
    TrigMatrix kinetic_integral() const;
    
    TrigMatrix potential_integral(const QVector<PointCharge> &C) const;
    TrigMatrix potential_integral(const QVector<PointCharge> &C, int m) const;

private:
    /** All of the orbital pairs */
    TrigArray2D<PP_GTO> orbs;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the vector from the center of the p orbital shell
    to the center of mass of the combined gaussian */
inline const Vector& PS_GTO::P_minus_A() const
{
    return p_minus_a;
}

/** Synonym for P_minus_A */
inline const Vector& PS_GTO::P_minus_B() const
{
    return PS_GTO::P_minus_A();
}

/** Synonym for P_minus_A */
inline const Vector& PS_GTO::Q_minus_C() const
{
    return PS_GTO::P_minus_A();
}

/** Synonym for P_minus_A */
inline const Vector& PS_GTO::Q_minus_D() const
{
    return PS_GTO::P_minus_A();
}

/** Return the additional scaling constant needed to normalise
    the integrals */
inline double PS_GTO::scale() const
{
    return norm_scl;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

//////////
////////// Integral functions involving P and S orbitals
//////////

Vector kinetic_integral(const PS_GTO &P);
Vector overlap_integral(const PS_GTO &P);

Vector potential_integral(const QVector<PointCharge> &C, const PS_GTO &P);
Vector potential_integral(const QVector<PointDipole> &C, const PS_GTO &P);

Vector potential_integral(const PointCharge &C, const PS_GTO &P);
Vector potential_integral(const PointDipole &C, const PS_GTO &P);

Vector potential_integral(const QVector<PointCharge> &C, const PS_GTO &P, int m);
Vector potential_integral(const QVector<PointDipole> &C, const PS_GTO &P, int m);

Vector potential_integral(const PointCharge &C, const PS_GTO &P, int m);
Vector potential_integral(const PointDipole &C, const PS_GTO &P, int m);

Matrix kinetic_integral(const PP_GTO &P);
Matrix overlap_integral(const PP_GTO &P);

Matrix potential_integral(const QVector<PointCharge> &C, const PP_GTO &P);
Matrix potential_integral(const QVector<PointDipole> &C, const PP_GTO &P);

Matrix potential_integral(const PointCharge &C, const PP_GTO &P);
Matrix potential_integral(const PointDipole &C, const PP_GTO &P);

Matrix potential_integral(const QVector<PointCharge> &C, const PP_GTO &P, int m);
Matrix potential_integral(const QVector<PointDipole> &C, const PP_GTO &P, int m);

Matrix potential_integral(const PointCharge &C, const PP_GTO &P, int m);
Matrix potential_integral(const PointDipole &C, const PP_GTO &P, int m);

///// Electron integrals involving S and P orbitals

Vector electron_integral(const PS_GTO &P, const SS_GTO &Q);
Vector electron_integral(const SS_GTO &P, const PS_GTO &Q);

Matrix electron_integral(const PS_GTO &P, const PS_GTO &Q);

Matrix electron_integral(const PP_GTO &P, const SS_GTO &Q);
Matrix electron_integral(const SS_GTO &P, const PP_GTO &Q);

SireBase::Array2D<Vector> electron_integral(const PP_GTO &P, const PS_GTO &Q);
SireBase::Array2D<Vector> electron_integral(const PS_GTO &P, const PP_GTO &Q);

void electron_integral(const PP_GTO &P, const PS_GTO &Q,
                       SireBase::Array2D<Vector> &matrix);
                       
void electron_integral(const PS_GTO &P, const PP_GTO &Q,
                       SireBase::Array2D<Vector> &matrix);

SireBase::Array2D<Matrix> electron_integral(const PP_GTO &P, const PP_GTO &Q);

void electron_integral(const PP_GTO &P, const PP_GTO &Q, 
                       SireBase::Array2D<Matrix> &matrix);

Vector electron_integral(const PS_GTO &P, const SS_GTO &Q, int m);
Vector electron_integral(const SS_GTO &P, const PS_GTO &Q, int m);

Matrix electron_integral(const PS_GTO &P, const PS_GTO &Q, int m);

Matrix electron_integral(const PP_GTO &P, const SS_GTO &Q, int m);
Matrix electron_integral(const SS_GTO &P, const PP_GTO &Q, int m);

SireBase::Array2D<Vector> 
electron_integral(const PP_GTO &P, const PS_GTO &Q, int m);
SireBase::Array2D<Vector> 
electron_integral(const PS_GTO &P, const PP_GTO &Q, int m);

void electron_integral(const PP_GTO &P, const PS_GTO &Q, int m,
                       SireBase::Array2D<Vector> &matrix);
                       
void electron_integral(const PS_GTO &P, const PP_GTO &Q, int m,
                       SireBase::Array2D<Vector> &matrix);

SireBase::Array2D<Matrix> 
electron_integral(const PP_GTO &P, const PP_GTO &Q, int m);

void electron_integral(const PP_GTO &P, const PP_GTO &Q, int m, 
                       SireBase::Array2D<Matrix> &matrix);

}

Q_DECLARE_METATYPE( Squire::P_GTO )
Q_DECLARE_METATYPE( Squire::PS_GTO )
Q_DECLARE_METATYPE( Squire::PP_GTO )

SIRE_EXPOSE_CLASS( Squire::P_GTO )
SIRE_EXPOSE_CLASS( Squire::PS_GTO )
SIRE_EXPOSE_CLASS( Squire::PP_GTO )
SIRE_EXPOSE_CLASS( Squire::PS_GTOs )
SIRE_EXPOSE_CLASS( Squire::PP_GTOs )

SIRE_EXPOSE_FUNCTION( Squire::kinetic_integral )
SIRE_EXPOSE_FUNCTION( Squire::overlap_integral )
SIRE_EXPOSE_FUNCTION( Squire::potential_integral )
SIRE_EXPOSE_FUNCTION( Squire::electron_integral )

SIRE_END_HEADER

#endif
