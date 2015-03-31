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

#include "pgto.h"
#include "pointcharge.h"
#include "pointdipole.h"

#include "SireMaths/boys.h"

#include "SireBase/array2d.hpp"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

//////////
////////// Implementation of P_GTO
//////////

static const RegisterMetaType<P_GTO> r_pgto;

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const P_GTO &p)
{
    writeHeader(ds, r_pgto, 1);
    
    ds << static_cast<const GTO&>(p);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, P_GTO &p)
{
    VersionID v = readHeader(ds, r_pgto);
    
    if (v == 1)
    {
        ds >> static_cast<GTO&>(p);
    }
    else
        throw version_error(v, "1", r_pgto, CODELOC);
        
    return ds;
}

/** Constructor */
P_GTO::P_GTO() : ConcreteProperty<P_GTO,GTO>()
{}

/** Construct with the specified alpha value and scaling factor */
P_GTO::P_GTO(double alpha, double scale)
      : ConcreteProperty<P_GTO,GTO>(alpha, scale)
{}

/** Copy constructor */
P_GTO::P_GTO(const P_GTO &other) : ConcreteProperty<P_GTO,GTO>(other)
{}

/** Destructor */
P_GTO::~P_GTO()
{}

const char* P_GTO::typeName()
{
    return QMetaType::typeName( qMetaTypeId<P_GTO>() );
}

/** Copy assignment operator */
P_GTO& P_GTO::operator=(const P_GTO &other)
{
    GTO::operator=(other);
    return *this;
}

/** Comparison operator */
bool P_GTO::operator==(const P_GTO &other) const
{
    return GTO::operator==(other);
}

/** Comparison operator */
bool P_GTO::operator!=(const P_GTO &other) const
{
    return GTO::operator!=(other);
}

/** Return a string representation of this orbital */
QString P_GTO::toString() const
{
	if (GTO::scale() == 1)
		return QObject::tr("P (alpha = %1)").arg( GTO::alpha() );
	else
		return QObject::tr("P (alpha = %1, scale = %2)")
        			.arg( GTO::alpha() ).arg(GTO::scale());
}

/** The angular momentum of P-GTOs is 1 */
int P_GTO::angularMomentum() const
{
    return 1;
}

/** There are 3 P-orbitals per shell (px, py and pz) */
int P_GTO::nOrbitals() const
{
    return 3;
}

//////////
////////// Implementation of PS_GTO
//////////

static const RegisterMetaType<PS_GTO> r_psgto;

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const PS_GTO &ps)
{
    writeHeader(ds, r_psgto, 1);
    
    ds << ps.p_minus_a << ps.norm_scl << static_cast<const GTOPair&>(ps);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, PS_GTO &ps)
{
    VersionID v = readHeader(ds, r_psgto);
    
    if (v == 1)
    {
        ds >> ps.p_minus_a >> ps.norm_scl >> static_cast<GTOPair&>(ps);
    }
    else
        throw version_error(v, "1", r_psgto, CODELOC);
        
    return ds;
}

/** Constructor */
PS_GTO::PS_GTO() 
       : ConcreteProperty<PS_GTO,GTOPair>(), norm_scl(0)
{}

static double getQ( const Matrix &e )
{
    const double *mat = e.constData();
    
    double max_e = 0;
    
    for (int i=0; i<9; ++i)
    {
    	max_e = qMax(max_e, mat[i]);
    }

	return std::sqrt(max_e);
}

/** Construct combining orbital 'a' at position 'A' with orbital 'b' at
    position 'B' */
PS_GTO::PS_GTO(const Vector &A, const P_GTO &a,
               const Vector &B, const S_GTO &b)
       : ConcreteProperty<PS_GTO,GTOPair>(A, a, B, b)
{
    p_minus_a = P() - A;
    norm_scl = std::sqrt(4 * a.alpha());
    
    GTOPair::setQ( ::getQ(electron_integral(*this,*this)) );
}

/** Construct combining orbital 'a' at position 'A' with orbital 'b' at
    position 'B' */
PS_GTO::PS_GTO(const Vector &A, const S_GTO &a,
               const Vector &B, const P_GTO &b)
       : ConcreteProperty<PS_GTO,GTOPair>(B, b, A, a)
{
    p_minus_a = P() - B;
    norm_scl = std::sqrt(4 * b.beta());
    
    GTOPair::setQ( ::getQ(electron_integral(*this,*this)) );
}
  
/** Copy constructor */     
PS_GTO::PS_GTO(const PS_GTO &other)
       : ConcreteProperty<PS_GTO,GTOPair>(other),
         p_minus_a(other.p_minus_a), norm_scl(other.norm_scl)
{}

/** Destructor */
PS_GTO::~PS_GTO()
{}

const char* PS_GTO::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PS_GTO>() );
}

/** Copy assignment operator */
PS_GTO& PS_GTO::operator=(const PS_GTO &other)
{
    if (this != &other)
    {
        GTOPair::operator=(other);
        p_minus_a = other.p_minus_a;
        norm_scl = other.norm_scl;
    }
    
    return *this;
}

/** Comparison operator */
bool PS_GTO::operator==(const PS_GTO &other) const
{
    return GTOPair::operator==(other) and p_minus_a == other.p_minus_a
             and norm_scl == other.norm_scl;
}

/** Comparison operator */
bool PS_GTO::operator!=(const PS_GTO &other) const
{
    return not this->operator==(other);
}

/** Return the angular momentum of the first GTO shell in this pair */
int PS_GTO::angularMomentum0() const
{
    return 1;
}

/** Return the angular momentum of the second GTO shell in this pair */
int PS_GTO::angularMomentum1() const
{
    return 0;
}

/** Return the number of orbitals in the first GTO shell in this pair */
int PS_GTO::nOrbitals0() const
{
    return 3;
}

/** Return the number of orbitals in the second GTO shell in this pair */
int PS_GTO::nOrbitals1() const
{
    return 1;
}

//////////
////////// Implementation of PP_GTO
//////////

static const RegisterMetaType<PP_GTO> r_ppgto;

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const PP_GTO &pp)
{
    writeHeader(ds, r_ppgto, 1);
    
    ds << pp.p_minus_a << pp.p_minus_b << pp.norm_scl
       << static_cast<const GTOPair&>(pp);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, PP_GTO &pp)
{
    VersionID v = readHeader(ds, r_ppgto);
    
    if (v == 1)
    {
        ds >> pp.p_minus_a >> pp.p_minus_b >> pp.norm_scl
           >> static_cast<GTOPair&>(pp);
    }
    else
        throw version_error(v, "1", r_ppgto, CODELOC);
        
    return ds;
}

/** Constructor */
PP_GTO::PP_GTO() : ConcreteProperty<PP_GTO,GTOPair>(), norm_scl(0)
{}

/** Construct combining orbital 'a' at position 'A' with orbital 'b' at
    position 'B' */
PP_GTO::PP_GTO(const Vector &A, const P_GTO &a,
               const Vector &B, const P_GTO &b)
       : ConcreteProperty<PP_GTO,GTOPair>(A, a, B, b)
{
    p_minus_a = P() - A;
    p_minus_b = P() - B;

    norm_scl = std::sqrt( 16 * a.alpha() * b.beta() );
    
    Array2D<Matrix> e = electron_integral(*this, *this);
    
    double max_e = 0;
    
    for (int i=0; i<3; ++i)
    {
    	for (int j=0; j<3; ++j)
        {
        	const double *mat = e.constData()[ e.offset(i,j) ].constData();
            
            for (int k=0; k<9; ++k)
            {
             	max_e = qMax(max_e, mat[k]);
            }
        }
    }
    
    GTOPair::setQ( std::sqrt(max_e) );
}
  
/** Copy constructor */     
PP_GTO::PP_GTO(const PP_GTO &other)
       : ConcreteProperty<PP_GTO,GTOPair>(other),
         p_minus_a(other.p_minus_a), p_minus_b(other.p_minus_b), 
         norm_scl(other.norm_scl)
{}

/** Destructor */
PP_GTO::~PP_GTO()
{}

const char* PP_GTO::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PP_GTO>() );
}

/** Copy assignment operator */
PP_GTO& PP_GTO::operator=(const PP_GTO &other)
{
    if (this != &other)
    {
        GTOPair::operator=(other);
        p_minus_a = other.p_minus_a;
        p_minus_b = other.p_minus_b;
        norm_scl = other.norm_scl;
    }
    
    return *this;
}

/** Comparison operator */
bool PP_GTO::operator==(const PP_GTO &other) const
{
    return GTOPair::operator==(other) and p_minus_a == other.p_minus_a and
           p_minus_b == other.p_minus_b and norm_scl == other.norm_scl;
}

/** Comparison operator */
bool PP_GTO::operator!=(const PP_GTO &other) const
{
    return not this->operator==(other);
}

/** Return the vector from the center of the first P orbital shell
    to the center of mass of the combined gaussian */
const Vector& PP_GTO::P_minus_A() const
{
    return p_minus_a;
}

/** Return the vector from the center of the second P orbital shell
    to the center of mass of the combined gaussian */
const Vector& PP_GTO::P_minus_B() const
{
    return p_minus_b;
}

/** Synonym for P_minus_A */
const Vector& PP_GTO::Q_minus_C() const
{
    return PP_GTO::P_minus_A();
}

/** Synonym for Q_minus_D */
const Vector& PP_GTO::Q_minus_D() const
{
    return PP_GTO::P_minus_B();
}

/** Return the additional scaling constant needed to normalise the
    integrals involving this shell-pair */
double PP_GTO::scale() const
{
    return norm_scl;
}

/** Return the angular momentum of the first GTO shell in this pair */
int PP_GTO::angularMomentum0() const
{
    return 1;
}

/** Return the angular momentum of the second GTO shell in this pair */
int PP_GTO::angularMomentum1() const
{
    return 1;
}

/** Return the number of orbitals in the first GTO shell in this pair */
int PP_GTO::nOrbitals0() const
{
    return 3;
}

/** Return the number of orbitals in the second GTO shell in this pair */
int PP_GTO::nOrbitals1() const
{
    return 3;
}

//////////
////////// Integrals
//////////

namespace Squire
{

// Integrals derived using Obara-Saika recursion - see Mathematica
// notebook (techdocs/obara_saika.nb) for details

//////
////// The overlap, kinetic and potential integrals
//////

/** Return the overlap integral (p||s) */
Vector SQUIRE_EXPORT overlap_integral(const PS_GTO &ps)
{
    // (pi||s) = (Pi - Ai)(s||s)
    return (ps.scale() * ps.ss()) * ps.P_minus_A();
}

/** Return the kinetic energy integral, (p|nabla|s) */
Vector SQUIRE_EXPORT kinetic_integral(const PS_GTO &ps)
{
    // (p|nabla|s) = ss xi (5 - 2 R^2 xi) (Pi-Ai)
    return ( ps.scale() * ps.ss() * ps.xi() * 
              (5 - 2*ps.R2()*ps.xi()) ) * ps.P_minus_A();
}

/** Return the potential energy integral with an array of point charges,
    (p|C|s) */
Vector SQUIRE_EXPORT potential_integral(const QVector<PointCharge> &C, const PS_GTO &P)
{
    // (p|C|s) = 2 ss Sqrt[zeta/pi] ( F0(U) (Pi-Ai) - F1(U) (Pi-Ci) )
    // U = zeta * (P-C)^2

    const int nc = C.count();
    const PointCharge *c = C.constData();

    Vector total;
    
    for (int i=0; i<nc; ++i)
    {
        const Vector P_minus_C = P.P() - c[i].center();
        const double U = P.zeta() * P_minus_C.length2();

        double boys[2];
        multi_boys_2(U, boys);

        total += c[i].charge() * (boys[0] * P.P_minus_A() - boys[1] * P_minus_C);
    }

    return (-2 * P.scale() * P.ss() * std::sqrt(P.zeta() * one_over_pi)) * total;
}

/** Return the potential energy integral with a single point charge, (p|C|s) */
Vector SQUIRE_EXPORT potential_integral(const PointCharge &C, const PS_GTO &ps)
{
    QVector<PointCharge> Cs(1, C);
    return potential_integral(Cs, ps);
}

/** Return the potential energy integral with an array of point charges,
    (p|C|s)^m */
Vector SQUIRE_EXPORT potential_integral(const QVector<PointCharge> &C, 
                                        const PS_GTO &P, int m)
{
    // (p|C|s)^m = 2 ss Sqrt[zeta/pi] ( Fm(U) (Pi-Ai) - Fm+1(U) (Pi-Ci) )
    // U = zeta * (P-C)^2

    const int nc = C.count();
    const PointCharge *c = C.constData();

    Vector total;
    
    for (int i=0; i<nc; ++i)
    {
        const Vector P_minus_C = P.P() - c[i].center();
        const double U = P.zeta() * P_minus_C.length2();

        double boys[2];
        multi_boys_2(U, boys, m);

        total += c[i].charge() * (boys[0] * P.P_minus_A() - boys[1] * P_minus_C);
    }

    return (-2 * P.scale() * P.ss() * std::sqrt(P.zeta() * one_over_pi)) * total;
}

/** Return the potential energy integral with a single point charge, (p|C|s)^m */
Vector SQUIRE_EXPORT potential_integral(const PointCharge &C, 
                                        const PS_GTO &ps, int m)
{
    QVector<PointCharge> Cs(1, C);
    return potential_integral(Cs, ps, m);
}

/** Return the potential energy integral with an array of point dipoles, (p|C|s) */
Vector SQUIRE_EXPORT potential_integral(const QVector<PointDipole> &C,
                                        const PS_GTO &ps)
{
    throw SireError::incomplete_code( "Not implemented", CODELOC );
    return Vector();
}

/** Return the potential energy integral with a single point dipole, (p|C|s) */
Vector SQUIRE_EXPORT potential_integral(const PointDipole &C, const PS_GTO &ps)
{
    QVector<PointDipole> Cs(1, C);
    return potential_integral(Cs, ps);
}

/** Return the potential energy integral with an array of point dipoles, (p|C|s)^m */
Vector SQUIRE_EXPORT potential_integral(const QVector<PointDipole> &C,
                                        const PS_GTO &ps, int m)
{
    throw SireError::incomplete_code( "Not implemented", CODELOC );
    return Vector();
}

/** Return the potential energy integral with a single point dipole, (p|C|s)^m */
Vector SQUIRE_EXPORT potential_integral(const PointDipole &C, 
                                        const PS_GTO &ps, int m)
{
    QVector<PointDipole> Cs(1, C);
    return potential_integral(Cs, ps, m);
}

/** Return the kinetic energy integral (p|nabla|p) */
Matrix SQUIRE_EXPORT kinetic_integral(const PP_GTO &pp)
{
    // (p|k|p) = -(1/2) (s||s) xi ( 2 (-7 + 2 R^2 xi) (Pi-Ai) (Pj-Bj) + 
    //                              delta_ij zeta (-5 + 2 R^2 xi) ) 

    const double prefac0 = -0.5 * pp.scale() * pp.ss() * pp.xi();
                            
    
    const double prefac1 = prefac0 * ( 2 * (-7 + 2*pp.R2()*pp.xi()) );
    
    Matrix mat;
    double *m = mat.data();

    const double *pa = pp.P_minus_A().constData();
    const double *pb = pp.P_minus_B().constData();
                                                                
    //do the off-diagonal elements
    for (int i=0; i<2; ++i)
    {
        const double pa_prefac = pa[i] * prefac1;
    
        for (int j=i+1; j<3; ++j)
        {
            const double val = pa_prefac * pb[j];
            m[ mat.offset(i,j) ] = val;
            m[ mat.offset(j,i) ] = val;
        }
    }
    
    //do the diagonal elements
    const double extra = prefac0 *
                           pp.zeta() * (-5 + 2 * pp.R2() * pp.xi());
    
    for (int i=0; i<3; ++i)
    {
        m[ mat.offset(i,i) ] = (prefac1 * pa[i] * pb[i]) + extra;
    }
    
    return mat;
}

/** Return the overlap integral (p||p) */
Matrix SQUIRE_EXPORT overlap_integral(const PP_GTO &pp)
{
    // (p||p) = ss (Pi-Ai) (Pj-Bj) + delta_ij (1/2 ss zeta) 
    Matrix mat;
    
    const double prefac = pp.scale() * pp.ss();
        
    const double *pa = pp.P_minus_A().constData();
    const double *pb = pp.P_minus_B().constData();
   
    double *m = mat.data();
    
    //do the off-diagonal elements
    for (int i=0; i<2; ++i)
    {
        const double prefac_pa = prefac * pa[i];
    
        for (int j=i+1; j<3; ++j)
        {
            const double val = prefac_pa * pb[j];
            m[ mat.offset(i,j) ] = val;
            m[ mat.offset(j,i) ] = val;
        }
    }

    //do the diagonal elements
    const double extra = 1 / (2 * pp.zeta());

    for (int i=0; i<3; ++i)
    {
        m[ mat.offset(i,i) ] = prefac * (pa[i]*pb[i] + extra);
    }

    return mat;
}

/** Return the potential integral with an array of point charges, (p|C|p) */
Matrix SQUIRE_EXPORT potential_integral(const QVector<PointCharge> &C, 
                                        const PP_GTO &P)
{
    // (p|C|p) = -2 * Q * ss * Sqrt[zeta/pi] * {
    //                             F0(U) (Pi-Ai)(Pj-Bj) +
    //                             F2(U) (Pi-Ci)(Pj-Cj) - 
    //                             F1(U)[ (Pj-Bj)(Pi-Ci) + (Pi-Ai)(Pj-Cj) ]
    //                                     + delta_ij( 0.5*zeta( F0(U) - F1(U) ) }

    // U = zeta * (P-C)^2

    const int nc = C.count();
    const PointCharge *c = C.constData();
    
    const double prefac = -2 * P.scale() * P.ss() * std::sqrt(P.zeta()*one_over_pi);
    
    Matrix mat(double(0));
    double *m = mat.data();
    
    const double *pa = P.P_minus_A().constData();
    const double *pb = P.P_minus_B().constData();
    
    for (int ic=0; ic<nc; ++ic)
    {
        const double U = P.zeta() * Vector::distance2(P.P(), c[ic].center());
        
        double boys[3];
        multi_boys_3(U, boys);
        
        const double q_prefac = c[ic].charge() * prefac;
        const Vector P_minus_C = P.P() - c[ic].center();
        const double *pc = P_minus_C.constData();
        
        //do the off-diagonal elements
        for (int i=0; i<2; ++i)
        {
            for (int j=i+1; j<3; ++j)
            {
                const double val = q_prefac * 
                                    ( boys[0]*pa[i]*pb[j] +
                                      boys[2]*pc[i]*pc[j] -
                                      boys[1]*(pb[j]*pc[i] + pa[i]*pc[j]) );
                
                m[ mat.offset(i,j) ] += val;
            }
        }
        
        //do the diagonal elements
        const double extra = 0.5 * P.zeta() * (boys[0] - boys[1]);
        
        for (int i=0; i<3; ++i)
        {
            m[ mat.offset(i,i) ] = q_prefac * 
                                    ( boys[0]*pa[i]*pb[i] +
                                      boys[2]*pc[i]*pc[i] -
                                      boys[1]*(pb[i]*pc[i] + pa[i]*pc[i]) + extra );
        }
    }

    //now copy mat[i,j] to mat[j,i]
    for (int i=0; i<2; ++i)
    {
        for (int j=i+1; j<3; ++j)
        {
            m[ mat.offset(j,i) ] = m[ mat.offset(i,j) ];
        }
    }
    
    return mat;
}

/** Return the potential integral with a single point charge (p|C|p) */
Matrix SQUIRE_EXPORT potential_integral(const PointCharge &C, const PP_GTO &P)
{
    QVector<PointCharge> Cs(1, C);
    return potential_integral(Cs, P);
}

/** Return the potential integral with an array of point charges, (p|C|p)^m */
Matrix SQUIRE_EXPORT potential_integral(const QVector<PointCharge> &C, 
                                        const PP_GTO &P, int M)
{
    // (p|C|p)^m = -2 * Q * ss * Sqrt[zeta/pi] * {
    //                             Fm(U) (Pi-Ai)(Pj-Bj) +
    //                             Fm+2(U) (Pi-Ci)(Pj-Cj) - 
    //                             Fm+1(U)[ (Pj-Bj)(Pi-Ci) + (Pi-Ai)(Pj-Cj) ]
    //                                       + delta_ij( 0.5*zeta( F0(U) - F1(U) ) }

    // U = zeta * (P-C)^2

    const int nc = C.count();
    const PointCharge *c = C.constData();
    
    const double prefac = -2 * P.scale() * P.ss() * std::sqrt(P.zeta()*one_over_pi);
    
    Matrix mat(double(0));
    double *m = mat.data();
    
    const double *pa = P.P_minus_A().constData();
    const double *pb = P.P_minus_B().constData();
    
    for (int ic=0; ic<nc; ++ic)
    {
        const double U = P.zeta() * Vector::distance2(P.P(), c[ic].center());
        
        double boys[3];
        multi_boys_3(U, boys, M);
        
        const double q_prefac = c[ic].charge() * prefac;
        const Vector P_minus_C = P.P() - c[ic].center();
        const double *pc = P_minus_C.constData();
        
        //do the off-diagonal elements
        for (int i=0; i<2; ++i)
        {
            for (int j=i+1; j<3; ++j)
            {
                const double val = q_prefac * 
                                    ( boys[0]*pa[i]*pb[j] +
                                      boys[2]*pc[i]*pc[j] -
                                      boys[1]*(pb[j]*pc[i] + pa[i]*pc[j]) );
                
                m[ mat.offset(i,j) ] += val;
            }
        }
        
        //do the diagonal elements
        const double extra = 0.5 * P.zeta() * (boys[0] - boys[1]);
        
        for (int i=0; i<3; ++i)
        {
            m[ mat.offset(i,i) ] = q_prefac * 
                                    ( boys[0]*pa[i]*pb[i] +
                                      boys[2]*pc[i]*pc[i] -
                                      boys[1]*(pb[i]*pc[i] + pa[i]*pc[i]) + extra );
        }
    }

    //now copy mat[i,j] to mat[j,i]
    for (int i=0; i<2; ++i)
    {
        for (int j=i+1; j<3; ++j)
        {
            m[ mat.offset(j,i) ] = m[ mat.offset(i,j) ];
        }
    }
    
    return mat;
}

/** Return the potential integral with a single point charge (p|C|p)^m */
Matrix SQUIRE_EXPORT potential_integral(const PointCharge &C, const PP_GTO &P, int m)
{
    QVector<PointCharge> Cs(1, C);
    return potential_integral(Cs, P, m);
}

/** Return the potential integral with an array of point dipoles, (p|C|p) */
Matrix SQUIRE_EXPORT potential_integral(const QVector<PointDipole> &C, const PP_GTO &P)
{
    throw SireError::incomplete_code( "Not implemented", CODELOC );
    return Matrix();
}

/** Return the potential integral with a single point dipole, (p|C|p) */
Matrix SQUIRE_EXPORT potential_integral(const PointDipole &C, const PP_GTO &P)
{
    QVector<PointDipole> Cs(1,C);
    return potential_integral(Cs, P);
}

/** Return the potential integral with an array of point dipoles, (p|C|p)^m */
Matrix SQUIRE_EXPORT potential_integral(const QVector<PointDipole> &C, 
                                        const PP_GTO &P, int m)
{
    throw SireError::incomplete_code( "Not implemented", CODELOC );
    return Matrix();
}

/** Return the potential integral with a single point dipole, (p|C|p)^m */
Matrix SQUIRE_EXPORT potential_integral(const PointDipole &C, 
                                        const PP_GTO &P, int m)
{
    QVector<PointDipole> Cs(1,C);
    return potential_integral(Cs, P, m);
}

/////////
///////// The electron repulsion integrals
/////////

static Vector my_electron_integral(const PS_GTO &P, const SS_GTO &Q, 
                                   const double boys[2])
{
    // (ps|ss) = preFac [ Fm(T) (Pi-Ai) + Fm+1(T) (Wi-Pi) ]
    const double prefac = P.scale() * GTOPair::preFac(P, Q);
    const Vector W_minus_P = GTOPair::W(P,Q) - P.P();

    return prefac * ( (P.P_minus_A() * boys[0]) + (W_minus_P * boys[1]) );
}

/** Return the electron repulsion integral (ps|ss) */
Vector SQUIRE_EXPORT electron_integral(const PS_GTO &P, const SS_GTO &Q)
{
    double boys[2];
    multi_boys_2(GTOPair::T(P,Q), boys);

    return my_electron_integral(P, Q, boys);
}

/** Return the electron repulsion integral (ss|ps) */
Vector SQUIRE_EXPORT electron_integral(const SS_GTO &P, const PS_GTO &Q)
{
    return electron_integral(Q, P);
}

/** Return the auxilliary electron repulsion integral (ps|ss)^m */
Vector SQUIRE_EXPORT electron_integral(const PS_GTO &P, const SS_GTO &Q, int m)
{
    double boys[2];
    multi_boys_2(GTOPair::T(P,Q), boys, m);

    return my_electron_integral(P, Q, boys);
}

/** Return the auxilliary electron repulsion integral (ss|ps)^m */
Vector SQUIRE_EXPORT electron_integral(const SS_GTO &P, const PS_GTO &Q, int m)
{
    return electron_integral(Q, P, m);
}

static Matrix my_electron_integral(const PS_GTO &P, const PS_GTO &Q,
                                   const double boys[3])
{
    // (ps|ps) = prefac * { F0(T) (Pi-Ai) (Qk-Ck) + 
    //                      F2(T) (Wi-Pi) (Wk-Qk) +
    //                      F1(T) [(Qk-Ck)(Wi-Pi) + (Pi-Ai)(Wk-Qk) +
    //                                (delta_ik / 2(zeta+eta)) ] }

    const double prefac = P.scale() * Q.scale() * GTOPair::preFac(P, Q);

    const Vector W = GTOPair::W(P, Q);
    const Vector W_minus_P = W - P.P();
    const Vector W_minus_Q = W - Q.Q();
    
    Matrix mat;
    double *m = mat.data();
    
    const double *pa = P.P_minus_A().constData();
    const double *qc = Q.Q_minus_C().constData();
    const double *wp = W_minus_P.constData();
    const double *wq = W_minus_Q.constData();
    
    for (int i=0; i<3; ++i)
    {
        for (int k=0; k<3; ++k)
        {
            const double val = prefac * (
                                 boys[0]*pa[i]*qc[k] + boys[2]*wp[i]*wq[k] + 
                                 boys[1]*(qc[k]*wp[i] + pa[i]*wq[k])
                                         );

            m[ mat.offset(i,k) ] = val;
        }
    }
    
    //do the diagonal i == k elements
    const double extra = 0.5 * boys[1] * prefac / (P.zeta()+Q.eta());
    
    for (int i=0; i<3; ++i)
    {
        m[ mat.offset(i,i) ] += extra;
    }
    
    return mat;
}

/** Return the electron repulsion integral (ps|ps) */
Matrix SQUIRE_EXPORT electron_integral(const PS_GTO &P, const PS_GTO &Q)
{
    double boys[3];
    multi_boys_3( GTOPair::T(P,Q), boys );
    return my_electron_integral(P, Q, boys);
}

/** Return the auxillary electron repulsion integral (ps|ps)^m */
Matrix SQUIRE_EXPORT electron_integral(const PS_GTO &P, const PS_GTO &Q, int m)
{
    double boys[3];
    multi_boys_3( GTOPair::T(P,Q), boys, m );
    return my_electron_integral(P, Q, boys);
}

static Matrix my_electron_integral(const PP_GTO &P, const SS_GTO &Q, 
                                   const double boys[3])
{
    // (pp|ss) = prefac * { F2(T) (Wi-Pi)(Wj-Pj) + 
    //                      F0(T) (Pi-Ai)(Pj-Bj) +
    //                      F1(T) [(Pj-Bj)(Wi-Pi) + (Pi-Ai)(Wj-Pj)] +
    //                      delta_ij/2zeta [F0(T) - rho/zeta F1(T)] }

    const double prefac = P.scale() * GTOPair::preFac(P,Q);
    const double rho = GTOPair::rho(P,Q);
    const Vector W_minus_P = GTOPair::W(P,Q) - P.P();

    Matrix mat;
    double *m = mat.data();
    
    const double *pa = P.P_minus_A().constData();
    const double *pb = P.P_minus_B().constData();
    const double *wp = W_minus_P.constData();
    
    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            const double val = prefac * (
                                   boys[0]*pa[i]*pb[j] + 
                                   boys[2]*wp[i]*wp[j] +
                                   boys[1]*(pb[j]*wp[i] + pa[i]*wp[j]) );
        
            m[ mat.offset(i,j) ] = val;
        }
    }
    
    //now do i == j
    const double one_over_zeta = 1 / P.zeta();
    
    const double extra = 0.5 * prefac * one_over_zeta * (boys[0] - 
                                                boys[1] * rho * one_over_zeta);
            
    for (int i=0; i<3; ++i)
    {
        m[ mat.offset(i,i) ] += extra;
    }
    
    return mat;
}

/** Return the electron repulsion integral (pp|ss) */
Matrix SQUIRE_EXPORT electron_integral(const PP_GTO &P, const SS_GTO &Q)
{
    double boys[3];
    multi_boys_3( GTOPair::T(P,Q), boys );
    return my_electron_integral(P, Q, boys);
}

/** Return the electron repulsion integral (ss|pp) */
Matrix SQUIRE_EXPORT electron_integral(const SS_GTO &P, const PP_GTO &Q)
{
    return electron_integral(Q,P);
}

/** Return the auxilliary electron repulsion integral (pp|ss)^m */
Matrix SQUIRE_EXPORT electron_integral(const PP_GTO &P, const SS_GTO &Q, int m)
{
    double boys[3];
    multi_boys_3( GTOPair::T(P,Q), boys, m );
    return my_electron_integral(P, Q, boys);
}

/** Return the auxilliary electron repulsion integral (ss|pp)^m */
Matrix SQUIRE_EXPORT electron_integral(const SS_GTO &P, const PP_GTO &Q, int m)
{
    return electron_integral(Q,P,m);
}

static void my_electron_integral(const PP_GTO &P, const PS_GTO &Q, const double boys[4],
                                 Array2D<Vector> &mat)
{
    // (pp|ps) = prefactor * { F3(T) (Wi-Pi)(Wj-Pj)(Wk-Qk) +
    // 
    //                         F0(T) (Pi-Ai)(Pj-Bj)(Qk-Ck) +
    //
    //                         F1(T) [(Pj-Bj)(Qk-Ck)(Wi-Pi) + 
    //                                (Pi-Ai)(Qk-Ck)(Wj-Pj) +
    //                                (Pi-Ai)(Pj-Bj)(Wk-Qk)] +
    //
    //                         F2(T) [(Qk-Ck)(Wi-Pi)(Wj-Pj) +
    //                                (Pj-Bj)(Wi-Pi)(Wk-Qk) +
    //                                (Pi-Ai)(Wj-Pj)(Wk-Qk)] +
    //
    //        delta_ij/(2zeta)[ F0(T)(Qk-Ck) + 
    //                          F1(T)[(Wk-Qj) - (rho/zeta)(Qk-Ck)] -
    //                          F2(T)[rho/zeta (Wk-Qk)] ] +
    //
    //        delta_ik/2(zeta+eta)[ F1(T)(Pj-Bj) + F2(T)(Wj-Pj) ] +
    //
    //        delta_jk/2(zeta+eta)[ F1(T)(Pi-Ai) + F2(T)(Wi-Pi) ] }

    const double prefac = P.scale() * Q.scale() * GTOPair::preFac(P,Q);
    const double rho = GTOPair::rho(P,Q);

    const Vector W = GTOPair::W(P,Q);
    const Vector W_minus_P = W - P.P();
    const Vector W_minus_Q = W - Q.Q();
    
    const double *pa = P.P_minus_A().constData();
    const double *pb = P.P_minus_B().constData();
    const double *qc = Q.Q_minus_C().constData();
    const double *wp = W_minus_P.constData();
    const double *wq = W_minus_Q.constData();

    const double inv_zeta = 0.5 * prefac / P.zeta();
    const double inv_zeta_eta = 0.5 / (P.zeta()+Q.eta());
    
    const Vector extra_ij = inv_zeta * 
                               (boys[0]*Q.Q_minus_C() + 
                                boys[1]*(W_minus_Q - (rho/P.zeta())*Q.Q_minus_C()) -
                                boys[2]*(rho/P.zeta())*W_minus_Q);
                                        
    const Vector extra_ik = inv_zeta_eta * (boys[1]*P.P_minus_B() + boys[2]*W_minus_P);
    const Vector extra_jk = inv_zeta_eta * (boys[1]*P.P_minus_A() + boys[2]*W_minus_P);

    const double *ik = extra_ik.constData();
    const double *jk = extra_jk.constData();

    mat.redimension(3,3);
    
    Vector *m = mat.data();
    
    //do the off-diagonal elements
    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            Vector v;
            double *vdata = v.data();
            
            //do the common parts of i,j,k
            for (int k=0; k<3; ++k)
            {
                vdata[k] = boys[3]*wp[i]*wp[j]*wq[k] + 
                           boys[0]*pa[i]*pb[j]*qc[k] +
                           boys[1]*(pb[j]*qc[k]*wp[i] + 
                                    pa[i]*qc[k]*wp[j] +
                                    pa[i]*pb[j]*wq[k]) +
                           boys[2]*(qc[k]*wp[i]*wp[j] +
                                    pb[j]*wp[i]*wq[k] +
                                    pa[i]*wp[j]*wq[k]);
            }
            
            //do i == k
            vdata[i] += ik[j];
            
            //do j == k
            vdata[j] += jk[i];
            
            m[ mat.offset(i,j) ] = prefac * v;
        }
    }

    //add the i==j parts onto the matrix
    for (int i=0; i<3; ++i)
    {
        m[ mat.offset(i,i) ] += extra_ij;
    }
}

/** Calculate the electron repulsion integral (pp|ps) into the passed matrix */
void SQUIRE_EXPORT electron_integral(const PP_GTO &P, const PS_GTO &Q,
                                     Array2D<Vector> &matrix)
{
    double boys[4];
    multi_boys_4( GTOPair::T(P,Q), boys );
    my_electron_integral(P, Q, boys, matrix);
    
    return;
}

/** Calculate the electron repulsion integral (ps|pp) into the passed matrix */
void SQUIRE_EXPORT electron_integral(const PS_GTO &P, const PP_GTO &Q,
                                     Array2D<Vector> &matrix)
{
    electron_integral(Q, P, matrix);
    return;
}

/** Return the electron repulsion integral (pp|ps) */
Array2D<Vector> SQUIRE_EXPORT electron_integral(const PP_GTO &P, const PS_GTO &Q)
{
    Array2D<Vector> matrix(3,3);
    electron_integral(P, Q, matrix);
    return matrix;
}

/** Return the electron repulsion integral (ps,pp) */
Array2D<Vector> SQUIRE_EXPORT electron_integral(const PS_GTO &P, const PP_GTO &Q)
{
    return electron_integral(Q,P);
}

/** Calculate the electron repulsion integral (pp|ps) into the passed matrix */
void SQUIRE_EXPORT electron_integral(const PP_GTO &P, const PS_GTO &Q, int m,
                                     Array2D<Vector> &matrix)
{
    double boys[4];
    multi_boys_4( GTOPair::T(P,Q), boys, m );
    my_electron_integral(P, Q, boys, matrix);
    
    return;
}

/** Calculate the electron repulsion integral (ps|pp) into the passed matrix */
void SQUIRE_EXPORT electron_integral(const PS_GTO &P, const PP_GTO &Q, int m,
                                     Array2D<Vector> &matrix)
{
    electron_integral(Q, P, m, matrix);
    return;
}

/** Return the auxilliary electron repulsion integral (pp|ps)^m */
Array2D<Vector> SQUIRE_EXPORT electron_integral(const PP_GTO &P, 
                                                    const PS_GTO &Q, int m)
{
    Array2D<Vector> matrix(3,3);
    electron_integral(P, Q, m, matrix);
    return matrix;
}

/** Return the auxilliary electron repulsion integral (ps,pp)^m */
Array2D<Vector> SQUIRE_EXPORT electron_integral(const PS_GTO &P, 
                                                const PP_GTO &Q, int m)
{
    return electron_integral(Q,P,m);
}

static void my_electron_integral(const PP_GTO &P, const PP_GTO &Q, const double boys[5],
                                 Array2D<Matrix> &mat)
{
    // (pp|pp) = prefac * {  F4(T)(WPi WPj WQk WQl) + 
    //                       
    //          F3(T)[ (QDl WPi WPj WQk) + (QCk WPi WPj WQl) +
    //                 (PBj WPi WQk WQl) + (PAi WPj WQk WQl) ] +
    //
    //          F0(T)[ PAi PBj QCk QDl ] +
    //
    //          F2(T)[ (QCk QDl WPi WPj) + (PBj QDl WPi WQk) +
    //                 (PAi QDl WPj WQk) + (PBj QCk WPi WQl) +
    //                 (PAi QCk WPj WQl) + (PAi PBj WQk WQl) ] +
    //
    //          F1(T)[ (PBj QCk QDl WPi) + (PAi QCk QDl WPj) + 
    //                 (PAi PBj QDl WQk) + (PAi PBj QCk WQl) ] +
    //
    //  delta_ij/2zeta{ F3(T)[-rho/zeta WQk WQl] + F0(T)[QCk QDl] +
    //                  F2(T)[-rho/zeta[QDl WQk + QCk WQl] + WQk WQl ] +
    //                  F1(T)[-rho/zeta[QCk QDl] + QDl WQk + QCk WQl ] } +
    //
    //  delta_kl/2eta { F3(T)[-rho/eta WPi WPj ] + F0(T)[PAi PBj] + 
    //                  F2(T)[-rho/eta [PBj WPi + PAi WPj] + WPi WPj ] +
    //                  F1(T)[-rho/eta [PAi PBj] + PBj WPi + PAi WPj ] } +
    //
    //  delta_ik/2(e+z){ F3(T)[WPj WQl] + F2(T)[QDl WPj + PBj WQl] + 
    //                   F1(T)[PBj QDl] } +
    //
    //  delta_il/2(e+z){ F3(T)[WPj WQk] + F2(T)[QCk WPj + PBj WQk] +
    //                   F1(T)[PBj QCk] } +
    //
    //  delta_jk/2(e+z){ F3(T)[WPi WQl] + F2(T)[QDl WPi + PAi WQl] + 
    //                   F1(T)[PAi QDl] } +
    //
    //  delta_jl/2(e+z){ F3(T)[WPi WQk] + F2(T)[QCk WPi + PAi WQk] + 
    //                   F1(T)[PAi QCk] } +
    //
    //  delta_il delta_jk/4(e+z)^2{ F2(T) } +
    //
    //  delta_ik delta_jl/4(e+z)^2{ F2(T) } }
    //
    //  delta_ij delta_kl/4ez { F0(T) + (rho^2/ez) F2(T) + 
    //                          F1(T)[-rho/eta - rho/zeta] } +

    const double prefac = P.scale() * Q.scale() * GTOPair::preFac(P,Q);
    
    const double rho = GTOPair::rho(P,Q);
    const Vector W = GTOPair::W(P,Q);
    
    const Vector W_minus_P = W - P.P();
    const Vector W_minus_Q = W - Q.Q();

    const double *pa = P.P_minus_A().constData();
    const double *pb = P.P_minus_B().constData();
    const double *qc = Q.Q_minus_C().constData();
    const double *qd = Q.Q_minus_D().constData();
    const double *wp = W_minus_P.constData();
    const double *wq = W_minus_Q.constData();
    
    mat.redimension(3,3);
    Matrix *m = mat.data();
    
    const double prefac_over_2zeta = 0.5 * prefac / P.zeta();
    const double prefac_over_2eta = 0.5 * prefac / Q.eta();
    const double prefac_over_2e_pl_z = 0.5 * prefac / (P.zeta()+Q.eta());
    const double prefac_over_4e_pl_z2 = 0.25 * prefac / pow_2(P.zeta()+Q.eta());
    const double prefac_over_4ez = 0.25 * prefac / (P.zeta()*Q.eta());
    
    const double rho_over_zeta = rho / P.zeta();
    const double rho_over_eta = rho / Q.eta();
    const double rho_over_ez = rho / (P.zeta()*Q.eta());
    
    //the delta_xy parts are 3x3 matricies
    Matrix delta_ij, delta_kl, delta_ik, delta_il, delta_jk, delta_jl;
    double *m_ij = delta_ij.data();
    double *m_kl = delta_kl.data();
    double *m_ik = delta_ik.data();
    double *m_il = delta_il.data();
    double *m_jk = delta_jk.data();
    double *m_jl = delta_jl.data();
    
    // i == j
    for (int k=0; k<3; ++k)
    {
        for (int l=0; l<3; ++l)
        {
            m_ij[ delta_ij.offset(k,l) ] = prefac_over_2zeta * 
                   ( boys[0]*qc[k]*qd[l] -
                     boys[3]*rho_over_zeta*wq[k]*wq[l] +
                     boys[2]*(wq[k]*wq[l] - rho_over_zeta*(qd[l]*wq[k] + qc[k]*wq[l])) +
                     boys[1]*(qd[l]*wq[k] + qc[k]*wq[l] - rho_over_zeta*qc[k]*qd[l]) );
        }
    }
    
    // k == l
    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            m_kl[ delta_kl.offset(i,j) ] = prefac_over_2eta * (
                    boys[0]*pa[i]*pb[j] -
                    boys[3]*rho_over_eta*wp[i]*wp[j] +
                    boys[2]*(wp[i]*wp[j] - rho_over_eta*(pb[j]*wp[i] + pa[i]*wp[j])) +
                    boys[1]*(pb[j]*wp[i] + pa[i]*wp[j] - rho_over_eta*pa[i]*pb[j]) );
        }
    }
    
    // i == k
    for (int j=0; j<3; ++j)
    {
        for (int l=0; l<3; ++l)
        {
            m_ik[ delta_ik.offset(j,l) ] = prefac_over_2e_pl_z * (
                        boys[1]*pb[j]*qd[l] + 
                        boys[2]*(qd[l]*wp[j] + pb[j]*wq[l]) +
                        boys[3]*wp[j]*wq[l] );
        }
    }
    
    // i == l
    for (int j=0; j<3; ++j)
    {
        for (int k=0; k<3; ++k)
        {
            m_il[ delta_il.offset(j,k) ] = prefac_over_2e_pl_z * (
                        boys[1]*pb[j]*qc[k] +
                        boys[2]*(qc[k]*wp[j] + pb[j]*wq[k]) +
                        boys[3]*wp[j]*wq[k] );
        }
    }
    
    // j == k
    for (int i=0; i<3; ++i)
    {
        for (int l=0; l<3; ++l)
        {
            m_jk[ delta_jk.offset(i,l) ] = prefac_over_2e_pl_z * (
                        boys[1]*pa[i]*qd[l] +
                        boys[2]*(qd[l]*wp[i] + pa[i]*wq[l]) +
                        boys[3]*wp[i]*wq[l] );
        }
    }
    
    // j == l
    for (int i=0; i<3; ++i)
    {
        for (int k=0; k<3; ++k)
        {
            m_jl[ delta_jl.offset(i,k) ] = prefac_over_2e_pl_z * (
                        boys[1]*pa[i]*qc[k] +
                        boys[2]*(qc[k]*wp[i] + pa[i]*wq[k]) +
                        boys[3]*wp[i]*wq[k] );
        }
    }
    
    // i == l and j == k
    for (int i=0; i<3; ++i)
    {
        m_il[ delta_il.offset(i,i) ] += prefac_over_4e_pl_z2 * boys[2];
    }

    // i == k and j == l
    for (int i=0; i<3; ++i)
    {
        m_ik[ delta_ik.offset(i,i) ] += prefac_over_4e_pl_z2 * boys[2];
    }
    
    // i == j and k == l
    for (int i=0; i<3; ++i)
    {
        m_ij[ delta_ij.offset(i,i) ] += prefac_over_4ez * (
                                boys[0] - 
                                boys[1]*(rho_over_eta + rho_over_zeta) + 
                                boys[2]*rho*rho_over_ez );
    }
    
    /////////
    ///////// Now everything is pre-computed, lets do the real work
    /////////
    
    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            Matrix &ij_mat = m[ mat.offset(i,j) ];
            double *ij_m = ij_mat.data();

            for (int k=0; k<3; ++k)
            {
                for (int l=0; l<3; ++l)
                {
                    //first the part common to all...
                    const double val = prefac * (
                                       boys[4]*wp[i]*wp[j]*wq[k]*wq[l] +
                                       boys[3]*(qd[l]*wp[i]*wp[j]*wq[k] + 
                                                qc[k]*wp[i]*wp[j]*wq[l] +
                                                pb[j]*wp[i]*wq[k]*wq[l] +
                                                pa[i]*wp[j]*wq[k]*wq[l]) +
                                       boys[0]*pa[i]*pb[j]*qc[k]*qd[l] +
                                       boys[2]*(qc[k]*qd[l]*wp[i]*wp[j] +
                                                pb[j]*qd[l]*wp[i]*wq[k] +
                                                pa[i]*qd[l]*wp[j]*wq[k] +
                                                pb[j]*qc[k]*wp[i]*wq[l] +
                                                pa[i]*qc[k]*wp[j]*wq[l] +
                                                pa[i]*pb[j]*wq[k]*wq[l]) +
                                       boys[1]*(pb[j]*qc[k]*qd[l]*wp[i] +
                                                pa[i]*qc[k]*qd[l]*wp[j] +
                                                pa[i]*pb[j]*qd[l]*wq[k] +
                                                pa[i]*pb[j]*qc[k]*wq[l]) );
                                                
                    ij_m[ ij_mat.offset(k,l) ] = val;
                }
            }

            
            //now add on the deltas
            
            // k == l
            for (int k=0; k<3; ++k)
            {
                ij_m[ ij_mat.offset(k,k) ] += m_kl[ delta_kl.offset(i,j) ];
            }

            // i == k, j == k
            for (int l=0; l<3; ++l)
            {
                ij_m[ ij_mat.offset(i,l) ] += m_ik[ delta_ik.offset(j,l) ];
                ij_m[ ij_mat.offset(j,l) ] += m_jk[ delta_jk.offset(i,l) ];
            }
            
            // i == l, j == l
            for (int k=0; k<3; ++k)
            {
                ij_m[ ij_mat.offset(k,i) ] += m_il[ delta_il.offset(j,k) ];
                ij_m[ ij_mat.offset(k,j) ] += m_jl[ delta_jl.offset(i,k) ]; 
            }
        }
    }

	//do the i == j diagonal
    for (int i=0; i<3; ++i)
    {        
        // i == j
        m[ mat.offset(i,i) ] += delta_ij;
    }
}

/** Calculate into the passed matrix the electron repulsion integral (pp|pp) */
void SQUIRE_EXPORT electron_integral(const PP_GTO &P, const PP_GTO &Q,
                                     Array2D<Matrix> &matrix)
{
    double boys[5];
    multi_boys_n( GTOPair::T(P,Q), boys, 5 );
    my_electron_integral(P, Q, boys, matrix);
    
    return;
}

/** Return the electron repulsion integral (pp|pp) */
Array2D<Matrix> SQUIRE_EXPORT electron_integral(const PP_GTO &P, const PP_GTO &Q)
{
    Array2D<Matrix> matrix(3,3);
    electron_integral(P, Q, matrix);
    return matrix;
}

/** Calculate into the passed matrix the electron repulsion integral (pp|pp)^m */
void SQUIRE_EXPORT electron_integral(const PP_GTO &P, const PP_GTO &Q, int m,
                                     Array2D<Matrix> &matrix)
{
    double boys[5];
    multi_boys_n( GTOPair::T(P,Q), boys, 5, m );
    my_electron_integral(P, Q, boys, matrix);
    
    return;
}

/** Return the auxilliary electron repulsion integral (pp|pp)^m */
Array2D<Matrix> SQUIRE_EXPORT electron_integral(const PP_GTO &P,
                                                const PP_GTO &Q, int m)
{
    Array2D<Matrix> matrix(3,3);
    electron_integral(P, Q, matrix);
    return matrix;
}

} // end of namespace Squire

///////////
/////////// Implementation of PS_GTOs
///////////

/** Constructor */
PS_GTOs::PS_GTOs()
{}

void PS_GTOs::_pvt_create(const QVector<S_GTO> &s_gtos,
                          const QVector<Vector> &s_centers,
                          const QVector<P_GTO> &p_gtos,
                          const QVector<Vector> &p_centers)
{
    if (s_gtos.count() != s_centers.count())
        throw SireError::invalid_arg( QObject::tr(
                "You must pass in the same number of centers (%1) as "
                "there are S orbitals (%2).")
                    .arg(s_centers.count()).arg(s_gtos.count()), CODELOC );

    if (p_gtos.count() != p_centers.count())
        throw SireError::invalid_arg( QObject::tr(
                "You must pass in the same number of centers (%1) as "
                "there are P orbitals (%2).")
                    .arg(p_centers.count()).arg(p_gtos.count()), CODELOC );
    
    if (s_gtos.isEmpty() or p_gtos.isEmpty())
        return;

    //build all of the orbital pairs
    orbs = Array2D<PS_GTO>(p_gtos.count(), s_gtos.count());
    
    const int np = p_gtos.count();
    const int ns = s_gtos.count();
    
    const P_GTO *p_gtos_array = p_gtos.constData();
    const Vector *p_centers_array = p_centers.constData();
    
    const S_GTO *s_gtos_array = s_gtos.constData();
    const Vector *s_centers_array = s_centers.constData();
    
    PS_GTO *orbs_array = orbs.data();
    
    //the PS array is row-major
    for (int i=0; i<np; ++i)
    {
        const P_GTO &p_gto = p_gtos_array[i];
        const Vector &p_center = p_centers_array[i];
    
        for (int j=0; j<ns; ++j)
        {
            *orbs_array = PS_GTO( p_center, p_gto,
                                  s_centers_array[j], s_gtos_array[j] );
                                  
            ++orbs_array;
        }
    }
}

/** Construct for the set of s orbitals s_gtos located at centers
    s_centers, and p orbitals p_gtos located at centers
    p_centers */
PS_GTOs::PS_GTOs(const QVector<S_GTO> &s_gtos,
                 const QVector<Vector> &s_centers,
                 const QVector<P_GTO> &p_gtos,
                 const QVector<Vector> &p_centers)
{
    this->_pvt_create(s_gtos, s_centers, p_gtos, p_centers);
}

/** Construct for the set of s orbitals s_gtos located at centers
    s_centers, and p orbitals p_gtos located at centers
    p_centers */
PS_GTOs::PS_GTOs(const QVector<P_GTO> &p_gtos,
                 const QVector<Vector> &p_centers,
                 const QVector<S_GTO> &s_gtos,
                 const QVector<Vector> &s_centers)
{
    this->_pvt_create(s_gtos, s_centers, p_gtos, p_centers);
}

/** Copy constructor */
PS_GTOs::PS_GTOs(const PS_GTOs &other)
        : orbs(other.orbs)
{}

/** Destructor */
PS_GTOs::~PS_GTOs()
{}

/** Copy assignment operator */
PS_GTOs& PS_GTOs::operator=(const PS_GTOs &other)
{
    orbs = other.orbs;
    return *this;
}

/** Comparison operator */
bool PS_GTOs::operator==(const PS_GTOs &other) const
{
    return orbs == other.orbs;
}

/** Comparison operator */
bool PS_GTOs::operator!=(const PS_GTOs &other) const
{
    return orbs != other.orbs;
}

/** Return the overlap integrals for each pair of orbitals */
NMatrix PS_GTOs::overlap_integral() const
{
    if (orbs.nRows() == 0 or orbs.nColumns() == 0)
        return NMatrix();

    const int nrows = orbs.nRows();
    const int ncols = orbs.nColumns();

    //3 p-orbitals per p shell
    NMatrix mat = NMatrix::createRowMajor(3 * nrows, ncols);
    
    double *m = mat.data();
    const PS_GTO *orbs_array = orbs.constData();
    
    const int sz = nrows * ncols;
    
    for (int i=0; i<sz; ++i)
    {
        const Vector v = Squire::overlap_integral( orbs_array[i] );

        m[ 3*i     ] = v.x();
        m[ 3*i + 1 ] = v.y();
        m[ 3*i + 2 ] = v.z();
    }
    
    return mat;
}

/** Return the kinetic energy integrals for each pair of orbitals */
NMatrix PS_GTOs::kinetic_integral() const
{
    if (orbs.nRows() == 0 or orbs.nColumns() == 0)
        return NMatrix();

    const int nrows = orbs.nRows();
    const int ncols = orbs.nColumns();

    //3 p-orbitals per p shell
    NMatrix mat = NMatrix::createRowMajor(3 * nrows, ncols);
    
    double *m = mat.data();
    const PS_GTO *orbs_array = orbs.constData();
    
    const int sz = nrows * ncols;
    
    for (int i=0; i<sz; ++i)
    {
        const Vector v = Squire::kinetic_integral( orbs_array[i] );

        m[ 3*i     ] = v.x();
        m[ 3*i + 1 ] = v.y();
        m[ 3*i + 2 ] = v.z();
    }
    
    return mat;
}

/** Return the potential energy integral for each pair of orbitals
    interacting with the point charges in 'C' */
NMatrix PS_GTOs::potential_integral(const QVector<PointCharge> &C) const
{
    if (orbs.nRows() == 0 or orbs.nColumns() == 0 or C.isEmpty())
        return NMatrix();

    const int nrows = orbs.nRows();
    const int ncols = orbs.nColumns();

    //3 p-orbitals per p shell
    NMatrix mat = NMatrix::createRowMajor(3 * nrows, ncols);
    
    double *m = mat.data();
    const PS_GTO *orbs_array = orbs.constData();
    
    const int sz = nrows * ncols;
    
    for (int i=0; i<sz; ++i)
    {
        const Vector v = Squire::potential_integral( C, orbs_array[i] );

        m[ 3*i     ] = v.x();
        m[ 3*i + 1 ] = v.y();
        m[ 3*i + 2 ] = v.z();
    }
    
    return mat;
}

/** Return the mth auxillary potential energy integral for each pair of orbitals
    interacting with the point charges in 'C' */
NMatrix PS_GTOs::potential_integral(const QVector<PointCharge> &C, int maux) const
{
    if (orbs.nRows() == 0 or orbs.nColumns() == 0 or C.isEmpty())
        return NMatrix();

    const int nrows = orbs.nRows();
    const int ncols = orbs.nColumns();

    //3 p-orbitals per p shell
    NMatrix mat = NMatrix::createRowMajor(3 * nrows, ncols);
    
    double *m = mat.data();
    const PS_GTO *orbs_array = orbs.constData();
    
    const int sz = nrows * ncols;
    
    for (int i=0; i<sz; ++i)
    {
        const Vector v = Squire::potential_integral( C, orbs_array[i], maux );

        m[ 3*i     ] = v.x();
        m[ 3*i + 1 ] = v.y();
        m[ 3*i + 2 ] = v.z();
    }
    
    return mat;
}

///////////
/////////// Implementation of PP_GTOs
///////////

/** Constructor */
PP_GTOs::PP_GTOs()
{}

/** Construct for the set of p orbitals p_gtos located at centers
    p_centers */
PP_GTOs::PP_GTOs(const QVector<P_GTO> &p_gtos,
                 const QVector<Vector> &p_centers)
{
    if (p_gtos.count() != p_centers.count())
        throw SireError::invalid_arg( QObject::tr(
                "You must pass in the same number of centers (%1) as "
                "there are P orbitals (%2).")
                    .arg(p_centers.count()).arg(p_gtos.count()), CODELOC );
    
    if (p_gtos.isEmpty())
        return;

    //build all of the orbital pairs
    const int n = p_gtos.count();

    orbs = TrigArray2D<PP_GTO>(n);
    
    const P_GTO *p = p_gtos.constData();
    const Vector *c = p_centers.constData();
    
    PP_GTO *orbs_data = orbs.data();
    
    for (int i=0; i<n; ++i)
    {
    	const P_GTO &pi = p[i];
        const Vector &ci = c[i];
    
    	for (int j=i; j<n; ++j)
        {
        	orbs_data[ orbs.offset(i,j) ] = PP_GTO(ci, pi, c[j], p[j]);
        }
    }
}

/** Copy constructor */
PP_GTOs::PP_GTOs(const PP_GTOs &other)
        : orbs(other.orbs)
{}

/** Destructor */
PP_GTOs::~PP_GTOs()
{}

/** Copy assignment operator */
PP_GTOs& PP_GTOs::operator=(const PP_GTOs &other)
{
    orbs = other.orbs;
    return *this;
}

/** Comparison operator */
bool PP_GTOs::operator==(const PP_GTOs &other) const
{
    return orbs == other.orbs;
}

/** Comparison operator */
bool PP_GTOs::operator!=(const PP_GTOs &other) const
{
    return orbs != other.orbs;
}

/** Return the overlap integrals for each pair of orbitals */
TrigMatrix PP_GTOs::overlap_integral() const
{
	const int n = orbs.count();
    
    if (n <= 0)
    	return TrigMatrix();
    
    //3 p orbitals per p shell
    const int nr = orbs.nRows();
    TrigMatrix mat( 3 * nr );

	const PP_GTO *orbs_data = orbs.constData();
    double *m = mat.data();
    
    for (int i=0; i<nr; ++i)
    {
        for (int j=i; j<nr; ++j)
        {
            const Matrix r = Squire::overlap_integral(orbs_data[i]);
    
            const int ioff = 3*i;
            const int joff = 3*j;
    
            for (int ii=0; ii<3; ++ii)
            {
                for (int jj=0; jj<3; ++jj)
                {
                    m[ orbs.offset( ioff + ii, joff + jj ) ] = r(ii,jj); 
                }
            }
        }
    }
    
    return mat;
}

/** Return the kinetic energy integrals for each pair of orbitals */
TrigMatrix PP_GTOs::kinetic_integral() const
{
	const int n = orbs.count();
    
    if (n <= 0)
    	return TrigMatrix();
    
    //3 p orbitals per p shell
    const int nr = orbs.nRows();
    TrigMatrix mat( 3 * nr );

	const PP_GTO *orbs_data = orbs.constData();
    double *m = mat.data();
    
    for (int i=0; i<nr; ++i)
    {
        for (int j=i; j<nr; ++j)
        {
            const Matrix r = Squire::kinetic_integral(orbs_data[i]);
    
            const int ioff = 3*i;
            const int joff = 3*j;
    
            for (int ii=0; ii<3; ++ii)
            {
                for (int jj=0; jj<3; ++jj)
                {
                    m[ orbs.offset( ioff + ii, joff + jj ) ] = r(ii,jj); 
                }
            }
        }
    }
    
    return mat;
}

/** Return the potential energy integral for each pair of orbitals
    interacting with the point charges in 'C' */
TrigMatrix PP_GTOs::potential_integral(const QVector<PointCharge> &C) const
{
	const int n = orbs.count();
    
    if (n <= 0 or C.isEmpty())
    	return TrigMatrix();
    
    //3 p orbitals per p shell
    const int nr = orbs.nRows();
    TrigMatrix mat( 3 * nr );

	const PP_GTO *orbs_data = orbs.constData();
    double *m = mat.data();
    
    for (int i=0; i<nr; ++i)
    {
        for (int j=i; j<nr; ++j)
        {
            const Matrix r = Squire::potential_integral(C, orbs_data[i]);
    
            const int ioff = 3*i;
            const int joff = 3*j;
    
            for (int ii=0; ii<3; ++ii)
            {
                for (int jj=0; jj<3; ++jj)
                {
                    m[ orbs.offset( ioff + ii, joff + jj ) ] = r(ii,jj); 
                }
            }
        }
    }
    
    return mat;
}

/** Return the mth auxillary potential energy integral for each pair of orbitals
    interacting with the point charges in 'C' */
TrigMatrix PP_GTOs::potential_integral(const QVector<PointCharge> &C, int maux) const
{
	const int n = orbs.count();
    
    if (n <= 0 or C.isEmpty())
    	return TrigMatrix();
    
    //3 p orbitals per p shell
    const int nr = orbs.nRows();
    TrigMatrix mat( 3 * nr );

	const PP_GTO *orbs_data = orbs.constData();
    double *m = mat.data();
    
    for (int i=0; i<nr; ++i)
    {
        for (int j=i; j<nr; ++j)
        {
            const Matrix r = Squire::potential_integral(C, orbs_data[i], maux);
    
            const int ioff = 3*i;
            const int joff = 3*j;
    
            for (int ii=0; ii<3; ++ii)
            {
                for (int jj=0; jj<3; ++jj)
                {
                    m[ orbs.offset( ioff + ii, joff + jj ) ] = r(ii,jj); 
                }
            }
        }
    }
    
    return mat;
}
