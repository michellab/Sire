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

#include "sgto.h"
#include "pointcharge.h"
#include "pointdipole.h"

#include "SireBase/trigarray2d.h"

#include "SireMaths/maths.h"
#include "SireMaths/boys.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

//////////
////////// Implementation of S_GTO
//////////

static const RegisterMetaType<S_GTO> r_sgto;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const S_GTO &sgto)
{
    writeHeader(ds, r_sgto, 1);

    ds << static_cast<const GTO&>(sgto);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, S_GTO &sgto)
{
    VersionID v = readHeader(ds, r_sgto);

    if (v == 1)
    {
        ds >> static_cast<GTO&>(sgto);
    }
    else
        throw version_error(v, "1", r_sgto, CODELOC);

    return ds;
}

/** Constructor */
S_GTO::S_GTO() : ConcreteProperty<S_GTO,GTO>()
{}

/** Construct with a specified value of alpha and (unnormalised) scale factor */
S_GTO::S_GTO(double alpha, double scale)
      : ConcreteProperty<S_GTO,GTO>(alpha, scale)
{}

/** Copy constructor */
S_GTO::S_GTO(const S_GTO &other) : ConcreteProperty<S_GTO,GTO>(other)
{}

/** Destructor */
S_GTO::~S_GTO()
{}

/** Copy assignment operator */
S_GTO& S_GTO::operator=(const S_GTO &other)
{
    GTO::operator=(other);
    return *this;
}

/** Comparison operator */
bool S_GTO::operator==(const S_GTO &other) const
{
    return GTO::operator==(other);
}

/** Comparison operator */
bool S_GTO::operator!=(const S_GTO &other) const
{
    return GTO::operator!=(other);
}

/** Return a string representation of this orbital */
QString S_GTO::toString() const
{
	if (GTO::scale() == 1)
		return QObject::tr("S (alpha = %1)").arg( GTO::alpha() );
	else
		return QObject::tr("S (alpha = %1, scale = %2)")
        			.arg( GTO::alpha() ).arg(GTO::scale());
}

/** Return the angular momentum of this orbital shell (l==0) */
int S_GTO::angularMomentum() const
{
    return 0;
}

/** Return the number of orbitals in this shell (1) */
int S_GTO::nOrbitals() const
{
    return 1;
}

const char* S_GTO::typeName()
{
    return QMetaType::typeName( qMetaTypeId<S_GTO>() );
}

//////////
////////// Implementation of SS_GTO
//////////

static const RegisterMetaType<SS_GTO> r_ssgto;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SS_GTO &ssgto)
{
    writeHeader(ds, r_ssgto, 1);

    ds << static_cast<const GTOPair&>(ssgto);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SS_GTO &ssgto)
{
    VersionID v = readHeader(ds, r_ssgto);

    if (v == 1)
    {
        ds >> static_cast<GTOPair&>(ssgto);
    }
    else
        throw version_error(v, "1", r_ssgto, CODELOC);

    return ds;
}

/** Constructor */
SS_GTO::SS_GTO() : ConcreteProperty<SS_GTO,GTOPair>()
{}

/** Construct between the passed two S GTOs at the specified points */
SS_GTO::SS_GTO(const Vector &A, const S_GTO &a,
               const Vector &B, const S_GTO &b)
       : ConcreteProperty<SS_GTO,GTOPair>(A, a, B, b)
{
	this->setQ( std::sqrt(electron_integral(*this, *this)) );
}

/** Copy constructor */
SS_GTO::SS_GTO(const SS_GTO &other) : ConcreteProperty<SS_GTO,GTOPair>(other)
{}

/** Destructor */
SS_GTO::~SS_GTO()
{}

/** Copy assignment operator */
SS_GTO& SS_GTO::operator=(const SS_GTO &other)
{
    GTOPair::operator=(other);
    return *this;
}

/** Comparison operator */
bool SS_GTO::operator==(const SS_GTO &other) const
{
    return GTOPair::operator==(other);
}

/** Comparison operator */
bool SS_GTO::operator!=(const SS_GTO &other) const
{
    return GTOPair::operator!=(other);
}

/** Return the angular momentum of the first GTO shell in this pair */
int SS_GTO::angularMomentum0() const
{
    return 0;
}

/** Return the angular momentum of the second GTO shell in this pair */
int SS_GTO::angularMomentum1() const
{
    return 0;
}

/** Return the number of orbitals in the first GTO shell in this pair */
int SS_GTO::nOrbitals0() const
{
    return 1;
}

/** Return the number of orbitals in the second GTO shell in this pair */
int SS_GTO::nOrbitals1() const
{
    return 1;
}

const char* SS_GTO::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SS_GTO>() );
}

////////////
//////////// Integrals involving just SS_GTO shell pairs
////////////
//////////// These are derived in;
//////////// "Efficient recursive computation of molecular integrals
////////////  over Cartesian Gaussian functions"
////////////
//////////// Obara and Saika
//////////// J. Chem. Phys., 84 (7), 3963-3974, 1986
////////////
//////////// I then used Mathematica to expand the recurrance relations
//////////// in Obara and Saika's paper. The Mathematica notebook
//////////// is in techdocs/obara_saika.nb
////////////

namespace Squire
{

/** Return the kinetic energy integral, (s|nabla|s) */
double kinetic_integral(const SS_GTO &P)
{
    // (s|nabla|s) = xi {3 - 2 xi (A-B)^2} (s||s)
    return P.xi() * (3 - 2*P.xi()*P.R2()) * P.ss();
}

/** Return the overlap integral, (s||s) */
double overlap_integral(const SS_GTO &P)
{
    // (s||s) = (pi/eta)^(3/2) exp{-eta(A-B)^2}
    return P.ss();
}

/** Return the potential energy integral with an array of
    point charges, (s|C|s) */
double potential_integral(const QVector<PointCharge> &C,
                                        const SS_GTO &P)
{
    // (s|C|s) = 2 ss Sqrt[zeta/pi] F0(U)
    // U = zeta * (P-C)^2

    const int nc = C.count();
    const PointCharge *c = C.constData();

    double total = 0;

    for (int i=0; i<nc; ++i)
    {
        const double U = P.zeta() * ((P.P() - c[i].center()).length2());
        total += c[i].charge() * boys_f0(U);
    }

    if (total != 0)
    {
        return -2 * std::sqrt( P.zeta() * one_over_pi ) * P.ss() * total;
    }
    else
        return 0;
}

/** Return the potential energy integral with an array of
    point charges, (s|C|s)^m */
double potential_integral(const QVector<PointCharge> &C,
                                        const SS_GTO &P, int m)
{
    const int nc = C.count();
    const PointCharge *c = C.constData();

    double total = 0;

    for (int i=0; i<nc; ++i)
    {
        const double U = P.zeta() * ((P.P() - c[i].center()).length2());
        total += c[i].charge() * boys(m,U);
    }

    if (total != 0)
    {
        return -2 * std::sqrt( P.zeta() * one_over_pi ) * P.ss() * total;
    }
    else
        return 0;
}

/** Return the potential energy integral with a single point charge,
    (s|C|s) */
double potential_integral(const PointCharge &C, const SS_GTO &P)
{
    QVector<PointCharge> Cs(1, C);
    return potential_integral(Cs, P);
}

/** Return the potential energy integral with a single point charge,
    (s|C|s)^m */
double potential_integral(const PointCharge &C, const SS_GTO &P, int m)
{
    QVector<PointCharge> Cs(1, C);
    return potential_integral(Cs, P, m);
}

/** Return the potential energy integral with an array of point dipoles,
    (s|C|s) */
double potential_integral(const QVector<PointDipole> &C, const SS_GTO &P)
{
    throw SireError::incomplete_code("Not implemented", CODELOC);
    return 0;
}

/** Return the potential energy integral with an array of point dipoles,
    (s|C|s)^m */
double potential_integral(const QVector<PointDipole> &C,
                                        const SS_GTO &P, int m)
{
    throw SireError::incomplete_code("Not implemented", CODELOC);
    return 0;
}

/** Return the potential energy integral with a single point dipole, (s|C|s) */
double potential_integral(const PointDipole &C, const SS_GTO &P)
{
    QVector<PointDipole> Cs(1, C);
    return potential_integral(Cs, P);
}

/** Return the potential energy integral with a single point dipole, (s|C|s)^m */
double potential_integral(const PointDipole &C, const SS_GTO &P, int m)
{
    QVector<PointDipole> Cs(1, C);
    return potential_integral(Cs, P, m);
}

/** Return the electron repulsion integral, (ss|ss) */
double electron_integral(const SS_GTO &P, const SS_GTO &Q)
{
    return GTOPair::preFac(P,Q) * boys_f0( GTOPair::T(P,Q) );
}

/** Return the electron repulsion integral, (ss|ss)^m */
double electron_integral(const SS_GTO &P, const SS_GTO &Q, int m)
{
    return GTOPair::preFac(P,Q) * boys_f1( GTOPair::T(P,Q) );
}

} // end of namespace Squire

//////////////
////////////// Implementation of SS_GTOs
//////////////

/** Null constructor */
SS_GTOs::SS_GTOs()
{}

/** Construct for the passed set of S orbitals (on the associated
    centers)

    \throw SireError::incompatible_error
*/
SS_GTOs::SS_GTOs(const QVector<S_GTO> &s_gtos,
        		 const QVector<Vector> &centers)
{
	const int n = s_gtos.count();

    if (n != centers.count())
    {
    	throw SireError::incompatible_error( QObject::tr(
        		"The number of S orbitals (%1) does not equal the number "
                "of centers (%2)!")
                	.arg(n).arg(centers.count()), CODELOC );
    }

	if (n <= 0)	return;

	orbs = TrigArray2D<SS_GTO>(n);

    SS_GTO *orbs_data = orbs.data();

	const S_GTO *s = s_gtos.constData();
    const Vector *c = centers.constData();

    for (int i=0; i<n; ++i)
    {
    	const S_GTO &si = s[i];
        const Vector &ci = c[i];

    	for (int j=i; j<n; ++j)
        {
        	orbs_data[ orbs.offset(i,j) ] = SS_GTO(ci, si, c[j], s[j]);
        }
    }
}

/** Copy constructor */
SS_GTOs::SS_GTOs(const SS_GTOs &other) : orbs(other.orbs)
{}

/** Destructor */
SS_GTOs::~SS_GTOs()
{}

/** Copy assignment operator */
SS_GTOs& SS_GTOs::operator=(const SS_GTOs &other)
{
	orbs = other.orbs;
    return *this;
}

/** Return the overlap integrals between all pairs of S orbitals in this set */
TrigMatrix SS_GTOs::overlap_integral() const
{
	const int n = orbs.count();

    if (n <= 0)
    	return TrigMatrix();

    TrigMatrix mat(orbs.nRows());

	const SS_GTO *orbs_data = orbs.constData();
    double *m = mat.data();

    for (int i=0; i<n; ++i)
    {
    	m[i] = Squire::overlap_integral(orbs_data[i]);
    }

    return mat;
}

/** Return the kinetic energy integral between all pairs of S orbitals in this set */
TrigMatrix SS_GTOs::kinetic_integral() const
{
	const int n = orbs.count();

    if (n <= 0) return TrigMatrix();

	TrigMatrix mat(orbs.nRows());

    const SS_GTO *orbs_data = orbs.constData();
    double *m = mat.data();

    for (int i=0; i<n; ++i)
    {
    	m[i] = Squire::kinetic_integral(orbs_data[i]);
    }

    return mat;
}

/** Return the potential energy integral of all pairs of S orbitals in this
    set with the array of point charges in 'C' */
TrigMatrix SS_GTOs::potential_integral(const QVector<PointCharge> &C) const
{
	const int n = orbs.count();

    if (n <= 0 or C.count() == 0) return TrigMatrix();

	TrigMatrix mat(orbs.nRows());

    const SS_GTO *orbs_data = orbs.constData();
    double *m = mat.data();

    for (int i=0; i<n; ++i)
    {
    	m[i] = Squire::potential_integral(C, orbs_data[i]);
    }

    return mat;
}

/** Return the mth auxilliary potential energy integral of all pairs of S
    orbitals in this set with the array of point charges in 'C' */
TrigMatrix SS_GTOs::potential_integral(const QVector<PointCharge> &C, int maux) const
{
	if (maux == 0)
    	return this->potential_integral(C);

	const int n = orbs.count();

    if (n <= 0 or C.count() == 0) return TrigMatrix();

	TrigMatrix mat(orbs.nRows());

    const SS_GTO *orbs_data = orbs.constData();
    double *m = mat.data();

    for (int i=0; i<n; ++i)
    {
    	m[i] = Squire::potential_integral(C, orbs_data[i], maux);
    }

    return mat;
}

/** Return the coulomb integral between this set of pair of
    S orbitals and the pairs other S orbitals in 'other' */
TrigMatrix SS_GTOs::coulomb_integral(const SS_GTOs &other) const
{
    // J_ij = D_ij Sum_kl < s_i s_j | s_k s_l >

    // (this does not include Dij)

    const int n = orbs.count();
    const int other_n = other.orbs.count();

    if (n <= 0 or other_n <= 0)
        return TrigMatrix();

    TrigMatrix mat(orbs.nRows());

    const SS_GTO *orbs_data = orbs.constData();
    const SS_GTO *other_orbs_data = other.orbs.constData();

    double *m = mat.data();

    for (int i=0; i<n; ++i)
    {
        const SS_GTO &orb = orbs_data[i];

        double j_sum = 0;

        for (int j=0; j<other_n; ++j)
        {
            j_sum += Squire::electron_integral( orb, other_orbs_data[j] );
        }

        m[i] = j_sum;
    }

    return mat;
}

/** Return the exchange integral between this set of pair of
    orbitals and the other set of pairs in 'other' */
TrigMatrix SS_GTOs::exchange_integral(const SS_GTOs &other) const
{
    // K_ij = D_ij Sum_kl < s_i s_k | s_j s_l >

    // (this does not include D_ij)

    const int n = orbs.count();
    const int other_n = other.orbs.count();

    if (n <= 0 or other_n <= 0) return TrigMatrix();

    TrigMatrix mat(orbs.nRows());

    //const SS_GTO *orbs_data = orbs.constData();
    //const SS_GTO *other_orbs_data = other.orbs.constData();

    //double *m = mat.data();

    throw SireError::incomplete_code( QObject::tr("NEED TO FINISH"), CODELOC );

    return TrigMatrix();
}
