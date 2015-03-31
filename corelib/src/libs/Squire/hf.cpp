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

#include "hf.h"
#include "sgto.h"
#include "pgto.h"
#include "pointcharge.h"
#include "pointdipole.h"

#include "SireMaths/boys.h"
#include "SireMaths/gamma.h"
#include "SireMaths/vector.h"
#include "SireMaths/maths.h"

#include "SireMaths/nmatrix.h"
#include "SireMaths/nvector.h"

#include "SireBase/array2d.hpp"

#include "SireError/errors.h"

#include <QDebug>

using namespace Squire;
using namespace SireMaths;
using namespace SireBase;            

/** Constructor */
HF::HF()
{}

/** Destructor */
HF::~HF()
{}

/** Add an orbital that does not to be located at a particular
    point in space */
void HF::add(const Orbital &orbital)
{
}

/** Add an orbital that needs to be placed at 'center' */
void HF::add(const Vector &center, const Orbital &orbital)
{
    if (orbital.isA<S_GTO>())
    {
        s_centers.append(center);
        s_orbs.append(orbital.asA<S_GTO>());
    }
    else if (orbital.isA<P_GTO>())
    {
        p_centers.append(center);
        p_orbs.append(orbital.asA<P_GTO>());
    }
    else
        throw SireError::unsupported( QObject::tr(
                "The HF program does not support orbitals of type %1 (%2)")
                    .arg(orbital.toString()).arg(orbital.what()), CODELOC );
}

void HF::add(const Vector &center, const SireUnits::Dimension::Charge &charge)
{
    chgs.append( PointCharge(center,charge) );
}

void HF::add(const Vector &center, const Vector &dipole)
{
    dipols.append( PointDipole(center, dipole) );
}

static NMatrix make_overlap_matrix(const Array2D<SS_GTO> &orbitals)
{
    const int norbitals = orbitals.nRows();
    
    NMatrix overlap_matrix(norbitals, norbitals);

    //calculate the overlap matrix - this is the overlap integral
    //of all pairs of orbitals
    for (int i=0; i<norbitals; ++i)
    {
        for (int j=0; j<norbitals; ++j)
        {
            overlap_matrix(i,j) = overlap_integral( orbitals(i,j) );
        }
    }

    return overlap_matrix;
}

static NMatrix make_core_fock_matrix(const Array2D<SS_GTO> &orbitals,
                                     const QVector<PointCharge> &charges)
{
    const int norbitals = orbitals.nRows();
    const int ncharges = charges.count();

    NMatrix fock_matrix(norbitals, norbitals);
    NMatrix kinetic_matrix(norbitals, norbitals);
    NMatrix nuclear_matrix(norbitals, norbitals, 0);
    
    //kinetic energy
    for (int i=0; i<norbitals; ++i)
    {
        for (int j=0; j<norbitals; ++j)
        {
            kinetic_matrix(i,j) = kinetic_integral( orbitals(i,j) ); 
        }
    }
    
    qDebug() << "KINETIC MATRIX\n" << kinetic_matrix.toString();
    
    //nuclear-electron energy
    for (int k=0; k<ncharges; ++k)
    {
        NMatrix my_nuclear_matrix(norbitals, norbitals);
    
        for (int i=0; i<norbitals; ++i)
        {
            for (int j=0; j<norbitals; ++j)
            {
                my_nuclear_matrix(i,j) = potential_integral(charges[k], orbitals(i,j));
            }
        }

        nuclear_matrix += my_nuclear_matrix;
        
        qDebug() << "NUCLEAR MATRIX" << k << "\n" << my_nuclear_matrix.toString();
    }

    qDebug() << "NUCLEAR MATRIX\n" << nuclear_matrix.toString();
    
    for (int i=0; i<norbitals; ++i)
    {
        for (int j=0; j<norbitals; ++j)
        {
            fock_matrix(i,j) = kinetic_matrix(i,j) + nuclear_matrix(i,j);
        }
    }
    
    for (int i=0; i<norbitals; ++i)
    {
        for (int j=0; j<norbitals; ++j)
        {
            for (int k=0; k<norbitals; ++k)
            {
                for (int l=0; l<norbitals; ++l)
                {
                    qDebug() << i+1 << j+1 << k+1 << l+1 
                             << electron_integral(orbitals(i,j), orbitals(k,l));
                }
            }
        }
    }
    
    return fock_matrix;
}

NMatrix make_G(const NMatrix &P, const Array2D<SS_GTO> &orbitals)
{
    const int norbitals = orbitals.nRows();
    
    NMatrix G(norbitals, norbitals);
    
    for (int i=0; i<norbitals; ++i)
    {
        for (int j=0; j<norbitals; ++j)
        {
            double v = 0;
            
            for (int k=0; k<norbitals; ++k)
            {
                for (int l=0; l<norbitals; ++l)
                {
                    v += P(k,l)*
                            ( electron_integral(orbitals(i,j), orbitals(k,l))
                              - 0.5*electron_integral(orbitals(i,l),orbitals(k,j)) );
                }
            }
            
            G(i,j) = v;
        }
    }
    
    return G;
}

double calc_e(const NMatrix &P, const NMatrix &H, const NMatrix &F)
{
    double e_elec = 0;
    
    const int norbitals = P.nRows();
    
    for (int i=0; i<norbitals; ++i)
    {
        for (int j=0; j<norbitals; ++j)
        {
            e_elec += 0.5*P(i,j)*(H(i,j) + F(i,j));
        }
    }
    
    return e_elec;
}

double calc_delta(const NMatrix &P, const NMatrix &NEW_P)
{
    double delta(0);
    
    const int norbitals = P.nRows();
    
    for (int i=0; i<norbitals; ++i)
    {
        for (int j=0; j<norbitals; ++j)
        {
            delta += pow_2( P(i,j) - NEW_P(i,j) );
        }
    }
    
    return std::sqrt(0.25 * delta);
}

double nuclear_energy(const QVector<PointCharge> &charges)
{
    double nrg = 0;

    for (int i=0; i<charges.count()-1; ++i)
    {
        for (int j=i+1; j<charges.count(); ++j)
        {
            double r = Vector::distance( charges[i].center(), charges[j].center() );
            
            nrg += charges[i].charge() * charges[j].charge() / r;
        }
    }
    
    return nrg;
}

void HF::solve()
{
    const int norbitals = s_orbs.count();

    if (norbitals == 0)
        return;

    //make matrix of all orbital shell pairs
    Array2D<SS_GTO> ss_orbitals(norbitals, norbitals);
    
    for (int i=0; i<norbitals; ++i)
    {
        for (int j=0; j<norbitals; ++j)
        {
            ss_orbitals(i,j) = SS_GTO( s_centers.at(i), s_orbs.at(i),
                                       s_centers.at(j), s_orbs.at(j) );
        }
    }

    //need to create the overlap matrix
    NMatrix S = make_overlap_matrix(ss_orbitals);
    
    //and the Fock matrix
    NMatrix H = make_core_fock_matrix(ss_orbitals, chgs);
    
    qDebug() << "S\n" << S.toString();
    qDebug() << "H\n" << H.toString();

    //make the orthonormalization matrix (canonical orthonormalization)
    std::pair<NVector,NMatrix> eig = S.diagonalise();

    NVector eigval = eig.first;
    NMatrix X = eig.second;
    
    for (int i=0; i<eigval.count(); ++i)
    {
        for (int j=0; j<eigval.count(); ++j)
        {
            X(i,j) /= std::sqrt(eigval[j]);
        }
    }
    
    qDebug() << "X\n" << X.toString();
    
    NMatrix P(norbitals, norbitals, 0);

    //build the new fock matrix from the current density matrix
    NMatrix G = make_G(P, ss_orbitals);
    NMatrix F = H + G;

    double e_elec = calc_e(P, H, F);

    qDebug() << "F\n" << F.toString();
    qDebug() << "G\n" << G.toString();
    qDebug() << "electronic energy" << e_elec;
       
    double delta;
       
    do
    {
        //transform the fock matrix -  F' = X F X^T
        NMatrix FPRIME = X.transpose() * F * X;
        
        qDebug() << "FPRIME\n" << FPRIME.toString();
        
        //diagonalise the transformed fock matrix
        std::pair<NVector,NMatrix> orbeig = FPRIME.diagonalise();

        NVector E = orbeig.first;
        NMatrix CPRIME = orbeig.second;
                    
        qDebug() << "E" << E.toString();

        qDebug() << "CPRIME\n" << CPRIME.toString();
        
        // C' == xfock
        NMatrix C = X * CPRIME;
        
        qDebug() << "C\n" << C.toString();
        
        //compute new density matrix
        NMatrix NEW_P(norbitals, norbitals, 0);

        for (int k=0; k<norbitals/2; ++k)
        {
            for (int i=0; i<norbitals; ++i)
            {
                for (int j=0; j<norbitals; ++j)
                {
                    NEW_P(i,j) = 2*C(i,k)*C(j,k);
                }
            }
        }

        //calculate the change in density matrix
        delta = calc_delta(P, NEW_P);
    
        P = NEW_P;

        //build the new fock matrix from the current density matrix
        G = make_G(P, ss_orbitals);
        F = H + G;

        qDebug() << "G\n" << G.toString();
        qDebug() << "F\n" << F.toString();

        e_elec = calc_e(P, H, F);

        qDebug() << "Convergence ==" << delta;

        qDebug() << "electronic energy" << e_elec;
        qDebug() << "P" << P.toString();
    }
    while (delta > 1e-4);
    
    double e_nuc = nuclear_energy(chgs);
    
    qDebug() << "Potential energy =" << e_nuc << "Hartrees";
    
    qDebug() << "Total energy =" << e_nuc+e_elec << "Hartrees";
    
    qDebug() << "The end!";
}
