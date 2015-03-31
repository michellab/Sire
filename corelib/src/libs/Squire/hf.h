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

#ifndef SQUIRE_HF_H
#define SQUIRE_HF_H

#include <QVector>

#include "SireMaths/vector.h"
#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace Squire
{

using SireMaths::Vector;

class Orbital;
class S_GTO;
class P_GTO;
class PointCharge;
class PointDipole;

/** This is my first attempt at a small HF program */
class SQUIRE_EXPORT HF
{
public:
    HF();
    
    ~HF();
    
    void solve();

    void add(const Orbital &orbital);
    void add(const Vector &point, const Orbital &orbital);
    void add(const Vector &point, const SireUnits::Dimension::Charge &charge);
    void add(const Vector &point, const Vector &dipole);

private:
    /** All of the coordinates of the s_orbital centers */
    QVector<Vector> s_centers;

    /** All of the S-basis functions */
    QVector<S_GTO> s_orbs;
    
    /** All of the coordinates of the p_orbital centers */
    QVector<Vector> p_centers;
    
    /** All of the P-basis functions */
    QVector<P_GTO> p_orbs;
    
    /** All of the point charges (nuclei) */
    QVector<PointCharge> chgs;
    
    /** All of the point dipoles */
    QVector<PointDipole> dipols;
};

}

SIRE_EXPOSE_CLASS( Squire::HF )

SIRE_END_HEADER

#endif
