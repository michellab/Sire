/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SQUIRE_LATTICECHARGES_H
#define SQUIRE_LATTICECHARGES_H

#include <QVector>

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Element;
}

namespace Squire
{

using SireMaths::Vector;

/** This small internal class is used to hold information
    about a single point lattice charge */
class LatticeCharge
{
public:
    LatticeCharge();
    LatticeCharge(const Vector &v, double charge);
    LatticeCharge(const Vector &v, double charge, const SireMol::Element &element);
    
    LatticeCharge(double x, double y, double z, double charge);
    LatticeCharge(double x, double y, double z, double charge,
                  const SireMol::Element &element);

    LatticeCharge(const LatticeCharge &other);
    
    ~LatticeCharge();

    LatticeCharge& operator=(const LatticeCharge &other);
    
    bool operator==(const LatticeCharge &other) const;
    bool operator!=(const LatticeCharge &other) const;
    
    double x() const
    {
        return d[0];
    }
    
    double y() const
    {
        return d[1];
    }
    
    double z() const
    {
        return d[2];
    }
    
    double charge() const
    {
        return d[3];
    }
    
    SireMol::Element element() const;
    
    quint32 nProtons() const
    {
        return nprotons;
    }
    
private:
    /** The lattice point charge data */
    double d[4];
    
    /** The proton number */
    quint32 nprotons;
};

/** This is a small internal class that is used to pass information about
    lattice charges to the QM programs
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT LatticeCharges
{
public:
    LatticeCharges();
    
    LatticeCharges(const LatticeCharges &other);
    
    ~LatticeCharges();

    void reserve(int n);

    void add(const LatticeCharge &charge);

    bool isEmpty() const;

    int nCharges() const;
    int count() const;
    
    void set(int i, const LatticeCharge &charge);
    
    void setElement(int i, const SireMol::Element &element);
    void setCoordinates(int i, const Vector &coords);
    void setCharge(int i, double charge);
    
    const LatticeCharge& operator[](int i) const;

    const LatticeCharge* constData() const;

private:
    /** All of the lattice charges */
    QVector<LatticeCharge> lattice_charges;
};

}

Q_DECLARE_TYPEINFO( Squire::LatticeCharge, Q_MOVABLE_TYPE );

SIRE_END_HEADER

#endif
