/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIREMM_CLJCALCULATOR_H
#define SIREMM_CLJCALCULATOR_H

#include "cljfunction.h"
#include "cljboxes.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJCalculator;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJCalculator&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJCalculator&);

namespace SireMM
{

/** This class is used to organise the calculation of CLJ energies
    between atoms in CLJ boxes. This is the class that contains
    all of the parallel, Intel Threaded Building Blocks magic ;-)
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJCalculator
{

friend QDataStream& ::operator<<(QDataStream&, const CLJCalculator&);
friend QDataStream& ::operator>>(QDataStream&, CLJCalculator&);

public:
    CLJCalculator(bool reproducible_sum = false);
    CLJCalculator(const CLJCalculator &other);
    ~CLJCalculator();
    
    CLJCalculator& operator=(const CLJCalculator &other);
    
    bool operator==(const CLJCalculator &other) const;
    bool operator!=(const CLJCalculator &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    boost::tuple<double,double> calculate(const CLJFunction &func,
                                          const CLJBoxes &boxes) const;
    
    boost::tuple< QVector<double>, QVector<double> >
            calculate( const QVector<CLJFunctionPtr> &funcs,
                       const CLJBoxes &boxes) const;

    boost::tuple<double,double> calculate(const CLJFunction &func,
                                          const CLJBoxes &boxes0,
                                          const CLJBoxes &boxes1) const;

    boost::tuple<double,double> calculate(const CLJFunction &func,
                                          const CLJAtoms &atoms0,
                                          const CLJBoxes &boxes1) const;

    boost::tuple< QVector<double>, QVector<double> >
            calculate( const QVector<CLJFunctionPtr> &funcs,
                       const CLJBoxes &boxes0, const CLJBoxes &boxes1) const;

    boost::tuple< QVector<double>, QVector<double> >
            calculate( const QVector<CLJFunctionPtr> &funcs,
                       const CLJAtoms &atoms0, const CLJBoxes &boxes1) const;

private:
    /** Whether or not the energy calculation should give the same
        result regardless of the order of summation (i.e. gives the same
        result even if different numbers of processors are used) */
    bool reproducible_sum;
};

}

Q_DECLARE_METATYPE( SireMM::CLJCalculator )

SIRE_EXPOSE_CLASS( SireMM::CLJCalculator )

SIRE_END_HEADER

#endif
