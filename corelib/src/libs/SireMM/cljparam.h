/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREMM_CLJPARAM_H
#define SIREMM_CLJPARAM_H

#include "SireUnits/dimensions.h"

#include "SireFF/atomicffparameters.hpp"

#include "ljparameter.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJParam;
}

QDataStream& operator<<(QDataStream&, const SireMM::CLJParam&);
QDataStream& operator>>(QDataStream&, SireMM::CLJParam&);

namespace SireMM
{

/** This is a simple internal class that holds the combined CLJ parameter
    as a reduced charge (charge / sqrt(4 pi epsilon_0)) and the LJ
    parameter as an ID number for the LJParameterDB 
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJParam
{
public:
    CLJParam();
    
    CLJParam(const SireUnits::Dimension::Charge &charge,
             const LJParameter &ljparam, bool auto_lock=true);
             
    CLJParam(const CLJParam &other);
    
    ~CLJParam();
    
    CLJParam& operator=(const CLJParam &other);
    
    bool operator==(const CLJParam &other) const;
    bool operator!=(const CLJParam &other) const;
    
    static const char* typeName();
    
    double q() const;
    SireUnits::Dimension::Charge charge() const;

    quint32 ljID() const;
    LJParameter lj() const;

private:
    /** The reduced charge (q / sqrt(4 pi epsilon_0)) */
    double reduced_chg;
    
    /** The ID of the LJ parameter in the global LJParameterDB */
    quint32 ljid;
};

/** Return the reduced charge parameter (charge / sqrt(4 pi epsilon_0)) */
inline double CLJParam::q() const
{
    return reduced_chg;
}

/** Return the LJ parameter */
inline quint32 CLJParam::ljID() const
{
    return ljid;
}

typedef SireFF::AtomicFFParameters<CLJParam> CLJParams;
typedef SireFF::AtomicFFParametersArray<CLJParam> CLJParamsArray;

}

Q_DECLARE_TYPEINFO( SireMM::CLJParam, Q_MOVABLE_TYPE );

Q_DECLARE_METATYPE( SireMM::CLJParams )
Q_DECLARE_METATYPE( SireMM::CLJParamsArray )

SIRE_END_HEADER

#endif
