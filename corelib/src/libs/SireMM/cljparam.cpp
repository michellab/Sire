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

#include "cljparam.h"
#include "ljparameterdb.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"

using namespace SireMM;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJParam &cljparam)
{
    ds << cljparam.charge().to(mod_electron) << cljparam.lj();
       
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJParam &cljparam)
{
    double q;
    LJParameter lj;
    
    ds >> q >> lj;
    
    cljparam = CLJParam( q * mod_electron, lj );
    
    return ds;
}

/** Null constructor */
CLJParam::CLJParam() : reduced_chg(0), ljid(0)
{}

/** Construct with the passed charge and LJ parameter. If 'auto_lock'
    is false then the LJParameterDB is assumed to be already locked,
    so is not locked */
CLJParam::CLJParam(const Charge &charge, const LJParameter &ljparam, bool auto_lock)
         : reduced_chg(0), ljid(0)
{
    static const double sqrt_one_over_4pieps0 
                            = std::sqrt(SireUnits::one_over_four_pi_eps0);

    reduced_chg = charge * sqrt_one_over_4pieps0;
    
    if (auto_lock)
    {
        ljid = LJParameterDB::addLJParameter(ljparam);
    }
    else
    {
        ljid = LJParameterDB::_locked_addLJParameter(ljparam);
    }
}
         
/** Copy constructor */
CLJParam::CLJParam(const CLJParam &other)   
         : reduced_chg(other.reduced_chg), ljid(other.ljid)
{}

/** Destructor */
CLJParam::~CLJParam()
{}

/** Copy assignment operator */
CLJParam& CLJParam::operator=(const CLJParam &other)
{
    reduced_chg = other.reduced_chg;
    ljid = other.ljid;
    return *this;
}

/** Comparison operator */
bool CLJParam::operator==(const CLJParam &other) const
{
    return ljid == other.ljid and reduced_chg == other.reduced_chg;
}

/** Comparison operator */
bool CLJParam::operator!=(const CLJParam &other) const
{
    return ljid != other.ljid or reduced_chg != other.reduced_chg;
}

const char* CLJParam::typeName()
{
    return "SireMM::CLJParam";
}

/** Return the charge parameter */
Charge CLJParam::charge() const
{
    static const double sqrt_4pieps0 = std::sqrt(SireUnits::four_pi_eps0);
    return Charge( reduced_chg * sqrt_4pieps0 );
}

/** Return the LJ parameter */
LJParameter CLJParam::lj() const
{
    return LJParameterDB::getLJParameter(ljid);
}
