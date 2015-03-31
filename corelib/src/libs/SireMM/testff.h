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

#ifndef SIREMM_TESTFF_H
#define SIREMM_TESTFF_H

#include "cljfunction.h"
#include "cljatoms.h"
#include "cljboxes.h"

#include "SireMol/molecules.h"

#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireMM
{

/** This is a simple forcefield that is designed to let me test new
    ways of calculating energies before rolling out the new design
    to the rest of Sire
    
    @author Christopher Woods
*/
class SIREMM_EXPORT TestFF
{
public:
    TestFF();
    TestFF(const TestFF &other);
    ~TestFF();
    
    TestFF& operator=(const TestFF &other);
    
    void add(const Molecules &molecules);
    void addFixedAtoms(const Molecules &molecules);
    
    void calculateEnergy();
    
    void setCutoff(Length coul_cutoff, Length lj_cutoff);
    
private:
    CLJAtoms atoms0;
    CLJAtoms atoms1;

    CLJBoxes cljboxes0;
    CLJBoxes cljboxes1;

    boost::shared_ptr<CLJFunction> cljfunc;
};

}

SIRE_EXPOSE_CLASS( SireMM::TestFF )

SIRE_END_HEADER

#endif
