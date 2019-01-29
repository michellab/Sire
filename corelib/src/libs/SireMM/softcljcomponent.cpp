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

#include "softcljcomponent.h"

#include "SireID/index.h"

#include "SireBase/quickcopy.hpp"

#include "SireStream/datastream.h"

#include "SireError/errors.h"

using namespace SireMM;
using namespace SireMM::detail;
using namespace SireCAS;
using namespace SireID;
using namespace SireBase;

using namespace SireStream;

///////////
/////////// Implementation of SoftCLJEnergy
///////////

/** Null constructor */
SoftCLJEnergy::SoftCLJEnergy()
{
    #ifdef SIRE_USE_SSE
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        nrgs[i] = _mm_set_pd(0,0);
    }
    #else
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        nrgs[i][0] = 0;
        nrgs[i][1] = 0;
    }
    #endif
}

/** Copy constructor */
SoftCLJEnergy::SoftCLJEnergy(const SoftCLJEnergy &other)
{
    quickCopy<double>( (double*)nrgs, (double*)(other.nrgs), MAX_ALPHA_VALUES*2 );
}

/** Destructor */
SoftCLJEnergy::~SoftCLJEnergy()
{}

/** Copy assignment operator */
SoftCLJEnergy& SoftCLJEnergy::operator=(const SoftCLJEnergy &other)
{
    if (this != &other)
    {
        quickCopy<double>( (double*)nrgs, (double*)(other.nrgs), MAX_ALPHA_VALUES*2 );
    }
    
    return *this;
}

/** Self-addition operator */
SoftCLJEnergy& SoftCLJEnergy::operator+=(const SoftCLJEnergy &other)
{
    #ifdef SIRE_USE_SSE
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        nrgs[i] = _mm_add_pd( nrgs[i], other.nrgs[i] );
    }
    #else
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        nrgs[i][0] += other.nrgs[i][0];
        nrgs[i][1] += other.nrgs[i][1];
    }
    #endif
    
    return *this;
}

/** Self-subtraction operator */
SoftCLJEnergy& SoftCLJEnergy::operator-=(const SoftCLJEnergy &other)
{
    #ifdef SIRE_USE_SSE
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        nrgs[i] = _mm_sub_pd( nrgs[i], other.nrgs[i] );
    }
    #else
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        nrgs[i][0] -= other.nrgs[i][0];
        nrgs[i][1] -= other.nrgs[i][1];
    }
    #endif
    
    return *this;
}

/** Addition operator */
SoftCLJEnergy SoftCLJEnergy::operator+(const SoftCLJEnergy &other) const
{
    SoftCLJEnergy ret(*this);
    ret += other;
    return ret;
}

/** Subtraction operator */
SoftCLJEnergy SoftCLJEnergy::operator-(const SoftCLJEnergy &other) const
{
    SoftCLJEnergy ret(*this);
    ret -= other;
    return ret;
}

/** Set the energy of the ith component - set the coulomb component
    to 'cnrg' and the LJ component to 'ljnrg'
    
    \throw SireError::invalid_index
*/
void SoftCLJEnergy::setEnergy(int i, double cnrg, double ljnrg)
{
    i = Index(i).map( MAX_ALPHA_VALUES );
    
    #ifdef SIRE_USE_SSE
    nrgs[i] = _mm_setr_pd(cnrg, ljnrg);
    #else
    nrgs[i][0] = cnrg;
    nrgs[i][1] = ljnrg;
    #endif
}

/** Return the sum of all of the coulomb components */
double SoftCLJEnergy::coulomb() const
{
    #ifdef SIRE_USE_SSE
    double sum = *((const double*)&(nrgs[0]));
    
    for (int i=1; i<MAX_ALPHA_VALUES; ++i)
    {
        sum += *((const double*)&(nrgs[i]));
    }
    
    return sum;
    
    #else
    double sum = nrgs[0][0];
    
    for (int i=1; i<MAX_ALPHA_VALUES; ++i)
    {
        sum += nrgs[i][0];
    }
    
    return sum;
    
    #endif
}

/** Return the sum of all of the LJ components */
double SoftCLJEnergy::lj() const
{
    #ifdef SIRE_USE_SSE
    double sum = *((const double*)&(nrgs[0]) + 1);
    
    for (int i=1; i<MAX_ALPHA_VALUES; ++i)
    {
        sum += *((const double*)&(nrgs[i]) + 1);
    }
    
    return sum;
    
    #else
    double sum = nrgs[0][1];
    
    for (int i=1; i<MAX_ALPHA_VALUES; ++i)
    {
        sum += nrgs[i][1];
    }
    
    return sum;
    
    #endif
}

/** Return the coulomb energy of the ith alpha value

    \throw SireError::invalid_index
*/
double SoftCLJEnergy::coulomb(int i) const
{
    i = SireID::Index(i).map( MAX_ALPHA_VALUES );

    #ifdef SIRE_USE_SSE
    return *((const double*)&(nrgs[i]));
    #else
    return nrgs[i][0];
    #endif
}

/** Return the LJ energy of the ith alpha value 

    \throw SireError::invalid_index
*/
double SoftCLJEnergy::lj(int i) const
{
    i = SireID::Index(i).map( MAX_ALPHA_VALUES );

    #ifdef SIRE_USE_SSE
    return *((const double*)&(nrgs[i]) + 1);
    #else
    return nrgs[i][1];
    #endif
}

/** Return the total energy of all of the alpha values */
double SoftCLJEnergy::total() const
{
    #ifdef SIRE_USE_SSE
    __m128d sum = nrgs[0];
    
    for (int i=1; i<MAX_ALPHA_VALUES; ++i)
    {
        sum = _mm_add_pd(sum, nrgs[i]);
    }
    
    return *((const double*)&sum) + *( ((const double*)&sum) + 1 );
    #else
    double sum[2] = { nrgs[0][0], nrgs[0][1] };
    
    for (int i=1; i<MAX_ALPHA_VALUES; ++i)
    {
        sum[0] += nrgs[i][0];
        sum[1] += nrgs[i][1];
    }
    
    return sum[0] + sum[1];
    
    #endif
}

/** Return the total energy of the ith alpha value

    \throw SireError::invalid_index
*/
double SoftCLJEnergy::total(int i) const
{
    i = SireID::Index(i).map(MAX_ALPHA_VALUES);
    
    #ifdef SIRE_USE_SSE
    return *((const double*)&(nrgs[i])) + *((const double*)&(nrgs[i]) + 1);
    #else
    return nrgs[i][0] + nrgs[i][1];
    #endif
}

/** Return the sum of all of the coulomb components */
double SoftCLJEnergy::component(const CoulombComponent&) const
{
    return this->coulomb();
}

/** Return the sum of all of the LJ components */
double SoftCLJEnergy::component(const LJComponent&) const
{
    return this->lj();
}

/** Return the coulomb energy of the ith alpha value

    \throw SireError::invalid_index
*/
double SoftCLJEnergy::component(const CoulombComponent&, int i) const
{
    return this->coulomb(i);
}

/** Return the LJ energy of the ith alpha value

    \throw SireError::invalid_index
*/
double SoftCLJEnergy::component(const LJComponent&, int i) const
{
    return this->lj(i);
}

///////////
/////////// Implementation of SoftCLJComponent
///////////

static const RegisterMetaType<SoftCLJComponent> r_softcljcomp;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                      const SoftCLJComponent &softcljcomp)
{
    writeHeader(ds, r_softcljcomp, 1);
    
    ds << static_cast<const CLJComponent&>(softcljcomp);
    
    return ds;
}

/** Internal function used to rebuild the components used to represent the 
    coulomb and LJ energies of each individual alpha value */
void SoftCLJComponent::rebuildComponents()
{
    FFName ffname = this->forceFieldName();
    
    alpha_components = QVector<CLJComponent>( MAX_ALPHA_VALUES );
    
    CLJComponent *alpha_components_array = alpha_components.data();
    
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        alpha_components_array[i] = CLJComponent(ffname, QString::number(i));
    }
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SoftCLJComponent &softcljcomp)
{
    VersionID v = readHeader(ds, r_softcljcomp);
    
    if (v == 1)
    {
        ds >> static_cast<CLJComponent&>(softcljcomp);
        softcljcomp.rebuildComponents();
    }
    else
        throw version_error(v, "1", r_softcljcomp, CODELOC);
        
    return ds;
}

/** Construct a null set of SoftCLJComponents */
SoftCLJComponent::SoftCLJComponent() : CLJComponent()
{}

/** Construct the SoftCLJComponents for the passed name */
SoftCLJComponent::SoftCLJComponent(const FFName &name)
                 : CLJComponent(name)
{
    this->rebuildComponents();
}

/** Construct the SoftCLJComponents from the passed symbol */
SoftCLJComponent::SoftCLJComponent(const SireCAS::Symbol &symbol)
                 : CLJComponent(symbol)
{
    this->rebuildComponents();
}

/** Copy constructor */
SoftCLJComponent::SoftCLJComponent(const SoftCLJComponent &other)
                 : CLJComponent(other), alpha_components(other.alpha_components)
{}

/** Destructor */
SoftCLJComponent::~SoftCLJComponent()
{}

/** Copy assignment operator */
SoftCLJComponent& SoftCLJComponent::operator=(const SoftCLJComponent &other)
{
    CLJComponent::operator=(other);
    alpha_components = other.alpha_components;
    
    return *this;
}

/** Return the component representing the total coulomb energy
    of all of the alpha values */
const CoulombComponent& SoftCLJComponent::coulomb() const
{
    return CLJComponent::coulomb();
}

/** Return the component representing the total LJ energy
    of all of the alpha values */
const LJComponent& SoftCLJComponent::lj() const
{
    return CLJComponent::lj();
}

/** Return the component representing the total energy of 
    all of the alpha values */
const CLJComponent& SoftCLJComponent::total() const
{
    return CLJComponent::total();
}

/** Return the component representing the coulomb energy of 
    the ith alpha component
    
    \throw SireError::invalid_index
*/
const CoulombComponent& SoftCLJComponent::coulomb(int i) const
{
    return alpha_components.at( Index(i).map(alpha_components.count()) )
                           .coulomb();
}

/** Return the component representing the LJ energy of 
    the ith alpha component
    
    \throw SireError::invalid_index
*/
const LJComponent& SoftCLJComponent::lj(int i) const
{
    return alpha_components.at( Index(i).map(alpha_components.count()) )
                           .lj();
}

/** Return the component representing the total energy of 
    the ith alpha component
    
    \throw SireError::invalid_index
*/
const CLJComponent& SoftCLJComponent::total(int i) const
{
    return alpha_components.at( Index(i).map(alpha_components.count()) )
                           .total();
}

void SoftCLJComponent::setEnergy(FF&, const CLJEnergy&) const
{
    throw SireError::program_bug( QObject::tr(
            "It is a mistake to call this function!"), CODELOC );
}

void SoftCLJComponent::changeEnergy(FF&, const CLJEnergy&) const
{
    throw SireError::program_bug( QObject::tr(
            "It is a mistake to call this function!"), CODELOC );
}

/** Set the energy in the forcefield 'ff' of all of these components to the 
    values held in 'value' */
void SoftCLJComponent::setEnergy(FF &ff, const SoftCLJEnergy &value) const
{
    double cnrg = value.coulomb();
    double ljnrg = value.lj();
    
    FFComponent::setEnergy(ff, this->total(), cnrg + ljnrg);
    FFComponent::setEnergy(ff, this->coulomb(), cnrg);
    FFComponent::setEnergy(ff, this->lj(), ljnrg);
    
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        const CLJComponent &clj = alpha_components.constData()[i];
    
        cnrg = value.coulomb(i);
        ljnrg = value.lj(i);
                
        FFComponent::setEnergy(ff, clj.total(), cnrg + ljnrg);
        FFComponent::setEnergy(ff, clj.coulomb(), cnrg);
        FFComponent::setEnergy(ff, clj.lj(), ljnrg);
    }
}

/** Change the energy in the forcefield 'ff' of all of these components by the 
    values held in 'delta' */
void SoftCLJComponent::changeEnergy(FF &ff, const SoftCLJEnergy &delta) const
{
    double cnrg = delta.coulomb();
    double ljnrg = delta.lj();
    
    FFComponent::changeEnergy(ff, this->total(), cnrg + ljnrg);
    FFComponent::changeEnergy(ff, this->coulomb(), cnrg);
    FFComponent::changeEnergy(ff, this->lj(), ljnrg);
    
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        const CLJComponent &clj = alpha_components.constData()[i];
    
        cnrg = delta.coulomb(i);
        ljnrg = delta.lj(i);
                
        FFComponent::changeEnergy(ff, clj.total(), cnrg + ljnrg);
        FFComponent::changeEnergy(ff, clj.coulomb(), cnrg);
        FFComponent::changeEnergy(ff, clj.lj(), ljnrg);
    }
}

/** Return all of the symbols associated with these components */
SireCAS::Symbols SoftCLJComponent::symbols() const
{
    Symbols symbls;
    
    symbls.reserve( 3 * (MAX_ALPHA_VALUES+1) );
    
    symbls.insert(this->total());
    symbls.insert(this->coulomb());
    symbls.insert(this->lj());

    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        const CLJComponent &clj = alpha_components.constData()[i];
    
        symbls.insert(clj.total());
        symbls.insert(clj.coulomb());
        symbls.insert(clj.lj());
    }
    
    return symbls;
}

const char* SoftCLJComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SoftCLJComponent>() );
}
