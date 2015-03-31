/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include <cmath>

#include "montecarlo.h"

#include "SireFF/forcefields.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMove;
using namespace SireUnits;
using namespace SireSystem;
using namespace SireFF;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<MonteCarlo> r_mc(MAGIC_ONLY, "SireMove::MonteCarlo");

/** Serialise to a binary data stream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const MonteCarlo &mc)
{
    writeHeader(ds, r_mc, 2);

    SharedDataStream sds(ds);

    sds << mc.ensmble << mc.rangenerator
        << mc.naccept << mc.nreject
        << mc.optimise_moves
        << static_cast<const Move&>(mc);

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, MonteCarlo &mc)
{
    VersionID v = readHeader(ds, r_mc);

    mc.optimise_moves = false;

    if (v == 2)
    {
        SharedDataStream sds(ds);
    
        sds >> mc.ensmble
            >> mc.rangenerator
            >> mc.naccept >> mc.nreject
            >> mc.optimise_moves
            >> static_cast<Move&>(mc);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
    
        sds >> mc.ensmble
            >> mc.rangenerator
            >> mc.naccept >> mc.nreject
            >> static_cast<Move&>(mc);
    }
    else
        throw version_error(v, "1,2", r_mc, CODELOC);

    return ds;
}

/** Construct using the supplied random number generator */
MonteCarlo::MonteCarlo(const PropertyMap &map)
           : Move(map), naccept(0), nreject(0), optimise_moves(false)
{}

/** Copy constructor */
MonteCarlo::MonteCarlo(const MonteCarlo &other)
           : Move(other), ensmble(other.ensmble),
             rangenerator(other.rangenerator),
             naccept(other.naccept), nreject(other.nreject),
             optimise_moves(other.optimise_moves)
{}

/** Destructor */
MonteCarlo::~MonteCarlo()
{}

/** Internal function called by derived classes to set the ensemble
    for this move */
void MonteCarlo::setEnsemble(const Ensemble &ensemble)
{
    ensmble = ensemble;
}

/** Return the ensemble for this move */
Ensemble MonteCarlo::ensemble() const
{
    return ensmble;
}

/** Copy assignment */
MonteCarlo& MonteCarlo::operator=(const MonteCarlo &other)
{
    ensmble = other.ensmble;
    rangenerator = other.rangenerator;
    naccept = other.naccept;
    nreject = other.nreject;
    optimise_moves = other.optimise_moves;

    Move::operator=(other);

    return *this;
}

/** Comparison operator */
bool MonteCarlo::operator==(const MonteCarlo &other) const
{
    return rangenerator == other.rangenerator and
           naccept == other.naccept and
           nreject == other.nreject and
           optimise_moves == other.optimise_moves;
}

/** Comparison operator */
bool MonteCarlo::operator!=(const MonteCarlo &other) const
{
    return not this->operator==(other);
}

/** Set the random number generator to use for these moves */
void MonteCarlo::setGenerator(const RanGenerator &generator)
{
    rangenerator = generator;
}

/** Return the random number generator used for these moves */
const RanGenerator& MonteCarlo::generator() const
{
    return rangenerator;
}

/** Return the number of attempted moves */
quint32 MonteCarlo::nAttempted() const
{
    return naccept + nreject;
}

/** Return the number of accepted moves */
quint32 MonteCarlo::nAccepted() const
{
    return naccept;
}

/** Return the number of rejected moves */
quint32 MonteCarlo::nRejected() const
{
    return nreject;
}

/** Return the acceptance ratio (ratio of
    accepted moves to attempted moves) */
double MonteCarlo::acceptanceRatio() const
{
    int ntotal = this->nAttempted();
    
    if (ntotal == 0)
        return 0;
    else
        return double(this->nAccepted()) / double(ntotal);
}

/** Return the total number of these moves that have been performed */
int MonteCarlo::nMoves() const
{
    return naccept + nreject;
}

/** Zero the move statistics */
void MonteCarlo::clearStatistics()
{
    naccept = 0;
    nreject = 0;
}

/** Turn on use of optimised MC moves. This turns on newer (and potentially more buggy)
    code that aims to speed up the memory allocation and energy calculation for 
    MC moves. */
void MonteCarlo::enableOptimisedMoves()
{
    setUseOptimisedMoves(true);
}

/** Turn off use of optimised MC moves. This uses slightly slower, but likely
    less buggy code, and is worth using if you suspect there are problems with
    the optimised code */
void MonteCarlo::disableOptimisedMoves()
{
    setUseOptimisedMoves(false);
}

/** Switch on or off use of the optimised MC code */
void MonteCarlo::setUseOptimisedMoves(bool on)
{
    optimise_moves = on;
}

/** Return whether or not the optimised MC code is being used */
bool MonteCarlo::usingOptimisedMoves() const
{
    return optimise_moves;
}

/** Perform the NVT Monte Carlo test, using the supplied change in energy
    and the supplied change in biasing probabilities */
bool MonteCarlo::test(double new_energy, double old_energy,
                      double new_bias, double old_bias)
{
    double beta = -1.0 / (k_boltz * ensmble.temperature().value());

    double x = (new_bias / old_bias) * std::exp( beta*(new_energy - old_energy) );

    if (x > 1 or x > rangenerator.rand())
    {
        ++naccept;
        return true;
    }
    else
    {
        ++nreject;
        return false;
    }
}

/** Perform the NVT Monte Carlo test, using the supplied change in energy
    (no change in biasing factor) */
bool MonteCarlo::test(double new_energy, double old_energy)
{
    if (new_energy <= old_energy)
    {
        ++naccept;
        return true;
    }

    double beta = -1.0 / (k_boltz * ensmble.temperature().value());

    double x = std::exp( beta*(new_energy - old_energy) );

    if (x > 1 or x > rangenerator.rand())
    {
        ++naccept;
        return true;
    }
    else
    {
        ++nreject;
        return false;
    }
}

/** Perform the NPT Monte Carlo test, using the supplied change in 
    energy and supplied change in volume (no change in biasing factor) */
bool MonteCarlo::test(double new_energy, double old_energy,
                      int nmolecules,
                      const Volume &new_volume, const Volume &old_volume)
{
    double p_deltav = this->pressure() * (new_volume - old_volume);

    double vratio = nmolecules * std::log(new_volume / old_volume);

    double beta = -1.0 / (k_boltz * ensmble.temperature().value());

    double x = std::exp( beta * (new_energy - old_energy + p_deltav) + vratio );

    if ( x > 1 or x >= rangenerator.rand() )
    {
        ++naccept;
        return true;
    }
    else
    {
        ++nreject;
        return false;
    }
}

/** Perform the NPT Monte Carlo test, using the supplied change in 
    energy and supplied change in volume (with a change in biasing factor) */
bool MonteCarlo::test(double new_energy, double old_energy,
                      int nmolecules,
                      const Volume &new_volume, const Volume &old_volume,
                      double new_bias, double old_bias)
{
    double p_deltav = this->pressure() * (new_volume - old_volume);

    double vratio = nmolecules * std::log(new_volume / old_volume);

    double beta = -1.0 / (k_boltz * ensmble.temperature().value());

    double x =  (new_bias / old_bias) * 
                    std::exp( beta * (new_energy - old_energy + p_deltav) + vratio );

    if ( x > 1 or x >= rangenerator.rand() )
    {
        ++naccept;
        return true;
    }
    else
    {
        ++nreject;
        return false;
    }
}
