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

#ifndef SIREMM_CLJGRID_H
#define SIREMM_CLJGRID_H

#include "cljfunction.h"
#include "cljatoms.h"
#include "cljboxes.h"

#include "SireVol/gridinfo.h"
#include "SireVol/aabox.h"

#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJGrid;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJGrid&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJGrid&);

namespace SireMM
{

using SireVol::AABox;
using SireVol::GridInfo;

/** This class holds a 3D grid of the coulomb potential
    at points in space created by a set of atoms, and calculates
    the coulomb and LJ energies of atoms with that grid
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJGrid
{

friend QDataStream& ::operator<<(QDataStream&, const CLJGrid&);
friend QDataStream& ::operator>>(QDataStream&, CLJGrid&);

public:
    CLJGrid();
    CLJGrid(const AABox &grid_dimensions);
    CLJGrid(const AABox &grid_dimensions, Length grid_spacing);
    CLJGrid(const GridInfo &grid);
    CLJGrid(const CLJFunction &cljfunc);
    CLJGrid(const CLJFunction &cljfunc, const AABox &grid_dimensions);
    CLJGrid(const CLJFunction &cljfunc, const AABox &grid_dimensions, Length grid_spacing);
    CLJGrid(const CLJFunction &cljfunc, const GridInfo &grid);
    
    CLJGrid(const CLJGrid &other);
    
    ~CLJGrid();
    
    CLJGrid& operator=(const CLJGrid &other);
    
    bool operator==(const CLJGrid &other) const;
    bool operator!=(const CLJGrid &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    CLJGrid* clone() const;
    
    QString toString() const;
    
    void addFixedAtoms(const CLJAtoms &atoms);
    void setFixedAtoms(const CLJAtoms &atoms);
    void setFixedAtoms(const CLJBoxes &atoms);
    
    CLJAtoms fixedAtoms() const;
    
    int nFixedAtoms() const;
    
    bool isEmpty() const;
    
    void setCLJFunction(const CLJFunction &function);
    
    const CLJFunction& cljFunction() const;
    
    void setGridSpacing(Length grid_spacing);
    Length gridSpacing() const;
    
    void setGridBuffer(Length grid_buffer);
    Length gridBuffer() const;
    
    void setGridDimensions(const AABox &grid_dimensions);

    void setGridDimensions(const CLJAtoms &atoms);
    void setGridDimensions(const CLJAtoms &atoms, Length grid_spacing);
    void setGridDimensions(const CLJAtoms &atoms, Length grid_spacing, Length buffer);

    AABox gridDimensions() const;
    
    void setGrid(const GridInfo &grid);

    GridInfo grid() const;
    
    void setUseGrid(bool on);
    
    void enableGrid();
    void disableGrid();
    
    bool usesGrid() const;
    bool functionSupportsGrid() const;
    
    void enableParallelCalculation();
    void disableParallelCalculation();
    void setUseParallelCalculation(bool on);
    bool usesParallelCalculation() const;
    
    void enableReproducibleCalculation();
    void disableReproducibleCalculation();
    void setUseReproducibleCalculation(bool on);
    bool usesReproducibleCalculation() const;
    
    void total(const CLJAtoms &atoms, double &cnrg, double &ljnrg) const;
    void total(const CLJBoxes &atoms, double &cnrg, double &ljnrg) const;
    
    boost::tuple<double,double> calculate(const CLJAtoms &atoms) const;
    boost::tuple<double,double> calculate(const CLJBoxes &atoms) const;
    
    double coulomb(const CLJAtoms &atoms) const;
    double coulomb(const CLJBoxes &atoms) const;
    
    double lj(const CLJAtoms &atoms) const;
    double lj(const CLJBoxes &atoms) const;
    
    static qint32 idOfFixedAtom();
    
private:
    void clearGrid();
    void calculateGrid();
    void checkIfGridSupported();

    /** Description of the grid */
    GridInfo grid_info;
    
    /** The buffer used when building the grid */
    float grid_buffer;
    
    /** The actual grid of potentials */
    QVector<float> grid_pots;
    
    /** The CLJFunction used to calculate the potential */
    CLJFunctionPtr cljfunc;
    
    /** The atoms whose potentials are contained in the grid */
    CLJBoxes cljboxes;
    
    /** The atoms from cljboxes that are within the LJ cutoff of any
        of the grid points */
    CLJBoxes close_atoms;
    
    /** Whether or not to use a grid */
    bool use_grid;
    
    /** Whether or not the CLJFunction supports use of a grid */
    bool cljfunc_supports_grid;
    
    /** Whether or not to use a parallel energy calculation */
    bool parallel_calc;
    
    /** Whether or not to use an exact summing algorithm to 
        ensure that parallel calculations always give the same energy */
    bool repro_sum;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return whether or not this grid is empty (has no fixed atoms) */
inline bool CLJGrid::isEmpty() const
{
    return cljboxes.isEmpty();
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::CLJGrid )

SIRE_EXPOSE_CLASS( SireMM::CLJGrid )

SIRE_END_HEADER

#endif
