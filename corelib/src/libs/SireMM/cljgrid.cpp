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

#include "cljgrid.h"
#include "cljcalculator.h"
#include "cljshiftfunction.h"

#include "SireVol/cartesian.h"

#include "SireMaths/multidouble.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QElapsedTimer>
#include <QMutex>

using namespace SireMM;
using namespace SireVol;
using namespace SireMaths;
using namespace SireUnits;
using namespace SireStream;

static const RegisterMetaType<CLJGrid> r_grid(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const CLJGrid &grid)
{
    writeHeader(ds, r_grid, 2);
    
    SharedDataStream sds(ds);
    
    sds << grid.grid_info << grid.grid_buffer
        << grid.cljfunc << grid.grid_pots << grid.cljboxes
        << grid.close_atoms << grid.use_grid
        << grid.parallel_calc << grid.repro_sum;
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJGrid &grid)
{
    VersionID v = readHeader(ds, r_grid);
    
    if (v <= 2)
    {
        SharedDataStream sds(ds);
        
        sds >> grid.grid_info >> grid.grid_buffer
            >> grid.cljfunc >> grid.grid_pots >> grid.cljboxes
            >> grid.close_atoms >> grid.use_grid;
        
        grid.checkIfGridSupported();
        
        if (grid.grid_info.nPoints() != grid.grid_pots.count())
        {
            qDebug() << "Discrepancy in grid dimensions!";
            qDebug() << grid.grid_info.nPoints() << grid.grid_pots.count();
            
            GridInfo old = grid.grid_info;
            grid.grid_info = GridInfo();
            grid.clearGrid();
            grid.setGrid(old);
        }
        
        if (v == 2)
        {
            sds >> grid.parallel_calc >> grid.repro_sum;
        }
        else
        {
            grid.parallel_calc = true;
            grid.repro_sum = false;
        }
    }
    else
        throw version_error(v, "1,2", r_grid, CODELOC);
    
    return ds;
}

static SireBase::PropPtr<CLJFunction> global_func( new CLJShiftFunction() );

static const Length global_grid_buffer = 2 * angstrom;
static const Length global_grid_spacing = 1 * angstrom;

/** Internal function used to check and then cache whether or
    not the CLJFunction supports a grid */
void CLJGrid::checkIfGridSupported()
{
    cljfunc_supports_grid = cljfunc.read().supportsGridCalculation() and
                            //(not cljfunc.read().isPeriodic()) and
                            cljfunc.read().hasCutoff() and
       cljfunc.read().coulombCutoff().value() > (cljfunc.read().ljCutoff().value() + 5);
}

/** Constructor */
CLJGrid::CLJGrid() : grid_buffer(global_grid_buffer),
                     cljfunc(global_func), use_grid(true),
                     parallel_calc(true), repro_sum(false)
{
    checkIfGridSupported();
}

/** Construct, specifying the dimensions of the grid */
CLJGrid::CLJGrid(const AABox &grid_dimensions)
        : grid_buffer(global_grid_buffer), cljfunc(global_func), use_grid(true),
          parallel_calc(true), repro_sum(false)
{
    setGrid( GridInfo(grid_dimensions, global_grid_spacing) );
    checkIfGridSupported();
}

/** Construct, specifying the dimensions and spacing for the grid */
CLJGrid::CLJGrid(const AABox &grid_dimensions, Length spacing)
        : grid_buffer(global_grid_buffer), cljfunc(global_func), use_grid(true),
          parallel_calc(true), repro_sum(false)
{
    setGrid( GridInfo(grid_dimensions,spacing) );
    checkIfGridSupported();
}

/** Construct, specifying the grid */
CLJGrid::CLJGrid(const GridInfo &grid)
        : grid_buffer(global_grid_buffer), cljfunc(global_func), use_grid(true),
          parallel_calc(true), repro_sum(false)
{
    setGrid(grid);
    checkIfGridSupported();
}

/** Construct, specifying the function to use to calculate the energy */
CLJGrid::CLJGrid(const CLJFunction &func)
        : grid_buffer(global_grid_buffer), cljfunc(func), use_grid(true),
          parallel_calc(true), repro_sum(false)
{
    checkIfGridSupported();
}

/** Construct, specifying the function to use to calculate the energy and 
    the grid dimensions */
CLJGrid::CLJGrid(const CLJFunction &func, const AABox &grid_dimensions)
        : grid_buffer(global_grid_buffer), cljfunc(func), use_grid(true),
          parallel_calc(true), repro_sum(false)
{
    setGrid( GridInfo(grid_dimensions,global_grid_spacing) );
    checkIfGridSupported();
}

/** Construct, specifying the function to use to calculate the energy and
    the grid dimensions and grid spacing */
CLJGrid::CLJGrid(const CLJFunction &func, const AABox &grid_dimensions, Length spacing)
        : grid_buffer(global_grid_buffer), cljfunc(func), use_grid(true),
          parallel_calc(true), repro_sum(false)
{
    setGrid( GridInfo(grid_dimensions,spacing) );
    checkIfGridSupported();
}

/** Construct, specifying the grid and the energy function */
CLJGrid::CLJGrid(const CLJFunction &func, const GridInfo &grid)
        : grid_buffer(global_grid_buffer), cljfunc(func), use_grid(true),
          parallel_calc(true), repro_sum(false)
{
    setGrid(grid);
    checkIfGridSupported();
}

/** Copy constructor */
CLJGrid::CLJGrid(const CLJGrid &other)
        : grid_info(other.grid_info),
          grid_buffer(other.grid_buffer), grid_pots(other.grid_pots),
          cljfunc(other.cljfunc), cljboxes(other.cljboxes), close_atoms(other.close_atoms),
          use_grid(other.use_grid), cljfunc_supports_grid(other.cljfunc_supports_grid),
          parallel_calc(other.parallel_calc), repro_sum(other.repro_sum)
{}

/** Destructor */
CLJGrid::~CLJGrid()
{}

/** Copy assignment operator */
CLJGrid& CLJGrid::operator=(const CLJGrid &other)
{
    if (this != &other)
    {
        grid_info = other.grid_info;
        grid_buffer = other.grid_buffer;
        grid_pots = other.grid_pots;
        cljfunc = other.cljfunc;
        cljboxes = other.cljboxes;
        close_atoms = other.close_atoms;
        use_grid = other.use_grid;
        cljfunc_supports_grid = other.cljfunc_supports_grid;
        parallel_calc = other.parallel_calc;
        repro_sum = other.repro_sum;
    }
    
    return *this;
}

/** Comparison operator */
bool CLJGrid::operator==(const CLJGrid &other) const
{
    return grid_info == other.grid_info and
           grid_buffer == other.grid_buffer and cljfunc == other.cljfunc and
           cljboxes == other.cljboxes and
           parallel_calc == other.parallel_calc and
           repro_sum == other.repro_sum;
}

/** Comparison operator */
bool CLJGrid::operator!=(const CLJGrid &other) const
{
    return not operator==(other);
}

const char* CLJGrid::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJGrid>() );
}

const char* CLJGrid::what() const
{
    return CLJGrid::typeName();
}

CLJGrid* CLJGrid::clone() const
{
    return new CLJGrid(*this);
}

/** Return the ID number of a fixed atom. All fixed atoms are given this ID, so that
    you can mask out interactions with them. This is a negative number and unlikely
    to be used by any other part of the code */
qint32 CLJGrid::idOfFixedAtom()
{
    return -439284;
}

/** Return the number of fixed atoms */
int CLJGrid::nFixedAtoms() const
{
    return cljboxes.nAtoms();
}

QString CLJGrid::toString() const
{
    return QObject::tr("CLJGrid( grid() == %1, "
                       "cljFunction() == %2, nFixedAtoms() == %3 )")
                            .arg(grid().toString())
                            .arg(cljFunction().toString())
                            .arg(nFixedAtoms());
}

void CLJGrid::clearGrid()
{
    grid_pots.clear();
    close_atoms = CLJBoxes();
}

/** Add the passed atoms onto the set of fixed atoms */
void CLJGrid::addFixedAtoms(const CLJAtoms &atoms)
{
    CLJAtoms atms(atoms);
    atms.setAllID( idOfFixedAtom() );
    cljboxes = CLJBoxes( (cljboxes.atoms() + atms).squeeze() );
    clearGrid();
}

/** Set the fixed atoms equal to the passed atoms */
void CLJGrid::setFixedAtoms(const CLJAtoms &atoms)
{
    CLJAtoms atms(atoms);
    atms.setAllID( idOfFixedAtom() );
    cljboxes = CLJBoxes(atms.squeeze());
    clearGrid();
}

/** Set the fixed atoms equal to the passed atoms */
void CLJGrid::setFixedAtoms(const CLJBoxes &atoms)
{
    CLJAtoms atms( atoms.atoms() );
    atms.setAllID( idOfFixedAtom() );
    cljboxes = CLJBoxes(atms.squeeze());
    clearGrid();
}

/** Return all of the fixed atoms */
CLJAtoms CLJGrid::fixedAtoms() const
{
    return cljboxes.atoms().squeeze();
}

/** Set the function used to calculate the coulomb and LJ energy */
void CLJGrid::setCLJFunction(const CLJFunction &function)
{
    cljfunc = function;
    checkIfGridSupported();
    clearGrid();
}

/** Return the function used to calculate the coulomb and LJ energy */
const CLJFunction& CLJGrid::cljFunction() const
{
    return cljfunc.read();
}

/** Enable use of a parallel multicore algorithm to calculate energies */
void CLJGrid::enableParallelCalculation()
{
    setUseParallelCalculation(true);
}

/** Disable use of the parallel algorithm for calculating energies */
void CLJGrid::disableParallelCalculation()
{
    setUseParallelCalculation(false);
}

/** Switch on or off use of a parallel multicore algorithm for calculating
    energies */
void CLJGrid::setUseParallelCalculation(bool on)
{
    parallel_calc = on;
}

/** Return whether or not a parallel algorithm is being used
    to calculate grid energies */
bool CLJGrid::usesParallelCalculation() const
{
    return parallel_calc;
}

/** Turn on an energy summing algorithm that guarantees the same energy
    regardless of whether a single core or multicore calculation is being
    performed (i.e. rounding errors in both cases will be identical) */
void CLJGrid::enableReproducibleCalculation()
{
    setUseReproducibleCalculation(true);
}

/** Turn off an energy summing algorithm that guarantees the same energy
    regardless of whether a single core or multicore calculation is being
    performed (i.e. rounding errors in both cases will not be identical) */
void CLJGrid::disableReproducibleCalculation()
{
    setUseReproducibleCalculation(false);
}

/** Switch on or off use of an energy summing algorithm that guarantees the 
    same energy regardless of whether a single core or multicore calculation 
    is being performed */
void CLJGrid::setUseReproducibleCalculation(bool on)
{
    repro_sum = on;
}

/** Return whether or not a reproducible energy summing algorithm is being
    used to accumulate the energies */
bool CLJGrid::usesReproducibleCalculation() const
{
    return repro_sum;
}

/** Set the grid spacing */
void CLJGrid::setGridSpacing(Length spacing)
{
    if (grid().spacing() != spacing)
    {
        setGrid( GridInfo(grid().dimensions(), spacing) );
    }
}

/** Return the grid spacing */
Length CLJGrid::gridSpacing() const
{
    return grid_info.spacing();
}

/** Set the grid buffer */
void CLJGrid::setGridBuffer(Length buffer)
{
    grid_buffer = buffer.value();
}

/** Return the grid buffer */
Length CLJGrid::gridBuffer() const
{
    return Length(grid_buffer);
}

/** Set the grid to be used for the optimised calculation */
void CLJGrid::setGrid(const GridInfo &grid)
{
    if (grid_info != grid)
    {
        grid_info = grid;
        clearGrid();
    }
}

/** Set the dimensions of the grid */
void CLJGrid::setGridDimensions(const AABox &dimensions)
{
    if (grid_info.dimensions() != dimensions)
    {
        setGrid( GridInfo(dimensions, grid_info.spacing()) );
    }
}

/** Set the dimensions of the grid so that it encompasses all of the atoms in 'atoms'
    (with a buffer of 'buffer' around all atoms) using a grid with spacing 'grid_spacing' */
void CLJGrid::setGridDimensions(const CLJAtoms &atoms, Length spacing, Length buffer)
{
    grid_buffer = buffer.value();
    
    if (atoms.x().isEmpty())
    {
        //build an empty grid
        setGrid( GridInfo( AABox::from( -Vector(buffer), Vector(buffer) ), spacing ) );
    }
    else
    {
        MultiFloat min_x( std::numeric_limits<float>::max() );
        MultiFloat max_x( -std::numeric_limits<float>::max() );
        
        MultiFloat min_y( min_x );
        MultiFloat max_y( max_x );
        
        MultiFloat min_z( min_x );
        MultiFloat max_z( max_x );
        
        const int nats = atoms.x().count();
        
        const MultiFloat *x = atoms.x().constData();
        const MultiFloat *y = atoms.y().constData();
        const MultiFloat *z = atoms.z().constData();
        const MultiInt *id = atoms.ID().constData();
        
        const qint32 dummy_id = CLJAtoms::idOfDummy()[0];
        
        for (int i=0; i<nats; ++i)
        {
            bool has_dummy = false;
        
            for (int j=0; j<MultiInt::count(); ++j)
            {
                if (id[i][j] == dummy_id)
                {
                    has_dummy = true;
                    break;
                }
            }
        
            if (has_dummy)
            {
                for (int j=0; j<MultiInt::count(); ++j)
                {
                    if (id[i][j] != dummy_id)
                    {
                        min_x.set(j, qMin(min_x[j], x[i][j]));
                        min_y.set(j, qMin(min_y[j], y[i][j]));
                        min_z.set(j, qMin(min_z[j], z[i][j]));
                        
                        max_x.set(j, qMax(max_x[j], x[i][j]));
                        max_y.set(j, qMax(max_y[j], y[i][j]));
                        max_z.set(j, qMax(max_z[j], z[i][j]));
                    }
                }
            }
            else
            {
                min_x = min_x.min(x[i]);
                max_x = max_x.max(x[i]);

                min_y = min_y.min(y[i]);
                max_y = max_y.max(y[i]);

                min_z = min_z.min(z[i]);
                max_z = max_z.max(z[i]);
            }
        }
        
        Vector mincoords( min_x.min(), min_y.min(), min_z.min() );
        Vector maxcoords( max_x.max(), max_y.max(), max_z.max() );
        
        mincoords -= Vector(grid_buffer);
        maxcoords += Vector(grid_buffer);
        
        setGrid( GridInfo( AABox::from(mincoords,maxcoords), spacing ) );
    }
}

/** Set the dimensions of the grid so that it encompasses all of the atoms in 'atoms'
    (plus a gridBuffer() buffer around all atoms) */
void CLJGrid::setGridDimensions(const CLJAtoms &atoms)
{
    setGridDimensions(atoms, grid_info.spacing(), Length(grid_buffer));
}

/** Set the dimensions of the grid so that it encompasses all of the atoms in 'atoms'
    (plus a gridBuffer() buffer around all atoms) using a grid with spacing 'grid_spacing' */
void CLJGrid::setGridDimensions(const CLJAtoms &atoms, Length spacing)
{
    setGridDimensions(atoms, spacing, Length(grid_buffer));
}

/** Return the dimensions of the grid */
AABox CLJGrid::gridDimensions() const
{
    return grid_info.dimensions();
}

/** Return the grid */
GridInfo CLJGrid::grid() const
{
    return grid_info;
}

/** Switch on or off the use of a grid - note that the grid may
    still be disabled if the CLJFunction is not compatible with
    use of a grid */
void CLJGrid::setUseGrid(bool on)
{
    use_grid = on;
    
    if (not use_grid)
        clearGrid();
}

/** Enable used of the grid - note that the grid may
    still be disabled if the CLJFunction is not compatible with
    use of a grid */
void CLJGrid::enableGrid()
{
    use_grid = true;
}

/** Disable use of the grid */
void CLJGrid::disableGrid()
{
    use_grid = false;
    clearGrid();
}

/** Return whether or not the grid will be used in the calculation */
bool CLJGrid::usesGrid() const
{
    return use_grid and cljfunc_supports_grid;
}

/** Return whether or not the CLJFunction supports use of the grid.
    To support a grid, the function must intrinsically support the grid,
    and the coulomb cutoff must be much greater than the LJ cutoff */
bool CLJGrid::functionSupportsGrid() const
{
    return cljfunc_supports_grid;
}

static QMutex const_update_mutex;

/** Function called to calculate the potential grid. Note that this is called by 
    a const function, so we must be thread-safe when we update the potential */
void CLJGrid::calculateGrid()
{
    //first, divide the atoms into two sets - those within the LJ cutoff of
    //the grid, and those outside
    CLJAtoms far_atoms;
    CLJAtoms near_atoms;
    {
        QVector<CLJAtom> near_atms;
        QVector<CLJAtom> far_atms;

        const Space &space = cljfunc.read().space();
        const float lj_cutoff = cljfunc->ljCutoff();
        const float coul_cutoff = cljfunc->coulombCutoff();

        for (CLJBoxes::const_iterator it = cljboxes.constBegin();
             it != cljboxes.constEnd();
             ++it)
        {
            const CLJAtoms &atoms = it->read().atoms();
            
            for (int i=0; i<atoms.count(); ++i)
            {
                CLJAtom atom = atoms[i];
                
                if (atom.ID() != 0)
                {
                    double mindist = space.minimumDistance(atom.coordinates(),
                                                           grid_info.dimensions());
                    
                    if (mindist < lj_cutoff)
                    {
                        near_atms.append(atom);
                    }
                    else if (mindist < coul_cutoff and atom.charge().value() != 0)
                    {
                        far_atms.append(atom);
                    }
                }
            }
        }
        
        near_atoms = CLJAtoms(near_atms);
        far_atoms = CLJAtoms(far_atms);
    }
    
    //now, go through any far atoms and add their potentials to the grid
    //if (parallel_calc)
    //{
    //    write a parallel algorithm for calculating the grid - divide the entire
    //    grid into boxes that can be evaluated in parallel
    //}
    //else
    //{
            QVector<float> pot = cljfunc.read().calculate(far_atoms, grid_info);
    //}
    
    //update the object - note that because this is called from a const function
    //we have to be doubly sure that this has not been called twice from two
    //different threads.
    QMutexLocker lkr(&const_update_mutex);
    
    if (grid_pots.isEmpty())
    {
        grid_pots = pot;
        close_atoms = CLJBoxes(near_atoms);
    }
}

/** Calculate the total energy of interaction between the passed atoms and
    all of the atoms added to the grid */
void CLJGrid::total(const CLJBoxes &atoms, double &cnrg, double &ljnrg) const
{
    cnrg = 0;
    ljnrg = 0;
    
    if (cljboxes.isEmpty() or atoms.isEmpty())
        return;
    
    if (use_grid and cljfunc_supports_grid)
    {
        //there is a big enough difference between the coulomb and LJ cutoffs that
        //a grid is worthwhile
        if (grid_pots.isEmpty())
        {
            //we haven't yet calculated the grid! We need to (just in time) calculate it.
            const_cast<CLJGrid*>(this)->calculateGrid();
        }

        //calculate the energy between the atoms and the close atoms
        tuple<double,double> nrgs;
        
        if (parallel_calc)
        {
            CLJCalculator cljcalc(repro_sum);
            nrgs = cljcalc.calculate(cljfunc.read(), close_atoms, atoms);
        }
        else
        {
            nrgs = cljfunc.read().calculate(atoms, close_atoms);
        }

        //now calculate the grid energy of each atom
        const float *gridpot_array = grid_pots.constData();

        bool all_within_grid = true;

        QVector<MultiInt> grid_corners;
        QVector<MultiFloat> grid_weights;

        MultiDouble grid_nrg(0);
    
        const qint32 dummy_id = CLJAtoms::idOfDummy()[0];
        const qint32 grid_id = idOfFixedAtom();

        const MultiInt m_dummy_id(dummy_id);
        const MultiInt m_grid_id(grid_id);

        for (CLJBoxes::const_iterator it = atoms.constBegin();
             it != atoms.constEnd();
             ++it)
        {
            const CLJAtoms &atms = it->read().atoms();
        
            const MultiFloat *x = atms.x().constData();
            const MultiFloat *y = atms.y().constData();
            const MultiFloat *z = atms.z().constData();
            const MultiFloat *q = atms.q().constData();
            const MultiInt *id = atms.ID().constData();
            
            const int nats = atms.x().count();

            for (int i=0; i<nats; ++i)
            {

                int n_in_grid = grid_info.pointToGridCorners(x[i], y[i], z[i],
                                                             grid_corners, grid_weights);

                if (n_in_grid != MultiFloat::count())
                {
                    // NEED TO CHECK IF ANY OF THE POINTS ARE DUMMIES. IF THEY ARE
                    // THEN WE NEED TO CALCULATE THE SUM MANUALLY
                    int ndummies = 0;
                    
                    for (int j=0; j<MultiFloat::count(); ++j)
                    {
                        if (id[i][j] == dummy_id)
                        {
                            ndummies += 1;
                            
                            //put the dummy into the first grid box
                            for (int k=0; k<8; ++k)
                            {
                                grid_corners[k].set(j, 0);
                                grid_weights[k].set(j, 0.0);
                            }
                        }
                    }
                    
                    if (ndummies + n_in_grid != MultiFloat::count())
                    {
                        //at least one of the grid points is outside of the grid
                        qDebug() << "POINT" << x[i].toString() << y[i].toString()
                                            << z[i].toString() << "LIES OUTSIDE OF"
                                 << "THE GRID?" << grid_info.toString();
                        
                        all_within_grid = false;
                        break;
                    }
                }

                MultiFloat phi(0);
                
                for (int j=0; j<8; ++j)
                {
                    phi += MultiFloat(gridpot_array, grid_corners.constData()[j]) *
                                      grid_weights.constData()[j];
                }


                //add the energy of these atoms onto the total, taking care to ignore
                //dummy atoms and to ignore atoms with IDs equal to the grid ID
                //(this allows atoms to be screened, e.g. for "intramolecular" calculations
                // where some of the molecule is fixed and on the grid, while the rest
                // is mobile)
                grid_nrg += (q[i] * phi).logicalAndNot( id[i].compareEqual(m_dummy_id) )
                                        .logicalAndNot( id[i].compareEqual(m_grid_id) );
            }
            
            if (not all_within_grid)
                break;
        }
        
        if (all_within_grid)
        {
            cnrg = nrgs.get<0>() + grid_nrg.sum();
            ljnrg = nrgs.get<1>();
            
            return;
        }
    }
    
    //either the grid is not used or something went wrong with the grid calculation
    //Do not use a grid
    if (parallel_calc)
    {
        CLJCalculator cljcalc(repro_sum);
        tuple<double,double> nrgs = cljcalc.calculate(cljfunc.read(), atoms, cljboxes);
        
        cnrg = nrgs.get<0>();
        ljnrg = nrgs.get<1>();
    }
    else
    {
        tuple<double,double> nrgs = cljfunc.read().calculate(atoms, cljboxes);
        
        cnrg = nrgs.get<0>();
        ljnrg = nrgs.get<1>();
    }
}

/** Calculate the total energy of interaction between the passed atoms and
    all of the atoms added to the grid */
void CLJGrid::total(const CLJAtoms &atoms, double &cnrg, double &ljnrg) const
{
    CLJBoxes boxes(atoms);
    this->total(boxes, cnrg, ljnrg);
}

/** Return the coulomb and LJ energies of the passed atoms with the fixed
    atoms added to this grid */
boost::tuple<double,double> CLJGrid::calculate(const CLJAtoms &atoms) const
{
    double cnrg, ljnrg;
    this->total(atoms, cnrg, ljnrg);
    return boost::tuple<double,double>(cnrg,ljnrg);
}

/** Return the coulomb and LJ energies of the passed atoms with the fixed
    atoms added to this grid */
boost::tuple<double,double> CLJGrid::calculate(const CLJBoxes &atoms) const
{
    double cnrg, ljnrg;
    this->total(atoms, cnrg, ljnrg);
    return boost::tuple<double,double>(cnrg,ljnrg);
}

/** Return the coulomb energy of the passed atoms interacting with 
    the fixed atoms on this grid */
double CLJGrid::coulomb(const CLJAtoms &atoms) const
{
    double cnrg, ljnrg;
    this->total(atoms, cnrg, ljnrg);
    return cnrg;
}

/** Return the coulomb energy of the passed atoms interacting with 
    the fixed atoms on this grid */
double CLJGrid::coulomb(const CLJBoxes &atoms) const
{
    double cnrg, ljnrg;
    this->total(atoms, cnrg, ljnrg);
    return cnrg;
}

/** Return the LJ energy of the passed atoms interacting with
    the fixed atoms on this grid */
double CLJGrid::lj(const CLJAtoms &atoms) const
{
    double cnrg, ljnrg;
    this->total(atoms, cnrg, ljnrg);
    return ljnrg;
}

/** Return the LJ energy of the passed atoms interacting with
    the fixed atoms on this grid */
double CLJGrid::lj(const CLJBoxes &atoms) const
{
    double cnrg, ljnrg;
    this->total(atoms, cnrg, ljnrg);
    return ljnrg;
}
