/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2012  Christopher Woods
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

#include "idassigner.h"

#include <QVarLengthArray>

#include "system.h"
#include "closemols.h"
#include "delta.h"

#include "SireSystem/errors.h"

#include "SireMaths/linearap.h"
#include "SireMaths/nmatrix.h"
#include "SireMaths/nvector.h"

#include "SireMol/molecules.h"
#include "SireMol/molecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/moleditor.h"
#include "SireMol/viewsofmol.h"

#include "SireVol/space.h"

#include "tostring.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QTime>
#include <QDebug>

using namespace SireSystem;
using namespace SireFF;
using namespace SireMol;
using namespace SireBase;
using namespace SireMaths;
using namespace SireVol;
using namespace SireStream;

static const RegisterMetaType<IDAssigner> r_idassigner;

QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds, const IDAssigner &idassigner)
{
    writeHeader(ds, r_idassigner, 2);
    
    SharedDataStream sds(ds);
    
    sds << idassigner.molgroup << idassigner.identity_points
        << idassigner.spce << idassigner.map;
    
    return ds;
}

QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, IDAssigner &idassigner)
{
    VersionID v = readHeader(ds, r_idassigner);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        MolGroupPtr molgroup;
        QVector<PointPtr> points;
        SpacePtr space;
        PropertyMap map;
        
        sds >> molgroup >> points >> space >> map;
        
        idassigner = IDAssigner(points, molgroup, space, map);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        MolGroupPtr molgroup;
        QVector<PointPtr> points;
        PropertyMap map;
        
        sds >> molgroup >> points >> map;
        
        idassigner = IDAssigner(points, molgroup, map);
    }
    else
        throw version_error(v, "1,2", r_idassigner, CODELOC);
        
    return ds;
}

void IDAssigner::validatePoints(const QVector<PointPtr> &points) const
{
    for (QVector<PointPtr>::const_iterator it = points.constBegin();
         it != points.constEnd();
         ++it)
    {
        if ( not (it->read().isIntraMoleculePoint() or
                  it->read().isExtraMoleculePoint()) )
        {
            throw SireError::incompatible_error( QObject::tr(
                    "An IDAssigner can only be used with "
                    "intramolecular or extramolecular identity points."),
                        CODELOC );
        }
    }
}

void IDAssigner::validateGroup(const MoleculeGroup &new_group) const
{
    //run through each molecule and ensure that the coordinate properties
    //of them all are compatible
    
    if (new_group.nMolecules() < 2)
        return;
    
    const PropertyName &coords_property = map["coordinates"];
    
    const ViewsOfMol &first_mol = new_group.moleculeAt(0);
    
    const AtomCoords &coords = first_mol.data().property(coords_property)
                                               .asA<AtomCoords>();
    
    for ( Molecules::const_iterator it = new_group.constBegin(); 
          it != new_group.constEnd(); 
          ++it )
    {
        if (not coords.isCompatibleWith(it->data().info()))
            throw SireError::incompatible_error( QObject::tr(
                    "The coordinates property (%1) of molecule number %2 "
                    "is not compatible with that of molecule number %3. "
                    "Only molecules with completely compatible coordinates "
                    "properties can be used in an IDAssigner.")
                        .arg(coords_property.toString())
                        .arg(first_mol.number())
                        .arg(it.key()), CODELOC );
    }
}


/** Constructor */
IDAssigner::IDAssigner() 
           : ConcreteProperty<IDAssigner,Property>(), distances_changed(false)
{}

/** Construct to find the identity of the molecules from
    'molgroup' that match the points in 'points' - 
    this creates a list of n molecules, where the ith molecule 
    is matched to the ith point */
IDAssigner::IDAssigner(const QVector<PointPtr> &points,
                       const MoleculeGroup &group,
                       const Space &space,
                       const PropertyMap &property_map)
           : ConcreteProperty<IDAssigner,Property>(),
             molgroup(group), identity_points(points), spce(space), map(property_map),
             distances_changed(false)
{
    this->validateGroup(group);
    this->validatePoints(points);
    identity_points.squeeze();

    if (not identity_points.isEmpty())
    {
        const int npoints = identity_points.count();
        
        //construct a CloseMols to record the closest 'npoints' molecules
        //to each of the identity points
        points_with_mols.reserve(npoints);
        
        for (int i=0; i<npoints; ++i)
        {
            points_with_mols.append( CloseMols(identity_points.at(i), 
                                               molgroup, spce, npoints, map) );
        }
        
        points_with_mols.squeeze();
        
        this->rebuildMolToMolNum();
        this->recalculateDistances();
        this->assignMoleculesToPoints();
    }
}

/** Construct to find the identity of the molecules from
    'molgroup' that match the points in 'points' - 
    this creates a list of n molecules, where the ith molecule 
    is matched to the ith point */
IDAssigner::IDAssigner(const QVector<PointPtr> &points,
                       const MoleculeGroup &group,
                       const PropertyMap &property_map)
           : ConcreteProperty<IDAssigner,Property>()
{
    this->operator=( IDAssigner(points, group, Cartesian(), property_map) );
}

/** Constructor */
IDAssigner::IDAssigner(const PointRef &point,
                       const MoleculeGroup &group,
                       const PropertyMap &property_map)
           : ConcreteProperty<IDAssigner,Property>()
{
    QVector<PointPtr> points;
    points.append( PointPtr(point) );
    this->operator=( IDAssigner(points,group,property_map) );
}

/** Constructor */
IDAssigner::IDAssigner(const PointRef &point,
                       const MoleculeGroup &group,
                       const Space &space,
                       const PropertyMap &property_map)
           : ConcreteProperty<IDAssigner,Property>()
{
    QVector<PointPtr> points;
    points.append( PointPtr(point) );
    this->operator=( IDAssigner(points,group,space,property_map) );
}
  
/** Copy constructor */                         
IDAssigner::IDAssigner(const IDAssigner &other)
           : ConcreteProperty<IDAssigner,Property>(other),
             molgroup(other.molgroup),
             identity_points(other.identity_points),
             spce(other.spce),
             map(other.map), 
             points_with_mols(other.points_with_mols),
             mol_to_molnum(other.mol_to_molnum),
             point_distances(other.point_distances),
             point_to_mol(other.point_to_mol),
             distances_changed(other.distances_changed)
{}

/** Destructor */
IDAssigner::~IDAssigner()
{}

/** Copy assignment operator */
IDAssigner& IDAssigner::operator=(const IDAssigner &other)
{
    if (this != &other)
    {
        molgroup = other.molgroup;
        identity_points = other.identity_points;
        spce = other.spce;
        map = other.map;
        points_with_mols = other.points_with_mols;
        mol_to_molnum = other.mol_to_molnum;
        point_distances = other.point_distances;
        point_to_mol = other.point_to_mol;
        distances_changed = other.distances_changed;
    }
    
    return *this;
}

/** Comparison function */
bool IDAssigner::operator==(const IDAssigner &other) const
{
    return molgroup == other.molgroup and identity_points == other.identity_points
              and map == other.map;
}

/** Comparison function */
bool IDAssigner::operator!=(const IDAssigner &other) const
{
    return not IDAssigner::operator==(other);
}

Q_GLOBAL_STATIC( MoleculeGroup, nullMoleculeGroup );

/** Return the molecule group operated on by this constraint */
const MoleculeGroup& IDAssigner::moleculeGroup() const
{
    if (molgroup.isNull())
        return *(nullMoleculeGroup());
    else
        return molgroup.read();
}

/** Return the points used to identify the molecules */
QVector<PointPtr> IDAssigner::points() const
{
    return identity_points;
}

/** Return the space used to calculate distances between the molecules
    and the identity points */
const Space& IDAssigner::space() const
{
    return spce.read();
}

/** Return the property map used to find the properties required
    by this constraint */
const PropertyMap& IDAssigner::propertyMap() const
{
    return map;
}

/** Update the space used to calculate the distances between
    the points and the molecules - this returns whether or 
    not the space has changed */
bool IDAssigner::updateSpace(const System &system)
{
    const Space &new_space = system.property(map["space"]).asA<Space>();
    
    if (spce != new_space)
    {
        spce = new_space;
        return true;
    }
    else
        return false;
}

/** Update the points in this constraint from the passed system - 
    this returns whether or not this changes any points */
bool IDAssigner::updatePoints(const System &system)
{
    int npoints = identity_points.count();
    
    const PointPtr *const_points_array = identity_points.constData();
    
    bool need_update = false;
    
    for (int i=0; i<npoints; ++i)
    {
        if (const_points_array[i].read().usesMoleculesIn(system))
        {
            need_update = true;
            break;
        }
    }
    
    if (need_update)
    {
        bool changed = false;
        
        PointPtr *points_array = identity_points.data();
        
        for (int i=0; i<npoints; ++i)
        {
            bool this_changed = points_array[i].edit().update(system);
            
            changed = changed or this_changed;
        }
        
        return changed;
    }
    else
        return false;
}

/** Update the molecule group whose molecules are constrained to equal
    the version from the passed system. */
bool IDAssigner::updateGroup(const System &system)
{
    const MoleculeGroup &old_group = molgroup.read();

    if (not system.contains(old_group.number()))
    {
        molgroup.edit().update(system.molecules());
        return true;
    }
        
    const MoleculeGroup &new_group = system[old_group.number()];
    
    if (new_group.version() == old_group.version())
        return false;
        
    else if (new_group.version().majorVersion() != old_group.version().majorVersion())
    {
        //the group's contents have changed - ensure that all
        //of the molecules are still compatible
        this->validateGroup(new_group);
    }
    
    molgroup = new_group;
    return true;
}

const char* IDAssigner::typeName()
{
    return QMetaType::typeName( qMetaTypeId<IDAssigner>() );
}

const char* IDAssigner::what() const
{
    return IDAssigner::typeName();
}

IDAssigner* IDAssigner::clone() const
{
    return new IDAssigner(*this);
}

/** Internal function used to recalculate all of the distances^2 
    between the molecules in 'mol_to_molnum' and all of the identity points */
void IDAssigner::recalculateDistances()
{
    point_distances = QHash<MolNum,NVector>();
    
    const Molecules &molecules = this->moleculeGroup().molecules();
    
    const int nmols = mol_to_molnum.count();
    const MolNum *mol_to_molnum_array = mol_to_molnum.constData();
    
    point_distances.reserve(nmols);
    
    const int npoints = identity_points.count();
    const PointPtr *const_points_array = identity_points.constData();
    
    QVector< Vector > points(npoints);
    Vector *points_array = points.data();
    
    for (int i=0; i<npoints; ++i)
    {
        points_array[i] = const_points_array[i].read().point();
    }
    
    const PropertyName &coords_property = map["coordinates"];
    
    for (int i=0; i<nmols; ++i)
    {
        const MolNum &molnum = mol_to_molnum_array[i];
    
        Molecules::const_iterator it = molecules.constFind(molnum);
        BOOST_ASSERT( it != molecules.constEnd() );
    
        const AtomCoords &coords = it->data().property(coords_property)
                                             .asA<AtomCoords>();
                                             
        Vector center = coords.array().aaBox().center();
        
        NVector distances(npoints);
        double *distances_array = distances.data();
        
        for (int j=0; j<npoints; ++j)
        {
            double dist2 = space().calcDist2(center, points_array[j]);
            distances_array[j] = dist2;
        }
        
        point_distances.insert( molnum, distances );
    }
    
    distances_changed = true;
}

/** Internal function used to get the identities of the candidate molecules
    for the points (the first 'npoints' molecules in the molecule group, plus
    the closest 'npoints' molecules to each identity point) and the mapping
    of candidate molecule index to molecule number (with molecule number
    matching the order in the molecule group) */
void IDAssigner::rebuildMolToMolNum()
{
    mol_to_molnum = QVector<MolNum>();
    
    const int npoints = points_with_mols.count();
    const QVector<MolNum> &molnums = molgroup.read().molNums();
    const int nmols = molnums.count();
    
    //reserve space - worst case is we have to record all npoints molecules
    //from all npoints points, together with the npoints assigned molecules
    mol_to_molnum.reserve( qMin(nmols, npoints*(npoints+1)) );
    
    //loop through the molecules in the order they appear in the 
    //molecule group and see if they are a candidate - the first 'npoints'
    //molecules are automatically candidates (as they are assigned)
    for (int i=0; i<nmols; ++i)
    {
        const MolNum &molnum = molnums.at(i);
    
        //is this molecule one of the close molecules to any 
        //of the points?
        for (int j=0; j<npoints; ++j)
        {
            if (points_with_mols.at(j).isClose(molnum))
            {
                mol_to_molnum.append(molnum);
                break;
            }
        }
    }

    mol_to_molnum.squeeze();
}

/** Internal function used to invert a mapping of points to molecules
    into a mapping of molecules to points */
static QVector<int> invert(const QVector<int> &point_to_mol)
{
    QVector<int> mol_to_point = point_to_mol;
    
    const int *point_to_mol_array = point_to_mol.constData();
    int *mol_to_point_array = mol_to_point.data();
    
    const int npoints = point_to_mol.count();
    
    for (int i=0; i<npoints; ++i)
    {
        mol_to_point_array[ point_to_mol_array[i] ] = i;
    }
    
    return mol_to_point;
}

/** This function uses the distances between all points and molecules
    stored in 'point_distances' to work out the optimum assignment
    of molecules to points such that the total distance between
    each molecule and each point is minimised */
void IDAssigner::assignMoleculesToPoints()
{
    if (not distances_changed)
        return;
    
    //the order of molecules may have changed - recalculate 
    //the correct order

    //use the 'mol_to_molnum' array as this holds the numbers
    //of all of the candidate molecules in the same order as
    //they appear in the molecule group
    const int nmols = mol_to_molnum.count();
    const MolNum *mol_to_molnum_array = mol_to_molnum.constData();
    
    const int npoints = identity_points.count();
    
    //construct the matrix that contains the distances between every
    //candidate molecule and every point - one molecule per row, one point
    //per column - this has to be a square matrix, so missing rows/columns
    //are given a value of 0
    NMatrix distmatrix;
    
    if (nmols == npoints)
    {
        distmatrix = NMatrix(nmols, nmols);
        distmatrix = distmatrix.transpose(); // change to row-major memory order

        for (int i=0; i<nmols; ++i)
        {
            const MolNum &molnum = mol_to_molnum_array[i];
            
            QHash<MolNum,NVector>::const_iterator 
                                       it = point_distances.constFind(molnum);
                                                
            BOOST_ASSERT( it != point_distances.constEnd() );
            
            distmatrix.setRow(i, it.value());
        }
    }
    else if (nmols > npoints)
    {
        //there are more molecules than points - we create some extra
        //points which have zero distance to all molecules
        distmatrix = NMatrix(nmols, nmols);
        distmatrix = distmatrix.transpose(); // change to row-major memory order
        
        const int nzeroes = nmols - npoints;
        NVector new_row(npoints + nzeroes);
        
        for (int i=0; i<nmols; ++i)
        {
            const MolNum &molnum = mol_to_molnum_array[i];
            
            QHash<MolNum,NVector>::const_iterator 
                                   it = point_distances.constFind(molnum);
                                                
            BOOST_ASSERT( it != point_distances.constEnd() );
            
            const NVector &distances = it.value();
            
            memcpy( new_row.data(), distances.constData(), npoints*sizeof(double) );
            
            distmatrix.setRow( i, new_row );
        }
    }
    else
    {
        //there are more points than molecules - we create some extra
        //molecules that are all equally a very long way from all of the points
        distmatrix = NMatrix(npoints, npoints, std::numeric_limits<double>::max());
        distmatrix = distmatrix.transpose(); // change to row-major memory order

        //copy the distances to the real molecules
        for (int i=0; i<nmols; ++i)
        {
            const MolNum &molnum = mol_to_molnum_array[i];
            
            QHash<MolNum,NVector>::const_iterator 
                                        it = point_distances.constFind(molnum);
                                                
            BOOST_ASSERT( it != point_distances.constEnd() );
            
            distmatrix.setRow(i, it.value());
        }
    }

    //find the smallest and second smallest absolute difference between distances
    // - this is used to see if there are any degeneracies
    //if (false)
    {
        double delta0 = std::numeric_limits<double>::max();
        double delta1 = delta0;
    
        const int nrows = distmatrix.nRows();
        const int ncolumns = qMin( distmatrix.nColumns(), npoints );
        
        const double *distmatrix_array = distmatrix.constData();
        
        for (int i=0; i<nrows; ++i)
        {
            const double *row = &(distmatrix_array[ distmatrix.offset(i,0) ]);
            
            for (int j=0; j<ncolumns-1; ++j)
            {
                for (int k=j+1; k<ncolumns; ++k)
                {
                    const double delta = std::abs(row[k] - row[j]);
                    
                    if (delta <= delta0)
                    {
                        if (delta < delta0)
                        {
                            delta1 = delta0;
                            delta0 = delta;
                        }
                    }
                    else if (delta < delta1)
                    {
                        delta1 = delta;
                    }
                }
            }
        }
    
        if (delta0 == 0)
        {
            //there are degeneracies - some molecules are an identical
            //distance from some points
            if (delta1 == std::numeric_limits<double>::max())
                //really degenerate!
                delta1 = 0.1;
    
            //add a small penalty function that make molecule 'i' prefer
            //point 'i', or the point as close to point 'i' as possible. This
            //is used to help remove degeneracies caused by using points with
            //the same coordinates. If multiple points have the same coordinates,
            //then the molecule with index closest to the index of the point
            //will be preferred - this penalty function has to be kept to 
            //about 10% of the smallest non-zero difference between distances,
            //so that it doesn't affect the assignment of points with no degeneracies

            const double scl = 0.1 * delta1 / (ncolumns * nrows);

            for (int i=0; i<nrows; ++i)
            {
                for (int j=0; j<ncolumns; ++j)
                {
                    distmatrix(i,j) += scl * (i-j) * (i-j);
                }
            }
        }
    }
    
    //now calculate optimum assignment of molecules to points that
    //minimises the total distance between each molecule and its
    //assigned point
    point_to_mol = solve_linear_assignment(distmatrix, true);
    
    point_to_mol = ::invert(point_to_mol);

    //if there are more molecules than points, then the last set of 
    //molecules are not associated with points. To ensure deterministic
    //mapping the extra molecules must be ordered so that they are
    //assigned to the null points in order
    if (nmols > npoints)
    {
        qSort( point_to_mol.data() + npoints, point_to_mol.data() + nmols );
    }
    
    point_to_mol.squeeze();
    distances_changed = false;
}

/** Return the number of identity points (and thus the number of 
    identified molecules) */
int IDAssigner::nPoints() const
{
    return identity_points.count();
}

/** Return a string representation of this assigner */
QString IDAssigner::toString() const
{
    QStringList lines;
    
    for (QVector<PointPtr>::const_iterator it = identity_points.constBegin();
         it != identity_points.constEnd();
         ++it)
    {
        lines.append( (*it)->toString() );
    }
    
    return QObject::tr("IDAssigner( points() => [ %1 ], moleculeGroup() = %2 )")
            .arg(lines.join(", "), molgroup->toString());
}

/** Update the assigner with the passed system.  */
void IDAssigner::update(const System &system)
{
    if (identity_points.isEmpty())
        return;

    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);
    bool new_space = this->updateSpace(system);

    if (new_group or new_points or new_space)
    {
        //update all of the points
        const int npoints = points_with_mols.count();

        bool closemols_changed = false;
    
        for (int i=0; i<npoints; ++i)
        {
            bool this_changed = points_with_mols[i].update(system);
            closemols_changed = closemols_changed or this_changed;
        }

        if (closemols_changed)
            this->rebuildMolToMolNum();

        this->recalculateDistances();
    }

    this->assignMoleculesToPoints();
}

/** Returns the list of identified molecules from the system, 
    which are returned in the same order as the list of identity points
*/
QVector<PartialMolecule> IDAssigner::identifiedMolecules() const
{
    //find each matching molecule in turn
    const Molecules molecules = this->moleculeGroup().molecules();
    
    int n_to_match = qMin( identity_points.count(), molecules.count() );
    
    const int *point_to_mol_array = point_to_mol.constData();

    QVector<PartialMolecule> matched_mols;
    matched_mols.reserve(point_to_mol.count());
    
    //the match uses mol_to_molnum, as only a subset of molecules
    //are candidate molecules
    const MolNum *mol_to_molnum_array = mol_to_molnum.constData();
    
    if (mol_to_molnum.count() < n_to_match)
        throw SireError::program_bug( QObject::tr(
                "WEIRD: Number of matched molecules is greater than the number of "
                "identified molecules? %1, %2, %3")
                    .arg(n_to_match).arg(mol_to_molnum.count())
                    .arg(Sire::toString(mol_to_molnum)), CODELOC );

    if (point_to_mol.count() < n_to_match)
        throw SireError::program_bug( QObject::tr(
                "WEIRD: Number of point_to_mol matched molecules is less than the number of "
                "identified molecules? %1, %2, %3")
                    .arg(n_to_match).arg(point_to_mol.count())
                    .arg(Sire::toString(point_to_mol)), CODELOC );
    
    for (int i=0; i<n_to_match; ++i)
    {
        const MolNum molnum = mol_to_molnum_array[ point_to_mol_array[i] ];

        matched_mols.append( PartialMolecule(molecules[molnum]) );
    }
    
    return matched_mols;
}
