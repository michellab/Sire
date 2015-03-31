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

#ifndef SIREMM_GRIDFF_H
#define SIREMM_GRIDFF_H

#include "intercljff.h"
#include "SireMaths/histogram.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class GridFF;
}

QDataStream& operator<<(QDataStream&, const SireMM::GridFF&);
QDataStream& operator>>(QDataStream&, SireMM::GridFF&);

namespace SireMM
{

/** This class calculates the coulomb and LJ energy between
    all molecules in group 1 and all molecules in group 2.
    The calculation is optimised, as the molecules in group 2
    are represented using a grid. This is ideal for situations
    where the molecules on group 2 move little, or not at all.
    
    @author Christopher Woods
*/
class SIREMM_EXPORT GridFF 
            : public SireBase::ConcreteProperty<GridFF,InterGroupCLJFF>
{

friend QDataStream& ::operator<<(QDataStream&, const GridFF&);
friend QDataStream& ::operator>>(QDataStream&, GridFF&);

public:
    GridFF();
    GridFF(const QString &name);
    
    GridFF(const GridFF &other);
    
    ~GridFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    GridFF& operator=(const GridFF &other);
    
    bool operator==(const GridFF &other) const;
    bool operator!=(const GridFF &other) const;
    
    GridFF* clone() const;

    void addFixedAtoms(const MoleculeView &fixed_atoms,
                       const SireBase::PropertyMap &map = SireBase::PropertyMap());
    void addFixedAtoms(const SireMol::Molecules &fixed_atoms,
                       const SireBase::PropertyMap &map = SireBase::PropertyMap());
    void addFixedAtoms(const SireMol::MoleculeGroup &fixed_atoms,
                       const SireBase::PropertyMap &map = SireBase::PropertyMap());

    void addFixedAtoms(const GridFF &other);

    void setBuffer(SireUnits::Dimension::Length buffer);
    void setGridSpacing(SireUnits::Dimension::Length spacing);

    void setCoulombCutoff(SireUnits::Dimension::Length cutoff);
    void setLJCutoff(SireUnits::Dimension::Length cutoff);

    bool setShiftElectrostatics(bool on);
    
    bool setUseReactionField(bool on);
    bool setReactionFieldDielectric(double dielectric);
    
    SireUnits::Dimension::Length buffer() const;
    SireUnits::Dimension::Length spacing() const;

    SireUnits::Dimension::Length coulombCutoff() const;
    SireUnits::Dimension::Length ljCutoff() const;

    void mustNowRecalculateFromScratch();    

protected:
    void recalculateEnergy();

    void _pvt_added(quint32 groupid, const PartialMolecule &mol, 
                    const PropertyMap &map);

    void _pvt_removed(quint32 groupid, const PartialMolecule &mol);

    void _pvt_changed(quint32 groupid, const SireMol::Molecule &molecule, bool auto_commit);
    void _pvt_changed(quint32 groupid, const QList<SireMol::Molecule> &molecules, bool auto_commit);
    
    void _pvt_removedAll(quint32 groupid);

private:
    void rebuildGrid();

    typedef InterGroupCLJFF::Parameters CLJParameters;
    typedef InterGroupCLJFF::Molecule CLJMolecule;
    typedef InterGroupCLJFF::Molecules CLJMolecules;

    void calculateEnergy(const SireVol::CoordGroup &coords,
                         const CLJParameters::Array &params,
                         double &cnrg, double &ljnrg);

    class Vector4
    {
    public:
        Vector4() : x(0), y(0), z(0), q(0)
        {}
        
        Vector4(const SireMaths::Vector &v, double q);

        Vector4(const Vector4 &other)
            : x(other.x), y(other.y), z(other.z), q(other.q)
        {}
        
        ~Vector4()
        {}
        
        Vector4& operator=(const Vector4 &other)
        {
            x = other.x;
            y = other.y;
            z = other.z;
            q = other.q;
            return *this;
        }
        
        double x;
        double y;
        double z;
        double q;
    };

    static void appendTo(QVector<Vector4> &coords_and_charges,
                         const Vector *coords, 
                         const detail::CLJParameter *params,
                         int nats);
                         
    void addToGrid(const QVector<Vector4> &coords_and_charges);

    /** The AABox that describes the grid */
    SireVol::AABox gridbox;

    /** The distance from the atoms and the edge of the grid.
        This prevents the grid from being re-evaluated whenever
        the atoms move */
    double buffer_size;

    /** The grid spacing */
    double grid_spacing;

    /** The coulomb cutoff (cutoff applies from the
        center of the grid) */
    double coul_cutoff;

    /** The LJ cutoff (cutoff applies from the
        center of the grid) */
    double lj_cutoff;

    /** The number of grid points in x, y, and z */
    quint32 dimx, dimy, dimz;

    /** The grid of coulomb potentials */
    QVector<double> gridpot;
    
    /** The set of coordinates and parameters for the fixed atoms.
        These are atoms which exist only in this GridFF, thereby
        allowing them to be present in the energy expression without
        needing to be present in the system */
    QVector<SireMaths::Vector> fixedatoms_coords;
    QVector<detail::CLJParameter> fixedatoms_params;
    
    /** The set of coordinates and parameters of group 2 
        molecules that are within the LJ cutoff of the center of the grid */
    QVector<SireMaths::Vector> closemols_coords;
    QVector<detail::CLJParameter> closemols_params;
    
    /** The old energy of each molecule */
    QHash<SireMol::MolNum,CLJEnergy> oldnrgs;
};

}

Q_DECLARE_METATYPE( SireMM::GridFF )

SIRE_EXPOSE_CLASS( SireMM::GridFF )

SIRE_END_HEADER

#endif
