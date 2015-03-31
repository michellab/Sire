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

#include <QFile>
#include <QTextStream>

#include "cube.h"

#include "SireMol/molecule.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"

#include "SireVol/grid.h"

#include "SireUnits/units.h"
#include "SireUnits/convert.h"

#include "SireError/errors.h"

using namespace SireIO;
using namespace SireMol;
using namespace SireFF;
using namespace SireVol;
using namespace SireBase;
using namespace SireUnits;

Cube::Cube() : cutoff(10000)
{}

Cube::Cube(SireUnits::Dimension::MolarEnergy c) 
     : cutoff( convertTo(c.value(), kcal_per_mol) )
{}

Cube::Cube(const Cube &other) : cutoff(other.cutoff)
{}

Cube::~Cube()
{}

Cube& Cube::operator=(const Cube &other)
{
    cutoff = other.cutoff;
    return *this;
}

bool Cube::operator==(const Cube &other) const
{
    return cutoff == other.cutoff;
}

bool Cube::operator!=(const Cube &other) const
{
    return not Cube::operator==(other);
}

static void assertValidGrid(const PotentialTable &table)
{
    if (table.nGrids() != 1)
        throw SireError::incompatible_error( QObject::tr(
                "Sire can only write Cube data files for potential tables "
                "that contain a single Grid."), CODELOC );

    const Grid &grid = table.gridData()[0].grid();
    
    if (not grid.isA<RegularGrid>())
        throw SireError::incompatible_error( QObject::tr(
                "Sire can only write Cube data files when using a RegularGrid. "
                "It cannot work with the Grid %1.")
                    .arg(grid.what()), CODELOC );
}

static void writeHeader(QTextStream &ts, int nats,  const Grid &grid)
{
    const RegularGrid &rgrid = grid.asA<RegularGrid>();

    ts << "SIRE CUBE FILE.\n";
    ts << "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n";
    ts.setFieldWidth(4);
    ts << nats;
    ts.setFieldWidth(12);
    ts.setRealNumberPrecision(6);
    ts << convertTo(grid.minCoords().x(), bohr_radii) 
       << convertTo(grid.minCoords().y(), bohr_radii) 
       << convertTo(grid.minCoords().z(), bohr_radii) << "\n";

    Vector basis = convertTo(rgrid.gridSpacing().value(), bohr_radii) 
                                        * rgrid.basis().column0();
    ts.setFieldWidth(4);
    ts << rgrid.dimX();
    ts.setFieldWidth(12);
    ts << basis.x() << basis.y() << basis.z() << "\n";
    
    basis = convertTo(rgrid.gridSpacing(), bohr_radii) * rgrid.basis().column1();
    ts.setFieldWidth(4);
    ts << rgrid.dimY();
    ts.setFieldWidth(12);
    ts << basis.x() << basis.y() << basis.z() << "\n";
    
    basis = convertTo(rgrid.gridSpacing(), bohr_radii) * rgrid.basis().column2();
    ts.setFieldWidth(4);
    ts << rgrid.dimZ();
    ts.setFieldWidth(12);
    ts << basis.x() << basis.y() << basis.z() << "\n";
}

void Cube::write(const PotentialTable &table,
                 const QString &filename, const PropertyMap &map) const
{
    assertValidGrid(table);
        
    QFile f(filename);
    
    if (not f.open(QIODevice::WriteOnly))
        throw SireError::file_error(f, CODELOC);

    QTextStream ts(&f);
    
    writeHeader(ts, 0, table.gridData()[0].grid());
        
    f.close();
}

template<class T>
static void writeTable(QTextStream &ts, const MolPotentialTable &moltable,
                       const T &molgroup, const PropertyMap &map, double cutoff)
{
    MolNum molnum = moltable.molNum();
    
    if (not molgroup.contains(molnum))
    {
    }
    else
    {
        Molecule mol = molgroup[molnum].molecule();

        const AtomCoords &coords = mol.property( map["coordinates"] )
                                      .asA<AtomCoords>();
                                          
        const AtomElements &elements = mol.property( map["element"] )
                                          .asA<AtomElements>();
    
        if (moltable.selectedAll())
        {
            for (CGIdx i(0); i<moltable.nCutGroups(); ++i)
            {
                const MolarEnergy *nrg = moltable.constData(i);
                const Vector *c = coords.constData(i);
                const Element *e = elements.constData(i);
                
                for (int j=0; j<coords.nAtoms(i); ++j)
                {
                    ts.setFieldWidth(4);
                    ts << e[j].nProtons();
                    
                    ts.setFieldWidth(12);

                    double val = convertTo( nrg[j], kcal_per_mol );
                    if (val > cutoff)
                        val = cutoff;
                    else if (val < -cutoff)
                        val = -cutoff;
                    
                    ts << SireUnits::convertTo(nrg[j], kcal_per_mol)
                       << SireUnits::convertTo(c[j].x(), bohr_radii)
                       << SireUnits::convertTo(c[j].y(), bohr_radii)
                       << SireUnits::convertTo(c[j].z(), bohr_radii) << "\n";
                }
            }
        }
    }
}

static void writeGrid(QTextStream &ts, const GridPotentialTable &gridtable,
                      double cutoff)
{
    const RegularGrid &grid = gridtable.grid().asA<RegularGrid>();
    
    ts.setFieldWidth(12);
    
    int count = 0;

    const MolarEnergy *data = gridtable.constData();
    
    for (int i=0; i<grid.dimX(); ++i)
    {
        for (int j=0; j<grid.dimY(); ++j)
        {
            for (int k=0; k<grid.dimZ(); ++k)
            {
                double val = convertTo( data->value(), kcal_per_mol );
            
                if (val > cutoff)
                    ts << cutoff << " ";
                else if (val < -cutoff)
                    ts << -cutoff << " ";
                else
                    ts << val << " ";
                    
                ++data;
                                
                ++count;
                
                if (count == 5)
                {
                    count = 0;
                    ts << "\n";
                }
            }
        }
    }
}

static int countAtoms(const PotentialTable &table)
{
    int nats = 0;
    
    for (int i=0; i<table.nMolecules(); ++i)
    {
        const MolPotentialTable &moltable = table.moleculeData()[i];
        nats += moltable.nValues();
    }
    
    return nats;
}

void Cube::write(const PotentialTable &table, const MoleculeGroup &molgroup,
                 const QString &filename, const PropertyMap &map) const
{
    if (table.nMolecules() == 0)
    {
        this->write(table, filename, map);
        return;
    }

    //get the number of atoms
    int nats = countAtoms(table);

    if (nats == 0)
    {
        this->write(table, filename, map);
        return;
    }

    assertValidGrid(table);

    QFile f(filename);
    
    if (not f.open(QIODevice::WriteOnly))
        throw SireError::file_error(f, CODELOC);
        
    QTextStream ts(&f);
    
    writeHeader(ts, nats, table.gridData()[0].grid());
    
    for (int i=0; i<table.nMolecules(); ++i)
    {
        writeTable( ts, table.moleculeData()[i], molgroup, map, cutoff );
    }

    writeGrid( ts, table.gridData()[0], cutoff );
        
    f.close();
}

void Cube::write(const PotentialTable &table, const MolGroupsBase &molgroups,
                 const QString &filename, const PropertyMap &map) const
{
    assertValidGrid(table);

    if (table.nMolecules() == 0)
    {
        this->write(table, filename, map);
        return;
    }

    //get the number of atoms
    int nats = countAtoms(table);

    if (nats == 0)
    {
        this->write(table, filename, map);
        return;
    }

    assertValidGrid(table);

    QFile f(filename);
    
    if (not f.open(QIODevice::WriteOnly))
        throw SireError::file_error(f, CODELOC);
        
    QTextStream ts(&f);
    
    writeHeader(ts, nats, table.gridData()[0].grid());
    
    for (int i=0; i<table.nMolecules(); ++i)
    {
        writeTable( ts, table.moleculeData()[i], molgroups, map, cutoff );
    }

    writeGrid( ts, table.gridData()[0], cutoff );
        
    f.close();
}
