/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#include "SireIO/dcd.h"
#include "SireIO/fortranfile.h"

#include "SireSystem/system.h"

#include "SireIO/amberformat.h"

#include "SireMol/mgname.h"
#include "SireMol/molidx.h"
#include "SireMol/molecule.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/atomforces.h"
#include "SireMol/moleditor.h"
#include "SireMol/core.h"

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/generalunitproperty.h"
#include "SireBase/timeproperty.h"
#include "SireBase/booleanproperty.h"
#include "SireBase/getinstalldir.h"
#include "SireBase/unittest.h"

#include "SireIO/errors.h"

#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <QFile>
#include <QDataStream>
#include <QDebug>

using namespace SireIO;
using namespace SireIO::detail;
using namespace SireMaths;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireVol;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<DCD> r_dcd;
const RegisterParser<DCD> register_dcd;

QDataStream &operator<<(QDataStream &ds, const DCD &dcd)
{
    writeHeader(ds, r_dcd, 1);

    SharedDataStream sds(ds);

    sds << dcd.ttle
        << dcd.current_time
        << dcd.coords << dcd.vels << dcd.frcs
        << dcd.box_dims << dcd.box_angs
        << dcd.parse_warnings
        << static_cast<const MoleculeParser&>(dcd);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, DCD &dcd)
{
    VersionID v = readHeader(ds, r_dcd);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> dcd.ttle
            >> dcd.current_time
            >> dcd.coords >> dcd.vels >> dcd.frcs
            >> dcd.box_dims >> dcd.box_angs
            >> dcd.parse_warnings
            >> static_cast<MoleculeParser&>(dcd);
    }
    else
        throw version_error(v, "1", r_dcd, CODELOC);

    return ds;
}

static Vector cubic_angs(90,90,90);

/** Constructor */
DCD::DCD() : ConcreteProperty<DCD,MoleculeParser>(), current_time(0)
{}

/** Return the format name that is used to identify this file format within Sire */
QString DCD::formatName() const
{
    return "DCD";
}

/** Return the suffixes that DCD files will typically use */
QStringList DCD::formatSuffix() const
{
    static const QStringList suffixes = { "dcd" };
    return suffixes;
}

/** This is not a text file */
bool DCD::isTextFile() const
{
    return false;
}

/** Return a description of the file format */
QString DCD::formatDescription() const
{
    return QObject::tr("DCD coordinate/velocity binary trajectory files "
                       "based on charmm / namd / x-plor format.");
}

SireIO::detail::DCDFile::DCDFile()
               : timestep(0), istart(0),
                 nsavc(0), nfixed(0), natoms(0),
                 nframes(0), first_frame_line(0),
                 CHARMM_FORMAT(false),
                 HAS_EXTRA_BLOCK(false),
                 HAS_FOUR_DIMS(false)
{}

SireIO::detail::DCDFile::~DCDFile()
{}

void SireIO::detail::DCDFile::readHeader(FortranFile &file)
{
    auto line = file[0];

    auto typ = line.readChar(4);

    if (typ != "CORD")
        throw SireIO::parse_error(QObject::tr(
            "This does not look like a DCD file, because it does "
            "not start with 'CORD'. Starts with %1.").arg(typ), CODELOC);

    auto ints = line.readInt32(9);

    nframes = ints[0];
    istart = ints[1];
    nsavc = ints[2];
    nfixed = ints[8];

    // now need to check the value at buffer[80] as, if this is non-zero,
    // then this is a CHARMM format DCD file
    CHARMM_FORMAT = line.readInt32At(80) != 0;

    // the value at buffer[44] says if there is an extra block
    HAS_EXTRA_BLOCK = line.readInt32At(44) != 0;

    // the value at buffer[48] says whether or not we have four dimensions
    HAS_FOUR_DIMS = line.readInt32At(48) != 0;

    timestep = 0;

    if (CHARMM_FORMAT)
    {
        timestep = line.readFloat32At(40);
    }
    else
    {
        timestep = line.readFloat64At(40);
    }

    line = file[1];

    int ntitle = line.readInt32(1)[0];

    for (int i=0; i<ntitle; ++i)
    {
        QString t = line.readChar(32).simplified().replace('\0', "");

        if (not t.isEmpty())
            title.append(t);
    }

    line = file[2];

    natoms = line.readInt32(1)[0];

    int linenum = 3;

    if (nfixed != 0)
    {
        line = file[linenum];
        linenum += 1;

        fixed_atoms = line.readInt32(nfixed);
    }

    first_frame_line = linenum;

    if (nfixed != 0)
    {
        // we have to read in the first set of coordinates, as these
        // hold the fixed atoms as well as the movable atoms
        if (CHARMM_FORMAT and HAS_EXTRA_BLOCK)
        {
            line = file[linenum];
            linenum += 1;

            line.readFloat64(6);
        }

        line = file[linenum];
        linenum += 1;
        auto x = line.readFloat32(natoms);

        line = file[linenum];
        linenum += 1;
        auto y = line.readFloat32(natoms);

        line = file[linenum];
        linenum += 1;
        auto z = line.readFloat32(natoms);

        first_frame = QVector<Vector>(natoms);

        for (int i=0; i<natoms; ++i)
        {
            first_frame[i] = Vector(x[i], y[i], z[i]);
        }
    }

    // now sanity check the rest of the file
    int num_lines_per_frame = 3;

    if (HAS_FOUR_DIMS)
    {
        num_lines_per_frame += 1;
    }

    if (CHARMM_FORMAT and HAS_EXTRA_BLOCK)
    {
        num_lines_per_frame += 1;
    }

    if (file.nRecords() != first_frame_line + (num_lines_per_frame * nframes))
    {
        throw SireIO::parse_error(QObject::tr(
            "Wrong number of records in the DCD file. Expect to have %1 "
            "for %2 frames, but actually have %3.")
                .arg(first_frame_line + (num_lines_per_frame*nframes))
                .arg(nframes).arg(file.nRecords()), CODELOC);
    }
}

QVector<double> SireIO::detail::DCDFile::readUnitCell(FortranFile &file, int frame)
{
    if (frame < 0 or frame >= nframes)
    {
        throw SireIO::parse_error(QObject::tr(
            "Trying to access an invalid frame (%1) from a DCD with %2 frames.")
                .arg(frame).arg(nframes), CODELOC);
    }

    // get the line number for this frame
    int num_lines_per_frame = 3;

    if (HAS_FOUR_DIMS)
        num_lines_per_frame += 1;

    if (CHARMM_FORMAT and HAS_EXTRA_BLOCK)
    {
        num_lines_per_frame += 1;

        int linenum = first_frame_line + (frame*num_lines_per_frame);

        auto line = file[linenum];

        return line.readFloat64(6);
    }

    return QVector<double>();
}

QVector<Vector> SireIO::detail::DCDFile::readCoordinates(FortranFile &file, int frame)
{
    if (frame < 0 or frame >= nframes)
    {
        throw SireIO::parse_error(QObject::tr(
            "Trying to access an invalid frame (%1) from a DCD with %2 frames.")
                .arg(frame).arg(nframes), CODELOC);
    }

    // get the line number for this frame
    int num_lines_per_frame = 3;
    int skip_unitcell = 0;

    if (HAS_FOUR_DIMS)
        num_lines_per_frame += 1;

    if (CHARMM_FORMAT and HAS_EXTRA_BLOCK)
    {
        num_lines_per_frame += 1;
        skip_unitcell = 1;
    }

    int linenum = first_frame_line + (frame*num_lines_per_frame) + skip_unitcell;

    if (nfixed == 0)
    {
        // read in the x, y, and z data
        auto line = file[linenum];
        auto x = line.readFloat32(natoms);
        line = file[linenum+1];
        auto y = line.readFloat32(natoms);
        line = file[linenum+2];
        auto z = line.readFloat32(natoms);

        QVector<Vector> coords(natoms);

        for (int i=0; i<natoms; ++i)
        {
            coords[i] = Vector(x[i], y[i], z[i]);
        }

        return coords;
    }
    else if (frame == 0)
    {
        return first_frame;
    }
    else
    {
        // read the coordinates and map them into a copy of first_frame
        QVector<Vector> frame = first_frame;

        return frame;
    }
}

/** Parse the data contained in the lines - this clears any pre-existing
    data in this object */
void DCD::parse(const QString &filename, const PropertyMap &map)
{
    FortranFile file(filename);

    SireIO::detail::DCDFile dcd;
    dcd.readHeader(file);

    ttle = dcd.title.join(" ");
    current_time = dcd.istart * dcd.timestep;

    coords = dcd.readCoordinates(file, 0);
    auto unitcell = dcd.readUnitCell(file, 0);

    //need to convert unit cell...

    //set the score, and save the warnings
    double score = 100.0 / (parse_warnings.count()+1);
    this->setScore(score);
}

/** Construct by parsing the passed file */
DCD::DCD(const QString &filename, const PropertyMap &map)
    : ConcreteProperty<DCD,MoleculeParser>(map), current_time(0)
{
    this->parse(filename, map);
}

/** Construct by parsing the data in the passed text lines */
DCD::DCD(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<DCD,MoleculeParser>(lines, map), current_time(0)
{
    throw SireIO::parse_error( QObject::tr(
            "You cannot create a binary DCD file from a set of text lines!"),
                CODELOC );
}

static QVector<Vector> getCoordinates(const Molecule &mol, const PropertyName &coords_property)
{
    if (not mol.hasProperty(coords_property))
    {
        return QVector<Vector>();
    }

    QVector<Vector> coords( mol.nAtoms() );

    const auto molcoords = mol.property(coords_property).asA<AtomCoords>();

    const auto molinfo = mol.info();

    for (int i=0; i<mol.nAtoms(); ++i)
    {
        //coords are already in angstroms :-)
        coords[i] = molcoords.at( molinfo.cgAtomIdx( AtomIdx(i) ) );
    }

    return coords;
}

static QVector<Vector> getVelocities(const Molecule &mol, const PropertyName &vels_property)
{
    if (not mol.hasProperty(vels_property))
    {
        return QVector<Vector>();
    }

    try
    {
        const auto molvels = mol.property(vels_property).asA<AtomVelocities>();
        const auto molinfo = mol.info();

        QVector<Vector> vels( mol.nAtoms() );

        const double units = 1.0 / (angstrom / (20.455*picosecond)).value();

        for (int i=0; i<mol.nAtoms(); ++i)
        {
            const auto atomvels = molvels.at( molinfo.cgAtomIdx( AtomIdx(i) ) );

            //need to convert the velocities into units of Angstroms / 20.455 picoseconds
            vels[i] = Vector( atomvels.x().value() * units,
                              atomvels.y().value() * units,
                              atomvels.z().value() * units );
        }

        return vels;
    }
    catch(...)
    {
        return QVector<Vector>();
    }
}

static QVector<Vector> getForces(const Molecule &mol, const PropertyName &forces_property)
{
    if (not mol.hasProperty(forces_property))
    {
        return QVector<Vector>();
    }

    try
    {
        const auto molforces = mol.property(forces_property).asA<AtomForces>();
        const auto molinfo = mol.info();

        QVector<Vector> forces( mol.nAtoms() );

        const double units = 1.0 / ( atomic_mass_constant * angstrom /
                                     (picosecond*picosecond) ).value();

        for (int i=0; i<mol.nAtoms(); ++i)
        {
            const auto atomforces = molforces.at( molinfo.cgAtomIdx( AtomIdx(i) ) );

            //need to convert the velocities into units of amu Angstroms / picosecond^2
            forces[i] = Vector( atomforces.x().value() * units,
                                atomforces.y().value() * units,
                                atomforces.z().value() * units );
        }

        return forces;
    }
    catch(...)
    {
        return QVector<Vector>();
    }
}

static bool hasData(const QVector< QVector<Vector> > &array)
{
    for (int i=0; i<array.count(); ++i)
    {
        if (not array[i].isEmpty())
            return true;
    }

    return false;
}

/** Construct by extracting the necessary data from the passed System */
DCD::DCD(const System &system, const PropertyMap &map)
    : ConcreteProperty<DCD,MoleculeParser>(), current_time(0)
{
    //get the MolNums of each molecule in the System - this returns the
    //numbers in MolIdx order
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    if (molnums.isEmpty())
    {
        //no molecules in the system
        this->operator=(DCD());
        return;
    }

    //get the coordinates (and velocities if available) for each molecule in the system
    {
        QVector< QVector<Vector> > all_coords(molnums.count());
        QVector< QVector<Vector> > all_vels(molnums.count());
        QVector< QVector<Vector> > all_forces(molnums.count());

        const auto coords_property = map["coordinates"];
        const auto vels_property = map["velocity"];
        const auto forces_property = map["force"];

        if (usesParallel())
        {
            tbb::parallel_for( tbb::blocked_range<int>(0,molnums.count()),
                               [&](const tbb::blocked_range<int> r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto mol = system[molnums[i]].molecule();

                    tbb::parallel_invoke(
                       [&](){ all_coords[i] = ::getCoordinates(mol, coords_property); },
                       [&](){ all_vels[i] = ::getVelocities(mol, vels_property); },
                       [&](){ all_forces[i] = ::getForces(mol, forces_property); }
                                        );
                }
            });
        }
        else
        {
            for (int i=0; i<molnums.count(); ++i)
            {
                const auto mol = system[molnums[i]].molecule();

                all_coords[i] = ::getCoordinates(mol, coords_property);
                all_vels[i] = ::getVelocities(mol, vels_property);
                all_forces[i] = ::getForces(mol, forces_property);
            }
        }

        coords.clear();
        vels.clear();
        frcs.clear();

        if (::hasData(all_coords))
        {
            coords = collapse(all_coords);
        }

        if (::hasData(all_vels))
        {
            vels = collapse(all_vels);
        }

        if (::hasData(all_forces))
        {
            frcs = collapse(all_forces);
        }
    }

    //extract the space of the system
    box_dims = Vector(0);
    box_angs = Vector(0);

    try
    {
        if (system.property(map["space"]).isA<PeriodicBox>())
        {
            box_dims = system.property(map["space"]).asA<PeriodicBox>().dimensions();
            box_angs = Vector(90,90,90);
        }
        else if (system.property(map["space"]).isA<TriclinicBox>())
        {
            Vector v0 = system.property(map["space"]).asA<TriclinicBox>().vector0();
            Vector v1 = system.property(map["space"]).asA<TriclinicBox>().vector1();
            Vector v2 = system.property(map["space"]).asA<TriclinicBox>().vector2();

            // Store the box magnitudes.
            box_dims = Vector(v0.magnitude(), v1.magnitude(), v2.magnitude());

            double alpha = system.property(map["space"]).asA<TriclinicBox>().alpha();
            double beta  = system.property(map["space"]).asA<TriclinicBox>().alpha();
            double gamma = system.property(map["space"]).asA<TriclinicBox>().gamma();

            // Store the angles between the box vectors.
            box_angs = Vector(alpha, beta, gamma);
        }
    }
    catch(...)
    {}

    //extract the current time for the system
    try
    {
        const Property &prop = system.property( map["time"] );

        Time time;

        if (prop.isA<TimeProperty>())
            time = prop.asA<TimeProperty>().value();
        else
            time = prop.asA<GeneralUnitProperty>();

        current_time = time.to(picosecond);
    }
    catch(...)
    {}

    //extract the title for the system
    ttle = system.name().value();
}

/** Copy constructor */
DCD::DCD(const DCD &other)
    : ConcreteProperty<DCD,MoleculeParser>(other),
      ttle(other.ttle), current_time(other.current_time),
      coords(other.coords), vels(other.vels), frcs(other.frcs),
      box_dims(other.box_dims), box_angs(other.box_angs),
      parse_warnings(other.parse_warnings)
{}

/** Destructor */
DCD::~DCD()
{}

DCD& DCD::operator=(const DCD &other)
{
    if (this != &other)
    {
        ttle = other.ttle;
        current_time = other.current_time;
        coords = other.coords;
        vels = other.vels;
        frcs = other.frcs;
        box_dims = other.box_dims;
        box_angs = other.box_angs;
        parse_warnings = other.parse_warnings;

        MoleculeParser::operator=(other);
    }

    return *this;
}

bool DCD::operator==(const DCD &other) const
{
    return MoleculeParser::operator==(other);
}

bool DCD::operator!=(const DCD &other) const
{
    return not DCD::operator==(other);
}

const char* DCD::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DCD>() );
}

const char* DCD::what() const
{
    return DCD::typeName();
}

QString DCD::toString() const
{
    if (nAtoms() == 0)
    {
        return QObject::tr("DCD( nAtoms() = 0 )");
    }
    else
    {
        return QObject::tr("DCD( title() = %1, nAtoms() = %2, "
                "hasCoordinates() = %3, hasVelocities() = %4, "
                "hasForces() = %5 )")
                .arg(title()).arg(nAtoms())
                .arg(hasCoordinates()).arg(hasVelocities()).arg(hasForces());
    }
}

/** Parse from the passed file */
DCD DCD::parse(const QString &filename)
{
    return DCD(filename);
}

/** Internal function used to add the data from this parser into the passed System */
void DCD::addToSystem(System &system, const PropertyMap &map) const
{
    const bool has_coords = hasCoordinates();
    const bool has_vels = hasVelocities();
    const bool has_forces = hasForces();

    if (not (has_coords or has_vels or has_forces))
    {
        //nothing to add...
        return;
    }

    //first, we are going to work with the group of all molecules, which should
    //be called "all". We have to assume that the molecules are ordered in "all"
    //in the same order as they are in this restart file, with the data
    //in MolIdx/AtomIdx order (this should be the default for all parsers!)
    MoleculeGroup allmols = system[MGName("all")];

    const int nmols = allmols.nMolecules();

    QVector<int> atom_pointers(nmols+1, -1);

    int natoms = 0;

    for (int i=0; i<nmols; ++i)
    {
        atom_pointers[i] = natoms;
        const int nats = allmols[MolIdx(i)].data().info().nAtoms();
        natoms += nats;
    }

    atom_pointers[nmols] = natoms;

    if (natoms != this->nAtoms())
        throw SireIO::parse_error( QObject::tr(
                "Incompatibility between the files, as this DCD file contains data "
                "for %1 atom(s), while the other file(s) have created a system with "
                "%2 atom(s)").arg(this->nAtoms()).arg(natoms), CODELOC );

    //next, copy the coordinates and optionally the velocities into the molecules
    QVector<Molecule> mols(nmols);
    Molecule *mols_array = mols.data();

    const PropertyName coords_property = map["coordinates"];
    const PropertyName vels_property = map["velocity"];
    const PropertyName forces_property = map["force"];

    const Vector *coords_array = this->coordinates().constData();
    const Vector *vels_array = this->velocities().constData();
    const Vector *forces_array = this->forces().constData();

    auto add_moldata = [&](int i)
    {
        const int atom_start_idx = atom_pointers.constData()[i];
        auto mol = allmols[MolIdx(i)].molecule();
        const auto molinfo = mol.data().info();

        auto moleditor = mol.edit();

        if (has_coords)
        {
            auto coords = QVector< QVector<Vector> >(molinfo.nCutGroups());

            for (int j=0; j<molinfo.nCutGroups(); ++j)
            {
                coords[j] = QVector<Vector>(molinfo.nAtoms(CGIdx(j)));
            }

            for (int j=0; j<mol.nAtoms(); ++j)
            {
                auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(j));

                const int atom_idx = atom_start_idx + j;

                coords[cgatomidx.cutGroup()][cgatomidx.atom()] = coords_array[atom_idx];
            }

            moleditor.setProperty( coords_property,AtomCoords(CoordGroupArray(coords)) );
        }

        if (has_vels)
        {
            auto vels = AtomVelocities(molinfo);

            for (int j=0; j<mol.nAtoms(); ++j)
            {
                auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(j));

                const int atom_idx = atom_start_idx + j;

                //velocity is Angstroms per 1/20.455 ps
                const auto vel_unit = (1.0 / 20.455) * angstrom / picosecond;

                const Vector &vel = vels_array[atom_idx];
                vels.set(cgatomidx, Velocity3D(vel.x() * vel_unit,
                                               vel.y() * vel_unit,
                                               vel.z() * vel_unit));
            }

            moleditor.setProperty( vels_property, vels );
        }

        if (has_forces)
        {
            auto forces = AtomForces(molinfo);

            for (int j=0; j<mol.nAtoms(); ++j)
            {
                auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(j));

                const int atom_idx = atom_start_idx + j;

                //force is amu Angstroms per ps^2
                const auto force_unit = atomic_mass_constant * angstrom / (picosecond*picosecond);

                const Vector &f = forces_array[atom_idx];
                forces.set(cgatomidx, Force3D(f.x() * force_unit,
                                              f.y() * force_unit,
                                              f.z() * force_unit));
            }

            moleditor.setProperty( forces_property, forces );
        }

        mols_array[i] = moleditor.commit();
    };

    if (usesParallel())
    {
        tbb::parallel_for( tbb::blocked_range<int>(0,nmols),
                           [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                add_moldata(i);
            }
        });
    }
    else
    {
        for (int i=0; i<nmols; ++i)
        {
            add_moldata(i);
        }
    }

    system.update( Molecules(mols) );

    PropertyName space_property = map["space"];
    if (space_property.hasValue())
    {
        system.setProperty("space", space_property.value());
    }
    else if ((not box_dims.isZero()) and space_property.hasSource())
    {
        // PeriodicBox.
        if (box_angs == cubic_angs)
        {
            system.setProperty( space_property.source(),
                                SireVol::PeriodicBox(box_dims) );
        }
        // TriclinicBox.
        else
        {
            system.setProperty( space_property.source(),
                                SireVol::TriclinicBox(box_dims.x(),
                                                      box_dims.y(),
                                                      box_dims.z(),
                                                      box_angs.x()*degrees,
                                                      box_angs.y()*degrees,
                                                      box_angs.z()*degrees) );
        }
    }

    //update the System fileformat property to record that it includes
    //data from this file format
    QString fileformat = this->formatName();

    PropertyName fileformat_property = map["fileformat"];

    try
    {
        QString last_format = system.property(fileformat_property).asA<StringProperty>().value();
        fileformat = QString("%1,%2").arg(last_format,fileformat);
    }
    catch(...)
    {}

    if (fileformat_property.hasSource())
    {
        system.setProperty(fileformat_property.source(), StringProperty(fileformat));
    }
    else
    {
        system.setProperty("fileformat", StringProperty(fileformat));
    }

    PropertyName time_property = map["time"];

    if (time_property.hasSource())
    {
        system.setProperty(time_property.source(), GeneralUnitProperty(current_time*picosecond));
    }
    else
    {
        system.setProperty("time", GeneralUnitProperty(current_time*picosecond));
    }
}

/** Return the title of the file */
QString DCD::title() const
{
    return ttle;
}

/** Return the current time of the simulation from which this restart
    file was written. Returns 0 if there is no time set. If there are
    multiple frames, then the time of the first frame is returned */
double DCD::time() const
{
    return current_time;
}

/** Return the number of atoms whose data are contained in this DCD file */
int DCD::nAtoms() const
{
    if (not coords.isEmpty())
        return coords.count();
    else if (not vels.isEmpty())
        return vels.count();
    else if (not frcs.isEmpty())
        return frcs.count();
    else
        return 0;
}

/** Return whether or not this DCD file provides coordinates */
bool DCD::hasCoordinates() const
{
    return not coords.isEmpty();
}

/** Return whether or not this DCD file also provides velocities */
bool DCD::hasVelocities() const
{
    return not vels.isEmpty();
}

/** Return whether or not this DCD file also provides forces */
bool DCD::hasForces() const
{
    return not frcs.isEmpty();
}

/** Return the parsed coordinate data. */
QVector<SireMaths::Vector> DCD::coordinates() const
{
    return coords;
}

/** Return the parsed velocity */
QVector<SireMaths::Vector> DCD::velocities() const
{
    return vels;
}

/** Return the parsed force data. */
QVector<SireMaths::Vector> DCD::forces() const
{
    return frcs;
}

/** Return the parsed box dimensions, or Vector(0) if there is no space information */
SireMaths::Vector DCD::boxDimensions() const
{
    return box_dims;
}

/** Return the parsed box angles, or Vector(0) if there is no space information. */
SireMaths::Vector DCD::boxAngles() const
{
    return box_angs;
}

/** Return any warnings that were triggered during parsing */
QStringList DCD::warnings() const
{
    return parse_warnings;
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr DCD::construct(const QString &filename,
                                 const PropertyMap &map) const
{
    //don't construct from a pointer as it could leak
    return MoleculeParserPtr( DCD(filename,map) );
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr DCD::construct(const QStringList &lines,
                                 const PropertyMap &map) const
{
    //don't construct from a pointer as it could leak
    return MoleculeParserPtr( DCD(lines,map) );
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr DCD::construct(const SireSystem::System &system,
                                 const PropertyMap &map) const
{
    //don't construct from a pointer as it could leak
    return MoleculeParserPtr( DCD(system,map) );
}

/** Write this DCD to a file called 'filename'. This will write out
    the data in this object to the DCD format */
void DCD::writeToFile(const QString &filename) const
{
    // (write all types to files, as needed)
}
