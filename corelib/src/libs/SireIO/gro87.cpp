/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#include "SireIO/gro87.h"
#include "SireIO/amberformat.h"

#include "SireSystem/system.h"

#include "SireMol/molecule.h"
#include "SireMol/moleculeinfo.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/mgname.h"
#include "SireMol/molidx.h"
#include "SireMol/moleditor.h"
#include "SireMol/core.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireBase/parallel.h"
#include "SireBase/timeproperty.h"
#include "SireBase/numberproperty.h"
#include "SireBase/stringproperty.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QRegularExpression>
#include <QDebug>
#include <QElapsedTimer>

using namespace SireIO;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireVol;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

const RegisterParser<Gro87> register_gro87;
static const RegisterMetaType<Gro87> r_gro87;

QDataStream &operator<<(QDataStream &ds, const Gro87 &gro87)
{
    writeHeader(ds, r_gro87, 1);

    SharedDataStream sds(ds);

    sds << gro87.ttle << gro87.current_time << gro87.coords
        << gro87.vels << gro87.box_v1 << gro87.box_v2 << gro87.box_v3
        << gro87.resnums << gro87.resnams << gro87.atmnums
        << gro87.atmnams << gro87.parse_warnings
        << static_cast<const MoleculeParser&>(gro87);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Gro87 &gro87)
{
    VersionID v = readHeader(ds, r_gro87);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> gro87.ttle >> gro87.current_time >> gro87.coords
            >> gro87.vels >> gro87.box_v1 >> gro87.box_v2 >> gro87.box_v3
            >> gro87.resnums >> gro87.resnams >> gro87.atmnums
            >> gro87.atmnams >> gro87.parse_warnings
            >> static_cast<MoleculeParser&>(gro87);
    }
    else
        throw version_error(v, "1", r_gro87, CODELOC);

    return ds;
}

/** Constructor */
Gro87::Gro87() : ConcreteProperty<Gro87,MoleculeParser>()
{}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
Gro87::Gro87(const QString &filename, const PropertyMap &map)
      : ConcreteProperty<Gro87,MoleculeParser>(filename,map)
{
    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
Gro87::Gro87(const QStringList &lines, const PropertyMap &map)
      : ConcreteProperty<Gro87,MoleculeParser>(lines,map)
{
    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

static QVector<QString> toLines(const QVector<QString> &atmnams,
                                const QVector<QString> &resnams,
                                const QVector<qint64> &resnums,
                                const QVector<Vector> &coords,
                                const QVector<Vector> &vels,
                                int precision,
                                bool uses_parallel, QStringList *errors)
{
    //do any of the molecules have velocities?
    const bool has_velocities = not vels.isEmpty();

    //calculate the total number of atoms
    const int nats = coords.count();

    if (nats == 0)
    {
        //there are no atoms!
        return QVector<QString>();
    }

    //reserve space for all of the lines
    QVector<QString> lines( nats );
    auto lines_data = lines.data();

    auto write_line = [&](int iatm)
    {
        //the atom number is iatm+1
        int atmnum = iatm+1;

        //however, it cannot be larger than 99999, so it should be capped at this value
        if (atmnum > 99999)
            atmnum = 99999;

        int resnum = resnums.constData()[iatm];

        //similarly, the residue number cannot be greater than 99999
        if (resnum > 99999)
            resnum = 99999;

        const auto resnam = resnams.constData()[iatm];
        const auto atmnam = atmnams.constData()[iatm];
        Vector coord = 0.1 * coords.constData()[iatm];  // convert to nanometers

        if (has_velocities)
        {
            Vector vel = 0.1 * vels.constData()[iatm]; // convert to nanometers per picosecond

            lines_data[iatm] = QString("%1%2%3%4%5%6%7%8%9%10")
                                    .arg(resnum, 5)
                                    .arg(resnam.left(5), -5)
                                    .arg(atmnam.left(5), 5)
                                    .arg(atmnum, 5)
                                    .arg(coord.x(), precision+5, 'f', precision)
                                    .arg(coord.y(), precision+5, 'f', precision)
                                    .arg(coord.z(), precision+5, 'f', precision)
                                    .arg(vel.x(), precision+5, 'f', precision+1)
                                    .arg(vel.y(), precision+5, 'f', precision+1)
                                    .arg(vel.z(), precision+5, 'f', precision+1);
        }
        else
        {
            lines_data[iatm] = QString("%1%2%3%4%5%6%7")
                                    .arg(resnum, 5)
                                    .arg(resnam.left(5), -5)
                                    .arg(atmnam.left(5), 5)
                                    .arg(atmnum, 5)
                                    .arg(coord.x(), precision+5, 'f', precision)
                                    .arg(coord.y(), precision+5, 'f', precision)
                                    .arg(coord.z(), precision+5, 'f', precision);
        }
    };

    if (uses_parallel)
    {
        tbb::parallel_for( tbb::blocked_range<int>(0, nats),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                write_line(i);
            }
        });
    }
    else
    {
        for (int i=0; i<nats; ++i)
        {
            write_line(i);
        }
    }

    return lines;
}

static std::tuple< QVector<QString>, QVector<qint64>, QVector<QString> >
getIDs(const MoleculeInfo &mol)
{
    const int nats = mol.nAtoms();

    //get the atoms out in AtomIdx order
    if (nats == 0)
    {
        return std::tuple< QVector<QString>, QVector<qint64>, QVector<QString> >();
    }

    QVector<QString> resnams(nats);
    QVector<qint64> resnums(nats);
    QVector<QString> atmnams(nats);

    for (int i=0; i<nats; ++i)
    {
        const AtomIdx idx(i);

        const auto residx = mol.parentResidue(idx);

        atmnams[i] = mol.name(idx).value();
        resnams[i] = mol.name(residx).value();
        resnums[i] = mol.number(residx).value();
    }

    return std::make_tuple(resnams, resnums, atmnams);
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

        const double units = 1.0 / (angstrom/picosecond).value();

        for (int i=0; i<mol.nAtoms(); ++i)
        {
            const auto atomvels = molvels.at( molinfo.cgAtomIdx( AtomIdx(i) ) );

            //need to convert the velocities into units of Angstroms / picoseconds
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

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
Gro87::Gro87(const SireSystem::System &system, const PropertyMap &map)
      : ConcreteProperty<Gro87,MoleculeParser>(map)
{
    //get the MolNums of each molecule in the System - this returns the
    //numbers in MolIdx order
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    if (molnums.isEmpty())
    {
        //no molecules in the system
        this->operator=(Gro87());
        return;
    }

    //get the names, numbers coordinates (and velocities if available)
    // for each molecule in the system
    QVector< QVector<Vector> > all_coords(molnums.count());
    QVector< QVector<Vector> > all_vels(molnums.count());
    QVector< QVector<QString> > all_resnams(molnums.count());
    QVector< QVector<qint64> > all_resnums(molnums.count());
    QVector< QVector<QString> > all_atmnams(molnums.count());

    const auto vels_property = map["velocity"];

    if (usesParallel())
    {
        tbb::parallel_for( tbb::blocked_range<int>(0,molnums.count()),
                           [&](const tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                const auto mol = system[molnums[i]].molecule();

                auto coords_property = map["coordinates"];

                bool is_perturbable = false;
                try
                {
                    is_perturbable = mol.property("is_perturbable").asABoolean();
                }
                catch (...)
                {}

                if (is_perturbable)
                {
                    // Allow the user to override the default.
                    if ((map["coordinates"] == "coordinates0") or
                        (map["coordinates"] == "coordinates1"))
                    {
                            coords_property = map["coordinates"];
                    }
                    else
                    {
                        // Default to lambda = 0.
                        if (mol.hasProperty("coordinates0"))
                            coords_property = "coordinates0";
                        else if (mol.hasProperty("coordinates1"))
                            coords_property = "coordinates1";
                        else
                            throw SireError::incompatible_error(QObject::tr("Missing coordinates for perturbable molecule!"));
                    }
                }

                tbb::parallel_invoke(
                    [&]()
                    {
                        const auto ids = ::getIDs(mol.info());
                        all_resnams[i] = std::get<0>(ids);
                        all_resnums[i] = std::get<1>(ids);
                        all_atmnams[i] = std::get<2>(ids);
                    },
                    [&](){ all_coords[i] = ::getCoordinates(mol, coords_property); },
                    [&](){ all_vels[i] = ::getVelocities(mol, vels_property); }
                                    );
            }
        });
    }
    else
    {
        for (int i=0; i<molnums.count(); ++i)
        {
            const auto mol = system[molnums[i]].molecule();

            auto coords_property = map["coordinates"];

            bool is_perturbable = false;
            try
            {
                is_perturbable = mol.property("is_perturbable").asABoolean();
            }
            catch (...)
            {}

            if (is_perturbable)
            {
                // Allow the user to override the default.
                if ((map["coordinates"] == "coordinates0") or
                    (map["coordinates"] == "coordinates1"))
                {
                    coords_property = map["coordinates"];
                }
                else
                {
                    // Default to lambda = 0.
                    if (mol.hasProperty("coordinates0"))
                        coords_property = "coordinates0";
                    else if (mol.hasProperty("coordinates1"))
                        coords_property = "coordinates1";
                    else
                        throw SireError::incompatible_error(QObject::tr("Missing coordinates for perturbable molecule!"));
                }
            }

            const auto ids = ::getIDs(mol.info());

            all_resnams[i] = std::get<0>(ids);
            all_resnums[i] = std::get<1>(ids);
            all_atmnams[i] = std::get<2>(ids);

            all_coords[i] = ::getCoordinates(mol, coords_property);
            all_vels[i] = ::getVelocities(mol, vels_property);
        }
    }

    QStringList errors;

    //extract the space of the system
    SpacePtr space;

    try
    {
        space = system.property( map["space"] ).asA<Space>();
    }
    catch(...)
    {}

    //extract the current time for the system
    Time time(-1);

    try
    {
        time = system.property( map["time"] ).asA<TimeProperty>().value();
    }
    catch(...)
    {}

    //what precision should be used - this can be set by the user
    int precision = 6;
    try
    {
        precision = map["precision"].value().asA<NumberProperty>().value();

        if (precision < 1)
        {
            precision = 1;
        }
        else if (precision > 16)
        {
            precision = 16;
        }
    }
    catch(...)
    {}

    //now convert these into text lines that can be written as the file
    QVector<QString> lines = ::toLines(SireIO::detail::collapse(all_atmnams),
                                       SireIO::detail::collapse(all_resnams),
                                       SireIO::detail::collapse(all_resnums),
                                       SireIO::detail::collapse(all_coords),
                                       SireIO::detail::collapse(all_vels),
                                       precision, this->usesParallel(),
                                       &errors);

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error( QObject::tr(
            "Errors converting the system to a Gromacs Gro87 format...\n%1")
                .arg(errors.join("\n")), CODELOC );
    }

    //we don't need the coords and vels data any more, so free the memory
    all_coords.clear();
    all_vels.clear();
    all_resnams.clear();
    all_resnums.clear();
    all_atmnams.clear();

    //add the title, number of atoms and time to the top of the lines
    //(the number of atoms is the current number of lines)
    const int natoms = lines.count();
    lines.prepend( QString("%1").arg(natoms, 5) );

    if (time.value() >= 0)
    {
        lines.prepend( QString("%1, t= %2").arg(system.name().value())
                                           .arg(time.to(picosecond), 8, 'f', 7) );
    }
    else
    {
        lines.prepend(system.name().value());
    }

    //finally add on the box dimensions
    if (space.read().isA<PeriodicBox>())
    {
        Vector dims = space.read().asA<PeriodicBox>().dimensions();

        lines += QString(" %1 %2 %3").arg(0.1 * dims.x(), 9, 'f', 5)
                                     .arg(0.1 * dims.y(), 9, 'f', 5)
                                     .arg(0.1 * dims.z(), 9, 'f', 5);
    }
    else if (space.read().isA<TriclinicBox>())
    {
        Matrix cell = space.read().asA<TriclinicBox>().cellMatrix();

        lines += QString(" %1 %2 %3 %4 %5 %6 %7 %8 %9")
                    .arg(0.1 * cell.column0().x(), 9, 'f', 5)  // XX
                    .arg(0.1 * cell.column1().y(), 9, 'f', 5)  // YY
                    .arg(0.1 * cell.column2().z(), 9, 'f', 5)  // ZZ
                    .arg(0.1 * cell.column0().y(), 9, 'f', 5)  // XY
                    .arg(0.1 * cell.column0().z(), 9, 'f', 5)  // XZ
                    .arg(0.1 * cell.column1().x(), 9, 'f', 5)  // YX
                    .arg(0.1 * cell.column1().z(), 9, 'f', 5)  // YZ
                    .arg(0.1 * cell.column2().x(), 9, 'f', 5)  // ZX
                    .arg(0.1 * cell.column2().y(), 9, 'f', 5); // ZY
    }
    else
    {
        //we have to provide some space for Gro87. Supply a null vector
        lines += QString("   0.00000   0.00000   0.00000");
    }

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error( QObject::tr(
            "Errors converting the system to a Gromacs Gro87 format...\n%1")
                .arg(errors.join("\n")), CODELOC );
    }

    //now generate this object by re-reading these lines
    Gro87 parsed(lines.toList());

    this->operator=(parsed);
}

/** Copy constructor */
Gro87::Gro87(const Gro87 &other)
      : ConcreteProperty<Gro87,MoleculeParser>(other),
        ttle(other.ttle), current_time(other.current_time),
        coords(other.coords), vels(other.vels),
        box_v1(other.box_v1), box_v2(other.box_v2), box_v3(other.box_v3),
        resnums(other.resnums), resnams(other.resnams),
        atmnams(other.atmnams), atmnums(other.atmnums), parse_warnings(other.parse_warnings)
{}

/** Destructor */
Gro87::~Gro87()
{}

/** Copy assignment operator */
Gro87& Gro87::operator=(const Gro87 &other)
{
    if (this != &other)
    {
        ttle = other.ttle;
        current_time = other.current_time;
        coords = other.coords;
        vels = other.vels;
        box_v1 = other.box_v1;
        box_v2 = other.box_v2;
        box_v3 = other.box_v3;
        resnums = other.resnums;
        resnams = other.resnams;
        atmnams = other.atmnams;
        atmnums = other.atmnums;
        parse_warnings = other.parse_warnings;

        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool Gro87::operator==(const Gro87 &other) const
{
    return ttle == other.ttle and current_time == other.current_time and
           coords == other.coords and vels == other.vels and
           box_v1 == other.box_v1 and box_v2 == other.box_v2 and box_v3 == other.box_v3 and
           resnums == other.resnums and resnams == other.resnams and
           atmnams == other.atmnams and atmnums == other.atmnums and
           parse_warnings == other.parse_warnings and
           MoleculeParser::operator==(other);
}

/** Comparison operator */
bool Gro87::operator!=(const Gro87 &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* Gro87::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Gro87>() );
}

/** Return the C++ name for this class */
const char* Gro87::what() const
{
    return Gro87::typeName();
}


/** Return the number of frames in the file */
int Gro87::count() const
{
    return this->nFrames();
}

/** Return the number of frames in the file */
int Gro87::size() const
{
    return this->nFrames();
}

/** Gro87 can be a lead parser as well as a follower */
bool Gro87::isLead() const
{
    return true;
}

/** Gro87 can be a lead parser as well as a follower */
bool Gro87::canFollow() const
{
    return true;
}

/** Return the Gro87 object that contains only the information for the ith
    frame. This allows you to extract and create a system for the ith frame
    from a trajectory */
Gro87 Gro87::operator[](int i) const
{
    i = Index(i).map( this->nFrames() );

    if (nFrames() == 1)
        return *this;

    Gro87 ret(*this);

    if (not coords.isEmpty())
    {
        ret.coords = { coords[i] };
    }

    if (not vels.isEmpty())
    {
        ret.vels = { vels[i] };
    }

    if (not current_time.isEmpty())
    {
        ret.current_time = { current_time[i] };
    }

    if (not box_v1.isEmpty())
    {
        ret.box_v1 = { box_v1[i] };
    }

    if (not box_v2.isEmpty())
    {
        ret.box_v2 = { box_v2[i] };
    }

    if (not box_v3.isEmpty())
    {
        ret.box_v3 = { box_v3[i] };
    }

    ret.assertSane();

    return ret;
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr Gro87::construct(const QString &filename,
                                   const PropertyMap &map) const
{
    return Gro87(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr Gro87::construct(const QStringList &lines,
                                   const PropertyMap &map) const
{
    return Gro87(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr Gro87::construct(const SireSystem::System &system,
                                  const PropertyMap &map) const
{
    return Gro87(system,map);
}

/** Return a string representation of this parser */
QString Gro87::toString() const
{
    if (nAtoms() == 0)
    {
        return QObject::tr("Gro87( nAtoms() = 0 )");
    }
    else
    {
        return QObject::tr("Gro87( title() = %1, nAtoms() = %2, nResidues() = %6, nFrames() = %5, "
                "hasCoordinates() = %3, hasVelocities() = %4 )")
                .arg(title()).arg(nAtoms())
                .arg(hasCoordinates()).arg(hasVelocities())
                .arg(nFrames()).arg(nResidues());
    }
}

/** Return the format name that is used to identify this file format within Sire */
QString Gro87::formatName() const
{
    return "Gro87";
}

/** Return a description of the file format */
QString Gro87::formatDescription() const
{
    return QObject::tr("Gromacs Gro87 structure format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList Gro87::formatSuffix() const
{
    static const QStringList suffixes = { "gro" };
    return suffixes;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void Gro87::assertSane() const
{
    //the number of coordinate and velocity frames should be identical
    QStringList errors;

    if (coords.count() != vels.count())
    {
        if (not (coords.isEmpty() or vels.isEmpty()))
        {
            errors.append( QObject::tr( "Error in Gro87 file as the number of "
               "coordinates frames (%1) read does not equal the number of "
               "velocity frames (2) read!")
                    .arg(coords.count()).arg(vels.count()) );
        }
    }

    //make sure that the number of atoms is consistent
    int nats = this->nAtoms();

    for (int i=0; i<coords.count(); ++i)
    {
        if (coords.at(i).count() != nats)
        {
            errors.append( QObject::tr( "Error: The number of atoms for coordinate "
               "frame %1 (%2) is not consistent with the number of atoms (%3).")
                    .arg(i).arg(coords.at(i).count()).arg(nats) );
        }
    }

    for (int i=0; i<vels.count(); ++i)
    {
        if (vels.at(i).count() != nats)
        {
            errors.append( QObject::tr( "Error: The number of atoms for velocity "
               "frame %1 (%2) is not consistent with the number of atoms (%3).")
                    .arg(i).arg(vels.at(i).count()).arg(nats) );
        }
    }

    if (atmnams.count() != nats)
    {
        errors.append( QObject::tr("Error: The number of atom names (%1) does not "
           "equal the number of atoms (%2)!")
                .arg(atmnams.count()).arg(nats) );
    }

    if (atmnums.count() != nats)
    {
        errors.append( QObject::tr("Error: The number of atom numbers (%1) does not "
           "equal the number of atoms (%2)!")
                .arg(atmnums.count()).arg(nats) );
    }

    if (resnams.count() != nats)
    {
        errors.append( QObject::tr("Error: The number of residue names (%1) does not "
           "equal the number of atoms (%2)!")
                .arg(resnams.count()).arg(nats) );
    }

    if (resnums.count() != nats)
    {
        errors.append( QObject::tr("Error: The number of residue numbers (%1) does not "
           "equal the number of atoms (%2)!")
                .arg(resnums.count()).arg(nats) );
    }

    if (box_v1.count() != box_v2.count() or box_v1.count() != box_v3.count())
    {
        errors.append( QObject::tr("Error: The number of frames of box dimension "
          "information is not consistent: %1 vs %2 vs %3.")
            .arg(box_v1.count()).arg(box_v2.count()).arg(box_v3.count()) );
    }

    if (not box_v1.isEmpty())
    {
        if (box_v1.count() != this->nFrames())
        {
            errors.append( QObject::tr("Error: The number of frames of box dimension "
               "information (%1) is not equal to the number of frames of trajectory (%2).")
                    .arg(box_v1.count()).arg(this->nFrames()) );
        }
    }

    if (not current_time.isEmpty())
    {
        if (current_time.count() != this->nFrames())
        {
            errors.append( QObject::tr("Error: The number of times from the trajectory "
              "(%1) does not equal the number of frames of trajectory (%2).")
                .arg(current_time.count()).arg(this->nFrames()) );
        }
    }

    //make sure that every atom with the same residue number has the same residue name
    if (not resnams.isEmpty())
    {
        QHash<qint64,QString> resnum_to_nam;
        resnum_to_nam.reserve(resnams.count());

        for (int i=0; i<resnams.count(); ++i)
        {
            const auto resnum = resnums[i];
            const auto resnam = resnams[i];

            if (resnum_to_nam.contains(resnum))
            {
                if (resnum_to_nam[resnum] != resnam)
                {
                    errors.append( QObject::tr("Error: Disagreement of the residue name "
                      "for residue number %1 for atom at index %2. It should be %3, "
                      "but for this atom it is %4.")
                        .arg(resnum).arg(i).arg(resnum_to_nam[resnum]).arg(resnam) );
                }
                else
                {
                    resnum_to_nam.insert(resnum,resnam);
                }
            }
        }
    }

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error( QObject::tr("There were errors reading the Gro87 format "
          "file:\n%1").arg(errors.join("\n\n")), CODELOC );
    }
}

/** Return the title of the file */
QString Gro87::title() const
{
    return ttle;
}

/** Return the current time of the simulation from which this coordinate
    file was written. Returns 0 if there is no time set. If there are
    multiple frames, then the time of the first frame is returned */
double Gro87::time() const
{
    if (current_time.isEmpty())
    {
        return 0.0;
    }
    else
    {
        return current_time[0];
    }
}

/** Return the time for the structure at the specified frame */
double Gro87::time(int frame) const
{
    return current_time[ Index(frame).map(current_time.count()) ];
}

/** Return the number of atoms whose data is contained in this file */
int Gro87::nAtoms() const
{
    if (coords.isEmpty())
    {
        if (vels.isEmpty())
        {
            return 0;
        }
        else
        {
            return vels.count();
        }
    }
    else
    {
        return coords[0].count();
    }
}

/** Return the number of unique residues in this file */
int Gro87::nResidues() const
{
    if (resnums.isEmpty())
    {
        return 0;
    }
    else
    {
        QHash<qint64,qint64> res_nats;
        res_nats.reserve(resnums.count());

        for (const auto resnum : resnums)
        {
            if (res_nats.contains(resnum))
            {
                res_nats[resnum] += 1;
            }
            else
            {
                res_nats.insert(resnum,1);
            }
        }

        return res_nats.count();
    }
}

/** Return whether or not this file contained coordinate data */
bool Gro87::hasCoordinates() const
{
    return not coords.isEmpty();
}

/** Return whether or not this file contained velocity data */
bool Gro87::hasVelocities() const
{
    return not vels.isEmpty();
}

/** Return the coordinates of the atoms for the first frame of the trajectory */
QVector<SireMaths::Vector> Gro87::coordinates() const
{
    if (coords.isEmpty())
    {
        return QVector<SireMaths::Vector>();
    }
    else
    {
        return coords[0];
    }
}

/** Return the velocities of the atoms for the first frame of the trajectory */
QVector<SireMaths::Vector> Gro87::velocities() const
{
    if (vels.isEmpty())
    {
        return QVector<SireMaths::Vector>();
    }
    else
    {
        return vels[0];
    }
}

/** Return the numbers of all of the atoms. These are in the same order
    as the coordinates */
QVector<qint64> Gro87::atomNumbers() const
{
    return atmnums;
}

/** Return the residue number for each atom (one per atom), in the same
    order as the coordinates */
QVector<qint64> Gro87::residueNumbers() const
{
    return resnums;
}

/** Return the names of all of the atoms, in the same order as the coordinates */
QVector<QString> Gro87::atomNames() const
{
    return atmnams;
}

/** Return the residue name for each atom (one per atom), in the same
    order as the coordinates */
QVector<QString> Gro87::residueNames() const
{
    return resnams;
}

/** Return the number of frames of the trajectory loaded from the file */
int Gro87::nFrames() const
{
    return coords.count();
}

/** Return the coordinates of the atoms at frame 'frame' */
QVector<SireMaths::Vector> Gro87::coordinates(int frame) const
{
    return coords[ Index(frame).map(coords.count()) ];
}

/** Return the velocities of the atoms at frame 'frame' */
QVector<SireMaths::Vector> Gro87::velocities(int frame) const
{
    return vels[ Index(frame).map(vels.count()) ];
}

/** Return the box V1 vector for the first frame */
SireMaths::Vector Gro87::boxV1() const
{
    if (box_v1.isEmpty())
    {
        return Vector(0);
    }
    else
    {
        return box_v1[0];
    }
}

/** Return the box V2 vector for the first frame */
SireMaths::Vector Gro87::boxV2() const
{
    if (box_v2.isEmpty())
    {
        return Vector(0);
    }
    else
    {
        return box_v2[0];
    }
}

/** Return the box V3 vector for the first frame */
SireMaths::Vector Gro87::boxV3() const
{
    if (box_v3.isEmpty())
    {
        return Vector(0);
    }
    else
    {
        return box_v3[0];
    }
}

/** Return the box V1 vector for the frame 'frame' */
SireMaths::Vector Gro87::boxV1(int frame) const
{
    return box_v1[ Index(frame).map(box_v1.count()) ];
}

/** Return the box V2 vector for the frame 'frame' */
SireMaths::Vector Gro87::boxV2(int frame) const
{
    return box_v2[ Index(frame).map(box_v2.count()) ];
}

/** Return the box V3 vector for the frame 'frame' */
SireMaths::Vector Gro87::boxV3(int frame) const
{
    return box_v3[ Index(frame).map(box_v3.count()) ];
}

/** Return the warnings encountered when parsing the file. This
    is empty if everything was ok */
QStringList Gro87::warnings() const
{
    return parse_warnings;
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void Gro87::parseLines(const PropertyMap &map)
{
    //file is described here - http://manual.gromacs.org/online/gro.html

    if (lines().count() < 2)
    {
        //there is nothing in the file
        return;
    }

    //the first line should be the title, with an optional timestep
    ttle = lines()[0].simplified();

    QRegularExpression re("t= ([-\\d\\.]+)$");

    //see if there is a timestep in the title
    //time is given as "title t= X.X"
    auto match = re.match(ttle);

    if (match.hasMatch())
    {
        const auto captured = match.captured(1);

        bool ok;
        double time = captured.toDouble(&ok);

        if (ok)
        {
            //we have converted this to a double - see if we need to clean up
            //any extra puncuation from the title
            ttle = ttle.remove(match.captured(0)).simplified();

            if (ttle.endsWith(","))
            {
                ttle = ttle.left( ttle.size()-1 );
            }

            //convert the time to picoseconds
            current_time.append( time );
        }
    }

    //the next line should be the number of atoms
    int nats = 0;
    {
        bool ok;
        nats = lines()[1].toInt(&ok);

        if (not ok)
        {
            throw SireIO::parse_error( QObject::tr(
                    "This does not look like a Gro87 file, as the second line should "
                    "contain just a free form integer which gives the number of atoms. "
                    "In this file, the second line is '%1'")
                        .arg(lines()[1]), CODELOC );
        }
    }

    //a trajectory is multiple copies of the file concatenated together. We should
    //now parse this file looking for lots of copies with this number of atoms
    int iframe = 0;
    int iline = 0;

    resnums.resize(nats);
    resnams.resize(nats);
    atmnums.resize(nats);
    atmnams.resize(nats);

    auto resnums_data = resnums.data();
    auto resnams_data = resnams.data();
    auto atmnums_data = atmnums.data();
    auto atmnams_data = atmnams.data();

    bool has_velocities = false;

    while (true)
    {
        //a complete set of information is 2 lines for the title, plus 1 line per atom,
        //plus one line for the box information
        if (iline + 2 + nats + 1 > lines().count())
        {
            //there is no more file to read
            break;
        }

        QVector<Vector> frame_coords( nats );
        QVector<Vector> frame_vels( nats );

        auto frame_coords_data = frame_coords.data();
        auto frame_vels_data = frame_vels.data();

        //internal function used to parse a single atom line in the file
        auto parse_atoms = [&](const QString &line, int iatm, bool *has_vels, QStringList &errors)
        {
            if (line.length() < 25)
            {
                errors.append( QObject::tr( "Cannot parse the data "
                   "for atom %1 as it does not match the format! '%2'")
                        .arg(iatm).arg(line) );

                return;
            }

            //residue number is the first 5 characters (integer)
            bool ok;
            qint64 resnum = line.midRef(0,5).toInt(&ok);

            if (not ok)
            {
                errors.append( QObject::tr( "Cannot extract the residue number "
                  "for atom %1 from the residue number part (%2) from line '%3'")
                    .arg(iatm).arg(line.mid(0,5)).arg(line) );
                return;
            }

            //the residue name is the next 5 characters
            QString resnam = line.mid(5,5).simplified();

            //the atom name is the next 5 characters
            QString atmnam = line.mid(10,5).simplified();

            //the atom number is the next 5 characters (integer)
            qint64 atmnum = line.midRef(15,5).toInt(&ok);

            if (not ok)
            {
                errors.append( QObject::tr( "Cannot extract the atom number "
                  "for atom %1 from the atom number part (%2) from line '%3'")
                    .arg(iatm).arg(line.mid(15,5)).arg(line) );
                return;
            }

            if (iframe == 0)
            {
                atmnams_data[iatm] = atmnam;
                atmnums_data[iatm] = atmnum;
                resnams_data[iatm] = resnam;
                resnums_data[iatm] = resnum;
            }
            else
            {
                //validate that this is the same information as for previous frames
                if (atmnam != atmnams_data[iatm] or
                    atmnum != atmnums_data[iatm] or
                    resnam != resnams_data[iatm] or
                    resnum != resnums_data[iatm])
                {
                    errors.append( QObject::tr("Disagreement in the ID of atom %1 "
                      "in frame %2 compared to the first frame. In the first frame "
                      "the atom is '%3 %4 %5 %6', but in this frame it is "
                      "'%7 %8 %9 %10'")
                        .arg(iatm).arg(iframe).arg(atmnams_data[iatm])
                        .arg(atmnums_data[iframe]).arg(resnams_data[iatm])
                        .arg(resnums_data[iframe])
                        .arg(atmnam).arg(atmnum).arg(resnam).arg(resnum) );
                    return;
                }
            }

            //now we need to read in the coordinate and velocity data. The format
            //is 3 columns of N+5 width (FN+5.N) coordinate data, and 3 columns
            //of N+5 width (FN+4.N+1) velocity data.

            //We must use the gaps between decimal points to work out the value of N
            const auto vals = line.mid(20);

            //find the indicies of all of the decimal points
            QVarLengthArray<int> point_idxs;

            int start = 0;

            while (true)
            {
                int idx = vals.indexOf('.', start);

                if (idx == -1)
                {
                    //no more decimal points
                    break;
                }

                point_idxs.append(idx);

                start = idx+1;
            }

            if (point_idxs.count() < 3)
            {
                errors.append( QObject::tr("Could not find at least three numbers "
                  "that can form the coordinates for atom %1 from '%2' as part of "
                  "the line '%3'").arg(iatm).arg(vals).arg(line) );
                return;
            }

            //there should be exactly 3 or 6 decimal points...
            if (point_idxs.count() != 3 and point_idxs.count() != 6)
            {
                errors.append( QObject::tr("There should be exactly 3 or 6 decimal "
                  "points for the coordinates (and optionally velocities) for atom %1 "
                  "from '%2' in line '%3'. We have found %4 decimal points!")
                    .arg(iatm).arg(vals).arg(line).arg(point_idxs.count()) );
            }

            //calculate the value of N using the coordinates
            int n = point_idxs[1] - point_idxs[0];

            if (point_idxs[2] - point_idxs[1] != n)
            {
                errors.append( QObject::tr("There coordinate data should be written with "
                  "a fixed format of consistent size. The coordinates for atom %1 "
                  "from '%2' in line '%3' has an inconsistent width! %4 versus %5")
                    .arg(iatm).arg(vals).arg(line).arg(n).arg(point_idxs[2]-point_idxs[1]) );
                return;
            }

            if (vals.length() < 3*n)
            {
                errors.append( QObject::tr("The coordinate line for atom %1 is not long "
                  "enough to contain the data. It should be %2 characters, but is really "
                  "%3 characters!").arg(3*n).arg(vals.length()) );
                return;
            }

            bool ok_x, ok_y, ok_z;
            double x = vals.midRef(0,n).toDouble(&ok_x);
            double y = vals.midRef(n,n).toDouble(&ok_y);
            double z = vals.midRef(2*n,n).toDouble(&ok_z);

            if (not (ok_x and ok_y and ok_z))
            {
                errors.append( QObject::tr("There was a problem reading the coordinate "
                  "values of x, y, and z for atom %1 from the data '%2' in line '%3'")
                    .arg(iatm).arg(vals.mid(0,3*n)).arg(line) );
                return;
            }

            //coordinates in file in nanometers - convert to angstroms
            frame_coords_data[iatm] = Vector( 10.0 * x, 10.0 * y, 10.0 * z );

            if (point_idxs.count() < 6)
            {
                return;
            }

            //now read in the velocities
            int vlen = 6*n;

            if (vals.length() < vlen)
            {
                errors.append( QObject::tr("The velocity line for atom %1 is not long "
                  "enough to contain the data. It should be %2 characters, but is really "
                  "%3 characters!").arg(iatm).arg(vlen).arg(vals.length()) );
                return;
            }

            x = vals.midRef(3*n, n).toDouble(&ok_x);
            y = vals.midRef(4*n, n).toDouble(&ok_y);
            z = vals.midRef(5*n, n).toDouble(&ok_z);

            if (not (ok_x and ok_y and ok_z))
            {
                errors.append( QObject::tr("There was a problem reading the velocity "
                  "values of x, y, and z for atom %1 from the data '%2' in line '%3'")
                    .arg(iatm).arg(vals.mid(3*n,3*n)).arg(line) );
                return;
            }

            *has_vels = true;

            //convert from nanometers per picosecond to angstroms per picosecond
            frame_vels_data[iatm] = Vector( 10.0*x, 10.0*y, 10.0*z );
        };

        if (usesParallel())
        {
            QMutex mutex;

            tbb::parallel_for( tbb::blocked_range<int>(0,nats),
                               [&](const tbb::blocked_range<int> &r)
            {
                QStringList local_errors;
                bool local_has_vels = false;

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    parse_atoms( lines().constData()[iline+2+i], i, &local_has_vels, local_errors );
                }

                if (local_has_vels)
                {
                    QMutexLocker lkr(&mutex);

                    if (not has_velocities)
                    {
                        has_velocities = true;
                    }
                }

                if (not local_errors.isEmpty())
                {
                    QMutexLocker lkr(&mutex);
                    parse_warnings += local_errors;
                }
            });
        }
        else
        {
            for (int i=0; i<nats; ++i)
            {
                parse_atoms( lines().constData()[iline+2+i], i, &has_velocities, parse_warnings );
            }
        }

        //save the data
        coords.append( frame_coords );

        if (has_velocities)
        {
            //if no previous frame has velocities, then set them to zero
            if (vels.count() < iframe)
            {
                QVector<Vector> zero_vels( frame_vels.count(), Vector(0) );

                while (vels.count() < iframe)
                {
                    vels.append(zero_vels);
                }
            }

            vels.append( frame_vels );
        }

        //next, read in the periodic box information. This is a single line
        //of 3 or 9 space separated real numbers containing
        //v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y) (the last six values can be omitted)
        {
            const auto boxline = lines()[iline + 2 + nats];
            const auto words = boxline.split(" ", Qt::SkipEmptyParts);

            bool all_ok = false;
            Vector v1(0), v2(0), v3(0);

            if ( words.count() == 3 )
            {
                bool ok_x, ok_y, ok_z;
                double x = words[0].toDouble(&ok_x);
                double y = words[1].toDouble(&ok_y);
                double z = words[2].toDouble(&ok_z);

                if (ok_x and ok_y and ok_z)
                {
                    v1 = Vector(x, 0, 0);
                    v2 = Vector(0, y, 0);
                    v3 = Vector(0, 0, z);
                    all_ok = true;
                }
            }
            else if (words.count() == 9)
            {
                all_ok = true;
                double v[9];

                for (int k=0; k<9; ++k)
                {
                    bool this_ok;
                    v[k] = words[k].toDouble( &this_ok );

                    if (not this_ok)
                    {
                        all_ok = false;
                        break;
                    }
                }

                if (all_ok)
                {
                    v1 = Vector(v[0], v[3], v[4]);
                    v2 = Vector(v[5], v[1], v[6]);
                    v3 = Vector(v[7], v[8], v[2]);
                }
            }

            if (all_ok)
            {
                while (box_v1.count() < iframe)
                {
                    box_v1.append( Vector(0) );
                    box_v2.append( Vector(0) );
                    box_v3.append( Vector(0) );
                }

                //remember to convert nanometers to angstroms
                box_v1.append( 10.0 * v1 );
                box_v2.append( 10.0 * v2 );
                box_v3.append( 10.0 * v3 );
            }
            else
            {
                parse_warnings.append( QObject::tr( "Cannot read the periodic box information "
                  "for frame %1 from the line '%2'. This should be a space-separated list "
                  "of three or nine numbers...")
                    .arg(iframe).arg(boxline) );

                //no need to break here as we can still make progress with this frame
            }
        }

        //finally, the coords and velocities were read, so see if there was a time value
        //for this frame (the t= X.X in the title line, which is lines()[iline])
        if (iframe != 0)
        {
            auto match = re.match( lines()[iline] );

            if (match.hasMatch())
            {
                const auto captured = match.captured(1);

                bool ok;
                double time = captured.toDouble(&ok);

                if (ok)
                {
                    while (current_time.count() < iframe)
                    {
                        current_time.append(0);
                    }

                    //convert the time to picoseconds
                    current_time.append( time );
                }
            }
        }

        //increment the number of read frames and the start of the next line
        iframe += 1;
        iline += 2 + nats + 1;

        //there must be enough extra lines to contain another frame of trajectory...
        while (true)
        {
            if (iline + 2 + nats + 1 > lines().count())
            {
                //there is no more file to read
                break;
            }

            //the next line should be the number of atoms. If not then there may be an
            //extra (incorrect) blank line, and we need to keep advancing through the file
            //until we find that extra line
            bool ok;
            int new_nats = lines()[iline+1].toInt(&ok);

            if ( ok and (new_nats == nats) )
            {
                //the number of atoms has been found and matches the expected value
                break;
            }

            iline += 1;
        }
    }

    this->setScore(nats);
}

/** This function finds the index of the atom called 'atomname' with number
    'atomnum' in residue 'resname' with residue number 'resnum'. The passed
    hint suggest where in the array to start looking */
int Gro87::findAtom(const MoleculeInfoData &molinfo, int atmidx, int hint,
                    bool *ids_match) const
{
    const int natoms = atmnams.count();

    if (hint < 0 or hint >= natoms)
        hint = 0;

    //get the atom/residue name/number for this atom
    const auto residx = molinfo.parentResidue( AtomIdx(atmidx) );
    const auto atmname = molinfo.name( AtomIdx(atmidx) ).value();
    const auto atmnum = molinfo.number( AtomIdx(atmidx) ).value();
    const auto resname = molinfo.name(residx).value();
    const auto resnum = molinfo.number(residx).value();

    QMap<int,int> partial_matches;

    auto is_match = [&](int i)
    {
        bool same_atmname = (atmname == atmnams.constData()[i]);
        bool same_atmnum = (atmnum == atmnums.constData()[i]);
        bool same_resname = (resname == resnams.constData()[i]);
        bool same_resnum = (resnum == resnums.constData()[i]);

        //make sure that we set the flag 'ids_match' to false if it
        //is not already false and any of the IDs don't match up
        if (ids_match)
        {
            if (*ids_match and not (same_atmname and same_atmnum and same_resname and same_resnum))
            {
                *ids_match = false;
            }
        }

        if (same_atmname and same_resname and i == hint)
        {
            return true;
        }
        else if (same_atmname and same_resname and same_atmnum and same_resnum)
        {
            return true;
        }
        else if (same_atmname)
        {
            //score a partial match - MUST match atom name, then ideally
            //residue name, then atom number than residue number
            partial_matches.insert( (10 * same_resname) +
                                    (5 * same_atmnum) +
                                    (2 * same_resnum) +
                                    (1 * (i == hint)), i );
        }

        return false;
    };

    //scan through the atoms and find a match (starting from 'hint')
    int match = -1;

    for (int i=hint; i<natoms; ++i)
    {
        if (is_match(i))
        {
            match = i;
            break;
        }
    }

    if (match == -1)
    {
        for (int i=0; i<hint; ++i)
        {
            if (is_match(i))
            {
                match = i;
                break;
            }
        }
    }

    if (match == -1)
    {
        //no match found - did we find any partial matches?
        if (not partial_matches.isEmpty())
        {
            //return the match with the highest score
            return partial_matches.last();
        }
        else
        {
            //we can't find the atom, so just go with the hint
            return hint;
        }
    }

    return match;
}

/** Internal function used to add the final properties to a system */
void Gro87::finaliseSystem(System &system, const PropertyMap &map) const
{
    PropertyName space_property = map["space"];
    if (space_property.hasValue())
    {
        system.setProperty("space", space_property.value());
    }
    else if ((not box_v1.isEmpty()) and space_property.hasSource())
    {
        // Check whether the box is cubic, i.e (x,0,0), (0,y,0), (0,0,z)
        double x = box_v1.at(0).x();
        double y = box_v2.at(0).y();
        double z = box_v3.at(0).z();

        if (box_v1.at(0).manhattanLength() != x or
            box_v2.at(0).manhattanLength() != y or
            box_v3.at(0).manhattanLength() != z)
        {
            // This is a triclinic space.
            system.setProperty( space_property.source(), SireVol::TriclinicBox(box_v1.at(0),
                                                                               box_v2.at(0),
                                                                               box_v3.at(0)) );
        }
        else if (x + y + z > 0)
        {
            system.setProperty( space_property.source(), SireVol::PeriodicBox(Vector(x,y,z)) );
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

    if (not current_time.isEmpty())
    {
        PropertyName time_property = map["time"];

        if (time_property.hasSource())
        {
            system.setProperty(time_property.source(), TimeProperty(current_time[0]*picosecond));
        }
        else
        {
            system.setProperty("time", TimeProperty(current_time[0]*picosecond));
        }
    }
}

/** Use the data contained in this parser to create a new System from scratch.
    This will be a one-molecule system as Gro87 files don't divide atoms up
    into molecules. */
System Gro87::startSystem(const PropertyMap &map) const
{
    const bool has_coords = hasCoordinates();
    const bool has_vels = hasVelocities();

    if (not (has_coords or has_vels))
    {
        //nothing to add...
        return System();
    }

    //loop through all atoms and add them to a single molecule
    Molecule mol;
    {
        MolStructureEditor moleditor;

        //add all of the atoms, creating residues as needed
        const auto atmnams = this->atomNames();
        const auto atmnums = this->atomNumbers();
        const auto resnums = this->residueNumbers();
        const auto resnams = this->residueNames();

        int ncg = 0;

        QSet<ResNum> completed_residues;

        for (int i=0; i<atmnams.count(); ++i)
        {
            auto atom = moleditor.add( AtomNum(atmnums[i]) );
            atom = atom.rename( AtomName(atmnams[i]) );

            const ResNum resnum(resnums[i]);

            if (completed_residues.contains(resnum))
            {
                auto res = moleditor.residue(resnum);

                if (res.name().value() != resnams[i])
                {
                    //different residue
                    res = moleditor.add(resnum);
                    res = res.rename( ResName(resnams[i]) );
                    ncg += 1;
                    moleditor.add( CGName(QString::number(ncg)) );
                }

                atom = atom.reparent(res.index());
                atom = atom.reparent(CGName(QString::number(ncg)));
            }
            else
            {
                auto res = moleditor.add(resnum);
                res = res.rename( ResName(resnams[i]) );
                atom = atom.reparent(res.index());

                ncg += 1;
                auto cg = moleditor.add( CGName(QString::number(ncg)) );
                atom = atom.reparent(cg.index());

                completed_residues.insert(resnum);
            }
        }

        //we have created the molecule - now add in the coordinates/velocities as needed
        mol = moleditor.commit();
    }

    //now add the coordinates and velocities
    {
        auto moleditor = mol.edit();
        const auto molinfo = mol.info();

        if (has_coords)
        {
            auto coords = QVector< QVector<Vector> >(molinfo.nCutGroups());
            const auto coords_array = this->coordinates().constData();

            for (int i=0; i<molinfo.nCutGroups(); ++i)
            {
                coords[i] = QVector<Vector>(molinfo.nAtoms(CGIdx(i)));
            }

            for (int i=0; i<mol.nAtoms(); ++i)
            {
                auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(i));

                coords[cgatomidx.cutGroup()][cgatomidx.atom()] = coords_array[i];
            }

            moleditor.setProperty( map["coordinates"], AtomCoords(CoordGroupArray(coords)) );
        }

        if (has_vels)
        {
            auto vels = AtomVelocities(molinfo);
            const auto vels_array = this->velocities().constData();

            for (int i=0; i<mol.nAtoms(); ++i)
            {
                auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(i));

                //velocity is Angstroms per 1/20.455 ps
                const auto vel_unit = (1.0 / 20.455) * angstrom / picosecond;

                const Vector &vel = vels_array[i];
                vels.set(cgatomidx, Velocity3D(vel.x() * vel_unit,
                                               vel.y() * vel_unit,
                                               vel.z() * vel_unit));
            }

            moleditor.setProperty( map["velocity"], vels );
        }

        mol = moleditor.commit();
    }

    //now that we have the molecule, add this to the System
    System system(this->title());

    MoleculeGroup all("all");
    all.add(mol);

    system.add(all);

    this->finaliseSystem(system,map);

    return system;
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in
    the passed System that are missing coordinate data. */
void Gro87::addToSystem(System &system, const PropertyMap &map) const
{
    const bool has_coords = hasCoordinates();
    const bool has_vels = hasVelocities();

    if (not (has_coords or has_vels))
    {
        //nothing to add...
        return;
    }

    //first, we are going to work with the group of all molecules, which should
    //be called "all". We have to assume that the molecules are ordered in "all"
    //in the same order as they are in this gro file, with the data
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
        throw SireError::incompatible_error( QObject::tr(
                "Incompatibility between the files, as this .gro file contains data "
                "for %1 atom(s), while the other file(s) have created a system with "
                "%2 atom(s)").arg(this->nAtoms()).arg(natoms), CODELOC );

    //next, copy the coordinates and optionally the velocities into the molecules
    QVector<Molecule> mols(nmols);
    Molecule *mols_array = mols.data();

    const PropertyName coords_property = map["coordinates"];
    const PropertyName vels_property = map["velocity"];

    const Vector *coords_array = this->coordinates().constData();
    const Vector *vels_array = this->velocities().constData();

    auto add_moldata = [&](int i)
    {
        const int atom_start_idx = atom_pointers.constData()[i];
        auto mol = allmols[MolIdx(i)].molecule();
        const auto molinfo = mol.data().info();

        //first, get the index of each atom in the gro file
        QVector<int> idx_in_gro( molinfo.nAtoms(), -1 );

        bool ids_match = true;

        for (int j=0; j<molinfo.nAtoms(); ++j)
        {
            idx_in_gro[j] = findAtom( molinfo, j, atom_start_idx + j, &ids_match );
        }

        // Convert the vector to a set to check for duplicates.
        auto unique_idx_in_gro = convert_to_qset(idx_in_gro);

        // Duplicate matches were found!
        if (idx_in_gro.count() != unique_idx_in_gro.count())
        {
            throw SireError::incompatible_error( QObject::tr(
                "Incompatibility between the files, could not match all atoms to the "
                "corresponding topology file. Check atom and residue names!"), CODELOC );
        }

        //if the IDs don't match, then we need to update the ID information
        //of the atoms and residues in the molecule
        if (not ids_match)
        {
            QHash<AtomNum,AtomNum> renumbered_atoms;
            QHash<ResNum,ResNum> renumbered_residues;

            for (int j=0; j<molinfo.nAtoms(); ++j)
            {
                int idx = idx_in_gro[j];

                auto oldnum = molinfo.number( AtomIdx(j) );
                AtomNum newnum(atmnums.constData()[idx]);

                if (oldnum != newnum)
                    renumbered_atoms.insert(oldnum, newnum);
            }

            for (int j=0; j<molinfo.nResidues(); ++j)
            {
                auto oldnum = molinfo.number( ResIdx(j) );

                auto atoms_in_res = molinfo.getAtomsIn(ResIdx(j));

                if (not atoms_in_res.isEmpty())
                {
                    int idx = idx_in_gro[ atoms_in_res.at(0).value() ];

                    ResNum newnum(resnums.constData()[idx]);

                    if (oldnum != newnum)
                        renumbered_residues.insert(oldnum, newnum);
                }
            }

            mol = mol.edit().renumber(renumbered_atoms, renumbered_residues).commit();
        }

        //now use this index to locate the correct coordinate and/or velocity
        //data to add to the molecules
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

                const int atom_idx = idx_in_gro.constData()[j];

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

                const int atom_idx = idx_in_gro.constData()[j];

                //velocity is Angstroms per 1/20.455 ps
                const auto vel_unit = (1.0 / 20.455) * angstrom / picosecond;

                const Vector &vel = vels_array[atom_idx];
                vels.set(cgatomidx, Velocity3D(vel.x() * vel_unit,
                                               vel.y() * vel_unit,
                                               vel.z() * vel_unit));
            }

            moleditor.setProperty( vels_property, vels );
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

    MoleculeGroup new_group("all");

    for (const auto &mol : mols)
    {
        new_group.add(mol);
    }

    system.remove(MGName("all"));
    system.add(new_group);

    this->finaliseSystem(system, map);
}
