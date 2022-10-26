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

#include "SireIO/amberrst.h"
#include "SireIO/amberformat.h"
#include "SireIO/netcdffile.h"

#include "SireSystem/system.h"

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
#include "SireBase/timeproperty.h"
#include "SireBase/booleanproperty.h"
#include "SireBase/getinstalldir.h"
#include "SireBase/unittest.h"

#include "SireIO/errors.h"

#include "SireStream/shareddatastream.h"

using namespace SireIO;
using namespace SireIO::detail;
using namespace SireMaths;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireVol;
using namespace SireUnits;
using namespace SireStream;

static const RegisterMetaType<AmberRst> r_rst;
const RegisterParser<AmberRst> register_amberrst;

QDataStream &operator<<(QDataStream &ds, const AmberRst &rst)
{
    writeHeader(ds, r_rst, 1);

    SharedDataStream sds(ds);

    sds << rst.ttle << rst.current_time
        << rst.coords << rst.vels << rst.frcs
        << rst.box_dims << rst.box_angs
        << rst.convention_version << rst.creator_app
        << rst.parse_warnings << rst.created_from_restart
        << static_cast<const MoleculeParser&>(rst);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AmberRst &rst)
{
    VersionID v = readHeader(ds, r_rst);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> rst.ttle >> rst.current_time
            >> rst.coords >> rst.vels >> rst.frcs
            >> rst.box_dims >> rst.box_angs
            >> rst.convention_version >> rst.creator_app
            >> rst.parse_warnings >> rst.created_from_restart
            >> static_cast<MoleculeParser&>(rst);
    }
    else
        throw version_error(v, "1", r_rst, CODELOC);

    return ds;
}

static Vector cubic_angs(90,90,90);

/** Constructor */
AmberRst::AmberRst()
         : ConcreteProperty<AmberRst,MoleculeParser>()
{}

/** Return the format name that is used to identify this file format within Sire */
QString AmberRst::formatName() const
{
    return "RST";
}

/** Return the suffixes that AmberRst files will typically use */
QStringList AmberRst::formatSuffix() const
{
    static const QStringList suffixes = { "rst", "crd", "trj", "traj" };
    return suffixes;
}

/** Internal function called to assert that this object is in a sane state */
void AmberRst::assertSane() const
{
    const int nframes = this->nFrames();
    const int natoms = this->nAtoms();

    if (not coords.isEmpty())
    {
        SireBase::assert_equal( coords.count(), nframes, CODELOC );

        for (int i=0; i<nframes; ++i)
        {
            SireBase::assert_equal( coords[i].count(), natoms, CODELOC );
        }
    }

    if (not vels.isEmpty())
    {
        SireBase::assert_equal( vels.count(), nframes, CODELOC );

        for (int i=0; i<nframes; ++i)
        {
            SireBase::assert_equal( vels[i].count(), natoms, CODELOC );
        }
    }

    if (not frcs.isEmpty())
    {
        SireBase::assert_equal( frcs.count(), nframes, CODELOC );

        for (int i=0; i<nframes; ++i)
        {
            SireBase::assert_equal( frcs[i].count(), natoms, CODELOC );
        }
    }

    if (not box_dims.isEmpty())
    {
        SireBase::assert_equal( box_dims.count(), nframes, CODELOC );
        SireBase::assert_equal( box_angs.count(), nframes, CODELOC );
    }
}

/** This is not a text file */
bool AmberRst::isTextFile() const
{
    return false;
}

/** Return a description of the file format */
QString AmberRst::formatDescription() const
{
    return QObject::tr("Amber coordinate/velocity binary (netcdf) restart/trajectory files "
                       "supported since Amber 9, now default since Amber 16.");
}

/** Parse the data contained in the lines - this clears any pre-existing
    data in this object */
void AmberRst::parse(const NetCDFFile &netcdf, const PropertyMap &map)
{
    //the netcdf conventions for Amber restart/trajectory files are
    //given here - http://ambermd.org/netcdf/nctraj.xhtml

    const bool uses_parallel = this->usesParallel();

    //collect all of the conventions
    QString conventions = netcdf.getStringAttribute("Conventions");

    //these must contain "AMBER" or "AMBERRESTART", else this may not be a valid file
    created_from_restart = false;

    if (conventions.contains("AMBERRESTART"))
    {
        created_from_restart = true;
    }
    else if (not conventions.contains("AMBER"))
    {
        //this is not an amber file
        throw SireIO::parse_error( QObject::tr(
                "This does not look like a valid NetCDF Amber restart/trajectory file. "
                "Such files have a 'Conventions' attribute that contains either "
                "'AMBER' or 'AMBERRESTART'. The 'Conventions' in this file equals '%1'")
                    .arg(conventions), CODELOC );
    }

    bool ok;
    convention_version = netcdf.getStringAttribute("ConventionVersion").toDouble(&ok);

    if (not ok)
    {
        throw SireIO::parse_error( QObject::tr(
                "This does not look like a valid NetCDF Amber restart/trajectory file. "
                "Such files have a 'ConventionVersion' attribute that gives the version "
                "number of the file format. The value in this file is '%1'")
                    .arg(netcdf.getStringAttribute("ConventionVersion")), CODELOC );
    }

    if (convention_version != 1.0)
    {
        parse_warnings.append( QObject::tr(
                "Reading in an Amber NetCDF CRD/TRJ file with an unsupported version: %1")
                    .arg(convention_version) );
    }

    QString program = netcdf.getStringAttribute("program");
    QString program_version = netcdf.getStringAttribute("programVersion");

    try
    {
        QString application = netcdf.getStringAttribute("application");

        creator_app = QString("%1 %2 version %3").arg(application)
                            .arg(program).arg(program_version);
    }
    catch(...)
    {
        creator_app = QString("%1 version %2").arg(program).arg(program_version);
    }

    try
    {
        ttle = netcdf.getStringAttribute("title");
    }
    catch(...)
    {
        ttle = "default";
    }

    //get information about all of the variables
    const auto varinfos = netcdf.getVariablesInfo();

    //get information about the number of atoms in the file
    const int natoms = netcdf.getDimensions().value("atom", -1);

    if (natoms < 0)
    {
        throw SireIO::parse_error( QObject::tr(
                "This does not look like a valid NetCDF Amber restart/trajectory file. "
                "Such files have a 'atom' dimension that gives the number of atoms. "
                "This dimension is missing from this file. Available dimensions are %1")
                    .arg( QStringList(netcdf.getDimensions().keys()).join(",") ),
                        CODELOC );
    }

    //now lets create some lambda functions that do the different orthogonal
    //work needed to parse the file. These functions will either be called serially
    //or in parallel depending on whether this is supported or requested
    QMutex warnings_mutex;

    //function to parse in all of the base data in the file
    auto parseBase = [&]()
    {
        //extract the current time of the simulation
        current_time.clear();

        if (varinfos.contains("time"))
        {
            const auto data = netcdf.read( varinfos["time"] );

            const auto atts = data.attributes();

            double scale_factor = 1.0;

            if (atts.contains("scale_factor"))
            {
                scale_factor = atts["scale_factor"].toDouble();
            }

            if (atts.contains("units"))
            {
                QString units = atts["units"].toString();

                if (units != "picosecond")
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Units are not the default of picoseconds. "
                        "They are actually %1!").arg(units) );
                }

                //we need to interpret different units and set the scale_factor accordingly.
                //We want to store the time here in picoseconds
            }

            if (data.nValues() < 1)
            {
                QMutexLocker lkr(&warnings_mutex);
                parse_warnings.append( QObject::tr("Could not interpret the time "
                                                   "as there is no value?") );
            }
            else
            {
                if (created_from_restart)
                {
                    current_time.append( scale_factor * data.toDoubleArray()[0] );
                }
                else
                {
                    const auto times = data.toFloatArray();
                    const int nframes = times.count();
                    current_time.reserve(nframes);

                    for (int i=0; i<nframes; ++i)
                    {
                        current_time.append( scale_factor * times[i] );
                    }
                }
            }
        }

        //now lets get any periodic box information
        if (varinfos.contains("cell_lengths"))
        {
            const auto data = netcdf.read( varinfos["cell_lengths"] );

            const auto atts = data.attributes();

            double scale_factor = 1.0;

            if (atts.contains("scale_factor"))
            {
                scale_factor = atts["scale_factor"].toDouble();
            }

            if (atts.contains("units"))
            {
                QString units = atts["units"].toString();

                if (units != "angstrom")
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Units are not the default of angstroms. "
                        "They are actually %1!").arg(units) );
                }

                //we need to interpret different units and set the scale_factor accordingly.
                //We want to store the lengths here in angstroms
            }

            if (data.nValues() < 3)
            {
                QMutexLocker lkr(&warnings_mutex);
                parse_warnings.append( QObject::tr("Could not interpret the box lengths "
                                                   "as there are less than 3 values (%1)?")
                                                        .arg(data.nValues()) );
            }
            else if (data.nValues() % 3 != 0)
            {
                QMutexLocker lkr(&warnings_mutex);
                parse_warnings.append( QObject::tr("Could not interpret the box lengths "
                                                   "as the number of values is not evenly "
                                                   "divisible by 3 (%1)")
                                                    .arg(data.nValues()) );
            }
            else
            {
                int nframes = data.nValues() / 3;
                const auto vals = data.toDoubleArray();

                for (int i=0; i<nframes; ++i)
                {
                    box_dims.append( scale_factor * Vector( vals[3*i + 0],
                                                            vals[3*i + 1],
                                                            vals[3*i + 2] ) );
                }
            }
        }

        //now lets get any periodic box angle information
        if (varinfos.contains("cell_angles"))
        {
            const auto data = netcdf.read( varinfos["cell_angles"] );

            const auto atts = data.attributes();

            double scale_factor = 1.0;

            if (atts.contains("scale_factor"))
            {
                scale_factor = atts["scale_factor"].toDouble();
            }

            if (atts.contains("units"))
            {
                QString units = atts["units"].toString();

                if (units != "degree")
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Units are not the default of degrees. "
                        "They are actually %1!").arg(units) );
                }

                //we need to interpret different units and set the scale_factor accordingly.
                //We want to store the angles here in degrees
            }

            if (data.nValues() < 3)
            {
                QMutexLocker lkr(&warnings_mutex);
                parse_warnings.append( QObject::tr("Could not interpret the box angles "
                                                   "as there are less than 3 values (%1)?")
                                                        .arg(data.nValues()) );
            }
            else if (data.nValues() % 3 != 0)
            {
                QMutexLocker lkr(&warnings_mutex);
                parse_warnings.append( QObject::tr("Could not interpret the box angles "
                                                   "as the number of values is not evenly "
                                                   "divisible by 3 (%1)")
                                                    .arg(data.nValues()) );
            }
            else
            {
                int nframes = data.nValues() / 3;
                const auto vals = data.toDoubleArray();

                for (int i=0; i<nframes; ++i)
                {
                    box_angs.append( scale_factor * Vector( vals[3*i + 0],
                                                            vals[3*i + 1],
                                                            vals[3*i + 2] ) );
                }
            }
        }
    };

    //function to read in all of the coordinates
    auto parseCoordinates = [&]()
    {
        coords = QVector< QVector<Vector> >();

        if (varinfos.contains("coordinates"))
        {
            const auto data = netcdf.read( varinfos["coordinates"] );

            const auto atts = data.attributes();

            double scale_factor = 1.0;

            if (atts.contains("scale_factor"))
            {
                bool ok;
                scale_factor = atts["scale_factor"].toDouble(&ok);

                if (not ok)
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Could not extract the scale factor "
                       "for the coordinates from the string '%1'")
                            .arg(atts["scale_factor"].toString()) );
                }
            }

            if (atts.contains("units"))
            {
                QString units = atts["units"].toString();

                if (units != "angstrom")
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Units are not the default of angstroms. "
                        "They are actually %1!").arg(units) );
                }

                //we need to interpret different units and set the scale_factor accordingly.
                //We want to store the lengths here in angstroms
            }

            if (data.nValues() % 3 != 0)
            {
                QMutexLocker lkr(&warnings_mutex);
                parse_warnings.append( QObject::tr("Could not interpret the coordinates "
                                                   "as they are not divisible by three (%1)?")
                                                        .arg(data.nValues()) );
                return;
            }
            else
            {
                const auto vals = data.toDoubleArray();

                const int natoms_times_nframes = vals.count() / 3;

                if (natoms_times_nframes % natoms != 0)
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Could not interpret the coordinates "
                                    "as the amount of data (%1) does not equal 3 x nframes(%2) "
                                    "x natoms(%3)")
                                        .arg(vals.count()).arg(natoms_times_nframes/natoms)
                                        .arg(natoms) );
                    return;
                }

                const int nframes = natoms_times_nframes / natoms;

                coords = QVector< QVector<Vector> >(nframes);

                if (uses_parallel)
                {
                    for (int i=0; i<nframes; ++i)
                    {
                        QVector<Vector> c(natoms);

                        tbb::parallel_for( tbb::blocked_range<int>(0,natoms),
                                           [&](const tbb::blocked_range<int> &r)
                        {
                            for (int j=r.begin(); j<r.end(); ++j)
                            {
                                const int idx = 3*natoms*i + 3*j;

                                c[j] = scale_factor * Vector( vals[idx],
                                                              vals[idx+1],
                                                              vals[idx+2] );
                            }
                        });

                        coords[i] = c;
                    }
                }
                else
                {
                    for (int i=0; i<nframes; ++i)
                    {
                        QVector<Vector> c(natoms);

                        for (int j=0; j<natoms; ++j)
                        {
                            const int idx = 3*natoms*i + 3*j;

                            c[j] = scale_factor * Vector( vals[idx], vals[idx+1], vals[idx+2] );
                        }

                        coords[i] = c;
                    }
                }
            }
        }
    };

    //function to parse all of the velocities
    auto parseVelocities = [&]()
    {
        vels = QVector< QVector<Vector> >();

        if (varinfos.contains("velocities"))
        {
            const auto data = netcdf.read( varinfos["velocities"] );

            const auto atts = data.attributes();

            double scale_factor = 1.0;

            if (atts.contains("scale_factor"))
            {
                bool ok;
                scale_factor = atts["scale_factor"].toDouble(&ok);

                if (not ok)
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Could not extract the scale factor "
                       "for the velocities from the string '%1'")
                            .arg(atts["scale_factor"].toString()) );
                }
            }

            if (atts.contains("units"))
            {
                QString units = atts["units"].toString();

                if (units != "angstrom/picosecond")
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Units are not the default of "
                        "angstrom/picosecond. They are actually %1!").arg(units) );
                }

                //we need to interpret different units and set the scale_factor accordingly.
                //We want to store the velocities here in angstroms/picosecond
            }

            //we want to store the velocities here in amber units, i.e.
            //angstroms / 1/20.455 picosecond
            scale_factor /= 20.455;

            if (data.nValues() % 3 != 0)
            {
                QMutexLocker lkr(&warnings_mutex);
                parse_warnings.append( QObject::tr("Could not interpret the velocities "
                                                   "as they are not divisible by three (%1)?")
                                                        .arg(data.nValues()) );
            }
            else
            {
                const auto vals = data.toDoubleArray();

                const int natoms_times_nframes = vals.count() / 3;

                if (natoms_times_nframes % natoms != 0)
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Could not interpret the velocities "
                                    "as the amount of data (%1) does not equal 3 x nframes(%2) "
                                    "x natoms(%3)")
                                        .arg(vals.count()).arg(natoms_times_nframes/natoms)
                                        .arg(natoms) );
                    return;
                }

                const int nframes = natoms_times_nframes / natoms;

                vels = QVector< QVector<Vector> >(nframes);

                if (uses_parallel)
                {
                    for (int i=0; i<nframes; ++i)
                    {
                        QVector<Vector> v(natoms);

                        tbb::parallel_for( tbb::blocked_range<int>(0,natoms),
                                           [&](const tbb::blocked_range<int> &r)
                        {
                            for (int j=r.begin(); j<r.end(); ++j)
                            {
                                const int idx = 3*natoms*i + 3*j;

                                v[j] = scale_factor * Vector( vals[idx],
                                                              vals[idx+1],
                                                              vals[idx+2] );
                            }
                        });

                        vels[i] = v;
                    }
                }
                else
                {
                    for (int i=0; i<nframes; ++i)
                    {
                        QVector<Vector> v(natoms);

                        for (int j=0; j<natoms; ++j)
                        {
                            const int idx = 3*natoms*i + 3*j;

                            v[j] = scale_factor * Vector( vals[idx], vals[idx+1], vals[idx+2] );
                        }

                        vels[i] = v;
                    }
                }
            }
        }
    };

    //function to parse all of the forces
    auto parseForces = [&]()
    {
        frcs = QVector< QVector<Vector> >();

        if (varinfos.contains("forces"))
        {
            const auto data = netcdf.read( varinfos["forces"] );

            const auto atts = data.attributes();

            double scale_factor = 1.0;

            if (atts.contains("scale_factor"))
            {
                bool ok;
                scale_factor = atts["scale_factor"].toDouble(&ok);

                if (not ok)
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Could not extract the scale factor "
                       "for the forces from the string '%1'")
                            .arg(atts["scale_factor"].toString()) );
                }
            }

            if (atts.contains("units"))
            {
                QString units = atts["units"].toString();

                if (units != "amu*angstrom/picosecond^2")
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Units are not the default of "
                        "amu*angstrom/picosecond^2. They are actually %1!").arg(units) );
                }

                //we need to interpret different units and set the scale_factor accordingly.
                //We want to store the forces here in amu*angstroms/picosecond^2
            }

            if (data.nValues() % 3 != 0)
            {
                QMutexLocker lkr(&warnings_mutex);
                parse_warnings.append( QObject::tr("Could not interpret the forces "
                                                   "as they are not divisible by three (%1)?")
                                                        .arg(data.nValues()) );
            }
            else
            {
                const auto vals = data.toDoubleArray();

                const int natoms_times_nframes = vals.count() / 3;

                if (natoms_times_nframes % natoms != 0)
                {
                    QMutexLocker lkr(&warnings_mutex);
                    parse_warnings.append( QObject::tr("Could not interpret the forces "
                                    "as the amount of data (%1) does not equal 3 x nframes(%2) "
                                    "x natoms(%3)")
                                        .arg(vals.count()).arg(natoms_times_nframes/natoms)
                                        .arg(natoms) );
                    return;
                }

                const int nframes = natoms_times_nframes / natoms;

                frcs = QVector< QVector<Vector> >(nframes);

                if (uses_parallel)
                {
                    for (int i=0; i<nframes; ++i)
                    {
                        QVector<Vector> f(natoms);

                        tbb::parallel_for( tbb::blocked_range<int>(0,natoms),
                                           [&](const tbb::blocked_range<int> &r)
                        {
                            for (int j=r.begin(); j<r.end(); ++j)
                            {
                                const int idx = 3*natoms*i + 3*j;

                                f[j] = scale_factor * Vector( vals[idx],
                                                              vals[idx+1],
                                                              vals[idx+2] );
                            }
                        });

                        frcs[i] = f;
                    }
                }
                else
                {
                    for (int i=0; i<nframes; ++i)
                    {
                        QVector<Vector> f(natoms);

                        for (int j=0; j<natoms; ++j)
                        {
                            const int idx = 3*natoms*i + 3*j;

                            f[j] = scale_factor * Vector( vals[idx], vals[idx+1], vals[idx+2] );
                        }

                        frcs[i] = f;
                    }
                }
            }
        }
    };

    //this is the list of all functions that need to be called to parse the file
    const QVector< std::function<void()> > parse_functions =
            { parseBase, parseCoordinates, parseVelocities, parseForces };

    //call the functions, either in parallel or in serial
    if (uses_parallel)
    {
        tbb::parallel_for( tbb::blocked_range<int>(0,parse_functions.count(),1),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int ifunc=r.begin(); ifunc<r.end(); ++ifunc)
            {
                parse_functions[ifunc]();
            }
        });
    }
    else
    {
        for (int ifunc=0; ifunc<parse_functions.count(); ++ifunc)
        {
            parse_functions[ifunc]();
        }
    }

    //set the score, and save the warnings
    double score = 100.0 / (parse_warnings.count()+1);
    this->setScore(score);

    //assert that this is sane
    this->assertSane();
}

/** Construct by parsing the passed file */
AmberRst::AmberRst(const QString &filename, const PropertyMap &map)
         : ConcreteProperty<AmberRst,MoleculeParser>(map),
           convention_version(0), created_from_restart(false)
{
    //open the NetCDF file and extract all of the data
    NetCDFFile netcdf(filename);
    this->parse(netcdf, map);
}

/** Construct by parsing the data in the passed text lines */
AmberRst::AmberRst(const QStringList &lines, const PropertyMap &map)
         : ConcreteProperty<AmberRst,MoleculeParser>(lines, map),
           convention_version(0), created_from_restart(false)
{
    throw SireIO::parse_error( QObject::tr(
            "You cannot create a binary Amber RST file from a set of text lines!"),
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
AmberRst::AmberRst(const System &system, const PropertyMap &map)
         : ConcreteProperty<AmberRst,MoleculeParser>(),
           convention_version(0), created_from_restart(false)
{
    //get the MolNums of each molecule in the System - this returns the
    //numbers in MolIdx order
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    if (molnums.isEmpty())
    {
        //no molecules in the system
        this->operator=(AmberRst());
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
            coords.append( collapse(all_coords) );
        }

        if (::hasData(all_vels))
        {
            vels.append( collapse(all_vels) );
        }

        if (::hasData(all_forces))
        {
            frcs.append( collapse(all_forces) );
        }
    }

    //extract the space of the system
    box_dims.clear();
    box_angs.clear();

    try
    {
        if (system.property(map["space"]).isA<PeriodicBox>())
        {
            box_dims.append(system.property(map["space"]).asA<PeriodicBox>().dimensions());
            box_angs.append(Vector(90,90,90));
        }
        else if (system.property(map["space"]).isA<TriclinicBox>())
        {
            Vector v0 = system.property(map["space"]).asA<TriclinicBox>().vector0();
            Vector v1 = system.property(map["space"]).asA<TriclinicBox>().vector1();
            Vector v2 = system.property(map["space"]).asA<TriclinicBox>().vector2();

            // Store the box magnitudes.
            box_dims.append(Vector(v0.magnitude(), v1.magnitude(), v2.magnitude()));

            double alpha = system.property(map["space"]).asA<TriclinicBox>().alpha();
            double beta  = system.property(map["space"]).asA<TriclinicBox>().alpha();
            double gamma = system.property(map["space"]).asA<TriclinicBox>().gamma();

            // Store the angles between the box vectors.
            box_angs.append(Vector(alpha, beta, gamma));
        }
    }
    catch(...)
    {}

    //extract the current time for the system
    try
    {
        const auto time = system.property( map["time"] ).asA<TimeProperty>().value();
        current_time.append( time.to(picosecond) );
    }
    catch(...)
    {}

    //extract whether or not this System was created from an Amber restart
    try
    {
        created_from_restart = system.property( map["created_from_restart"] )
                                     .asA<BooleanProperty>().value();
    }
    catch(...)
    {}

    //extract the title for the system
    ttle = system.name().value();

    //this is created by Sire
    creator_app = QString("Sire AmberRst %1").arg(SireBase::getReleaseVersion());

    //assert this is sane
    this->assertSane();
}

/** Copy constructor */
AmberRst::AmberRst(const AmberRst &other)
         : ConcreteProperty<AmberRst,MoleculeParser>(other),
           ttle(other.ttle), current_time(other.current_time),
           coords(other.coords), vels(other.vels), frcs(other.frcs),
           box_dims(other.box_dims), box_angs(other.box_angs),
           convention_version(other.convention_version),
           creator_app(other.creator_app),
           parse_warnings(other.parse_warnings),
           created_from_restart(other.created_from_restart)
{}

/** Destructor */
AmberRst::~AmberRst()
{}

/** Return the number of frames in the file */
int AmberRst::count() const
{
    return this->nFrames();
}

/** Return the number of frames in the file */
int AmberRst::size() const
{
    return this->nFrames();
}

/** Return the AmberRst object that contains only the information for the ith
    frame. This allows you to extract and create a system for the ith frame
    from a trajectory */
AmberRst AmberRst::operator[](int i) const
{
    i = Index(i).map( this->nFrames() );

    if (nFrames() == 1)
        return *this;

    AmberRst ret(*this);

    if (not coords.isEmpty())
    {
        ret.coords = { coords[i] };
    }

    if (not vels.isEmpty())
    {
        ret.vels = { vels[i] };
    }

    if (not frcs.isEmpty())
    {
        ret.frcs = { frcs[i] };
    }

    if (not current_time.isEmpty())
    {
        ret.current_time = { current_time[i] };
    }

    if (not box_dims.isEmpty())
    {
        ret.box_dims = { box_dims[i] };
    }

    if (not box_angs.isEmpty())
    {
        ret.box_angs = { box_angs[i] };
    }

    ret.assertSane();

    return ret;
}

AmberRst& AmberRst::operator=(const AmberRst &other)
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
        convention_version = other.convention_version;
        creator_app = other.creator_app;
        parse_warnings = other.parse_warnings;
        created_from_restart = other.created_from_restart;

        MoleculeParser::operator=(other);
    }

    return *this;
}

bool AmberRst::operator==(const AmberRst &other) const
{
    return MoleculeParser::operator==(other);
}

bool AmberRst::operator!=(const AmberRst &other) const
{
    return not AmberRst::operator==(other);
}

const char* AmberRst::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberRst>() );
}

const char* AmberRst::what() const
{
    return AmberRst::typeName();
}

QString AmberRst::toString() const
{
    if (nAtoms() == 0)
    {
        return QObject::tr("AmberRst( nAtoms() = 0 )");
    }
    else
    {
        return QObject::tr("AmberRst( title() = %1, nAtoms() = %2, nFrames() = %6, "
                "hasCoordinates() = %3, hasVelocities() = %4, "
                "hasForces() = %5 )")
                .arg(title()).arg(nAtoms())
                .arg(hasCoordinates()).arg(hasVelocities()).arg(hasForces())
                .arg(nFrames());
    }
}

/** Return the number of frames that have been loaded from the file */
int AmberRst::nFrames() const
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

/** Parse from the passed file */
AmberRst AmberRst::parse(const QString &filename)
{
    return AmberRst(filename);
}

/** Internal function used to add the data from this parser into the passed System */
void AmberRst::addToSystem(System &system, const PropertyMap &map) const
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
                "Incompatibility between the files, as this restart file contains data "
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
    else if ((not box_dims.isEmpty()) and space_property.hasSource())
    {
        // PeriodicBox.
        if (box_angs[0] == cubic_angs)
        {
            system.setProperty( space_property.source(),
                                SireVol::PeriodicBox(box_dims[0]) );
        }
        // TriclinicBox.
        else
        {
            system.setProperty( space_property.source(),
                                SireVol::TriclinicBox(box_dims[0].x(),
                                                      box_dims[0].y(),
                                                      box_dims[0].z(),
                                                      box_angs[0].x()*degrees,
                                                      box_angs[0].y()*degrees,
                                                      box_angs[0].z()*degrees) );
        }
    }

    PropertyName restart_property = map["created_from_restart"];
    if (restart_property.hasSource())
    {
        system.setProperty( restart_property.source(), BooleanProperty(created_from_restart) );
    }
    else
    {
        system.setProperty("created_from_restart", BooleanProperty(created_from_restart));
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

/** Return the title of the file */
QString AmberRst::title() const
{
    return ttle;
}

/** Return the current time of the simulation from which this restart
    file was written. Returns 0 if there is no time set. If there are
    multiple frames, then the time of the first frame is returned */
double AmberRst::time() const
{
    if (current_time.isEmpty())
        return 0;
    else
        return current_time[0];
}

/** Return the time of the 'ith' frame from the file */
double AmberRst::time(int i) const
{
    return current_time[ Index(i).map(current_time.count()) ];
}

/** Return the number of atoms whose data are contained in this restart file */
int AmberRst::nAtoms() const
{
    if (not coords.isEmpty())
        return coords[0].count();
    else if (not vels.isEmpty())
        return vels[0].count();
    else if (not frcs.isEmpty())
        return frcs[0].count();
    else
        return 0;
}

/** Return whether or not this restart file provides coordinates */
bool AmberRst::hasCoordinates() const
{
    return not coords.isEmpty();
}

/** Return whether or not this restart file also provides velocities */
bool AmberRst::hasVelocities() const
{
    return not vels.isEmpty();
}

/** Return whether or not this restart file also provides forces */
bool AmberRst::hasForces() const
{
    return not frcs.isEmpty();
}

/** Return the parsed coordinate data. If there are multiple frames,
    then only the first frame is returned */
QVector<SireMaths::Vector> AmberRst::coordinates() const
{
    if (coords.isEmpty())
        return QVector<Vector>();
    else
        return coords[0];
}

/** Return the coordinates of the 'ith' frame from the file */
QVector<SireMaths::Vector> AmberRst::coordinates(int i) const
{
    return coords[ Index(i).map(coords.count()) ];
}

/** Return the parsed coordinate data. If there are multiple frames,
    then only the first frame is returned */
QVector<SireMaths::Vector> AmberRst::velocities() const
{
    if (vels.isEmpty())
        return QVector<Vector>();
    else
        return vels[0];
}

/** Return the velocities of the 'ith' frame from the file */
QVector<SireMaths::Vector> AmberRst::velocities(int i) const
{
    return vels[ Index(i).map(vels.count()) ];
}

/** Return the parsed force data. If there are multiple frames,
    then only the first frame is returned */
QVector<SireMaths::Vector> AmberRst::forces() const
{
    if (frcs.isEmpty())
        return QVector<Vector>();
    else
        return frcs[0];
}

/** Return the forces of the 'ith' frame from the file */
QVector<SireMaths::Vector> AmberRst::forces(int i) const
{
    return frcs[ Index(i).map(frcs.count()) ];
}

/** Return the parsed box dimensions, or Vector(0) if there is no space information.
    If there are multiple frames, then only the first frame is returned */
SireMaths::Vector AmberRst::boxDimensions() const
{
    if (box_dims.isEmpty())
        return Vector(0);
    else
        return box_dims[0];
}

/** Return the box dimensions of the 'ith
' frame from the file */
SireMaths::Vector AmberRst::boxDimensions(int i) const
{
    return box_dims[ Index(i).map(box_dims.count()) ];
}

/** Return the parsed box angles, or Vector(0) if there is no space information.
    If there are multiple frames, then only the first frame is returned */
SireMaths::Vector AmberRst::boxAngles() const
{
    if (box_angs.isEmpty())
        return Vector(0);
    else
        return box_angs[0];
}

/** Return the box angles of the 'ith' frame from the file */
SireMaths::Vector AmberRst::boxAngles(int i) const
{
    return box_angs[ Index(i).map(box_angs.count()) ];
}

/** Return any warnings that were triggered during parsing */
QStringList AmberRst::warnings() const
{
    return parse_warnings;
}

/** Return the application that created the file that has been parsed */
QString AmberRst::creatorApplication() const
{
    return creator_app;
}

/** Return the version of the file format that was parsed */
double AmberRst::formatVersion() const
{
    return convention_version;
}

/** Return whether or not this was created from a restart (.rst) file */
bool AmberRst::createdFromRestart() const
{
    return created_from_restart;
}

/** Return whether or not this was created from a trajectory (.trj) file */
bool AmberRst::createdFromTrajectory() const
{
    return not createdFromRestart();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr AmberRst::construct(const QString &filename,
                                      const PropertyMap &map) const
{
    //don't construct from a pointer as it could leak
    return MoleculeParserPtr( AmberRst(filename,map) );
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr AmberRst::construct(const QStringList &lines,
                                      const PropertyMap &map) const
{
    //don't construct from a pointer as it could leak
    return MoleculeParserPtr( AmberRst(lines,map) );
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr AmberRst::construct(const SireSystem::System &system,
                                      const PropertyMap &map) const
{
    //don't construct from a pointer as it could leak
    return MoleculeParserPtr( AmberRst(system,map) );
}

/** Write this AmberRst to a file called 'filename'. This will write out
    the data in this object to the Amber NetCDF format */
void AmberRst::writeToFile(const QString &filename) const
{
    //create the hash of all global data
    QHash<QString,QString> globals;

    if (created_from_restart)
    {
        globals.insert( "Conventions", "AMBERRESTART" );
    }
    else
    {
        globals.insert( "Conventions", "AMBER" );
    }

    globals.insert( "ConventionVersion", "1.0" );

    globals.insert( "application", "Sire" );

    globals.insert( "program", "AmberRst" );

    globals.insert( "programVersion", SireBase::getReleaseVersion() );

    if (ttle.count() > 80)
    {
        globals.insert( "title", ttle.mid(0,80) );
    }
    else if (not ttle.isEmpty())
    {
        globals.insert( "title", ttle );
    }

    //convert all of the data into NetCDFData objects
    QHash<QString,NetCDFData> data;

    //start off with the label variables
    {
        QStringList dimensions = { "spatial" };
        QList<int> dimension_sizes = { 3 };
        QVector<char> values = { 'x', 'y', 'z' };
        data.insert( "spatial", NetCDFData("spatial", values, dimensions, dimension_sizes) );

        QStringList dimensions2 = { "cell_spatial" };
        dimension_sizes = { 3 };
        values = { 'a', 'b', 'c' };
        data.insert( "cell_spatial", NetCDFData("cell_spatial", values,
                                                dimensions2, dimension_sizes) );

        QStringList dimensions3 = { "cell_angular", "label" };
        dimension_sizes = { 3, 5 };
        values = { 'a', 'l', 'p', 'h', 'a',
                   'b', 'e', 't', 'a', ' ',
                   'g', 'a', 'm', 'm', 'a' };
        data.insert( "cell_angular", NetCDFData("cell_angular", values,
                                                dimensions3, dimension_sizes) );
    }

    //now the time
    if (not current_time.isEmpty())
    {
        QHash<QString,QVariant> attributes;

        attributes.insert( "units", QString("picosecond") );

        if (created_from_restart)
        {
            QVector<double> values = { current_time[0] };
            data.insert( "time", NetCDFData("time", values, QStringList(),
                                            QList<int>(), attributes ) );
        }
        else
        {
            QVector<float> values(current_time.count());

            for (int i=0; i<current_time.count(); ++i)
            {
                values[i] = current_time[i];
            }

            QStringList dimensions = { "frame" };
            QList<int> dimension_sizes = { current_time.count() };

            data.insert( "time", NetCDFData("time", values, dimensions,
                                            dimension_sizes, attributes ) );
        }
    }

    QHash<QString,QVariant> angstrom_units;
    angstrom_units.insert("units", "angstrom");

    //coordinates
    if (not coords.isEmpty())
    {
        if (created_from_restart)
        {
            QStringList dimensions = { "atom", "spatial" };
            QList<int> dimension_sizes = { coords[0].count(), 3 };
            QVector<double> values( coords[0].count() * 3 );

            for (int i=0; i<coords[0].count(); ++i)
            {
                const Vector &c = coords[0].constData()[i];

                values[3*i + 0] = c.x();
                values[3*i + 1] = c.y();
                values[3*i + 2] = c.z();
            }

            data.insert("coordinates", NetCDFData("coordinates", values,
                                                  dimensions, dimension_sizes,
                                                  angstrom_units));
        }
        else
        {
            QStringList dimensions = { "frame", "atom", "spatial" };
            QList<int> dimension_sizes = { coords.count(), coords[0].count(), 3 };
            QVector<float> values( coords.count() * coords[0].count() * 3 );

            for (int i=0; i<coords.count(); ++i)
            {
                for (int j=0; j<coords[0].count(); ++j)
                {
                    const Vector &c = coords[i].constData()[j];

                    const int idx = 3*coords[0].count()*i + 3*j;

                    values[idx + 0] = c.x();
                    values[idx + 1] = c.y();
                    values[idx + 2] = c.z();
                }
            }

            data.insert("coordinates", NetCDFData("coordinates", values,
                                                  dimensions, dimension_sizes,
                                                  angstrom_units));
        }
    }

    //cell lengths and angles
    if (not box_dims.isEmpty())
    {
        QHash<QString,QVariant> degree_units;
        degree_units.insert("units", "degree");

        if (created_from_restart)
        {
            QStringList dimensions = { "cell_spatial" };
            QList<int> dimension_sizes = { 3 };
            QVector<double> values = { box_dims[0].x(), box_dims[0].y(), box_dims[0].z() };

            data.insert( "cell_lengths", NetCDFData("cell_lengths", values,
                                                    dimensions, dimension_sizes,
                                                    angstrom_units) );

            QStringList dimensions2 = { "cell_angular" };
            values = { box_angs[0].x(), box_angs[0].y(), box_angs[0].z() };

            data.insert( "cell_angles", NetCDFData("cell_angles", values,
                                                    dimensions2, dimension_sizes,
                                                    degree_units) );
        }
        else
        {
            QStringList dimensions = { "frame", "cell_spatial" };
            QList<int> dimension_sizes = { box_dims.count(), 3 };
            QVector<float> values( 3 * box_dims.count() );

            for (int i=0; i<box_dims.count(); ++i)
            {
                values[3*i + 0] = box_dims[i].x();
                values[3*i + 1] = box_dims[i].y();
                values[3*i + 2] = box_dims[i].z();
            }

            data.insert( "cell_lengths", NetCDFData("cell_lengths", values,
                                                    dimensions, dimension_sizes,
                                                    angstrom_units) );

            dimensions[1] = "cell_angular";

            for (int i=0; i<box_dims.count(); ++i)
            {
                values[3*i + 0] = box_angs[i].x();
                values[3*i + 1] = box_angs[i].y();
                values[3*i + 2] = box_angs[i].z();
            }

            data.insert( "cell_angles", NetCDFData("cell_angles", values,
                                                    dimensions, dimension_sizes,
                                                    degree_units) );
        }
    }

    //velocities
    if (not vels.isEmpty())
    {
        QHash<QString,QVariant> vel_attributes;
        vel_attributes.insert("units", "angstrom/picosecond");

        if (created_from_restart)
        {
            vel_attributes.insert("scale_factor", double(20.455));

            QStringList dimensions = { "atom", "spatial" };
            QList<int> dimension_sizes = { vels[0].count(), 3 };
            QVector<double> values( vels[0].count() * 3 );

            for (int i=0; i<vels[0].count(); ++i)
            {
                const Vector &v = vels[0].constData()[i];

                values[3*i + 0] = v.x();
                values[3*i + 1] = v.y();
                values[3*i + 2] = v.z();
            }

            data.insert("velocities", NetCDFData("velocities", values,
                                                  dimensions, dimension_sizes,
                                                  vel_attributes));
        }
        else
        {
            vel_attributes.insert("scale_factor", float(20.455));
            QStringList dimensions = { "frame", "atom", "spatial" };
            QList<int> dimension_sizes = { vels.count(), vels[0].count(), 3 };
            QVector<float> values( vels.count() * vels[0].count() * 3 );

            for (int i=0; i<vels.count(); ++i)
            {
                for (int j=0; j<vels[0].count(); ++j)
                {
                    const Vector &v = vels[i].constData()[j];

                    const int idx = 3*i*vels[0].count() + 3*j;
                    values[idx + 0] = v.x();
                    values[idx + 1] = v.y();
                    values[idx + 2] = v.z();
                }
            }

            data.insert("velocities", NetCDFData("velocities", values,
                                                  dimensions, dimension_sizes,
                                                  vel_attributes));
        }
    }

    //forces
    if (not frcs.isEmpty())
    {
        QHash<QString,QVariant> force_units;
        force_units.insert("units", "amu*angstrom/picosecond^2");

        if (created_from_restart)
        {
            QStringList dimensions = { "atom", "spatial" };
            QList<int> dimension_sizes = { frcs[0].count(), 3 };
            QVector<double> values( frcs[0].count() * 3 );

            for (int i=0; i<frcs[0].count(); ++i)
            {
                const Vector &f = frcs[0].constData()[i];

                values[3*i + 0] = f.x();
                values[3*i + 1] = f.y();
                values[3*i + 2] = f.z();
            }

            data.insert("forces", NetCDFData("forces", values,
                                             dimensions, dimension_sizes,
                                             force_units));
        }
        else
        {
            QStringList dimensions = { "frame", "atom", "spatial" };
            QList<int> dimension_sizes = { frcs.count(), frcs[0].count(), 3 };
            QVector<float> values( frcs.count() * frcs[0].count() * 3 );

            for (int i=0; i<frcs.count(); ++i)
            {
                for (int j=0; j<frcs[0].count(); ++j)
                {
                    const Vector &f = frcs[i].constData()[j];

                    const int idx = 3*frcs[0].count()*i + 3*j;

                    values[idx + 0] = f.x();
                    values[idx + 1] = f.y();
                    values[idx + 2] = f.z();
                }
            }

            data.insert("forces", NetCDFData("forces", values,
                                             dimensions, dimension_sizes,
                                             force_units));
        }
    }

    NetCDFFile::write(filename, globals, data);
}
