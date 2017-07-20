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

#include "SireVol/periodicbox.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/timeproperty.h"

#include "SireIO/errors.h"

#include "SireStream/shareddatastream.h"

using namespace SireIO;
using namespace SireIO::detail;
using namespace SireMaths;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireUnits;
using namespace SireStream;

static const RegisterMetaType<AmberRst> r_rst;
const RegisterParser<AmberRst> register_amberrst;

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const AmberRst &rst)
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

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, AmberRst &rst)
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
         : ConcreteProperty<AmberRst,MoleculeParser>(),
           current_time(-1), box_dims(0), box_angs(cubic_angs)
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

    //now lets create some lambda functions that do the different orthogonal
    //work needed to parse the file. These functions will either be called serially
    //or in parallel depending on whether this is supported or requested
    QMutex warnings_mutex;

    //function to parse in all of the base data in the file
    auto parseBase = [&]()
    {
        //extract the current time of the simulation
        current_time = -1;
        
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
                current_time = scale_factor * data.toDoubleArray()[0];
            }
        }

        //now lets get any periodic box information
        box_dims = Vector(0);
        
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
            else
            {
                const auto vals = data.toDoubleArray();
                box_dims = scale_factor * Vector( vals[0], vals[1], vals[2] );
            }
        }

        //now lets get any periodic box angle information
        box_angs = Vector(0);
        
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
            else
            {
                const auto vals = data.toDoubleArray();
                box_angs = scale_factor * Vector( vals[0], vals[1], vals[2] );
            }
        }
    };
    
    //function to read in all of the coordinates
    auto parseCoordinates = [&]()
    {
        coords = QVector<Vector>();
        
        if (varinfos.contains("coordinates"))
        {
            const auto data = netcdf.read( varinfos["coordinates"] );
            
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
                    parse_warnings.append( QObject::tr("Units are not the default of angstroms. "
                        "They are actually %1!").arg(units) );
                }
                
                //we need to interpret different units and set the scale_factor accordingly.
                //We want to store the lengths here in angstroms
            }

            if (data.nValues() % 3 != 0)
            {
                parse_warnings.append( QObject::tr("Could not interpret the coordinates "
                                                   "as they are not divisible by three (%1)?")
                                                        .arg(data.nValues()) );
            }
            else
            {
                const auto vals = data.toDoubleArray();
                
                int nats = vals.count() / 3;
                
                coords = QVector<Vector>(nats);
                
                if (uses_parallel)
                {
                    tbb::parallel_for( tbb::blocked_range<int>(0,nats),
                                       [&](const tbb::blocked_range<int> &r)
                    {
                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            int idx = 3*i;
                            
                            coords[i] = scale_factor * Vector( vals[idx],
                                                               vals[idx+1],
                                                               vals[idx+2] );
                        }
                    });
                }
                else
                {
                    for (int i=0; i<nats; ++i)
                    {
                        int idx = 3 * i;
                    
                        coords[i] = scale_factor * Vector( vals[idx], vals[idx+1], vals[idx+2] );
                    }
                }
            }
        }
    };

    //function to parse all of the velocities
    auto parseVelocities = [&]()
    {
        vels = QVector<Vector>();
        
        if (varinfos.contains("velocities"))
        {
            const auto data = netcdf.read( varinfos["velocities"] );
            
            const auto atts = data.attributes();
            
            double scale_factor = 1.0;
            
            if (atts.contains("scale_factor"))
            {
                scale_factor = atts["scale_factor"].toDouble();
            }
            
            if (atts.contains("units"))
            {
                QString units = atts["units"].toString();
                
                if (units != "angstrom/picosecond")
                {
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
                parse_warnings.append( QObject::tr("Could not interpret the velocities "
                                                   "as they are not divisible by three (%1)?")
                                                        .arg(data.nValues()) );
            }
            else
            {
                const auto vals = data.toDoubleArray();
                
                int nats = vals.count() / 3;
                
                vels = QVector<Vector>(nats);
                
                if (uses_parallel)
                {
                    tbb::parallel_for( tbb::blocked_range<int>(0,nats),
                                       [&](const tbb::blocked_range<int> &r)
                    {
                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            int idx = 3*i;
                            vels[i] = scale_factor * Vector( vals[idx],
                                                             vals[idx+1],
                                                             vals[idx+2] );
                        }
                    });
                }
                else
                {
                    for (int i=0; i<nats; ++i)
                    {
                        int idx = 3 * i;
                        vels[i] = scale_factor * Vector( vals[idx], vals[idx+1], vals[idx+2] );
                    }
                }
            }
        }
    };

    //function to parse all of the forces
    auto parseForces = [&]()
    {
        frcs = QVector<Vector>();
        
        if (varinfos.contains("forces"))
        {
            const auto data = netcdf.read( varinfos["forces"] );
            
            const auto atts = data.attributes();
            
            double scale_factor = 1.0;
            
            if (atts.contains("scale_factor"))
            {
                scale_factor = atts["scale_factor"].toDouble();
            }
            
            if (atts.contains("units"))
            {
                QString units = atts["units"].toString();
                
                if (units != "amu*angstrom/picosecond^2")
                {
                    parse_warnings.append( QObject::tr("Units are not the default of "
                        "amu*angstrom/picosecond^2. They are actually %1!").arg(units) );
                }
                
                //we need to interpret different units and set the scale_factor accordingly.
                //We want to store the forces here in amu*angstroms/picosecond^2
            }

            if (data.nValues() % 3 != 0)
            {
                parse_warnings.append( QObject::tr("Could not interpret the forces "
                                                   "as they are not divisible by three (%1)?")
                                                        .arg(data.nValues()) );
            }
            else
            {
                const auto vals = data.toDoubleArray();
                
                int nats = vals.count() / 3;
                
                frcs = QVector<Vector>(nats);
                
                if (uses_parallel)
                {
                    tbb::parallel_for( tbb::blocked_range<int>(0,nats),
                                       [&](const tbb::blocked_range<int> &r)
                    {
                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            int idx = 3*i;
                            frcs[i] = scale_factor * Vector( vals[idx],
                                                             vals[idx+1],
                                                             vals[idx+2] );
                        }
                    });
                }
                else
                {
                    for (int i=0; i<nats; ++i)
                    {
                        int idx = 3 * i;
                        frcs[i] = scale_factor * Vector( vals[idx], vals[idx+1], vals[idx+2] );
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
            for (int i=r.begin(); i<r.end(); ++i)
            {
                parse_functions[i]();
            }
        });
    }
    else
    {
        for (int i=0; i<parse_functions.count(); ++i)
        {
            parse_functions[i]();
        }
    }

    //set the score, and save the warnings
    double score = 100.0 / (parse_warnings.count()+1);
    this->setScore(score);
}

/** Construct by parsing the passed file */
AmberRst::AmberRst(const QString &filename, const PropertyMap &map)
         : ConcreteProperty<AmberRst,MoleculeParser>(map),
           current_time(0), box_dims(0), box_angs(cubic_angs),
           convention_version(0), created_from_restart(false)
{
    //open the NetCDF file and extract all of the data
    NetCDFFile netcdf(filename);
    this->parse(netcdf, map);
}

/** Construct by parsing the data in the passed text lines */
AmberRst::AmberRst(const QStringList &lines, const PropertyMap &map)
         : ConcreteProperty<AmberRst,MoleculeParser>(lines, map),
           current_time(0), box_dims(0), box_angs(cubic_angs),
           convention_version(0), created_from_restart(false)
{
    throw SireIO::parse_error( QObject::tr(
            "You cannot create a binary Amber RST file from a set of text lines!"),
                CODELOC );
}

/** Construct by extracting the necessary data from the passed System */
AmberRst::AmberRst(const System &system, const PropertyMap &map)
         : ConcreteProperty<AmberRst,MoleculeParser>(),
           current_time(0), box_dims(0), box_angs(cubic_angs),
           convention_version(0), created_from_restart(false)
{
    //DO SOMETHING
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
        return QObject::tr("AmberRst( title() = %1, nAtoms() = %2, "
                "hasCoordinates() = %3, hasVelocities() = %4, "
                "hasForces() = %5 )")
                .arg(title()).arg(nAtoms())
                .arg(hasCoordinates()).arg(hasVelocities()).arg(hasForces());
    }
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
    else if (box_dims != Vector(0) and space_property.hasSource())
    {
        if (box_angs != cubic_angs)
        {
            throw SireIO::parse_error( QObject::tr(
                    "Sire cannot currently support a non-cubic periodic box! %1")
                        .arg(box_angs.toString()), CODELOC );
        }
        
        
        system.setProperty( space_property.source(), SireVol::PeriodicBox(box_dims) );
    }
    
    //update the System fileformat property to record that it includes
    //data from this file format
    QString fileformat = this->formatName();
    
    try
    {
        QString last_format = system.property("fileformat").asA<StringProperty>();
        fileformat = QString("%1,%2").arg(last_format,fileformat);
    }
    catch(...)
    {}
    
    system.setProperty( "fileformat", StringProperty(fileformat) );
    
    if (current_time >= 0)
    {
        system.setProperty( map["time"].source(), TimeProperty(current_time*picosecond) );
    }
}

/** Return the title of the file */
QString AmberRst::title() const
{
    return ttle;
}

/** Return the current time of the simulation from which this restart
    file was written */
double AmberRst::time() const
{
    return current_time;
}

/** Return the number of atoms whose data are contained in this restart file */
int AmberRst::nAtoms() const
{
    return qMax( coords.count(), qMax( vels.count(), frcs.count() ) );
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

/** Return the parsed coordinate data */
QVector<SireMaths::Vector> AmberRst::coordinates() const
{
    return coords;
}

/** Return the parsed coordinate data */
QVector<SireMaths::Vector> AmberRst::velocities() const
{
    return vels;
}

/** Return the parsed force data */
QVector<SireMaths::Vector> AmberRst::forces() const
{
    return frcs;
}

/** Return the parsed box dimensions */
SireMaths::Vector AmberRst::boxDimensions() const
{
    return box_dims;
}

/** Return the parsed box angles */
SireMaths::Vector AmberRst::boxAngles() const
{
    return box_angs;
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
