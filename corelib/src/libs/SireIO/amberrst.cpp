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
#include "SireMol/moleditor.h"

#include "SireVol/periodicbox.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"

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
        << rst.coords << rst.vels
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
            >> rst.coords >> rst.vels
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
           current_time(0), box_dims(0), box_angs(cubic_angs)
{}

/** Return the format name that is used to identify this file format within Sire */
QString AmberRst::formatName() const
{
    return "RST";
}

/** Return the suffixes that AmberRst files will typically use */
QStringList AmberRst::formatSuffix() const
{
    static const QStringList suffixes = { "rst", "crd" };
    return suffixes;
}

/** Return a description of the file format */
QString AmberRst::formatDescription() const
{
    return QObject::tr("Amber coordinate/velocity text (binary, netcdf) restart files "
                       "supported since Amber 9, now default since Amber 16.");
}

/** Parse the data contained in the lines - this clears any pre-existing
    data in this object */
void AmberRst::parse(const NetCDFFile &netcdf, const PropertyMap &map)
{
    //the netcdf conventions for Amber restart/trajectory files are
    //given here - http://ambermd.org/netcdf/nctraj.xhtml

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
        throw SireError::io_error( QObject::tr(
                "This does not look like a valid NetCDF Amber restart/trajectory file. "
                "Such files have a 'Conventions' attribute that contains either "
                "'AMBER' or 'AMBERRESTART'. The 'Conventions' in this file equals '%1'")
                    .arg(conventions), CODELOC );
    }

    bool ok;
    convention_version = netcdf.getStringAttribute("ConventionVersion").toDouble(&ok);

    if (not ok)
    {
        throw SireError::io_error( QObject::tr(
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

    //extract the current time of the simulation
    current_time = 0;
    
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
                parse_warnings.append( QObject::tr("Units are not the default of picoseconds. "
                    "They are actually %1!").arg(units) );
            }
            
            //we need to interpret different units and set the scale_factor accordingly.
            //We want to store the time here in picoseconds
        }

        if (data.nValues() < 1)
        {
            parse_warnings.append( QObject::tr("Could not interpret the time "
                                               "as there is no value?") );
        }
        else
        {
            current_time = scale_factor * data.toDoubleArray()[0];
        }
    }
    
    //now lets read in the coordinates
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
            
            for (int i=0; i<nats; ++i)
            {
                int idx = 3 * i;
            
                coords[i] = scale_factor * Vector( vals[idx], vals[idx+1], vals[idx+2] );
            }
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
                parse_warnings.append( QObject::tr("Units are not the default of angstroms. "
                    "They are actually %1!").arg(units) );
            }
            
            //we need to interpret different units and set the scale_factor accordingly.
            //We want to store the lengths here in angstroms
        }

        if (data.nValues() < 3)
        {
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
                parse_warnings.append( QObject::tr("Units are not the default of degrees. "
                    "They are actually %1!").arg(units) );
            }
            
            //we need to interpret different units and set the scale_factor accordingly.
            //We want to store the angles here in degrees
        }

        if (data.nValues() < 3)
        {
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
    throw SireError::io_error( QObject::tr(
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
           coords(other.coords), vels(other.vels),
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
    if (coords.isEmpty())
    {
        return QObject::tr("AmberRst::null");
    }
    else if (vels.isEmpty())
    {
        return QObject::tr("AmberRst( title() = %1, nAtoms() = %2, hasVelocities() = false )")
                .arg(title()).arg(nAtoms());
    }
    else
    {
        return QObject::tr("AmberRst( title() = %1, nAtoms() = %2, hasVelocities() = true )")
                .arg(title()).arg(nAtoms());
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
    const Vector *coords_array = this->coordinates().constData();
    
    const PropertyName coords_property = map["coordinates"];
    
    if (this->hasVelocities())
    {
        const PropertyName vels_property = map["velocity"];
        const Vector *vels_array = this->velocities().constData();
    
        auto add_coords_and_vels = [&](int i)
        {
            const int atom_start_idx = atom_pointers.constData()[i];
            auto mol = allmols[MolIdx(i)].molecule();
            const auto molinfo = mol.data().info();
            
            //create space for the coordinates and velocities
            auto coords = QVector< QVector<Vector> >(molinfo.nCutGroups());
            auto vels = AtomVelocities(molinfo);

            for (int j=0; j<molinfo.nCutGroups(); ++j)
            {
                coords[j] = QVector<Vector>(molinfo.nAtoms(CGIdx(j)));
            }
            
            for (int j=0; j<mol.nAtoms(); ++j)
            {
                auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(j));
                
                const int atom_idx = atom_start_idx + j;
                
                coords[cgatomidx.cutGroup()][cgatomidx.atom()] = coords_array[atom_idx];

                //velocity is Angstroms per 1/20.455 ps
                const auto vel_unit = (1.0 / 20.455) * angstrom / picosecond;
                
                const Vector &vel = vels_array[atom_idx];
                vels.set(cgatomidx, Velocity3D(vel.x() * vel_unit,
                                               vel.y() * vel_unit,
                                               vel.z() * vel_unit));
            }
            
            mols_array[i] = mol.edit()
                            .setProperty(vels_property, vels)
                            .setProperty(coords_property,
                                         AtomCoords(CoordGroupArray(coords)))
                            .commit();
        };
    
        if (usesParallel())
        {
            tbb::parallel_for( tbb::blocked_range<int>(0,nmols),
                               [&](tbb::blocked_range<int> r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    add_coords_and_vels(i);
                }
            });
        }
        else
        {
            for (int i=0; i<nmols; ++i)
            {
                add_coords_and_vels(i);
            }
        }
    }
    else
    {
        auto add_coords = [&](int i)
        {
            const int atom_start_idx = atom_pointers.constData()[i];
            auto mol = system[MolIdx(i)].molecule();
            const auto molinfo = mol.data().info();
            
            //create space for the coordinates
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
            
            mols_array[i] = mol.edit()
                            .setProperty(coords_property,
                                         AtomCoords(CoordGroupArray(coords)))
                            .commit();
        };
    
        if (usesParallel())
        {
            tbb::parallel_for( tbb::blocked_range<int>(0,nmols),
                               [&](tbb::blocked_range<int> r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    add_coords(i);
                }
            });
        }
        else
        {
            for (int i=0; i<nmols; ++i)
            {
                add_coords(i);
            }
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

/** Return the number of atoms whose coordinates are contained in this restart file */
int AmberRst::nAtoms() const
{
    return coords.count();
}

/** Return whether or not this restart file also provides velocities */
bool AmberRst::hasVelocities() const
{
    return not vels.isEmpty();
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
