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
#include "SireMol/trajectory.h"
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
#include <QFileInfo>
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


//// Thanks to the MDAnalysis parser
//// (https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/lib/formats/include/readdcd.h)
//// which really helped with the reverse engineering of the DCD fileformat

static const RegisterMetaType<DCD> r_dcd;
const RegisterParser<DCD> register_dcd;

QDataStream &operator<<(QDataStream &ds, const DCD &dcd)
{
    writeHeader(ds, r_dcd, 1);

    SharedDataStream sds(ds);

    sds << dcd.coords
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

        sds >> dcd.coords
            >> dcd.parse_warnings
            >> static_cast<MoleculeParser&>(dcd);

        try
        {
            dcd.dcd = DCDFile();
            FortranFile file(dcd.filename());
            dcd.dcd.readHeader(file);
        }
        catch(SireError::exception &e)
        {
            qDebug() << "WARNING: Failed to reload DCD file" << dcd.filename();
            qDebug() << e.what() << ":" << e.error();
            dcd.dcd = DCDFile();
        }
    }
    else
        throw version_error(v, "1", r_dcd, CODELOC);

    return ds;
}

static Vector cubic_angs(90,90,90);

/** Constructor */
DCD::DCD() : ConcreteProperty<DCD,MoleculeParser>()
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

SireIO::detail::DCDFile::DCDFile(const QString &filename)
               : timestep(0), istart(0),
                 nsavc(0), nfixed(0), natoms(0),
                 nframes(0), first_frame_line(0),
                 CHARMM_FORMAT(false),
                 HAS_EXTRA_BLOCK(false),
                 HAS_FOUR_DIMS(false)
{
    FortranFile file(filename);
    this->readHeader(file);
}

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

    // read the timestep between frames (it is assumed to be in picoseconds)
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

    // now read in the space
    spc = this->readSpace(file, 0);

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

    if (nframes != 0)
    {
        if (file.nRecords() != first_frame_line + (num_lines_per_frame * nframes))
        {
            throw SireIO::parse_error(QObject::tr(
                "Wrong number of records in the DCD file. Expect to have %1 "
                "for %2 frames, but actually have %3.")
                    .arg(first_frame_line + (num_lines_per_frame*nframes))
                    .arg(nframes).arg(file.nRecords()), CODELOC);
        }
    }
    else
    {
        // we need to calculate nframes
        nframes = (file.nRecords() - first_frame_line) / num_lines_per_frame;
    }
}

double SireIO::detail::DCDFile::getTimeAtFrame(int frame) const
{
    return (istart*timestep) + (frame*timestep);
}

double SireIO::detail::DCDFile::getCurrentTime() const
{
    return getTimeAtFrame(0);
}

void SireIO::detail::DCDFile::setCurrentTime(double time)
{
    if (timestep != 0)
    {
        istart = int(time / timestep);
    }
    else
    {
        timestep = time;
        istart = 1;
    }
}

void SireIO::detail::DCDFile::setSpace(const Space &s)
{
    spc = s;
}

const Space& SireIO::detail::DCDFile::getSpace() const
{
    return *spc;
}

SpacePtr SireIO::detail::DCDFile::readSpace(FortranFile &file, int frame) const
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

        auto boxinfo = line.readFloat64(6);

        // construct from the above boxinfo
        return SpacePtr();
    }

    return SpacePtr();
}

QVector<Vector> SireIO::detail::DCDFile::readCoordinates(FortranFile &file, int frame) const
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

Frame SireIO::detail::DCDFile::readFrame(FortranFile &file, int frame) const
{
    frame = SireID::Index(frame).map(nframes);

    auto space = this->readSpace(file, frame);
    auto coords = this->readCoordinates(file, frame);

    return Frame(coords, space, getTimeAtFrame(frame)*picosecond);
}

QString SireIO::detail::DCDFile::getTitle() const
{
    return title.join("");
}

void SireIO::detail::DCDFile::setTitle(QString t)
{
    // need to split into blocks of 32 characters
    title.clear();

    while(t.length() > 32)
    {
        title.append(t.mid(0, 32));
        t.remove(0, 32);
    }

    if (t.length() > 0)
    {
        title.append(t);
    }
}

double SireIO::detail::DCDFile::getTimeStep() const
{
    return timestep;
}

qint64 SireIO::detail::DCDFile::getFrameStart() const
{
    return istart;
}

qint64 SireIO::detail::DCDFile::getFrameDelta() const
{
    return nsavc;
}

qint64 SireIO::detail::DCDFile::nAtoms() const
{
    return natoms;
}

qint64 SireIO::detail::DCDFile::nFrames() const
{
    return nframes;
}

/** Parse the data contained in the lines - this clears any pre-existing
    data in this object */
void DCD::parse(const QString &filename, const PropertyMap &map)
{
    FortranFile file(filename);

    dcd = DCDFile();
    dcd.readHeader(file);

    coords = dcd.readCoordinates(file, 0);

    //need to convert unit cell...

    //set the score, and save the warnings
    double score = 100.0 / (parse_warnings.count()+1);
    this->setScore(score);
}

/** Construct by parsing the passed file */
DCD::DCD(const QString &fname, const PropertyMap &map)
    : ConcreteProperty<DCD,MoleculeParser>(map)
{
    MoleculeParser::setFilename(fname);
    this->parse(MoleculeParser::filename(), map);
}

/** Construct by parsing the data in the passed text lines */
DCD::DCD(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<DCD,MoleculeParser>(lines, map)
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
    : ConcreteProperty<DCD,MoleculeParser>()
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

        const auto coords_property = map["coordinates"];

        if (usesParallel())
        {
            tbb::parallel_for( tbb::blocked_range<int>(0,molnums.count()),
                               [&](const tbb::blocked_range<int> r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto mol = system[molnums[i]].molecule();

                    all_coords[i] = ::getCoordinates(mol, coords_property);
                }
            });
        }
        else
        {
            for (int i=0; i<molnums.count(); ++i)
            {
                const auto mol = system[molnums[i]].molecule();
                all_coords[i] = ::getCoordinates(mol, coords_property);
            }
        }

        coords.clear();

        if (::hasData(all_coords))
        {
            coords = collapse(all_coords);
        }
    }

    //extract the space of the system
    SpacePtr space;

    if (system.containsProperty(map["space"]))
    {
        space = system.property(map["space"]).asA<Space>();
    }

    //extract the current time for the system
    double current_time = 0;

    if (system.containsProperty(map["time"]))
    {
        const Property &prop = system.property( map["time"] );

        Time time;

        if (prop.isA<TimeProperty>())
            time = prop.asA<TimeProperty>().value();
        else
            time = prop.asA<GeneralUnitProperty>();

        current_time = time.to(picosecond);
    }

    dcd.setCurrentTime(current_time);

    //extract the title for the system
    dcd.setTitle(system.name().value());
}

/** Copy constructor */
DCD::DCD(const DCD &other)
    : ConcreteProperty<DCD,MoleculeParser>(other),
      coords(other.coords), dcd(other.dcd),
      parse_warnings(other.parse_warnings)
{}

/** Destructor */
DCD::~DCD()
{}

DCD& DCD::operator=(const DCD &other)
{
    if (this != &other)
    {
        coords = other.coords;
        dcd = other.dcd;
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
        return QObject::tr("DCD( title() = %1, nAtoms() = %2, nFrames() = %3 )")
                .arg(title()).arg(nAtoms()).arg(nFrames());
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
    if (coords.isEmpty())
        return;

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

    const Vector *coords_array = this->coordinates().constData();

    auto add_moldata = [&](int i)
    {
        const int atom_start_idx = atom_pointers.constData()[i];
        auto mol = allmols[MolIdx(i)].molecule();
        const auto molinfo = mol.data().info();

        auto moleditor = mol.edit();

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

        moleditor.setProperty(coords_property, AtomCoords(CoordGroupArray(coords)));
        mols_array[i] = moleditor.commit();
    };

    if (coords_property.hasSource())
    {
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
    }

    PropertyName space_property = map["space"];
    if (space_property.hasSource())
    {
        system.setProperty(space_property.source(), dcd.getSpace());
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

    PropertyName time_property = map["time"];

    if (time_property.hasSource())
    {
        system.setProperty(time_property.source(), GeneralUnitProperty(dcd.getCurrentTime()*picosecond));
    }
}

/** Return the title of the file */
QString DCD::title() const
{
    return dcd.getTitle();
}

/** Return the current time of the simulation from which this restart
    file was written. Returns 0 if there is no time set. If there are
    multiple frames, then the time of the first frame is returned */
SireUnits::Dimension::Time DCD::time() const
{
    return dcd.getCurrentTime() * picosecond;
}

/** Return the number of atoms whose data are contained in this DCD file */
int DCD::nAtoms() const
{
    return coords.count();
}

bool DCD::isFrame() const
{
    return true;
}

/** Return the number of frames in this DCD file */
int DCD::nFrames() const
{
    return dcd.nFrames();
}

/** Return the ith frame */
Frame DCD::getFrame(int i) const
{
    //will eventually look to see if we should cache this?
    QString f = this->filename();

    if (f.isEmpty())
        throw SireIO::parse_error(QObject::tr(
            "Cannot get the frame as the DCD filename has not been specified."),
                CODELOC);

    FortranFile file(this->filename());
    return dcd.readFrame(file, i);
}

/** Return the parsed coordinate data. */
QVector<SireMaths::Vector> DCD::coordinates() const
{
    return coords;
}

/** Return the parsed space */
const Space& DCD::space() const
{
    return dcd.getSpace();
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
