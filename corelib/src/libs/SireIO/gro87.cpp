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

#include "SireSystem/system.h"

#include "SireBase/parallel.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QRegularExpression>
#include <QDebug>

using namespace SireIO;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireStream;

const RegisterParser<Gro87> register_gro87;
static const RegisterMetaType<Gro87> r_gro87;

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Gro87 &gro87)
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

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Gro87 &gro87)
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

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
Gro87::Gro87(const SireSystem::System &system, const PropertyMap &map)
      : ConcreteProperty<Gro87,MoleculeParser>(map)
{
    //look through the system and extract out the information
    //that is needed to go into the file, based on the properties
    //found in 'map'. Do this to create the set of text lines that
    //will make up the file
    
    QStringList lines;  // = code used to generate the lines
    
    //now that you have the lines, reparse them back into a Gro87 object,
    //so that the information is consistent, and you have validated that
    //the lines you have written are able to be correctly read in by
    //this parser. This will also implicitly call 'assertSane()'
    Gro87 parsed( lines, map );
    
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
            }
            
            //coordinates in file in nanometers - convert to angstroms
            frame_coords_data[iatm] = Vector( 10.0 * x, 10.0 * y, 10.0 * z );
            
            if (point_idxs.count() < 6)
            {
                return;
            }
            
            //now read in the velocities
            int vlen = (3*n)+(3*(n+1));
            
            if (vals.length() < vlen)
            {
                errors.append( QObject::tr("The velocity line for atom %1 is not long "
                  "enough to contain the data. It should be %2 characters, but is really "
                  "%3 characters!").arg(iatm).arg(vlen).arg(vals.length()) );
                return;
            }
            
            x = vals.midRef(3*n, n+1).toDouble(&ok_x);
            y = vals.midRef(3*n + n + 1, n+1).toDouble(&ok_y);
            z = vals.midRef(3*n + 2*(n+1), n+1).toDouble(&ok_z);
            
            if (not (ok_x and ok_y and ok_z))
            {
                errors.append( QObject::tr("There was a problem reading the velocity "
                  "values of x, y, and z for atom %1 from the data '%2' in line '%3'")
                    .arg(iatm).arg(vals.mid(3*n,3*(n+1))).arg(line) );
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
        
        //increment the number of read frames and the start of the next line
        iframe += 1;
        iline += 2 + nats + 1;
    }

    if (not parse_warnings.isEmpty())
    {
        qDebug() << parse_warnings.join("\n");
    }

    this->setScore(0);
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in 
    the passed System that are missing coordinate data. */
void Gro87::addToSystem(System &system, const PropertyMap &map) const
{
    //you should loop through each molecule in the system and work out
    //which ones are described in the file, and then add data from the file
    //to these molecules
}
