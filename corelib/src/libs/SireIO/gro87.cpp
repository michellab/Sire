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

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

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
        << gro87.vels << gro87.box_dims << gro87.box_angs
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
            >> gro87.vels >> gro87.box_dims >> gro87.box_angs
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
        coords(other.coords), vels(other.vels), box_dims(other.box_dims),
        box_angs(other.box_angs), resnums(other.resnums), resnams(other.resnams),
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
        box_dims = other.box_dims;
        box_angs = other.box_angs;
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
           box_dims == other.box_dims and box_angs == other.box_angs and
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
        return QObject::tr("Gro87( title() = %1, nAtoms() = %2, nFrames() = %5, "
                "hasCoordinates() = %3, hasVelocities() = %4 )")
                .arg(title()).arg(nAtoms())
                .arg(hasCoordinates()).arg(hasVelocities())
                .arg(nFrames());
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
        errors.append( QObject::tr("Error: The number of atoms names (%1) does not "
           "equal the number of atoms (%2)!")
                .arg(atmnams.count()).arg(nats) );
    }
    
    if (atmnums.count() != nats)
    {
        errors.append( QObject::tr("Error: The number of atoms numbers (%1) does not "
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
    
    if (box_dims.count() != box_angs.count())
    {
        errors.append( QObject::tr("Error: The number of frames of box dimension "
          "information (%1) is not equal to the number of frames of box angle "
          "information (%2).").arg(box_dims.count()).arg(box_angs.count()) );
    }
    
    if (not box_dims.isEmpty())
    {
        if (box_dims.count() != this->nFrames())
        {
            errors.append( QObject::tr("Error: The number of frames of box dimension "
               "information (%1) is not equal to the number of frames of trajectory (%2).")
                    .arg(box_dims.count()).arg(this->nFrames()) );
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

/** Return the dimensions of the box for the first frame of the trajectory */
SireMaths::Vector Gro87::boxDimensions() const
{
    if (box_dims.isEmpty())
    {
        return Vector(0);
    }
    else
    {
        return box_dims[0];
    }
}

/** Return the angles of the box for the first frame of the trajectory */
SireMaths::Vector Gro87::boxAngles() const
{
    if (box_angs.isEmpty())
    {
        return Vector(0);
    }
    else
    {
        return box_angs[0];
    }
}

/** Return the dimensions of the box for the frame 'frame' */
SireMaths::Vector Gro87::boxDimensions(int frame) const
{
    return box_dims[ Index(frame).map(box_dims.count()) ];
}

/** Return the angles of the box for the frame 'frame' */
SireMaths::Vector Gro87::boxAngles(int frame) const
{
    return box_angs[ Index(frame).map(box_angs.count()) ];
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
