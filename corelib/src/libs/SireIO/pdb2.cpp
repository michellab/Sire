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

#include "SireIO/pdb2.h"

#include "SireMM/pdbparams.h"

#include "SireSystem/system.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireIO;
using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireStream;

const RegisterParser<PDB2> register_pdb;
static const RegisterMetaType<PDB2> r_PDB2;

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDB2 &PDB2)
{
    writeHeader(ds, r_PDB2, 1);
    
    ds << static_cast<const MoleculeParser&>(PDB2);
    
    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDB2 &PDB2)
{
    VersionID v = readHeader(ds, r_PDB2);
    
    if (v == 1)
    {
        ds >> static_cast<MoleculeParser&>(PDB2);
    }
    else
        throw version_error(v, "1", r_PDB2, CODELOC);

    return ds;
}

/** Constructor */
PDB2::PDB2() : ConcreteProperty<PDB2,MoleculeParser>()
{}

/** Construct to read in the data from the file called 'filename'. The 
    passed property map can be used to pass extra parameters to control
    the parsing */
PDB2::PDB2(const QString &filename, const PropertyMap &map)
     : ConcreteProperty<PDB2,MoleculeParser>(filename,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.
    
    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    //parse the data in the parse function
    this->parseLines(map);
    
    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
PDB2::PDB2(const QStringList &lines, const PropertyMap &map)
     : ConcreteProperty<PDB2,MoleculeParser>(lines,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.
    
    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    //parse the data in the parse function
    this->parseLines(map);
    
    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
PDB2::PDB2(const SireSystem::System &system, const PropertyMap &map)
     : ConcreteProperty<PDB2,MoleculeParser>(map)
{
    //look through the system and extract out the information
    //that is needed to go into the file, based on the properties
    //found in 'map'. Do this to create the set of text lines that
    //will make up the file
    
    QStringList lines;  // = code used to generate the lines
    
    //now that you have the lines, reparse them back into a PDB2 object,
    //so that the information is consistent, and you have validated that
    //the lines you have written are able to be correctly read in by
    //this parser. This will also implicitly call 'assertSane()'
    PDB2 parsed( lines, map );
    
    this->operator=(parsed);
}

/** Copy constructor */
PDB2::PDB2(const PDB2 &other)
     : ConcreteProperty<PDB2,MoleculeParser>(other)
{}

/** Destructor */
PDB2::~PDB2()
{}

/** Copy assignment operator */
PDB2& PDB2::operator=(const PDB2 &other)
{
    if (this != &other)
    {
        MoleculeParser::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool PDB2::operator==(const PDB2 &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool PDB2::operator!=(const PDB2 &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* PDB2::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDB2>() );
}

/** Return the C++ name for this class */
const char* PDB2::what() const
{
    return PDB2::typeName();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr PDB2::construct(const QString &filename,
                                  const PropertyMap &map) const
{
    return PDB2(filename,map);
}

/** Return the parser that has been constructed by reading in the passed  
    text lines using the passed properties */
MoleculeParserPtr PDB2::construct(const QStringList &lines,
                                  const PropertyMap &map) const
{
    return PDB2(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr PDB2::construct(const SireSystem::System &system,
                                  const PropertyMap &map) const
{
    return PDB2(system,map);
}

/** Return a string representation of this parser */
QString PDB2::toString() const
{
    return QObject::tr("PDB2::null");
}

/** Return the format name that is used to identify this file format within Sire */
QString PDB2::formatName() const
{
    return "PDB";
}

/** Return a description of the file format */
QString PDB2::formatDescription() const
{
    return QObject::tr("Protein Data Bank (PDB) format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList PDB2::formatSuffix() const
{
    static const QStringList suffixes = { "PDB" };
    return suffixes;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void PDB2::assertSane() const
{
    //check state, raise SireError::program_bug if we are in an invalid state
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void PDB2::parseLines(const PropertyMap &map)
{
    //you should write the code to parse the file here. The aim is to
    //parse the information into an intermediate format, e.g. arrays
    //of coordinates, arrays of atom names, etc. etc.
    //If there is any PDB-file specific information then create
    //specific PDBParams classes using SireMM::PDBParams
    
    //this information will then be interpreted by the "startSystem" and
    //"addToSystem" functions to extract actual Molecule and System objects
    //via a later function call.
    
    //You should set a score from parsing. The higher the score, the better
    //the parsing. Set a score of 0 if the text does not contain a viable
    //file. If it is really bad, you can also throw a SireIO::parse_error
    //exception

    this->setScore(0);
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map' */
System PDB2::startSystem(const PropertyMap &map) const
{
    //you should now take the data that you have parsed into the intermediate
    //format and use it to start a new System and create all molecules,
    //which are added to this System
    
    System system;
    
    
    return system;
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in 
    the passed System that are missing coordinate data. */
void PDB2::addToSystem(System &system, const PropertyMap &map) const
{
    //you should loop through each molecule in the system and work out
    //which ones are described in the file, and then add data from the file
    //to thise molecules.
}
