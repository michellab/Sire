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


#include "SireIO/sdf.h"

#include "SireSystem/system.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/errors.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"

#include "SireUnits/units.h"

#include <QFile>
#include <QtMath>

using namespace SireBase;
using namespace SireIO;
using namespace SireMol;
using namespace SireStream;
using namespace SireSystem;
using namespace SireUnits;

const RegisterParser<SDF> register_sdf;
static const RegisterMetaType<SDF> r_sdf;

QDataStream &operator<<(QDataStream &ds, const SDF &sdf)
{
    writeHeader(ds, r_sdf, 1);

    SharedDataStream sds(ds);

    sds << sdf.parse_warnings << static_cast<const MoleculeParser&>(sdf);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SDF &sdf)
{
    VersionID v = readHeader(ds, r_sdf);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> sdf.parse_warnings >> static_cast<MoleculeParser&>(sdf);
    }
    else
        throw version_error(v, "1", r_sdf, CODELOC);

    return ds;
}

/** Constructor */
SDF::SDF() : ConcreteProperty<SDF,MoleculeParser>()
{
}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
SDF::SDF(const QString &filename, const PropertyMap &map) :
    ConcreteProperty<SDF,MoleculeParser>(filename,map)
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
SDF::SDF(const QStringList &lines, const PropertyMap &map) :
    ConcreteProperty<SDF,MoleculeParser>(lines,map)
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
SDF::SDF(const SireSystem::System &system, const PropertyMap &map) :
    ConcreteProperty<SDF,MoleculeParser>(map)
{

    // write the file into lines, then parse them back for self-consistency
    QStringList lines;

    // Reparse the lines as a self-consistency check.
    SDF parsed(lines, map);

    this->operator=(parsed);
}

/** Copy constructor */
SDF::SDF(const SDF &other) :
    ConcreteProperty<SDF,MoleculeParser>(other),
    parse_warnings(other.parse_warnings)
{}

/** Destructor */
SDF::~SDF()
{}

/** Copy assignment operator */
SDF& SDF::operator=(const SDF &other)
{
    if (this != &other)
    {
        this->parse_warnings = other.parse_warnings;

        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool SDF::operator==(const SDF &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool SDF::operator!=(const SDF &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* SDF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SDF>() );
}

/** Return the C++ name for this class */
const char* SDF::what() const
{
    return SDF::typeName();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr SDF::construct(const QString &filename,
                                 const PropertyMap &map) const
{
    return SDF(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr SDF::construct(const QStringList &lines,
                                 const PropertyMap &map) const
{
    return SDF(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr SDF::construct(const SireSystem::System &system,
                                 const PropertyMap &map) const
{
    return SDF(system,map);
}

/** Return a string representation of this parser */
QString SDF::toString() const
{
    if (lines().isEmpty())
        return QObject::tr("SDF::null");
    else
    {
        return QObject::tr("SDF()");
    }
}

/** Return the format name that is used to identify this file format within Sire */
QString SDF::formatName() const
{
    return "SDF";
}

/** Return a description of the file format */
QString SDF::formatDescription() const
{
    return QObject::tr("Structure Data File (SDF) format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList SDF::formatSuffix() const
{
    static const QStringList suffixes = { "SDF" };
    return suffixes;
}

/** Return whether or not this is a lead parser. The lead parser is responsible
    for starting the process of turning the parsed file into the System. There
    must be one and one-only lead parser in a set of parsers creating a System */
bool SDF::isLead() const
{
    return true;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void SDF::assertSane() const
{
    QStringList errors;

    // NEED TO DO SOME CHECKS HERE!

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error(QObject::tr("There were errors reading the SDF format "
          "file:\n%1").arg(errors.join("\n\n")), CODELOC);
    }
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void SDF::parseLines(const PropertyMap &map)
{
    /* File format is decribed here:
        https://www.herongyang.com/Molecule/SDF-Format-Specification.html
        http://www.nonlinear.com/progenesis/sdf-studio/v0.9/faq/sdf-file-format-guidance.aspx
     */


    this->setScore(0);
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map'. */
System SDF::startSystem(const PropertyMap &map) const
{
    return System();
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in
    the passed System that are missing coordinate data. */
void SDF::addToSystem(System &system, const PropertyMap &map) const
{
    /* Here we add information to an existing system that has been
       generated by a lead parser.
     */
    return;

  	// Update the System fileformat property to record that it includes
    // data from this file format.
    QString fileformat = this->formatName();

    PropertyName fileformat_property = map["fileformat"];

    try
    {
        QString last_format = system.property(fileformat_property)
                                    .asA<StringProperty>().value();

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
}
