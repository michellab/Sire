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

#include "SireIO/supplementary.h"

#include "SireError/errors.h"

#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireIO;
using namespace SireStream;
using namespace SireSystem;

const RegisterParser<Supplementary> register_supp;
static const RegisterMetaType<Supplementary> r_supp;

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Supplementary &supp)
{
    writeHeader(ds, r_supp, 1);

    ds << static_cast<const MoleculeParser&>(supp);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Supplementary &supp)
{
    VersionID v = readHeader(ds, r_supp);

    if (v == 1)
    {
        ds >> static_cast<MoleculeParser&>(supp);
    }
    else
        throw version_error(v, "1", r_supp, CODELOC);

    return ds;
}

/** Constructor */
Supplementary::Supplementary() : ConcreteProperty<Supplementary,MoleculeParser>()
{}

/** Construct to read in the data from the file called 'filename'. The
passed property map can be used to pass extra parameters to control
the parsing */
Supplementary::Supplementary(const QString &filename, const PropertyMap &map)
     : ConcreteProperty<Supplementary,MoleculeParser>(filename,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    // Store the name of the input file.
    this->filename = filename;

    // Set score to 2. (Need a low, non-zero score)
    // This ensures that not all files will be flagged as supplementary,
    // i.e. they will be scored more highly by their native parser.
    // We can't set the score to 1 since certain supplementary NAMD data
    // records, e.g. box information, are given a score of 1 by the AmberRst7
    // parser due apparent to similarities in the record formatting.
    // TODO: Work out why AmberRst7 thinks xst/xsc box records are atoms.
    this->setScore(2);
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
Supplementary::Supplementary(const QStringList &lines, const PropertyMap &map)
     : ConcreteProperty<Supplementary,MoleculeParser>(lines,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    // Store the name of the input file.
    this->filename = filename;

    // Set score to one. (Need a low, non-zero score)
    // This ensures that not all files will be flagged as supplementary,
    // i.e. they will be scored more highly by their native parser.
    this->setScore(1);
}

/** Copy constructor */
Supplementary::Supplementary(const Supplementary &other) :
    ConcreteProperty<Supplementary,MoleculeParser>(other),
    filename(other.filename)
{
}

/** Destructor */
Supplementary::~Supplementary()
{}

/** Copy assignment operator */
Supplementary& Supplementary::operator=(const Supplementary &other)
{
    if (this != &other)
    {
        filename = other.filename;

        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool Supplementary::operator==(const Supplementary &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool Supplementary::operator!=(const Supplementary &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* Supplementary::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Supplementary>() );
}

/** Return the C++ name for this class */
const char* Supplementary::what() const
{
    return Supplementary::typeName();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr Supplementary::construct(const QString &filename,
                                  const PropertyMap &map) const
{
    return Supplementary(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr Supplementary::construct(const QStringList &lines,
                                    const PropertyMap &map) const
{
    return Supplementary(lines,map);
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr Supplementary::construct(const SireSystem::System &system,
                                    const PropertyMap &map) const
{
    throw SireError::program_bug( QObject::tr(
            "The Supplementary parser cannot construct from a SireSystem::System!"), CODELOC );

    return MoleculeParserPtr();
}

/** Return a string representation of this parser */
QString Supplementary::toString() const
{
    if (lines().isEmpty())
        return QObject::tr("Supplementary::null");
    else
    {
        return QObject::tr("Supplementary( filename = %1 )").arg(filename);
    }
}

/** Return the format name that is used to identify this file format within Sire */
QString Supplementary::formatName() const
{
    return "SUPPLEMENTARY";
}

/** Return a description of the file format */
QString Supplementary::formatDescription() const
{
    return QObject::tr("Files that are supplementary to a lead parser.");
}

/** Return the suffixes that these files are normally associated with */
QStringList Supplementary::formatSuffix() const
{
    static const QStringList suffixes = { "*" };
    return suffixes;
}

/** Return whether or not this parser can follow another lead parser, and add
    data to an existing molecular system. The Supplementary parser cannot follow. */
bool Supplementary::canFollow() const
{
    return false;
}
