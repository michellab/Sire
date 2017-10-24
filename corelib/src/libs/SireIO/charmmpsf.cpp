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

#include "SireIO/charmmpsf.h"

#include "SireSystem/system.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

using namespace SireIO;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireStream;

const RegisterParser<CharmmPSF> register_psf;
static const RegisterMetaType<CharmmPSF> r_psf;
static const RegisterMetaType<PSFAtom> r_psfatom(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PSFAtom &psfatom)
{
    writeHeader(ds, r_psfatom, 1);

    SharedDataStream sds(ds);

    sds << psfatom.number << psfatom.segment << psfatom.res_num << psfatom.res_name
        << psfatom.name << psfatom.type << psfatom.charge << psfatom.mass;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PSFAtom &psfatom)
{
    VersionID v = readHeader(ds, r_psfatom);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> psfatom.number >> psfatom.segment >> psfatom.res_num >> psfatom.res_name
            >> psfatom.name >> psfatom.type >> psfatom.charge >> psfatom.mass;
    }
    else
        throw version_error(v, "1", r_psfatom, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const CharmmPSF &psf)
{
    writeHeader(ds, r_psf, 1);

    ds << psf.atoms << psf.bonds << psf.angles << psf.dihedrals << psf.impropers
       << psf.cross_terms << psf.coords << static_cast<const MoleculeParser&>(psf);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, CharmmPSF &psf)
{
    VersionID v = readHeader(ds, r_psf);

    if (v == 1)
    {
        ds >> psf.atoms >> psf.bonds << psf.angles >> psf.dihedrals << psf.impropers
           >> psf.cross_terms >> psf.coords >> static_cast<MoleculeParser&>(psf);
    }
    else
        throw version_error(v, "1", r_psf, CODELOC);

    return ds;
}

/** Default constructor. */
PSFAtom::PSFAtom()
{

}

/** Constructor. */
PSFAtom::PSFAtom(const QString &line, QStringList &errors)
{

}

/** Constructor. */
PSFAtom::PSFAtom(const SireMol::Atom &atom, bool is_ter, QStringList &errors)
{

}

/** Generate a PDB record from the atom data. */
QString PSFAtom::toPSFRecord() const
{

}

const char* PSFAtom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PSFAtom>() );
}

/** Constructor */
CharmmPSF::CharmmPSF() : ConcreteProperty<CharmmPSF,MoleculeParser>()
{}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
CharmmPSF::CharmmPSF(const QString &filename, const PropertyMap &map)
     : ConcreteProperty<CharmmPSF,MoleculeParser>(filename,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    // Store the name of the input file.
    this->filename = filename;

    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
CharmmPSF::CharmmPSF(const QStringList &lines, const PropertyMap &map)
     : ConcreteProperty<CharmmPSF,MoleculeParser>(lines,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    // Store the name of the input file.
    this->filename = filename;

    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
CharmmPSF::CharmmPSF(const SireSystem::System &system, const PropertyMap &map)
     : ConcreteProperty<CharmmPSF,MoleculeParser>(map)
{
    //this->operator=(parsed);
}

/** Copy constructor */
CharmmPSF::CharmmPSF(const CharmmPSF &other) :
    ConcreteProperty<CharmmPSF,MoleculeParser>(other),
    parse_warnings(other.parse_warnings)
{}

/** Destructor */
CharmmPSF::~CharmmPSF()
{}

/** Copy assignment operator */
CharmmPSF& CharmmPSF::operator=(const CharmmPSF &other)
{
    if (this != &other)
    {
        parse_warnings = other.parse_warnings;

        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool CharmmPSF::operator==(const CharmmPSF &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool CharmmPSF::operator!=(const CharmmPSF &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* CharmmPSF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CharmmPSF>() );
}

/** Return the C++ name for this class */
const char* CharmmPSF::what() const
{
    return CharmmPSF::typeName();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr CharmmPSF::construct(const QString &filename,
                                  const PropertyMap &map) const
{
    return CharmmPSF(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr CharmmPSF::construct(const QStringList &lines,
                                  const PropertyMap &map) const
{
    return CharmmPSF(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr CharmmPSF::construct(const SireSystem::System &system,
                                  const PropertyMap &map) const
{
    return CharmmPSF(system,map);
}

/** Return a string representation of this parser */
QString CharmmPSF::toString() const
{
    if (lines().isEmpty())
        return QObject::tr("CharmmPSF::null");
    else
    {
        /*return QObject::tr("CharmmPSF( nMolecules() = %1, "
            "nResidues() = %2, nAtoms() = %3 )")
            .arg(nMolecules()).arg(nSubstructures()).arg(nAtoms());*/
    }
}

/** Return the format name that is used to identify this file format within Sire */
QString CharmmPSF::formatName() const
{
    return "CHARMMPSF";
}

/** Return a description of the file format */
QString CharmmPSF::formatDescription() const
{
    return QObject::tr("Charmm PSF format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList CharmmPSF::formatSuffix() const
{
    static const QStringList suffixes = { "psf" };
    return suffixes;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void CharmmPSF::assertSane() const
{
    //check state, raise SireError::program_bug if we are in an invalid state
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void CharmmPSF::parseLines(const PropertyMap &map)
{
    //this->setScore(nAtoms());
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map' */
System CharmmPSF::startSystem(const PropertyMap &map) const
{
    /*const int nmols = nMolecules();

    if (nmols == 0)
        return System();

    QVector<Molecule> mols(nmols);
    Molecule *mols_array = mols.data();

    if (usesParallel())
    {
        tbb::parallel_for( tbb::blocked_range<int>(0, nmols),
                           [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                mols_array[i] = this->getMolecule(i, map);
            }
        });
    }
    else
    {
        for (int i=0; i<nmols; ++i)
        {
            mols_array[i] = this->getMolecule(i, map);
        }
    }

    MoleculeGroup molgroup("all");

    for (auto mol : mols)
    {
        molgroup.add(mol);
    }

    System system;
    system.add(molgroup);
    system.setProperty(map["filename"].source(), StringProperty(filename));
    system.setProperty(map["fileformat"].source(), StringProperty(this->formatName()));*/

    //return system;
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in
    the passed System that are missing coordinate data. */
void CharmmPSF::addToSystem(System &system, const PropertyMap &map) const
{
    //you should loop through each molecule in the system and work out
    //which ones are described in the file, and then add data from the file
    //to thise molecules.
}

/** Internal function used to get the molecule structure for molecule 'imol'. */
MolStructureEditor CharmmPSF::getMolStructure(int imol, const PropertyName &cutting) const
{
    MolStructureEditor mol;

    return mol;
}

/** Internal function used to get the molecule structure for molecule 'imol'. */
MolEditor CharmmPSF::getMolecule(int imol, const PropertyMap &map) const
{
    /*return mol.setProperty(map["coordinates"], coords)
              .setProperty(map["charge"], charges)
              .setProperty(map["element"], elements)
              .commit();*/
}
