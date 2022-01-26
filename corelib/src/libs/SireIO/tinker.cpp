/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include <QString>
#include <QFile>
#include <QRegExp>
#include <QStringList>
#include <QTextStream>

#include "tinker.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/molecules.h"
#include "SireMol/element.h"

#include "SireMaths/vector.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireIO;
using namespace SireMol;
using namespace SireMaths;
using namespace SireStream;

/////////////
///////////// Implementation of TinkerParameters
/////////////

/** Constructor */
TinkerParameters::TinkerParameters() : IOParametersBase()
{}

/** Destructor */
TinkerParameters::~TinkerParameters()
{}

/////////////
///////////// Implementation of everything needed to get the
///////////// Tinker reader/writer working
/////////////

/** This is an internal class used to store all of the information
    available on a tinker XYZ atom line

    e.g.
     2  CT1  -22.301380   12.051195   -1.495658    24     1     3     8     9

    First entry is atom number, then atom name, then x,y,z coordinates,
    then the numbers of the atoms this atom is bonded to.

    Tinker XYZ does not separate the file into molecules - this is
    implied by the connectivity information
*/
class TinkerXYZAtom
{
public:
    TinkerXYZAtom();

    ~TinkerXYZAtom();

    static TinkerXYZAtom readFromLine(const QString &line, int linenum);

    QString writeToLine() const;

    Element guessElement() const;

    Vector coordinates() const;

    int linenum;

    int atomnum;
    QString atomname;

    double x, y, z;

    int param;

    QList<int> bonded_to;
};

TinkerXYZAtom::TinkerXYZAtom()
{}

TinkerXYZAtom::~TinkerXYZAtom()
{}

TinkerXYZAtom TinkerXYZAtom::readFromLine(const QString &line, int linenum)
{
    QStringList words = line.split(" ", Qt::SkipEmptyParts);

    if (words.count() < 6)
        throw SireIO::parse_error( QObject::tr(
                "Line %1 does not look like a valid Tinker XYZ line. A Tinker "
                "XYZ line should have at least 6 space separated words, but the "
                "number of space separated words on this line is just %2.\n%3")
                    .arg(linenum).arg(words.count()).arg(line), CODELOC );

    TinkerXYZAtom atom;

    atom.linenum = linenum;

    bool ok = true;

    atom.atomnum = words[0].toInt(&ok);

    if (not ok)
        throw SireIO::parse_error( QObject::tr(
                "Line %1 does not look like a valid Tinker XYZ line. The first "
                "space separated word in the line should be an integer, but "
                "the first word (%2) does not look like an integer.\n%3")
                    .arg(linenum).arg(words[0]).arg(line), CODELOC );

    atom.atomname = words[1];

    atom.x = words[2].toDouble(&ok);

    if (not ok)
        throw SireIO::parse_error( QObject::tr(
                "Line %1 does not look like a valid Tinker XYZ line. The third "
                "space separated word in the line should be the X coordinate, but "
                "the third word (%2) does not look like a number.\n%3")
                    .arg(linenum).arg(words[2]).arg(line), CODELOC );

    atom.y = words[3].toDouble(&ok);

    if (not ok)
        throw SireIO::parse_error( QObject::tr(
                "Line %1 does not look like a valid Tinker XYZ line. The fourth "
                "space separated word in the line should be the Y coordinate, but "
                "the fourth word (%2) does not look like a number.\n%3")
                    .arg(linenum).arg(words[3]).arg(line), CODELOC );

    atom.z = words[4].toDouble(&ok);

    if (not ok)
        throw SireIO::parse_error( QObject::tr(
                "Line %1 does not look like a valid Tinker XYZ line. The fifth "
                "space separated word in the line should be the Z coordinate, but "
                "the fifth word (%2) does not look like a number.\n%3")
                    .arg(linenum).arg(words[4]).arg(line), CODELOC );

    //now read in the parameter number
    atom.param = words[5].toInt(&ok);

    if (not ok)
        throw SireIO::parse_error( QObject::tr(
                "Line %1 does not look like a valid Tinker XYZ line. The sixth "
                "space separated word in the line should be the parameter number, "
                "but the size word (%2) does look like an integer.\n%3")
                    .arg(linenum).arg(words[5]).arg(line), CODELOC );

    //now read in the connectivity
    for (int i=6; i<words.count(); ++i)
    {
        atom.bonded_to.append( words[i].toInt(&ok) );

        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                    "Line %1 does not look like a valid Tinker XYZ line. The last "
                    "space separated words in the line should be integers giving the "
                    "connectivity of the molecule. Word %2 (%3) is not an integer.\n%4")
                        .arg(linenum).arg(i).arg(words[i]).arg(line), CODELOC );
    }

    return atom;
}

QString TinkerXYZAtom::writeToLine() const
{
    QString line;

    QTextStream ts(&line);

    ts.setRealNumberNotation( QTextStream::FixedNotation );
    ts.setRealNumberPrecision( 6 );

    ts.setFieldWidth(6);
    ts << atomnum;

    ts.setFieldWidth(1);
    ts << " ";

    ts.setFieldWidth(4);
    ts << atomname;

    ts.setFieldWidth(1);
    ts << " ";

    ts.setFieldWidth(11);
    ts << x;

    ts.setFieldWidth(1);
    ts << " ";

    ts.setFieldWidth(11);
    ts << y;

    ts.setFieldWidth(1);
    ts << " ";

    ts.setFieldWidth(11);
    ts << z;

    ts.setFieldWidth(1);

    ts << " ";
    ts.setFieldWidth(5);
    ts << param;

    ts.setFieldWidth(1);

    foreach (int bond, bonded_to)
    {
        ts << " ";

        ts.setFieldWidth(5);
        ts << bond;

        ts.setFieldWidth(1);
    }

    return line;
}

/** This is a Tinker XYZ molecule - this is used just to hold all of the
    atoms, and to help resolve the connectivity (so only a set of
    connected atoms will be added to the molecule) */
class TinkerXYZMolecule
{
public:
    TinkerXYZMolecule();
    ~TinkerXYZMolecule();

    bool contains(const TinkerXYZAtom &atom) const;

    void add(const TinkerXYZAtom &atom);

    QList<TinkerXYZAtom> atoms;

    QSet<int> contained_atoms;
};

TinkerXYZMolecule::TinkerXYZMolecule()
{}

TinkerXYZMolecule::~TinkerXYZMolecule()
{}

void TinkerXYZMolecule::add(const TinkerXYZAtom &atom)
{
    atoms.append(atom);

    contained_atoms.insert(atom.atomnum);

    foreach (int bond, atom.bonded_to)
    {
        contained_atoms.insert(bond);
    }
}

bool TinkerXYZMolecule::contains(const TinkerXYZAtom &atom) const
{
    if (contained_atoms.contains(atom.atomnum))
        return true;
    else
    {
        foreach (int bond, atom.bonded_to)
        {
            if (contained_atoms.contains(bond))
                return true;
        }

        return false;
    }
}

/** This class holds a tinker atom parameter */
/*class TinkerAtomParam
{
public:
    TinkerAtomParam();

    ~TinkerAtomParam();

    static bool matches(const QString &line);

    static TinkerAtomParam readFromLine(const QString &line, int linenum);

    int id_num;
    int lj_type;
    QString atm_type;
    QString description;
    Element element;
    double mass;
    int bond_order;
};

TinkerAtomParam::TinkerAtomParam() : id_num(0), lj_type(0), mass(0), bond_order(0)
{}

TinkerAtomParam::~TinkerAtomParam()
{}

bool TinkerAtomParam::matches(const QString &line)
{
    return line.startsWith("atom ", Qt::CaseInsensitive);
}

TinkerAtomParam TinkerAtomParam::readFromLine(const QString &line, int linenum)
{
    // The passed regexp must have the form
    //atom      1     1    HA      "Nonpolar Hydrogen"         1      1.008     1
    QRegExp atom_regexp("atom\\s+(\\d+)\\s+(\\d+)\\s+(\\w+)\\s+"
                        "\"([\\w\\s\\d]+)\"\\s+(\\d+)\\s+([\\d\\.]+)\\s+(\\d+)");

    if (atom_regexp.indexIn(line) == -1)
        throw SireIO::parse_error( QObject::tr(
                "Line number %1 does not look like a valid Tinker atom parameter "
                "line, despite it starting with \"atom\"\n%2")
                    .arg(linenum).arg(line), CODELOC );

    TinkerAtomParam param;

    param.id_num = atom_regexp.cap(1).toInt();
    param.lj_typ = atom_regexp.cap(2).toInt();
    param.atm_typ = atom_regexp.cap(3);
    param.description = atom_regexp.cap(4);
    param.element = Element( atom_regexp.cap(5).toInt() );
    param.mass = atom_regexp.cap(6).toDouble();
    param.bond_order = atom_regexp.cap(7).toInt();

    return param;
}
*/
/////////////
///////////// Implementation of Tinker
/////////////

static const RegisterMetaType<Tinker> r_tinker;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Tinker &tinker)
{
    writeHeader(ds, r_tinker, 1);

    SharedDataStream sds(ds);

    sds //<<
        << static_cast<const IOBase&>(tinker);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Tinker &tinker)
{
    VersionID v = readHeader(ds, r_tinker);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds //>>
            >> static_cast<IOBase&>(tinker);
    }
    else
        throw version_error(v, "1", r_tinker, CODELOC);

    return ds;
}

/** Constructor */
Tinker::Tinker() : ConcreteProperty<Tinker,IOBase>()
{}

/** Copy constructor */
Tinker::Tinker(const Tinker &other) : ConcreteProperty<Tinker,IOBase>(other)
{}

/** Destructor */
Tinker::~Tinker()
{}

const char* Tinker::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Tinker>() );
}

/** Copy assignment operator */
Tinker& Tinker::operator=(const Tinker &other)
{
    IOBase::operator=(other);
    return *this;
}

/** Comparison operator */
bool Tinker::operator==(const Tinker &other) const
{
    return IOBase::operator==(other);
}

/** Comparison operator */
bool Tinker::operator!=(const Tinker &other) const
{
    return not Tinker::operator==(other);
}

Q_GLOBAL_STATIC( TinkerParameters, tinkerParameters )

/** Return the default parameters used by this IO object */
const TinkerParameters& Tinker::parameters()
{
    return *(tinkerParameters());
}

/** Load the parameter file 'prmfile' into this reader. This will allow
    the reader to parameterise the molecules as they are being read */
void Tinker::loadParameters(const QString &prmfile)
{
    QFile f(prmfile);

    if (not f.open(QIODevice::ReadOnly))
        throw SireError::file_error(f, CODELOC);

    //QRegExp objects to match the parts of the file we are interested in

    //atom      1     1    HA      "Nonpolar Hydrogen"         1      1.008     1
    QRegExp atom_regexp("atom\\s+(\\d+)\\s+(\\d+)\\s+(\\w+)\\s+"
                        "\"([\\w\\s\\d]+)\"\\s+(\\d+)\\s+([\\d\\.]+)\\s+(\\d+)");

    QTextStream ts(&f);

    int linenum = -1;

    while (not ts.atEnd())
    {
        ++linenum;

        QString line = ts.readLine().simplified();

        if (line.isEmpty() or line.startsWith("#"))
            continue;

        //if (TinkerAtomParam::matches(line))
        //{
        //    TinkerAtomParam = TinkerAtomParam::readFromLine(line, linenum);
       // }
    }
}

/** Read the molecule group from 'data', using the parameters in 'map' */
MoleculeGroup Tinker::readMols(const QByteArray &data, const PropertyMap &map) const
{
    //connect a text stream to this bytearray so that we can read it as text
    QTextStream ts(data, QIODevice::ReadOnly | QIODevice::Text);

    int linenum = -1;

    QList<TinkerXYZMolecule> mols;

    //read all of the molecules from the file - place them into
    //the intermediate format 'mols'
    while (not ts.atEnd())
    {
        ++linenum;

        //read a line from the file
        QString line = ts.readLine().simplified();

        if (line.isEmpty())
            continue;

        //if this is not the header line, then read the atoms
        if (linenum > 0)
        {
            TinkerXYZAtom atom = TinkerXYZAtom::readFromLine(line, linenum);

            bool added_atom = false;

            for (QList<TinkerXYZMolecule>::iterator it = mols.begin();
                 it != mols.end();
                 ++it)
            {
                if (it->contains(atom))
                {
                    it->add(atom);
                    added_atom = true;
                    break;
                }
            }

            if (not added_atom)
            {
                TinkerXYZMolecule mol;
                mol.add(atom);
                mols.append(mol);
            }
        }
    }

    //now run through each molecule and create a SireMol::Molecule object
    //to hold it
    /*MoleculeGroup molgroup;

    foreach (TinkerXYZMolecule mol, mols)
    {
        Molecule molecule;
    }*/

    return MoleculeGroup();
}

/** Write the molecule group in 'molgroup' to an array, using the parameters
    in map */
QByteArray Tinker::writeMols(const MoleculeGroup &molgroup, const PropertyMap &map) const
{
    return QByteArray();
}

/** Write the molecules in 'molecules' to an  array, using the parameters
    in map */
QByteArray Tinker::writeMols(const Molecules &molecules, const PropertyMap &map) const
{
    return QByteArray();
}
