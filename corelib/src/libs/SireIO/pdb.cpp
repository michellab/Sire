/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include <QTextStream>
#include <QByteArray>
#include <QString>
#include <QStringList>
#include <QMap>

#include "pdb.h"

#include "SireMol/element.h"

#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/atomcharges.h"
#include "SireMol/connectivity.h"

#include "SireMol/mover.hpp"
#include "SireMol/selector.hpp"

#include "SireMol/molecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/molidx.h"
#include "SireMol/molidentifier.h"
#include "SireMol/moleditor.h"
#include "SireMol/segeditor.h"
#include "SireMol/chaineditor.h"
#include "SireMol/reseditor.h"
#include "SireMol/cgeditor.h"
#include "SireMol/atomeditor.h"
#include "SireMol/residue.h"
#include "SireMol/chain.h"
#include "SireMol/segment.h"

#include "SireMol/cuttingfunction.h"

#include "SireBase/stringmangler.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireStream/sharestrings.h"

#include "tostring.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireIO;
using namespace SireUnits;
using namespace SireStream;

///////////
/////////// Implementation of PDBParameters
///////////

PDBParameters::PDBParameters() : IOParametersBase()
{}

PDBParameters::~PDBParameters()
{}

PropertyName PDBParameters::animation_property("animation_frames");
PropertyName PDBParameters::alternatives_property("alternative");
PropertyName PDBParameters::icode_property("icode");
PropertyName PDBParameters::bfactor_property("b-factor");
PropertyName PDBParameters::formalcharge_property("formal-charge");

PropertyName PDBParameters::pdbatomnames_property("PDB-atom-name");
PropertyName PDBParameters::pdbresnames_property("PDB-residue-name");
PropertyName PDBParameters::pdbchainnames_property("PDB-chain-name");
PropertyName PDBParameters::pdbsegnames_property("PDB-segment-name");

PropertyName PDBParameters::frame_selector( "animation-frame-selector",
                                            Property::null() );
PropertyName PDBParameters::atomname_mangler( "atom-name-mangler",
                                              TrimString() );
PropertyName PDBParameters::resname_mangler( "residue-name-mangler",
                                             TrimString() );
PropertyName PDBParameters::chainname_mangler( "chain-name-mangler", 
                                               TrimString() );
PropertyName PDBParameters::segname_mangler( "segment-name-mangler",
                                             TrimString() );

///////////
/////////// Implementation of everything to get the PDB reader/writer working
///////////

/** This internal class is used to store all of the data
    that is used in a PDB ATOM line 

    // PBB ATOM line
    //
    // Format as described in PDB 2.3 format guide
    // (see link above)
    //
    // COLUMNS             DATA TYPE        FIELD      DEFINITION 
    // ------------------------------------------------------------- 
    //  1 -   6            Record name      "ATOM    " 
    //
    //  7 - 11             Integer          serial     Atom serial number. 
    //
    // 13 - 16             Atom             name       Atom name. 
    //
    // 17                  Character        altLoc     Alternate location 
    //                                                 indicator. 
    //
    // 18 - 20             Residue name     resName    Residue name. 
    //
    // 22                  Character        chainID    Chain identifier. 
    //
    // 23 - 26             Integer          resSeq     Residue sequence number. 
    //
    // 27                  AChar            iCode      Code for insertion of 
    //                                                 residues. 
    // 
    // 31 - 38             Real(8.3)        x          Orthogonal coordinates for 
    //                                                 X in Angstroms 
    //
    // 39 - 46             Real(8.3)        y          Orthogonal coordinates for 
    //                                                 Y in Angstroms 
    //
    // 47 - 54             Real(8.3)        z          Orthogonal coordinates for 
    //                                                 Z in Angstroms 
    //
    // 55 - 60             Real(6.2)        occupancy  Occupancy. 
    // 
    // 61 - 66             Real(6.2)        tempFactor Temperature factor. 
    //
    // 73 - 76             LString(4)       segid      Segment ID, left-justified,
    //                                                 (not part of 2.1 spec, but used
    //                                                  widely)
    //
    // 77 - 78             LString(2)       element    Element symbol, 
    //                                                 right-justified. 
    //
    // 79 - 80             LString(2)       charge     Charge on the atom. 

    // 000000000011111111112222222222333333333344444444445555555555666666666677777777778
    // 012345678901234567890123456789012345678901234567890123456789012345678901234567890 
    // ATOM      1  N   MET     2     -15.160  18.227  60.039  1.00 13.87      MA  N
    // ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92          N 
    // ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85          C 
    // ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34          C  
    // ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65          O 
    // ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88          C 
    // ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41          C 
    // ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64          C 
    // ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11          C 
    // ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58          C 
    // ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25          C 
    // ATOM    155 CA   CAL A  26      -5.000   3.234 -10.034  1.00  1.00          CA2+

    @author Christopher Woods
*/
class PDBAtom
{
public:
    PDBAtom();

    ~PDBAtom();

    static PDBAtom readFromLine(const QString &line, int linenum);
    
    QString toString() const;
    
    QString writeToLine() const;
    
    bool isATOM() const;
    bool isHETATM() const;
    
    Element guessElement() const;
    Vector coordinates() const;
    
    QString   record_name;   // ATOM or HETATM      - line[0:5]
    int       serial;        // ATOM ID number      - line[6:11]
    QString   name;          // atom name           - line[12:15]
    QString   altloc;        // atom altloc         - line[16]
    QString   resname;       // atom residue name   - line[17:19] 
    QString   chainid;       // atom chain name     - line[21]
    int       resseq;        // atom residue number - line[22:25]
    QString   icode;         // atom icode          - line[26]
    double    x;             // atom x coords       - line[30:37]
    double    y;             // atom y coords       - line[38:45]
    double    z;             // atom z coords       - line[46:53]
    double    occupancy;     // atom occupancy      - line[54:59]
    double    tempfactor;    // bfactor             - line[60:65]
    QString   segid;         // atom segment name   - line[72:75]
    QString   element;       // atom element        - line[76:77]
    int       charge;        // atom formal charge  - line[78:79]
  
    /** This is an ID string that is used 
        to help index an atom */
    QString idstring;
};

/** The collection of all information read about 
    a molecule from a PDB file */
class PDBMolecule
{
public:
    PDBMolecule();
    ~PDBMolecule();
    
    void addAtom(const PDBAtom &atom);
    
    QList<int> availableFrames() const;
    
    const QList<PDBAtom>& atoms(int frame=1) const;
    
    bool hasMultipleOccupancy(int frame=1) const;
    bool hasMultipleFrames() const;
    
    void openFrame(int framenum);
    void closeFrame();
    
    bool isEmpty() const;

private:    
    /** Extra sets of coordinates (in case of multiple models,
        e.g. from NMR or from a trajectory) */
    QMap< int, QList<PDBAtom> > frames;
    
    int current_frame;
};

/** The collection of all molecules read from a PDB file */
class PDBMolecules
{
public:
    PDBMolecules();
    ~PDBMolecules();
    
    PDBMolecule& nextMolecule();
    void closeMolecule();
    
    int openFrame(int frame_number);
    void closeFrame();
    
    bool frameOpened() const;
    bool moleculeOpened() const;

    bool isEmpty() const;

    const QList<PDBMolecule>& molecules() const;

private:
    QList<PDBMolecule> mols;
    
    int last_frame;
    int current_molecule;
    
    bool frame_opened;
    bool molecule_opened;
};

////////
//////// Implementation of PDBAtom
////////

PDBAtom::PDBAtom() : record_name("ATOM"), serial(0), name("UNK"),
                     altloc(" "), resname("UNK"), chainid(" "),
                     resseq(-1), icode(" "), x(0), y(0), z(0),
                     occupancy(0), tempfactor(0), segid("  "),
                     element("  "), charge(0)
{}

PDBAtom::~PDBAtom()
{}

QString PDBAtom::toString() const
{
    return QString("%1 %2 %3 %4 %5 %6 %7 %8 %9 %10 %11 %12 %13 %14 %15 %16")
                .arg(record_name).arg(serial).arg(name)
                .arg(altloc).arg(resname).arg(chainid).arg(resseq)
                .arg(icode).arg(x).arg(y).arg(z).arg(occupancy)
                .arg(tempfactor).arg(segid).arg(element).arg(charge);
}

PDBAtom PDBAtom::readFromLine(const QString &line, int linenum)
{
    if (line.length() < 54)
        //the line is too short to be a PDB atom line!
        throw SireIO::parse_error( QObject::tr(
            "Line %1 does not look like a valid PDB line. It only contains "
            "%2 characters, when it should contain at least 54.\n%3.")
                .arg(linenum).arg(line.length()).arg(line), CODELOC );

    PDBAtom atom;

    bool ok = true;

    atom.record_name = SireStream::shareString( line.mid(0,6).trimmed() );
    
    atom.serial = line.mid(6,5).toInt(&ok);
    
    if (not ok)
        throw SireIO::parse_error( QObject::tr(
            "Line %1 does not have a valid PDB atom number (%2).\n%3")
                .arg(linenum).arg(line.mid(6,5).trimmed(), line), CODELOC );

    atom.name = SireStream::shareString( line.mid(12,4) );
    atom.altloc = line.mid(16,1).trimmed();
    atom.resname = SireStream::shareString( line.mid(17,3) );
    atom.chainid = line.mid(21,1).trimmed();
    atom.resseq = line.mid(22,4).toInt(&ok);
    
    if (not ok)
        throw SireIO::parse_error( QObject::tr(
            "Line %1 does not have a valid PDB residue number (%2).\n%3")
                .arg(linenum).arg(line.mid(22,4).trimmed(), line), CODELOC );
                
    atom.icode = line.mid(26,1).trimmed();
    
    atom.x = line.mid(30,8).toDouble(&ok);
    
    if (not ok)
        throw SireIO::parse_error( QObject::tr(
            "Line %1 does not have a valid PDB x coordinate (%2)\n%3")
                .arg(linenum).arg(line.mid(30,8).trimmed(), line), CODELOC );
    
    atom.y = line.mid(38,8).toDouble(&ok);
    
    if (not ok)
        throw SireIO::parse_error( QObject::tr(
            "Line %1 does not have a valid PDB y coordinate (%2)\n%3")
                .arg(linenum).arg(line.mid(38,8).trimmed(), line), CODELOC );
    
    atom.z = line.mid(46,8).toDouble(&ok);
    
    if (not ok)
        throw SireIO::parse_error( QObject::tr(
            "Line %1 does not have a valid PDB z coordinate (%2)\n%3")
                .arg(linenum).arg(line.mid(46,8).trimmed(), line), CODELOC );

    if (line.length() >= 55)
    {
        //we can read the occupancy
        atom.occupancy = line.mid(54,6).toDouble(&ok);
        
        if (not ok)
            //something went wrong reading the occupancy - ignore the error
            atom.occupancy = 1;
    }
    else
        atom.occupancy = 1;
    
    if (line.length() >= 61)
    {
        //try to read the temperature factor
        atom.tempfactor = line.mid(60,6).toDouble(&ok);
        
        if (not ok)
            //something went wrong
            atom.tempfactor = 0;
    }
    else
        atom.tempfactor = 0;
    
    if (line.length() >= 73)
    {
        atom.segid = SireStream::shareString( line.mid(72,4).trimmed() );
    }
    else
        atom.segid = QString::null;
    
    if (line.length() >= 77)
    {
        atom.element = SireStream::shareString( line.mid(76,2).trimmed() );
    }
    else
        atom.element = QString::null;
    
    if (line.length() >= 79)
    {
        QString chgstring = line.mid(78,2).trimmed();
        
        //format should be 2+ or 1- or something like that
        int factor = 1;
        
        if (chgstring.contains("-"))
            factor = -1;
            
        chgstring.remove("-").remove("+");
        
        atom.charge = factor * chgstring.toInt(&ok);
        
        if (not ok)
            //something went wrong - ignore the charge
            atom.charge = 0;
    }
    else
        atom.charge = 0;
        
    char idline[31];
    idline[30] = '\0';
    
    qsnprintf(idline, 30, "%4s %3s %1s %4d%1s %4s",
              qPrintable(atom.name), qPrintable(atom.resname),
              qPrintable(atom.chainid), atom.resseq,
              qPrintable(atom.icode), qPrintable(atom.segid));
        
    atom.idstring = QString::fromLocal8Bit(idline);
        
    return atom;
}

QString PDBAtom::writeToLine() const
{
    QString chg;
    
    if (charge > 0)
        chg = QString("%1+").arg(charge);
    else if (charge < 0)
        chg = QString("%1-").arg(charge);

    char line[83];
    line[82] = '\0';
    
    int num = serial;
    
    if (num < 0)
        num = 0;
        
    while (num >= 100000)
    {
        num -= 100000;
    }

    int resnum = resseq;
    
    if (resnum < 0)
        resnum = 0;

    while (resnum >= 10000)
    {
        resnum -= 10000;
    }
    
    qsnprintf(line, 82,
     "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s",
            qPrintable(record_name), num, qPrintable(name),
            qPrintable(altloc), qPrintable(resname), qPrintable(chainid),
            resnum, qPrintable(icode), x, y, z, occupancy, tempfactor,
            qPrintable(segid), qPrintable(element), qPrintable(chg));

    return QString::fromLocal8Bit(line).trimmed();
}

bool PDBAtom::isATOM() const
{
    return QString::compare(record_name, "ATOM", Qt::CaseInsensitive) == 0;
}

bool PDBAtom::isHETATM() const
{
    return QString::compare(record_name, "HETATM", Qt::CaseInsensitive) == 0;
}

Element PDBAtom::guessElement() const
{
    if (element.isEmpty())
    {
        //we need to guess the element from the atom name - the
        //element is the first two characters of the name, right
        //justified
        return Element( element.mid(0,2).trimmed() );
    }
    else
        return Element(element);
}

Vector PDBAtom::coordinates() const
{
    return Vector(x,y,z);
}

////////
//////// Implementation of PDBMolecule
////////

void PDBMolecule::openFrame(int framenum)
{
    frames.insert(framenum, QList<PDBAtom>());
    current_frame = framenum;
}

void PDBMolecule::closeFrame()
{
    if (current_frame != -1)
    {
        if (frames.constFind(current_frame)->isEmpty())
        {
            frames.remove(current_frame);
            current_frame = -1;
        }
    }
}

PDBMolecule::PDBMolecule() : current_frame(-1)
{}

PDBMolecule::~PDBMolecule()
{}

void PDBMolecule::addAtom(const PDBAtom &atom)
{
    if (current_frame == -1)
        //open the default frame (1)
        this->openFrame(1);

    frames[current_frame].append(atom);
}

QList<int> PDBMolecule::availableFrames() const
{
    return frames.keys();
}

bool PDBMolecule::isEmpty() const
{
    return frames.isEmpty();
}

const QList<PDBAtom>& PDBMolecule::atoms(int frame) const
{
    QMap< int,QList<PDBAtom> >::const_iterator it = frames.constFind(frame);
    
    if (it == frames.constEnd())
        throw SireError::invalid_index( QObject::tr(    
            "There is no frame at index %1. Available frames are %2.")
                .arg(frame).arg( Sire::toString(frames.keys()) ), CODELOC );

    return *it;
}

bool PDBMolecule::hasMultipleOccupancy(int frame) const
{
    foreach (const PDBAtom &atom, frames.value(frame))
    {
        if (atom.occupancy < 1)
            //as the sum of occupancies for a site must equal 1,
            //this implies that we do have multiple occupancy in this frame
            return true;
    }
    
    return false;
}

bool PDBMolecule::hasMultipleFrames() const
{
    return frames.count() > 1;
}

////////
//////// Implementation of PDBMolecules
////////

PDBMolecules::PDBMolecules() : last_frame(0), current_molecule(-1), 
                               frame_opened(false), molecule_opened(false)
{}

PDBMolecules::~PDBMolecules()
{}

void PDBMolecules::closeMolecule()
{
    if (molecule_opened)
    {
        PDBMolecule &mol = mols[current_molecule];
        mol.closeFrame();
    
        if (mol.isEmpty())
        {
            //remove the empty molecule
            mols.removeAt(current_molecule);
            
            --current_molecule;
        }
    }
    
    molecule_opened = false;
}

PDBMolecule& PDBMolecules::nextMolecule()
{
    if (molecule_opened)
        this->closeMolecule();

    ++current_molecule;

    if (current_molecule >= mols.count())
    {
        mols.append(PDBMolecule());
        
        if (last_frame > 0)
            mols[current_molecule].openFrame(last_frame);
    }
        
    molecule_opened = true;
        
    return mols[current_molecule];
}

int PDBMolecules::openFrame(int frame_number)
{
    if (frame_number <= 0)
        frame_number = last_frame + 1;

    current_molecule = -1;
    
    for (QList<PDBMolecule>::iterator it = mols.begin();
         it != mols.end();
         ++it)
    {
        it->openFrame(frame_number);
    }
    
    frame_opened = true;
    last_frame = frame_number;

    return frame_number;
}

void PDBMolecules::closeFrame()
{
    for (QList<PDBMolecule>::iterator it = mols.begin();
         it != mols.end();
         ++it)
    {
        it->closeFrame();
    }
    
    frame_opened = false;
}

bool PDBMolecules::frameOpened() const
{
    return frame_opened;
}

bool PDBMolecules::moleculeOpened() const
{
    return molecule_opened;
}

const QList<PDBMolecule>& PDBMolecules::molecules() const
{
    return mols;
}

bool PDBMolecules::isEmpty() const
{
    return mols.isEmpty();
}

////////
//////// Implementation of PDB
////////

static const RegisterMetaType<PDB> r_pdb;

/** Serialise to a binary datastream */
QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDB &pdb)
{
    writeHeader(ds, r_pdb, 1);
    
    ds << static_cast<const IOBase&>(pdb);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDB &pdb)
{
    VersionID v = readHeader(ds, r_pdb);
    
    if (v == 1)
    {
        ds >> static_cast<IOBase&>(pdb);
    }
    else
        throw version_error(v, "1", r_pdb, CODELOC);
        
    return ds;
}

PDBParameters PDB::pdbparams;

/** Constructor */
PDB::PDB() : ConcreteProperty<PDB,IOBase>()
{}

/** Copy constructor */
PDB::PDB(const PDB &other) : ConcreteProperty<PDB,IOBase>(other)
{}

/** Destructor */
PDB::~PDB()
{}

/** Copy assignment operator */
PDB& PDB::operator=(const PDB &other)
{
    IOBase::operator=(other);
    return *this;
}

/** Comparison operator */
bool PDB::operator==(const PDB&) const
{
    return true;
}

/** Comparison operator */
bool PDB::operator!=(const PDB&) const
{
    return false;
}

class PDBResidue
{
public:
    PDBResidue()
    {}

    PDBResidue(const PDBAtom &atom)
         : resname(atom.resname), resnum(atom.resseq),
           icode(atom.icode), chainname(atom.chainid)
    {}

    ~PDBResidue()
    {}
    
    bool operator==(const PDBResidue &other) const
    {
        return resname == other.resname and resnum == other.resnum and
               icode == other.icode and chainname == other.chainname;
    }
    
    bool operator!=(const PDBResidue &other) const
    {
        return resname != other.resname or resnum != other.resnum or
               icode != other.icode or chainname != other.chainname;
    }
    
    bool isEmpty() const
    {
        return resname.isEmpty() and resnum.isNull() and icode.isEmpty()
                  and chainname.isEmpty();
    }
    
    ResName resname;
    ResNum resnum;
    QString icode;
    ChainName chainname;
};

/** Convert a list of PDBAtoms into a molecule */
static Molecule convert(const QList<PDBAtom> &pdbatoms,
                        const PropertyMap &map)
{
    if (pdbatoms.isEmpty())
        return Molecule();

    //get the manglers for the atom, residue, chain and segment names
    StringManglerPtr atommangler, resmangler, chainmangler, segmangler;

    PropertyName atommangler_property = map[PDB::parameters().atomNameMangler()];
    PropertyName resmangler_property = map[PDB::parameters().residueNameMangler()];
    PropertyName chainmangler_property = map[PDB::parameters().chainNameMangler()];
    PropertyName segmangler_property = map[PDB::parameters().segmentNameMangler()];
    
    if (atommangler_property.hasValue())
        atommangler = atommangler_property.value().asA<StringMangler>();
    
    if (resmangler_property.hasValue())
        resmangler = resmangler_property.value().asA<StringMangler>();
    
    if (chainmangler_property.hasValue())
        chainmangler = chainmangler_property.value().asA<StringMangler>();
    
    if (segmangler_property.hasValue())
        segmangler = segmangler_property.value().asA<StringMangler>();

    //editor for the molecule
    MolStructureEditor moleditor;

    //create the array that holds the location of each atom in the
    //molecule. This is because some of these atoms may be alternates,
    //so are not included directly in the molecule. This means that
    //the ith atom in pdbatoms may not be the ith atom in the molecule.
    //The 'atomlocations' array maps the ith atom in pdbatoms to 
    //the jth atom in the molecule (or -1 if this is an alternate)
    int natoms = pdbatoms.count();
    
    QVector<int> atomlocations(natoms);
    int *atomlocations_array = atomlocations.data();
    
    //create an array that holds the index of the residue to which
    //the ith atom in pdbatoms belongs (or -1 if this atom is an 
    //alternate or does not belong in any residue)
    QVector<int> resparent(natoms);
    int *resparent_array = resparent.data();

    //this is a list of all of the residues in the order they
    //appear in pdbatoms
    QList<PDBResidue> pdbresidues;

    //this contains the index of the primary atom of a set
    //of alternates (each alternate must have an identical idstring)
    QHash<QString,int> atom_idstrings;

    int atomidx = 0;
    int residx = 0;
    int current_res = -1;

    //create the 'old_res' residue - a new residue is created
    //whenever the residue information for the current atom is
    //different to that of the previous atom
    PDBResidue old_res(pdbatoms.at(0));
    
    if (not old_res.isEmpty())
    {
        pdbresidues.append(old_res);
        current_res = residx;
        
        ++residx;
    }

    //read through all of the atoms, and as a first step,
    //just assign each one to a place in the molecule and to a residue
    for (int i=0; i<natoms; ++i)
    {
        const PDBAtom &pdbatom = pdbatoms.at(i);
    
        if (not pdbatom.altloc.isEmpty())
        {
            //this atom has an alternate location ID
            if (atom_idstrings.contains(pdbatom.idstring))
            {
                //one of the other alternates of this atom has 
                //already been read, so we don't need to include
                //this atom in the molecule directly
                
                //though we do want to find the atom with the largest
                //occupancy...
                int alt_idx = atom_idstrings.value(pdbatom.idstring);
                const PDBAtom &alt_atom = pdbatoms.at(alt_idx);
                
                if (alt_atom.occupancy < pdbatom.occupancy)
                {
                    //this atom has a higher occupancy - make it the
                    //dominant alternative atom
                    atomlocations_array[i] = atomlocations_array[alt_idx];
                    atomlocations_array[alt_idx] = -1;
                    
                    resparent_array[i] = resparent_array[alt_idx];
                    resparent_array[alt_idx] = -1;
                    
                    atom_idstrings[pdbatom.idstring] = i;
                }
                else
                {
                    atomlocations_array[i] = -1;
                    resparent_array[i] = -1; 
                }
                
                //we've now finished with this atom
                continue;
            }
            else
            {
                //this is the first of a (possible) set of alternate
                //atoms. Record its location
                atom_idstrings.insert(pdbatom.idstring, i);
            }
        }

        //add this new atom to the molecule
        atomlocations_array[i] = atomidx;
        ++atomidx;
                
        PDBResidue pdbres(pdbatom);
                
        if (pdbres != old_res)
        {
            //this is a new residue as well
            if (not pdbres.isEmpty())
            {
                pdbresidues.append(pdbres);
                current_res = residx;
                ++residx;
            }
            else
                current_res = -1;

            old_res = pdbres;
        }
                
        resparent_array[i] = current_res;
    }
    
    //we now know the location of every residue in the molecule,
    //so lets create them all!
    QSet<ChainName> created_chains;
    QList<QString> pdbchains;
    
    foreach (const PDBResidue &pdbresidue, pdbresidues)
    {
        ResStructureEditor reseditor = moleditor.add(pdbresidue.resnum);
        
        if (not pdbresidue.resname.isEmpty())
            reseditor.rename( ResName(resmangler->mangle(pdbresidue.resname)) );
    
        if (not pdbresidue.chainname.isEmpty())
        {
            //this residue is part of a chain - add the residue to the chain
            if (not created_chains.contains(pdbresidue.chainname))
            {
                moleditor.add( ChainName(chainmangler->mangle(pdbresidue.chainname)) );
                created_chains.insert(pdbresidue.chainname);
                pdbchains.append(pdbresidue.chainname);
            }
        
            reseditor.reparent( ChainName(chainmangler->mangle(pdbresidue.chainname)) );
       }
    }
    
    //now we have created all of the residues, we can now create
    //all of the atoms
    atomidx = 0;
    QSet<QString> created_segments;
    QList<QString> pdbsegs;
    
    for (int i=0; i<natoms; ++i)
    {
        //is this an alternative atom?
        if (atomlocations_array[i] == -1)
        {
            //yes it is - we don't add these directly to the molecule
            continue;
        }
        
        const PDBAtom &pdbatom = pdbatoms.at(i);
        
        //add the atom to the molecule
        BOOST_ASSERT(atomlocations_array[i] == atomidx);
        ++atomidx;
        
        AtomStructureEditor atomeditor = moleditor.add( AtomNum(pdbatom.serial) );
        
        if (not pdbatom.name.isEmpty())
            atomeditor.rename( AtomName(atommangler->mangle(pdbatom.name)) );
            
        //place this atom into a residue (if it belongs in one)
        if (resparent_array[i] != -1)
            atomeditor.reparent( ResIdx(resparent_array[i]) );
            
        //place this atom into a segment (if it belongs in one)
        if (not pdbatom.segid.isEmpty())
        {
            if (not created_segments.contains(pdbatom.segid))
            {
                moleditor.add( SegName(segmangler->mangle(pdbatom.segid)) );
                created_segments.insert(pdbatom.segid);
                pdbsegs.append(pdbatom.segid);
            }
            
            atomeditor.reparent( SegName(segmangler->mangle(pdbatom.segid)) );
        }
    }
    
    //now we have the layout, use the supplied function to 
    //break this molecule down into CutGroups...
    CutFuncPtr cutfunc;
    
    PropertyName cutfunc_property = map[PDB::parameters().cuttingFunction()];
    
    if (cutfunc_property.hasValue())
        cutfunc = cutfunc_property.value().asA<CuttingFunction>();
        
    moleditor = (*cutfunc)(moleditor);

    //ok, we've now built the structure of the molecule, so commit it!
    Molecule molecule = moleditor.commit();

    //shall we save the coordinates?
    PropertyName coords_property = map[PDB::parameters().coordinates()];
    
    if (coords_property.hasSource())
    {
        //now build the coordinates - it is inefficient to do this
        //atom by atom, as the bounding sphere for each CutGroup would
        //be calculated after setting the coordinates of each atom. We
        //therefore have to place the coordinates all into a large
        //array, and then convert them into CoordGroups in one go
        QVector< QVector<Vector> > atomcoords( molecule.nCutGroups() );
        QVector<Vector> *atomcoords_array = atomcoords.data();
    
        for (CGIdx i(0); i<molecule.nCutGroups(); ++i)
        {
            atomcoords_array[i] = QVector<Vector>( molecule.data().info().nAtoms(i) );
        }

        for (int i=0; i<natoms; ++i)
        {
            //is this an alternative atom?
            if (atomlocations_array[i] == -1)
            {
                //yes it is - we need to add this to the alternative
                //atoms property...
                continue;
            }
        
            const PDBAtom &pdbatom = pdbatoms.at(i);
        
            const CGAtomIdx &cgatomidx = molecule.data().info()
                                           .cgAtomIdx(AtomIdx(atomlocations_array[i]));
        
            //set the coordinates
            atomcoords_array[cgatomidx.cutGroup()][cgatomidx.atom()]
                                = Vector(pdbatom.x, pdbatom.y, pdbatom.z);
        }
    
        molecule = molecule.edit()
                       .setProperty(coords_property.source(), AtomCoords(atomcoords))
                       .commit();
    }
    
    //shall we save the element types of the atom?
    PropertyName elements_property = map[PDB::parameters().element()];
    
    if (elements_property.hasSource())
    {
        AtomElements elements(molecule.data().info());
        
        for (int i=0; i<natoms; ++i)
        {
            if (atomlocations_array[i] == -1)
                continue;
                
            const PDBAtom &pdbatom = pdbatoms.at(i);
            
            const CGAtomIdx &cgatomidx = molecule.data().info()
                                           .cgAtomIdx(AtomIdx(atomlocations_array[i]));

            //save the element
            if (not pdbatom.element.isEmpty())
                elements.set(cgatomidx, Element(pdbatom.element));
            else
                //get the element from the first two characters of the 
                //atom name
                elements.set(cgatomidx, Element(pdbatom.name.mid(0,2).trimmed()));
        }
        
        molecule = molecule.edit()
                           .setProperty(elements_property.source(), elements)
                           .commit();
    }

    //shall we save the temperature (b) factors?
    PropertyName bfactor_property = map[PDB::parameters().bFactor()];
    
    if (bfactor_property.hasSource())
    {
        AtomFloatProperty bfactors( molecule.data().info() );
        
        for (int i=0; i<natoms; ++i)
        {
            if (atomlocations_array[i] == -1)
                continue;
                
            const PDBAtom &pdbatom = pdbatoms.at(i);
            
            const CGAtomIdx &cgatomidx = molecule.data().info()
                                           .cgAtomIdx(AtomIdx(atomlocations_array[i]));

            //save the b-factor
            bfactors.set(cgatomidx, pdbatom.tempfactor);
        }
        
        molecule = molecule.edit()
                           .setProperty(bfactor_property.source(), bfactors)
                           .commit();
    }
    
    //shall we save the atom's formal charge?
    PropertyName formalcharge_property = map[PDB::parameters().formalCharge()];
    
    if (formalcharge_property.hasSource())
    {
        AtomCharges charges( molecule.data().info() );
        
        for (int i=0; i<natoms; ++i)
        {
            if (atomlocations_array[i] == -1)
                continue;
                
            const PDBAtom &pdbatom = pdbatoms.at(i);
            
            const CGAtomIdx &cgatomidx = molecule.data().info()
                                           .cgAtomIdx(AtomIdx(atomlocations_array[i]));

            //save the charge
            charges.set(cgatomidx, pdbatom.charge * SireUnits::mod_electron);
        }
        
        molecule = molecule.edit()
                           .setProperty(formalcharge_property.source(), charges)
                           .commit();
    }
    
    //shall we save the original PDB atom name?
    PropertyName pdbatomname_property = map[PDB::parameters().pdbAtomName()];
    
    if (pdbatomname_property.hasSource())
    {     
        AtomStringProperty atomnames( molecule.data().info() );
        
        for (int i=0; i<natoms; ++i)
        {
            if (atomlocations_array[i] == -1)
                continue;
                
            const PDBAtom &pdbatom = pdbatoms.at(i);
            
            const CGAtomIdx &cgatomidx = molecule.data().info()
                                           .cgAtomIdx(AtomIdx(atomlocations_array[i]));

            //save the original atom name
            atomnames.set(cgatomidx, SireMol::cacheName(pdbatom.name));
        }
        
        molecule = molecule.edit()
                           .setProperty(pdbatomname_property.source(), atomnames)
                           .commit();
    }
    
    //shall we save the residue insertion code?
    PropertyName icode_property = map[PDB::parameters().iCode()];
    
    if (icode_property.hasSource())
    {
        ResStringProperty icodes( molecule.data().info() );
        
        for (int i=0; i<pdbresidues.count(); ++i)
        {
            const PDBResidue &pdbres = pdbresidues.at(i);
            
            icodes.set( ResIdx(i), SireMol::cacheName(pdbres.icode) );
        }
        
        molecule = molecule.edit()
                           .setProperty(icode_property.source(), icodes)
                           .commit();
    }
    
    //shall we save the original PDB residue names?
    PropertyName pdbresnames_property = map[PDB::parameters().pdbResidueName()];
    
    if (pdbresnames_property.hasSource())
    {
        ResStringProperty resnames( molecule.data().info() );
        
        for (int i=0; i<pdbresidues.count(); ++i)
        {
            const PDBResidue &pdbres = pdbresidues.at(i);
            
            resnames.set( ResIdx(i), SireMol::cacheName(pdbres.resname) );
        }
        
        molecule = molecule.edit()
                           .setProperty(pdbresnames_property.source(), resnames)
                           .commit();
    }
    
    //shall we save the original PDB chain names?
    PropertyName pdbchainnames_property = map[PDB::parameters().pdbChainName()];
    
    if (pdbchainnames_property.hasSource())
    {
        ChainStringProperty chainnames( molecule.data().info() );
        
        for (int i=0; i<pdbchains.count(); ++i)
        {
            chainnames.set( ChainIdx(i), SireMol::cacheName(pdbchains.at(i)) );
        }
        
        molecule = molecule.edit()
                           .setProperty(pdbchainnames_property.source(), chainnames)
                           .commit();
    }
    
    //shall we save the original PDB segment names?
    PropertyName pdbsegnames_property = map[PDB::parameters().pdbSegmentName()];
    
    if (pdbsegnames_property.hasSource())
    {
        SegStringProperty segnames( molecule.data().info() );
        
        for (int i=0; i<pdbsegs.count(); ++i)
        {
            segnames.set( SegIdx(i), SireMol::cacheName(pdbsegs.at(i)) );
        }
        
        molecule = molecule.edit()
                           .setProperty(pdbsegnames_property.source(), segnames)
                           .commit();
    }
    
    return molecule;
}

/** Convert a PDBMolecule to a Molecule */
static Molecule convert(const PDBMolecule &pdbmol, int first_frame,
                        int last_frame,
                        const PropertyMap &map)
{
    BOOST_ASSERT( first_frame <= last_frame );

    QList<int> available_frames = pdbmol.availableFrames();
    
    if (available_frames.isEmpty())
        return Molecule();
        
    //create a molecule from the first frame
    int first_available_frame = available_frames.takeFirst();
    
    BOOST_ASSERT( first_available_frame >= first_frame and
                  first_available_frame <= last_frame );
    
    Molecule molecule = convert(pdbmol.atoms(first_available_frame), map);
    
    if (not available_frames.isEmpty())
    {
        ///there are more frames to read...
        foreach (int framenum, available_frames)
        {
            Molecule next_frame = convert(pdbmol.atoms(framenum), map);
            
            //add this frame onto the molecule
        }
    }
    else if (first_available_frame != first_frame or
             first_frame != last_frame)
    {
        //tell the molecule that it exists only at a single
        //frame within an animation
    }
    
    return molecule;
}

/** Read a group of molecules from the data */
MoleculeGroup PDB::readMols(const QByteArray &data,
                            const PropertyMap &map) const
{
    // See the PDB format description for details on the format
    //
    // Version 3.1 <http://www.wwpdb.org/documentation/format3.1-20080211.pdf>
    // Version 2.3 <http://www.wwpdb.org/documentation/format2.3-0108-a4.pdf>
    //
    
    //connect a text stream to this bytearray so that we can read it as text
    QTextStream ts(data, QIODevice::ReadOnly | QIODevice::Text);

    PDBMolecules pdbmols;
    int linenum = -1;

    PDBMolecule *current_molecule = &(pdbmols.nextMolecule());
    
    QList<int> opened_frames;
    
    while (not ts.atEnd())
    {
        ++linenum;
    
        //read a line from the file
        QString line = ts.readLine();
        
        if (line.startsWith("END"))
        {
            QStringList words = line.split(" ", QString::SkipEmptyParts);

            if (words[0].length() == 3)
                //this is the end of the PDB file
                break;
        }
        
        if (line.startsWith("ATOM", Qt::CaseInsensitive) or
            line.startsWith("HETATM", Qt::CaseInsensitive))
        {
            if (not pdbmols.moleculeOpened())
                current_molecule = &(pdbmols.nextMolecule());
        
            current_molecule->addAtom( PDBAtom::readFromLine(line, linenum) );
        }
        else if (line.startsWith("TER", Qt::CaseInsensitive))
        {
            pdbmols.closeMolecule();
        }
        else if (line.startsWith("MODEL", Qt::CaseInsensitive))
        {
            //we've started a new model
            if (pdbmols.frameOpened())
            {
                throw SireIO::parse_error( QObject::tr(
                    "Encountered a 'MODEL' entry at line %1 "
                    "before the previous model has been terminated "
                    "(via a valid ENDMDL line).")
                        .arg(linenum), CODELOC );
            }
                
            bool ok;
            int frame_number = line.mid(10,4).toInt(&ok);
            
            if (not ok)
                frame_number = -1;
            
            opened_frames.append( pdbmols.openFrame(frame_number) );
        }
        else if (line.startsWith("ENDMDL", Qt::CaseInsensitive))
        {
            pdbmols.closeFrame();
        }
    }

    if (pdbmols.moleculeOpened())
    {
        pdbmols.closeMolecule();
    }
    
    if (pdbmols.frameOpened())
    {
        pdbmols.closeFrame();
    }

    if (pdbmols.isEmpty())
        //no molecules have been loaded!
        return MoleculeGroup();

    int first_frame, last_frame;

    //get the range of frames loaded
    if (not opened_frames.isEmpty())
    {
        qSort(opened_frames);
    
        first_frame = opened_frames.first();
        last_frame = opened_frames.last();
    }
    else
    {
        //this wasn't an animation
        first_frame = 1;
        last_frame = 1;
    }
    
    //we now need to convert each PDBMolecule into a real molecule
    MoleculeGroup molgroup;
    
    foreach (const PDBMolecule &pdbmol, pdbmols.molecules())
    {
        Molecule new_molecule = convert(pdbmol, first_frame, last_frame, map);
        
        if (new_molecule.nAtoms() != 0)
        {
            molgroup.add(new_molecule);
        }
    }

    return molgroup;
}

int PDB::writeMolecule(QTextStream &ts, const MoleculeView &molview,
                       int atomnum, const PropertyMap &map) const
{
    //get the manglers for the atom, residue, chain and segment names
    StringManglerPtr atommangler, resmangler, chainmangler, segmangler;

    PropertyName atommangler_property = map[PDB::parameters().atomNameMangler()];
    PropertyName resmangler_property = map[PDB::parameters().residueNameMangler()];
    PropertyName chainmangler_property = map[PDB::parameters().chainNameMangler()];
    PropertyName segmangler_property = map[PDB::parameters().segmentNameMangler()];
    
    PropertyName connectivity_property = map["connectivity"];
    
    if (atommangler_property.hasValue())
        atommangler = atommangler_property.value().asA<StringMangler>();
    
    if (resmangler_property.hasValue())
        resmangler = resmangler_property.value().asA<StringMangler>();
    
    if (chainmangler_property.hasValue())
        chainmangler = chainmangler_property.value().asA<StringMangler>();
    
    if (segmangler_property.hasValue())
        segmangler = segmangler_property.value().asA<StringMangler>();
    
    AtomSelection selected_atoms = molview.selection();
    
    Molecule mol(molview);
    
    //map of AtomIdx to PDB atomnum
    QHash<AtomIdx,int> atomidx_to_pdbnum;
    
    const AtomCoords &coords = mol.property(map[PDB::parameters().coordinates()])
                                  .asA<AtomCoords>();
        
    AtomElements elements;
    
    if (mol.hasProperty(map[PDB::parameters().element()]))
    {
        elements = mol.property(map[PDB::parameters().element()])
                      .asA<AtomElements>();
    }
    
    AtomFloatProperty bfactor;    
    
    if (mol.hasProperty(map[PDB::parameters().bFactor()]))
    {
        bfactor = mol.property(map[PDB::parameters().bFactor()])
                     .asA<AtomFloatProperty>();
    }
    
    AtomCharges charges;
    
    if (mol.hasProperty(map[PDB::parameters().formalCharge()]))
    {
        charges = mol.property(map[PDB::parameters().formalCharge()])
                     .asA<AtomCharges>();
    }
    
    AtomStringProperty pdbatomname;
    
    if (mol.hasProperty(map[PDB::parameters().pdbAtomName()]))
    {
        pdbatomname = mol.property(map[PDB::parameters().pdbAtomName()])
                         .asA<AtomStringProperty>();
    }
    
    ResStringProperty icode;
    
    if (mol.hasProperty(map[PDB::parameters().iCode()]))
    {
        icode = mol.property(map[PDB::parameters().iCode()])
                   .asA<ResStringProperty>();
    }
    
    ResStringProperty pdbresname;
    
    if (mol.hasProperty(map[PDB::parameters().pdbResidueName()]))
    {
        pdbresname = mol.property(PDB::parameters().pdbResidueName())
                        .asA<ResStringProperty>();
    }
    
    ChainStringProperty pdbchainname;
    
    if (mol.hasProperty(map[PDB::parameters().pdbChainName()]))
    {
        pdbchainname = mol.property(map[PDB::parameters().pdbChainName()])
                          .asA<ChainStringProperty>();
    }
    
    SegStringProperty pdbsegname;
    
    if (mol.hasProperty(map[PDB::parameters().pdbSegmentName()]))
    {
        pdbsegname = mol.property(map[PDB::parameters().pdbSegmentName()])
                        .asA<SegStringProperty>();
    }

    int natoms = mol.nAtoms();
    
    for (AtomIdx i(0); i<natoms; ++i)
    {
        if (not selected_atoms.selected(i))
            continue;
            
        Atom atom = mol.atom(i);
            
        PDBAtom pdbatom;
        
        ++atomnum;
        
        pdbatom.serial = atomnum;
        
        atomidx_to_pdbnum.insert(atom.index(), atomnum);
        
        if (pdbatomname.isEmpty() or pdbatomname[ atom.cgAtomIdx() ].isEmpty())
        {
            pdbatom.name = atommangler->mangle( atom.name() );
        }
        else
            pdbatom.name = pdbatomname[ atom.cgAtomIdx() ];

        if (not elements.isEmpty())
            pdbatom.element = elements[ atom.cgAtomIdx() ].symbol();
            
        const Vector &c = coords[ atom.cgAtomIdx() ];
        
        pdbatom.x = c.x();
        pdbatom.y = c.y();
        pdbatom.z = c.z();
        
        if (not charges.isEmpty())
            pdbatom.charge = charges[ atom.cgAtomIdx() ];
            
        if (not bfactor.isEmpty())
            pdbatom.tempfactor = bfactor[ atom.cgAtomIdx() ];
        else
            pdbatom.tempfactor = 1;
            
        if (atom.isWithinResidue())
        {
            Residue residue = atom.residue();
            
            if (pdbresname.isEmpty() or pdbresname[ residue.index() ].isEmpty())
                pdbatom.resname = resmangler->mangle(residue.name());
                
            else
                pdbatom.resname = pdbresname[ residue.index() ];

	    // JM May 11. PDB file format does not support more than 4 chars
	    if ( residue.number() > 9999 )
	      pdbatom.resseq = 9999;
	    else
	      pdbatom.resseq = residue.number();

            if (not icode.isEmpty())
                pdbatom.icode = icode[ residue.index() ];
                
            if (residue.isWithinChain())
            {
                Chain chain = residue.chain();
                
                if (pdbchainname.isEmpty() or pdbchainname[ chain.index() ].isEmpty())
                    pdbatom.chainid = chainmangler->mangle(chain.name());
                    
                else
                    pdbatom.chainid = pdbchainname[ chain.index() ];
            }
        }
        
        if (atom.isWithinSegment())
        {
            Segment segment = atom.segment();
            
            if (pdbsegname.isEmpty() or pdbsegname[ segment.index() ].isEmpty())
                pdbatom.segid = segmangler->mangle(segment.name());
                
            else
                pdbatom.segid = pdbsegname[ segment.index() ];
        }
        
        //write the atom to the file
        ts << pdbatom.writeToLine() << "\n";
    }
    
    if (mol.hasProperty("axis"))
    {
        const Property &axes_prop = mol.property("axis");
        
        if (axes_prop.isA<VariantProperty>())
        {
            AxisSet axes = axes_prop.asA<VariantProperty>().convertTo<AxisSet>();
            
            Vector v = axes.fromIdentity( Vector(0,0,0) );
            
            ++atomnum;
            PDBAtom pdbatom;
            pdbatom.serial = atomnum;
            pdbatom.name = "ORI";
            pdbatom.element = "Xx";
            pdbatom.x = v.x();
            pdbatom.y = v.y();
            pdbatom.z = v.z();
            
            ts << pdbatom.writeToLine() << "\n";
            
            v = axes.fromIdentity( Vector(1,0,0) );
            
            ++atomnum;
            pdbatom.serial = atomnum;
            pdbatom.name = "X";
            pdbatom.x = v.x();
            pdbatom.y = v.y();
            pdbatom.z = v.z();

            ts << pdbatom.writeToLine() << "\n";
            
            v = axes.fromIdentity( Vector(0,1,0) );
            
            ++atomnum;
            pdbatom.serial = atomnum;
            pdbatom.name = "Y";
            pdbatom.x = v.x();
            pdbatom.y = v.y();
            pdbatom.z = v.z();

            ts << pdbatom.writeToLine() << "\n";
            
            v = axes.fromIdentity( Vector(0,0,1) );
            
            ++atomnum;
            pdbatom.serial = atomnum;
            pdbatom.name = "Z";
            pdbatom.x = v.x();
            pdbatom.y = v.y();
            pdbatom.z = v.z();

            ts << pdbatom.writeToLine() << "\n";
        }
    }
    
    if (mol.hasProperty(connectivity_property))
    {
        const Connectivity &connectivity = mol.property(connectivity_property)
                                              .asA<Connectivity>();
                                              
        for (AtomIdx i(0); i<mol.nAtoms(); ++i)
        {
            QSet<AtomIdx> bonded_atoms = connectivity.connectionsTo(i);
            
            if (not bonded_atoms.isEmpty())
            {
                ts << "CONECT";
                ts.setFieldWidth(5);
                
                ts << atomidx_to_pdbnum.value(i);
                
                foreach (AtomIdx bonded_atom, bonded_atoms)
                {
                    ts << atomidx_to_pdbnum.value(bonded_atom);
                }
                
                ts << "\n";
            }
        }
    } 
    
    return atomnum;
}

//reserve at least 16 MB of space
static const int RESERVE_SIZE = 16 * 1024 * 1024;

/** Write a group of molecules to a bytearray */
QByteArray PDB::writeMols(const MoleculeGroup &molgroup,
                          const PropertyMap &map) const
{
    QByteArray data;
    data.reserve(RESERVE_SIZE);
    
    if (data.capacity() != RESERVE_SIZE)
    {
        qDebug() << CODELOC;
        qWarning() << "Could not reserve enough space!" << data.capacity();
    }
    
    QTextStream ts(&data, QIODevice::WriteOnly | QIODevice::Text);

    int nmols = molgroup.nMolecules();
    
    int atomid = 0;
    
    for (MolIdx i(0); i<nmols; ++i)
    {
        atomid = this->writeMolecule(ts, molgroup[i], atomid, map);
        
        if (i < nmols - 1)
            ts << "TER\n";
    }

    return data;
}

/** Write a group of molecules to a bytearray */
QByteArray PDB::writeMols(const Molecules &molecules,
                          const PropertyMap &map) const
{
    QByteArray data;
    data.reserve(RESERVE_SIZE);
    
    if (data.capacity() != RESERVE_SIZE)
    {
        qDebug() << CODELOC;
        qWarning() << "Could not reserve enough space!" << data.capacity();
    }

    QTextStream ts(&data, QIODevice::WriteOnly | QIODevice::Text);

    int atomid = 0;
    
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        atomid = this->writeMolecule(ts, *it, atomid, map);
        
        Molecules::const_iterator it2 = it;
        ++it2;
        
        if ( it2 != molecules.constEnd())
            ts << "TER\n";
    }

    return data;
}

const char* PDB::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDB>() );
}
