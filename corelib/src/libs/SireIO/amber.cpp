/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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


#include <QFile>
#include <QTextStream>
#include <QHash>

#include "amber.h"

#include "SireMol/element.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomelements.h"
#include "SireMol/connectivity.h"
#include "SireMol/selector.hpp"

#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/reseditor.h"
#include "SireMol/atomeditor.h"
#include "SireMol/cgatomidx.h"
#include "SireMol/residuecutting.h"
#include "SireMol/atomcutting.h"
#include "SireMol/molidx.h"
#include "SireMol/atomidx.h"

#include "SireMol/amberparameters.h"

#include "SireCAS/trigfuncs.h"

#include "SireMM/ljparameter.h"
#include "SireMM/atomljs.h"
#include "SireMM/internalff.h"
#include "SireMM/cljnbpairs.h"

#include "SireVol/cartesian.h"
#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireMaths/maths.h"
#include "SireUnits/units.h"

#include "SireBase/tempdir.h"
#include "SireBase/findexe.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMove/internalmove.h"
#include "SireMove/flexibility.h"

using namespace SireIO;
using namespace SireMol;
using namespace SireMM;
using namespace SireMaths;
using namespace SireCAS;
using namespace SireMove;
using namespace SireUnits;
using namespace SireVol;
using namespace SireStream;
using namespace SireBase;

/** FortranFormat is an internal class that holds the information about a
    format used in a top file entry
*/
class FortranFormat
{
public:
    FortranFormat() : repeat(0), type("I"), size(0), decimal(0)
    {}

    FortranFormat(int r, const QString &t, int s, int d)
        : repeat(r), type(t), size(s), decimal(d)
    {}

    ~FortranFormat()
    {}

    /// 10I8 'repeat' 'type' 'size', 'decimal' is 0
    // 5E16.8 'repeat' is 5, 'type' is 'E' 'size' is 16 'decimal' is 8
    int repeat;
    QString type;
    int size;
    int decimal;
};

/** enumerates the FLAGS in a TOP file*/
enum { UNKNOWN = 0, //a flag that is not known, e.g. in newer format top files
       TITLE = 1, //a TITLE flag in a top file
       POINTERS = 2, // a POINTERS flag in a top file
       ATOM_NAME = 3, //a ATOM_NAME flag in a top file
       CHARGE = 4, //a CHARGE flag in a top file
       MASS = 5, // a MASS flag
       ATOM_TYPE_INDEX = 6, //
       NUMBER_EXCLUDED_ATOMS = 7,
       NONBONDED_PARM_INDEX = 8,
       RESIDUE_LABEL = 9,
       RESIDUE_POINTER = 10,
       BOND_FORCE_CONSTANT = 11,
       BOND_EQUIL_VALUE = 12,
       ANGLE_FORCE_CONSTANT = 13,
       ANGLE_EQUIL_VALUE = 14,
       DIHEDRAL_FORCE_CONSTANT = 15,
       DIHEDRAL_PERIODICITY = 16,
       DIHEDRAL_PHASE = 17,
       SOLTY = 18,
       LENNARD_JONES_ACOEF = 19,
       LENNARD_JONES_BCOEF = 20,
       BONDS_INC_HYDROGEN = 21,
       BONDS_WITHOUT_HYDROGEN = 22,
       ANGLES_INC_HYDROGEN = 23,
       ANGLES_WITHOUT_HYDROGEN = 24,
       DIHEDRALS_INC_HYDROGEN = 25,
       DIHEDRALS_WITHOUT_HYDROGEN = 26,
       EXCLUDED_ATOMS_LIST = 27,
       HBOND_ACOEF = 28,
       HBOND_BCOEF = 29,
       HBCUT = 30,
       AMBER_ATOM_TYPE = 31,
       TREE_CHAIN_CLASSIFICATION = 32,
       JOIN_ARRAY = 33,
       IROTAT = 34,
       SOLVENT_POINTERS = 35,
       ATOMS_PER_MOLECULE = 36,
       BOX_DIMENSIONS = 37,
       RADIUS_SET = 38,
       RADII = 39,
       SCREEN = 40,
       ATOMIC_NUMBER = 41,
       SCEE_SCALE_FACTOR = 42,
       SCNB_SCALE_FACTOR = 43
     };

/** enumerates the POINTERS in a TOP file */
enum { NATOM = 0, // total number of atoms
       NTYPES = 1, // total number of distinct atom types
       NBONH = 2, // number of bonds containing hydrogen
       MBONA = 3, // number of bonds not containing hydrogen
       NTHETH = 4, // number of angles containing hydrogen
       MTHETA = 5, // number of angles not containing hydrogen
       NPHIH = 6, // number of dihedrals containing hydrogen
       MPHIA = 7, // number of dihedrals not containing hydrogen
       NHPARM = 8, // currently not used
       NPARM = 9, // currently  not used
       NEXT = 10, // number of excluded atoms
       NRES = 11, // number of residues
       NBONA = 12, // MBONA + number of constraint bonds
       NTHETA = 13, // MTHETA + number of constraint angles
       NPHIA = 14, // MPHIA + number of constraint dihedrals
       NUMBND = 15, // number of unique bond types
       NUMANG = 16, // number of unique angle types
       NPTRA = 17, // number of unique dihedral types
       NATYP = 18, // number of atom types in parameter file, see SOLTY below
       NPHB = 19, // number of distinct 10-12 hydrogen bond pair types
       IFPERT = 20, // set to 1 if perturbation info is to be read in
       NBPER = 21, // number of bonds to be perturbed
       NGPER = 22, // number of angles to be perturbed
       NDPER = 23, // number of dihedrals to be perturbed
       MBPER = 24, // number of bonds with atoms completely in perturbed group
       MGPER = 25, // number of agnles with atoms completely in perturbed group
       MDPER = 26, // number of dihedrals with atoms completely in perturbed groups
       IFBOX = 27, // set to 1 if standard periodic box, 2 when truncated octahedral
       NMXRS = 28, // number of atoms in the largest residue
       IFCAP = 29 // set to 1 if the CAP option from edit was specified
     };

/** enumerate SOLVENT pointers*/
enum { IPTRES = 0, //final residue that is considered part of the solute,
       //reset in sander and gibbs
       NSPM = 1,   //total number of molecules
       NSPOL = 2   //the first solvent "molecule"
     };

/** Enumerate cutting scheme **/
enum {PERRESIDUE=0,
      PERATOM=1
     };


// The partial charges in the top file are not in electrons
static const double AMBERCHARGECONV = 18.2223;

static const double AMBER14COUL = 1.0 / 1.2 ;
static const double AMBER14LJ = 0.50 ;

/** Processes a line that starts with %FLAG*/
static void processFlagLine(const QStringList &words, int &flag)
{
    if (words.count() < 2)
    {
        qDebug() << "WARNING: No flag given for this section of the topology file?";
        flag = UNKNOWN;
        return;
    }

    //%FLAG TITLE
    if (words[1] == "TITLE")
        flag = TITLE;
    else if (words[1] == "POINTERS")
        flag = POINTERS;
    else if (words[1] == "ATOM_NAME")
        flag = ATOM_NAME;
    else if (words[1] == "CHARGE")
        flag = CHARGE;
    else if (words[1] == "MASS")
        flag = MASS;
    else if (words[1] == "ATOM_TYPE_INDEX")
        flag = ATOM_TYPE_INDEX;
    else if (words[1] == "NUMBER_EXCLUDED_ATOMS")
        flag = NUMBER_EXCLUDED_ATOMS;
    else if (words[1] == "NONBONDED_PARM_INDEX")
        flag = NONBONDED_PARM_INDEX;
    else if (words[1] == "RESIDUE_LABEL")
        flag = RESIDUE_LABEL;
    else if (words[1] == "RESIDUE_POINTER")
        flag = RESIDUE_POINTER;
    else if (words[1] == "BOND_FORCE_CONSTANT")
        flag = BOND_FORCE_CONSTANT;
    else if (words[1] == "BOND_EQUIL_VALUE")
        flag = BOND_EQUIL_VALUE;
    else if (words[1] == "ANGLE_FORCE_CONSTANT")
        flag = ANGLE_FORCE_CONSTANT;
    else if (words[1] == "ANGLE_EQUIL_VALUE")
        flag = ANGLE_EQUIL_VALUE;
    else if (words[1] == "DIHEDRAL_FORCE_CONSTANT")
        flag = DIHEDRAL_FORCE_CONSTANT;
    else if (words[1] == "DIHEDRAL_PERIODICITY")
        flag = DIHEDRAL_PERIODICITY;
    else if (words[1] == "DIHEDRAL_PHASE")
        flag = DIHEDRAL_PHASE;
    else if (words[1] == "SOLTY")
        flag = SOLTY;
    else if (words[1] == "LENNARD_JONES_ACOEF")
        flag = LENNARD_JONES_ACOEF;
    else if (words[1] == "LENNARD_JONES_BCOEF")
        flag = LENNARD_JONES_BCOEF;
    else if (words[1] == "BONDS_INC_HYDROGEN")
        flag = BONDS_INC_HYDROGEN;
    else if (words[1] == "BONDS_WITHOUT_HYDROGEN")
        flag = BONDS_WITHOUT_HYDROGEN;
    else if (words[1] == "ANGLES_INC_HYDROGEN")
        flag = ANGLES_INC_HYDROGEN;
    else if (words[1] == "ANGLES_WITHOUT_HYDROGEN")
        flag = ANGLES_WITHOUT_HYDROGEN;
    else if (words[1] == "DIHEDRALS_INC_HYDROGEN")
        flag = DIHEDRALS_INC_HYDROGEN;
    else if (words[1] == "DIHEDRALS_WITHOUT_HYDROGEN")
        flag = DIHEDRALS_WITHOUT_HYDROGEN;
    else if (words[1] == "EXCLUDED_ATOMS_LIST")
        flag = EXCLUDED_ATOMS_LIST;
    else if (words[1] == "HBOND_ACOEF")
        flag = HBOND_ACOEF;
    else if (words[1] == "HBOND_BCOEF")
        flag = HBOND_BCOEF;
    else if (words[1] == "HBCUT")
        flag = HBCUT;
    else if (words[1] == "AMBER_ATOM_TYPE")
        flag = AMBER_ATOM_TYPE;
    else if (words[1] == "TREE_CHAIN_CLASSIFICATION")
        flag = TREE_CHAIN_CLASSIFICATION;
    else if (words[1] == "JOIN_ARRAY")
        flag = JOIN_ARRAY;
    else if (words[1] == "IROTAT")
        flag = IROTAT;
    else if (words[1] == "SOLVENT_POINTERS")
        flag = SOLVENT_POINTERS;
    else if (words[1] == "ATOMS_PER_MOLECULE")
        flag = ATOMS_PER_MOLECULE;
    else if (words[1] == "BOX_DIMENSIONS")
        flag = BOX_DIMENSIONS;
    else if (words[1] == "RADIUS_SET")
        flag = RADIUS_SET;
    else if (words[1] == "RADII")
        flag = RADII;
    else if (words[1] == "SCREEN")
        flag = SCREEN;
    else if (words[1] == "ATOMIC_NUMBER")
        flag = ATOMIC_NUMBER;
    else if (words[1] == "SCEE_SCALE_FACTOR")
        flag = SCEE_SCALE_FACTOR;
    else if (words[1] == "SCNB_SCALE_FACTOR")
        flag = SCNB_SCALE_FACTOR;
    else if (words[1] == "IPOL")
        flag = UNKNOWN;
    else
    {
        qDebug() << "WARNING: Skipping section" << words[1] << "in the topology file.";
        flag = UNKNOWN;
    }
}

/** Processes a line that starts with %FORMAT*/
static void processFormatLine(const QStringList &words, FortranFormat &format)
{
    //%FORMAT(20a4)
    QString tmp = words[0];

    tmp = tmp.remove( QString("%FORMAT(") );
    tmp = tmp.remove( QChar(')') ).toUpper();

    // QRegExp must match either (a,I,E) (case insensitive)
    QString match = " ";
    QRegExp rx("[AIE]");

    if (tmp.toUpper().contains(rx))
        match = rx.cap(0);
    else
    {
        qDebug() << "The format line" << words[0] << "cannot be read";
        throw SireError::program_bug( QObject::tr(
                            "The format line '%1' does not contain a supported format type")
                                      .arg(tmp), CODELOC);
    }

    QStringList elements = tmp.split(match);

    bool ok;
    int repeat = elements[0].toInt(&ok);

    QStringList elem2;
    int size;
    int decimal;

    if (elements[1].contains("."))
    {
        elem2 =  elements[1].split(".");
        size = elem2[0].toInt(&ok);
        decimal = elem2[1].toInt(&ok);
    }
    else
    {
        size = elements[1].toInt(&ok);
        decimal = 0;
    }

    if (match == "A")
        match = "a";

    format.repeat = repeat;
    format.type = match;
    format.size = size;
    format.decimal = decimal;
}

/** Processes the line that contains the VERSION of the top file*/
/*static void processVersionLine(const QStringList &words, QString &version)
{
    //%VERSION VERSION_STAMP = V0001.000 DATE = 05/22/06  12:10:21
    version = words[3];
}*/

/** Processes the title of the top file*/
static void processTitleLine(const QString &line, const FortranFormat &format,
                             QString &title)
{
    //FORMAT(20a4) reads up to 20 times 4 characters and treat them as a string
    // However, FORMAT is ignored in the TOP file, so we merely set title to line
    //qDebug() << "TITLE LINE IS " << line;
    title = line;
}

/** Processes a line of INTEGER in the top file */
static void processIntegerLine(const QString &line, const FortranFormat &format,
                               QList<int> &intarray)
{
    //FORMAT(10i8) reads up to 10 times 8 characters and treat them as an integer
    if ( format.type != "I")
        throw SireError::program_bug( QObject::tr(
                    "Format '%1' is not supported, should be I").arg(format.type), CODELOC);

    int count = 0;

    for (int i= 0 ; i < line.size() ; i = i + format.size )
    {
        if ( ++count > format.repeat || (i + format.size) > line.size() )
            break;

        QString element = line.mid(i, format.size );
        bool ok;
        int data = element.toInt(&ok);
        intarray.append(data);
    }
}

/** Processes a line of STRINGS in the top file */
static void processStringLine(const QString &line, const FortranFormat &format,
                              QStringList &stringarray)
{
    //FORMAT(20a4)
    if ( format.type != "a")
        throw SireError::program_bug( QObject::tr(
                "Format '%1' is not supported, should be a").arg(format.type), CODELOC);

    int count = 0;

    for (int i= 0 ; i < line.size() ; i = i + format.size )
    {
        if ( ++count > format.repeat || (i + format.size) > line.size() )
            break;

        QString str = line.mid(i, format.size );
        // Should we remove trailing spaces?
        str = str.trimmed();
        // Terminate if empty string

        if (str.size() == 0)
            break;

        stringarray.append(str);
    }
}

/** Processes a line of type DOUBLE in the top file*/
static void processDoubleLine(const QString &line, const FortranFormat &format,
                              QList<double> &doublearray)
{
    //FORMAT(5E16.8)
    if ( format.type != "E")
        throw SireError::program_bug( QObject::tr(
                "Format '%1' is not supported, should be E").arg(format.type), CODELOC);

    int count = 0;

    for (int i= 0 ; i < line.size() ; i = i + format.size )
    {
        if ( ++count > format.repeat || (i + format.size) > line.size() )
            break;

        QString element = line.mid(i, format.size );
        bool ok;
        double data = element.toDouble(&ok);

        //qDebug() << "element " << element << "Data " << data;
        doublearray.append(data);
    }
}

static void setAtomParameters(AtomEditor &editatom, MolEditor &editmol,
                              const QList<double> &crdCoords,
                              const PropertyName &coords_property,
                              const QList<int> &element,
                              const PropertyName &element_property,
                              const QList<double> &charge,
                              const PropertyName &charge_property,
                              const QList<double> &mass,
                              const PropertyName &mass_property,
                              const QList<int> &atom_type_index,
                              const QList<int> &nb_parm_index,
                              const QList<double> &lj_a_coeff,
                              const QList<double> &lj_b_coeff,
                              const PropertyName &lj_property,
                              const QStringList &amber_type,
                              const PropertyName &ambertype_property,
                              const QList<int> &pointers)
{
    //AtomEditor editatom = editmol.atom(AtomIdx(atomIndex));

    int atomNumber = editatom.number().value();

    //qDebug() << " Coordinates for number..." << atomNumber;
    // set the coordinates. index into coordinates array is idx = 3 * ( atom.number().value() -1)
    int array_index = 3 * ( atomNumber - 1 );

    //qDebug() << " array_index " << array_index;
    //qDebug() << " crdCoords.size() " << crdCoords.size();
    //qDebug() << " crdCoords[ array_index ] " << crdCoords[array_index];
    //qDebug() << " crdCoords[ array_index + 1 ] " << crdCoords[array_index + 1 ];
    //qDebug() << " crdCoords[ array_index + 2 ] " << crdCoords[array_index + 2 ];

    Vector coords = Vector(crdCoords[array_index],
                           crdCoords[array_index+1],
                           crdCoords[array_index+2]);

    editatom.setProperty( coords_property.source(), coords);

    //qDebug() << " Coords " << coords.x() << coords.y() << coords.z() ;

    //qDebug() << " Charges..." << charge[ atomNumber - 1 ] << " " <<  charge[atomNumber - 1] / AMBERCHARGECONV  ;
    // set the charges
    SireUnits::Dimension::Charge chg = ( charge[atomNumber - 1] / AMBERCHARGECONV )
                                       * mod_electron;
    editatom.setProperty( charge_property, chg);

    //qDebug() << " Masses...";
    // set the masses
    SireUnits::Dimension::MolarMass ms = mass[atomNumber - 1] * g_per_mol;
    editatom.setProperty( mass_property.source(), ms);

    // set the element (if it exists)
    if (not element.isEmpty())
    {
        Element elem;

        // Infer element from mass if atomic number is negative, which
        // can be the case for files generated by Acellara's "parameterize".
        if (element[atomNumber-1] < 0 and mass[atomNumber-1] > 0)
        {
            elem = Element::elementWithMass(mass[atomNumber-1] * g_per_mol);
        }
        else
        {
            elem = Element(element[atomNumber-1]);
        }

        editatom.setProperty( element_property.source(), elem );
    }

    // set the LJ parameters
    // For atom 'i' first get the 'itype' from atom_type_index
    // Then lookup the index 'inbparams' for an interaction of 'itype' with 'itype'
    // in nb_parm_index which is calculated as
    // 'inbparams' = pointers[NTYPES]*(itype-1)+itype
    // Then lookup the values of A & B in lj_a_coeff and lj_b_coeff
    // iAcoef = lj_a_coeff[inbparams] iBcoef =  lj_b_coeff[inbparams]

    //qDebug() << " LJ...";
    int itype = atom_type_index[ atomNumber - 1 ];
    int inbparams = nb_parm_index[ pointers[NTYPES] * (itype - 1) + itype - 1 ];
    double iAcoef = lj_a_coeff[ inbparams - 1 ];
    double iBcoef = lj_b_coeff[ inbparams - 1 ];
    double sigma, epsilon, rstar;

    // is there a SMALL?
    if (iAcoef < 1E-10)
    {
        sigma = 0.0;
        epsilon = 0.0;
        rstar = 0.0;
    }
    else
    {
        // and convert lj_a_coeff & lj_b_coeff into angstroms and kcal/mol-1
        sigma = std::pow( iAcoef / iBcoef ,  1/6. );
        epsilon = pow_2( iBcoef ) / (4*iAcoef);
        rstar = (sigma/2.)* std::pow(2.0, 1/6. ) ;
    }

    //qDebug() << " Atom " << atomNumber - 1 << " itype " << itype << " inbparams "
    //         << inbparams << " iAcoef " << iAcoef << " iBcoef " << iBcoef << " sigma "
    //         << sigma << " epsilon " << epsilon << " rstar " << rstar ;

    // Note that the amber par files give rstar=(sigma/2)*(2**(1/6.)) instead of sigma
    LJParameter lj( sigma * angstrom, epsilon * kcal_per_mol);
    editatom.setProperty( lj_property.source(), lj);

    //qDebug() << " Type...";

    // set the Amber atom type
    QString ambertype = amber_type[ atomNumber - 1 ];
    editatom.setProperty( ambertype_property.source(), ambertype);

    editmol = editatom.molecule();
}

/** Set the connectivity property of molecule editmol*/
static int setConnectivity(MolEditor &editmol, int pointer,
                           const QList<int> &bondsArray,
                           ConnectivityEditor &connectivity,
                           const PropertyName &connectivity_property,
                           int last_idx = -1)
{
    int atomStart = editmol.atoms()(0).number().value();
    int atomEnd = editmol.atoms()(-1).number().value();
    int start_idx = last_idx;
    bool got_to_end = true;

    /*if (last_idx < 0)
        start_idx = 0;
    else if (last_idx >= pointer)
        return last_idx;*/

    start_idx = 0;

    for (int i=start_idx; i < pointer ; ++i)
    {
        int index0 = bondsArray[ 3*i ] / 3 + 1 ;
        int index1 = bondsArray[ 3*i + 1 ] / 3 + 1 ;

        if ( index0 < atomStart or index1 < atomStart )
        {
            //we haven't yet reached the indices for this molecule
            continue;
        }
        else if ( index0 > atomEnd or index1 > atomEnd )
        {
            //the indices are sequential, so we have finished with this molecule
            //(assumes that the indicies refer to increasing molecule, i.e. sorted)
            last_idx = i+1;
            got_to_end = false;
            continue;
        }
        else
        {
            AtomNum number0 = AtomNum( index0 );
            AtomNum number1 = AtomNum( index1 );
            AtomIdx atom0 = editmol.select( number0 ).index();
            AtomIdx atom1 = editmol.select( number1 ).index();
            connectivity.connect( atom0, atom1 );
        }
    }

    if (got_to_end)
        last_idx = pointer;

    editmol.setProperty( connectivity_property.source(), connectivity.commit() );

    return last_idx;
}

/** Set the property bonds for molecule editmol*/
static int setBonds(MolEditor &editmol, int pointer,
                    const QList<int> &bondsArray,
                    const QList<double> &bond_force_constant,
                    const QList<double> &bond_equil_value,
                    TwoAtomFunctions &bondfuncs,
                    const PropertyName &bond_property,
                    AmberParameters &amberparams,
                    const PropertyName &amberparameters_property,
                    int last_idx=-1)
{
    int atomStart = editmol.atoms()(0).number().value();
    int atomEnd = editmol.atoms()(-1).number().value();

    int start_idx = last_idx;
    /*if (last_idx < 0)
        start_idx = 0;
    else if (last_idx >= pointer)
        return last_idx;*/

    bool got_to_end = true;

    start_idx = 0;

    for (int i=start_idx; i<pointer; ++i)
    {
        int index0 = bondsArray[ i*3 ] / 3 + 1 ;
        int index1 = bondsArray[ i*3 + 1] / 3 + 1 ;

        if ( index0 < atomStart or index1 < atomStart )
        {
            //we haven't yet reached the indices for this molecule
            continue;
        }
        else if ( index0 > atomEnd or index1 > atomEnd )
        {
            //the indices are sequential, so we have finished with this molecule
            //(assumes that the indicies refer to increasing molecule, i.e. sorted)
            last_idx = i+1;
            got_to_end = false;
            continue;
        }
        else
        {
            int paramIndex = bondsArray[ i*3 + 2 ];

            AtomNum number0 = AtomNum( index0 );
            AtomNum number1 = AtomNum( index1 );
            AtomIdx atom0 = editmol.select( number0 ).index();
            AtomIdx atom1 = editmol.select( number1 ).index();

            Symbol r = InternalPotential::symbols().bond().r();

            double k = bond_force_constant[ paramIndex -1 ];
            double r0 = bond_equil_value[ paramIndex - 1 ];

            Expression bondfunc = k * SireMaths::pow_2(r - r0);
            bondfuncs.set( atom0, atom1, bondfunc );

            BondID bond = BondID(atom0, atom1);
            amberparams.add( bond, k, r0);
        }
    }

    if (got_to_end)
        last_idx = pointer;

    editmol.setProperty( bond_property.source(), bondfuncs );

    editmol.setProperty( amberparameters_property.source(), amberparams);

    return last_idx;
}

static int setAngles(MolEditor &editmol, int pointer,
                     const QList<int> &anglesArray,
                     const QList<double> &ang_force_constant,
                     const QList<double> &ang_equil_value,
                     ThreeAtomFunctions &anglefuncs,
                     const PropertyName &angle_property,
                     AmberParameters &amberparams,
                     const PropertyName &amberparameters_property,
                     int last_idx=-1)
{
    //QSet<AtomNum> moleculeAtomNumbers = _pvt_selectAtomsbyNumber(editmol);
    int atomStart = editmol.atoms()(0).number().value();
    int atomEnd = editmol.atoms()(-1).number().value();

    int start_idx = last_idx;
    /*if (last_idx < 0)
        start_idx = 0;
    else if (last_idx >= pointer)
        return last_idx;*/

    bool got_to_end = true;

    start_idx = 0;

    for (int i=start_idx; i<pointer; ++i)
    {
        int index0 = anglesArray[ i*4 ] / 3 + 1 ;
        int index1 = anglesArray[ i*4 + 1] / 3 + 1 ;
        int index2 = anglesArray[ i*4 + 2] / 3 + 1 ;

        if ( index0 < atomStart or index1 < atomStart or index2 < atomStart )
        {
            //we haven't yet reached the indices for this molecule
            continue;
        }
        else if ( index0 > atomEnd or index1 > atomEnd or index2 > atomEnd )
        {
            //the indices are sequential, so we have finished with this molecule
            //(assumes that the indicies refer to increasing molecule, i.e. sorted)
            last_idx = i+1;
            got_to_end = false;
            continue;
        }
        else
        {
            AtomNum number0 = AtomNum( index0 );
            AtomNum number1 = AtomNum( index1 );
            AtomNum number2 = AtomNum( index2 );

            int paramIndex = anglesArray[ i*4 + 3 ];
            AtomIdx atom0 = editmol.select( number0 ).index();
            AtomIdx atom1 = editmol.select( number1 ).index();
            AtomIdx atom2 = editmol.select( number2 ).index();

            Symbol theta = InternalPotential::symbols().angle().theta();

            double k = ang_force_constant[ paramIndex - 1 ];
            double theta0 = ang_equil_value[ paramIndex - 1 ];// radians

            Expression anglefunc = k * SireMaths::pow_2(theta - theta0);
            anglefuncs.set( atom0, atom1, atom2, anglefunc );

            AngleID angle = AngleID(atom0, atom1, atom2);
            amberparams.add( angle, k, theta0);
        }
    }

    if (got_to_end)
        last_idx = pointer;

    editmol.setProperty( angle_property.source(), anglefuncs );

    editmol.setProperty( amberparameters_property.source(), amberparams);

    return last_idx;
}

static int setDihedrals(MolEditor &editmol, int pointer,
                         const QList<int> &dihedralsArray,
                         const QList<double> &dih_force_constant,
                         const QList<double> &dih_periodicity,
                         const QList<double> &dih_phase,
                         const QList<double> &sceefactor,
                         const QList<double> &scnbfactor,
                         FourAtomFunctions &dihedralfuncs,
                         const PropertyName &dihedral_property,
                         FourAtomFunctions &improperfuncs,
                         const PropertyName &improper_property,
                         QHash<AtomNum, QList<AtomNum> > &atoms14,
                         QHash<AtomNum, QHash<AtomNum, double> > &atoms14sclee,
                         QHash<AtomNum, QHash<AtomNum, double> > &atoms14sclnb,
                         double coul_14scl, double lj_14scl,
                         AmberParameters &amberparams,
                         const PropertyName &amberparameters_property,
                         int last_idx)
{
    //QSet<AtomNum> moleculeAtomNumbers = _pvt_selectAtomsbyNumber(editmol);
    int atomStart = editmol.atoms()(0).number().value();
    int atomEnd = editmol.atoms()(-1).number().value();

    QHash<DofID,Expression> improper_hash;
    QHash<DofID,Expression> dihedral_hash;

    int start_idx = last_idx;
    /*if (last_idx < 0)
        start_idx = 0;
    else if (last_idx >= pointer)
        return last_idx;*/

    bool got_to_end = true;

    start_idx = 0;

    for (int i=start_idx; i<pointer; ++i)
    {
        bool ignored = false;
        bool improper = false;

        //qDebug() << " RAW IDXs " << dihedralsArray[ i ] << dihedralsArray[ i +1 ]
        //         << dihedralsArray[ i + 2 ] << dihedralsArray[ i + 3]
        //         << dihedralsArray[ i + 4];

        int index0 = dihedralsArray[ i*5 ] ;
        index0 = index0 / 3 + 1 ;

        int index1 = dihedralsArray[ i*5 + 1 ] ;
        index1 = index1 / 3 + 1 ;

        int index2 = dihedralsArray[ i*5 + 2 ] ;

        // Note that if index2 is negative it indicates that end group interactions
        // re ignored (= non bonded) this could be because this quad of atoms has
        // already been set in the top file (for multi term dihedrals)
        // or because the 1,4 interaction has already been counted (for ring systems).
        if (index2 < 0)
        {
            ignored = true;
            index2 = -index2 ;
        }

        index2 = index2 / 3 + 1 ;

        int index3 = dihedralsArray[ i*5 + 3 ] ;

        // Note that if index3 is negative, it indicates that the dihedral is an improper
        if (index3 < 0)
        {
            improper = true;
            index3 = -index3;
        }

        index3 = index3 / 3 + 1 ;

        if ( index0 < atomStart or index1 < atomStart or index2 < atomStart or index3 < atomStart )
        {
            //we haven't yet reached the indices for this molecule
            continue;
        }
        else if ( index0 > atomEnd or index1 > atomEnd or index2 > atomEnd or index3 > atomEnd )
        {
            //the indices are sequential, so we have finished with this molecule
            //(assumes that the indicies refer to increasing molecule, i.e. sorted)
            last_idx = i+1;
            got_to_end = false;
            continue;
        }
        else
        {
            int paramIndex = dihedralsArray[ i*5 + 4 ];

            AtomNum number0 = AtomNum( index0 );
            AtomNum number1 = AtomNum( index1 );
            AtomNum number2 = AtomNum( index2 );
            AtomNum number3 = AtomNum( index3 );

            Symbol phi = InternalPotential::symbols().dihedral().phi();

            double k = dih_force_constant[ paramIndex - 1 ];// kcal_per_mol
            double periodicity = dih_periodicity[ paramIndex - 1] * radians.to(radians);
            double phase = dih_phase[ paramIndex - 1];

            // Assume default values for 14 scaling factors
            // Note that these are NOT inversed after reading from input
            double sclee14 = 1/coul_14scl;
            double sclnb14 = 1/lj_14scl;
            if (sceefactor.size() > 0)
                sclee14 = sceefactor[ paramIndex - 1 ];
            if (scnbfactor.size() > 0)
                sclnb14 = scnbfactor[ paramIndex - 1 ];

            Expression dihedral_func = k * ( 1 + Cos( periodicity * ( phi - 0 ) - phase ) );

            Atom atom0 = editmol.select( number0 );
            Atom atom1 = editmol.select( number1 );
            Atom atom2 = editmol.select( number2 );
            Atom atom3 = editmol.select( number3 );

            if (improper)
            {
                ImproperID dih = ImproperID( atom0.index(), atom1.index(), atom2.index(), atom3.index() );
                amberparams.add( dih, k, periodicity, phase);
            }
            else
            {
                DihedralID dih = DihedralID( atom0.index(), atom1.index(), atom2.index(), atom3.index() );
                amberparams.add( dih, k, periodicity, phase);
            }

            // Actually, we just save the terms in an array of atom indices
            // because some dihedrals may have multi-term
            // JM Feb 13. Maybe better to create improper/dihedral functions with null terms in case they are changed by a perturbation?
            //if (improper and k > 0.00001)
            if (improper)
            {
                DofID improperid = DofID( atom0.index(), atom1.index(),
                                          atom2.index(), atom3.index() );

                if ( improper_hash.contains(improperid) )
                    improper_hash[improperid] += dihedral_func;
                else
                    improper_hash.insert(improperid, dihedral_func);
            }
            //else if ( k > 0.00001)
            else
            {
                DofID dihid = DofID( atom0.index(), atom1.index(),
                                     atom2.index(), atom3.index() );

                if ( dihedral_hash.contains(dihid) )
                    dihedral_hash[dihid] += dihedral_func;
                else
                    dihedral_hash.insert( dihid, dihedral_func);
            }

            if (not ignored and not improper)
            {
	      //if (sclee14 < 0.00001)
              //  {
	      //	  /*         throw SireError::program_bug( QObject::tr(
	      //" A 1,4 pair has a coulombic scaling factor of 0.0 in the top file which would mean an infinite energy ! "),
	      //CODELOC );*/
	      //	}
              //  if (sclnb14 < 0.00001)
              //  {
	      //	  /*throw SireError::program_bug( QObject::tr(
	      //" A 1,4 pair has a LJ scaling factor of 0.0 in the top file which would mean an infinite energy ! "),
	      //CODELOC );*/
	      //	}

                /**Save this in 14 array */
                if (not atoms14.contains(atom0.number()))
                {
                    QList<AtomNum> list;
                    atoms14.insert( atom0.number(), list);
                }

                /** Not sure this can happens but to be safe.. */
                if ( not atoms14[atom0.number()].contains(atom3.number()) )
                {
                    atoms14[atom0.number()].append(atom3.number());
                    /* JM 07/14 Save scale factor for this pair*/
                    //atoms14sclee[atom0.number()][atom3.number()] = 1/sclee14;
                    //atoms14sclnb[atom0.number()][atom3.number()] = 1/sclnb14;
		    if (sclee14 < 0.00001)
		      atoms14sclee[atom0.number()][atom3.number()] = 0.0;
		    else
		      atoms14sclee[atom0.number()][atom3.number()] = 1/sclee14;
		    if (sclnb14 < 0.00001)
		      atoms14sclnb[atom0.number()][atom3.number()] = 0.0;
		    else
		      atoms14sclnb[atom0.number()][atom3.number()] = 1/sclnb14;
                    // Add pair (atom0,atom3) = (1/sclee14, 1/sclnb14) to amber parameters object
                    BondID pair = BondID(atom0.index() , atom3.index() );
                    amberparams.add14Pair( pair, 1/sclee14, 1/sclnb14 );
                }

                if (not atoms14.contains(atom3.number()))
                {
                    QList<AtomNum> list;
                    atoms14.insert( atom3.number(), list);
                }

                if ( not atoms14[atom3.number()].contains(atom0.number()) )
                {
                    atoms14[atom3.number()].append(atom0.number());
		    if (sclee14 < 0.00001)
		      atoms14sclee[atom3.number()][atom0.number()] = 0.0;
		    else
		      atoms14sclee[atom3.number()][atom0.number()] = 1/sclee14;

		    if (sclee14 < 0.00001)
		      atoms14sclnb[atom3.number()][atom0.number()] = 0.0;
		    else
		      atoms14sclnb[atom3.number()][atom0.number()] = 1/sclnb14;
		    /* JM 07/14 Save scale factor for this pair*/
		    //atoms14sclee[atom3.number()][atom0.number()] = 1/sclee14;
		    //atoms14sclnb[atom3.number()][atom0.number()] = 1/sclnb14;
                }
            }
        }
    }

    if (got_to_end)
        last_idx = pointer;

    /** We can now create the appropriate dihedrals*/
    for (QHash<DofID,Expression>::const_iterator it = dihedral_hash.constBegin();
            it != dihedral_hash.constEnd();
            ++it)
    {
        dihedralfuncs.set( it.key().atom0(), it.key().atom1(), it.key().atom2(),
                           it.key().atom3(), it.value() );
    }

    for (QHash<DofID,Expression>::const_iterator it = improper_hash.constBegin();
            it != improper_hash.constEnd();
            ++it)
    {
        improperfuncs.set( it.key().atom0(), it.key().atom1(), it.key().atom2(),
                           it.key().atom3(), it.value() );
    }

    editmol.setProperty( dihedral_property.source(), dihedralfuncs );
    editmol.setProperty( improper_property.source(), improperfuncs );

    editmol.setProperty( amberparameters_property.source(), amberparams);

    return last_idx;
}

static void setNonBondedPairs(MolEditor &editmol, int pointer,
                              const QList<int> &num_excluded_atoms,
                              const QList<int> &exc_atom_list,
                              CLJNBPairs &nbpairs,
                              const PropertyName &nb_property,
                              const QHash<AtomNum, QList<AtomNum> > &atoms14,
			      const QHash<AtomNum, QHash<AtomNum, double> > &atoms14sclee,
			      const QHash<AtomNum, QHash<AtomNum, double> > &atoms14sclnb)
//                              double coul_14scl, double lj_14scl)
{
    // For each pair of atoms within a molecule
    // --> if 1,2 or 1,3 CLJScaleFactor(0,0)
    // --> if 1,4 CLJScaleFactor( 1/ 1.2 , 1 / 2. )
    // the TOP file stores this information in num_excluded_atoms,
    // exc_atom_list, pointers[NEXT]

    // however it does not discriminates 1,4 from 1,2 or 1,3. So the
    // solution is to get setDihedrals to return a list of 1,4 atoms.

    // --> if 1,5 or more CLJScaleFactor ( 1, 1 )

    // if the number of atoms is less than or equal to 3, then everything is bonded
    if (editmol.nAtoms() <= 3)
    {
        nbpairs = CLJNBPairs(editmol.data().info(), CLJScaleFactor(0,0));
        editmol.setProperty( nb_property, nbpairs );
        return;
    }

    // this is the default situation
    nbpairs = CLJNBPairs(editmol.data().info(), CLJScaleFactor(1.0,1.0));
    int natoms = editmol.nAtoms();

    // find the excluded atoms of i
    // IEXCL = SUM(NUMEX(j), j=1,i-1)
    // then excluded atoms are NATEX(IEXCL) to NATEX(IEXCL+NUMEX(i))
    // if excluded and in 1,4  CLJScaleFactor( 1/ 1.2 , 1 / 2. )
    // if excluded and not in 1,4 CLJScaleFactor(0,0)
    // if not excluded CLJScaleFactor ( 1, 1 )
    for (int i = 0; i < natoms ; i ++ )
    {
        Atom atom0 = editmol.atom(AtomIdx(i));

        // skip itself
        nbpairs.set( atom0.cgAtomIdx(), atom0.cgAtomIdx(),
                     CLJScaleFactor(0.0, 0.0) );

        // Excluded atoms of atom0?
        int iexcl = 0;
        int atomNum = atom0.number();

        // TODO OPT 1: save previous values so no need to recompute from scratch
        for ( int j = 0 ; j < ( atomNum - 1 ); j++ )
        {
            iexcl += num_excluded_atoms[j];
        }

        QList<Atom> excludedAtoms;
        //qDebug() << " Looking at ATOM " << atomNum << atom0.toString();
        //qDebug() << " iexcl is " << iexcl << " num_excluded_atoms[ i ] " <<  num_excluded_atoms[ atomNum - 1 ];

        for ( int j = iexcl ; j < iexcl + num_excluded_atoms[ atomNum - 1 ] ; j++ )
        {
            int jnumber = exc_atom_list[ j ];
            //qDebug() << " jnumber " << jnumber;
            if (jnumber > 0)
            {
                Atom excludedAtom = editmol.atom( AtomNum( jnumber  ) );
                //qDebug() << " EXCLUDED ATOM OF " << atomNum << " : " << jnumber;
                excludedAtoms.append(excludedAtom);
            }
        }

        int nexcluded = excludedAtoms.size();

        if ( nexcluded == 0 )
            continue;

        for ( int j = 0 ; j < nexcluded ; j++ )
        {
            double cscl;
            double ljscl;

            Atom atom1 = excludedAtoms.at(j);

            if ( atoms14[ atom0.number() ].contains( atom1.number() ) )
            {
                //qDebug() << " ATOMS " << atom0.number() << " and " << atom1.number() << " are 14";

		//cscl = coul_14scl;
		//ljscl = lj_14scl;

		if ( not atoms14sclee.contains( atom0.number() )  )
		  {
		    throw SireError::program_bug( QObject::tr(
                "It should not happen that atoms14sclee does not have an entry for atom0"),
                                             CODELOC );
		  }

		if ( not atoms14sclee[atom0.number()].contains( atom1.number() ) )
		  {
		    throw SireError::program_bug( QObject::tr(
                "It should not happen that atoms14sclee does not have an entry for [atom0]]atom1]"),
                                             CODELOC );
		  }

		cscl = atoms14sclee[ atom0.number() ][ atom1.number() ];

		if ( not atoms14sclnb.contains( atom0.number() ) )
		  {
		    throw SireError::program_bug( QObject::tr(
                "It should not happen that atoms14sclnb does not have an entry for atom0"),
                                             CODELOC );
		  }

		if ( not atoms14sclnb[atom0.number()].contains( atom1.number() ) )
		  {
		    throw SireError::program_bug( QObject::tr(
                "It should not happen that atoms14sclnb does not have an entry for [atom0]]atom1]"),
                                             CODELOC );
		  }

		ljscl = atoms14sclnb[ atom0.number() ][ atom1.number() ];

		//qDebug() << " cscl " << cscl << " ljscl " << ljscl << "\n";
            }
            else
            {
                //qDebug() << " ATOMS " << atom0.number() << " and "
                //         << atom1.number() << " are 12 or 13";
                cscl = 0.0;
                ljscl = 0.0;
            }

            CGAtomIdx atom0cgidx = atom0.cgAtomIdx();
            CGAtomIdx atom1cgidx = atom1.cgAtomIdx();

            nbpairs.set( atom0.cgAtomIdx(), atom1.cgAtomIdx(),
                         CLJScaleFactor(cscl, ljscl) );

            // setting atom0/atom1 automatically sets atom1/atom0
            //nbpairs.set( atom1.cgAtomIdx(), atom0.cgAtomIdx(),
            //             CLJScaleFactor(cscl, ljscl) );
        }
    }

    editmol.setProperty( nb_property, nbpairs );
}

static void walk(int atom,
                 const QMultiHash<int, int> &atoms12,
                 int totalMolecules,
                 QHash<int, int> &atIsInMol,
                 QMultiHash<int, int> &atomsInMolecule)
{
    QList<int> neighbors = atoms12.values( atom );
    //qDebug() << " neighbors of " << atom << " are : " << neighbors;

    foreach ( int neighbor, neighbors )
    {
        QList<int> nn = atomsInMolecule.values( totalMolecules );

        if ( not nn.contains(neighbor) )
        {
            atIsInMol[ neighbor ] = totalMolecules;
            atomsInMolecule.insert( totalMolecules, neighbor );
            walk( neighbor, atoms12, totalMolecules, atIsInMol, atomsInMolecule);
        }
    }
}

static void calcNumberMolecules(int &totalMolecules,
                                QList<int> &atoms_per_mol,
                                const QList<int> &bond_inc_h,
                                const QList<int> &bonds_exc_h,
                                int natoms, int nbondsh, int nbondsa)
{
    QHash<int, int> atIsInMol;
    QMultiHash<int, int> atoms12 ;

    /** First construct a list of atoms bonded to each atom*/
    for ( int i = 1 ; i <= natoms ; i++ )
    {
        atIsInMol[ i ] = -1;
    }

    for ( int j = 0 ; j < 3 * nbondsa ; j = j + 3 )
    {
        int atom0 = bonds_exc_h[ j ] / 3 + 1 ;
        int atom1 = bonds_exc_h[ j + 1 ] / 3 + 1 ;
        atoms12.insert(atom0, atom1);
        atoms12.insert(atom1, atom0);
    }

    for ( int j = 0 ; j < 3 * nbondsh ; j = j + 3 )
    {
        int atom0 = bond_inc_h[ j ] / 3 + 1 ;
        int atom1 = bond_inc_h[ j + 1 ] / 3 + 1 ;
        atoms12.insert(atom0, atom1);
        atoms12.insert(atom1, atom0);
    }

    // Then recursively walk along each atom find all the atoms that
    // are in the same molecule
    totalMolecules = 0;
    QMultiHash<int, int> atomsInMolecule;

    for ( int i = 1 ; i <= natoms ; i++ )
    {
        if (atIsInMol[ i ] == -1)
        {
            totalMolecules++;
            atIsInMol[ i ] = totalMolecules;
            atomsInMolecule.insert( totalMolecules, i );

            // Recursive walk
            //qDebug() << " Calling walk ";
            walk(i, atoms12, totalMolecules, atIsInMol, atomsInMolecule);

            //qDebug() << " Molecule " << totalMolecules << " has "
            //         << atomsInMolecule.values(totalMolecules).size() << " atoms ";
            atoms_per_mol.append( atomsInMolecule.values(totalMolecules).size() );
        }
    }
}

/** Internal function to create a set of atom numbers in an editmol*/
/*static QSet<AtomNum> selectAtomsbyNumber(const MolEditor &editmol)
{
    QSet<AtomNum> atomnums;

    Selector<Atom> atoms = editmol.selectAllAtoms();

    for (int i=0; i<atoms.count(); ++i)
    {
        atomnums.insert( atoms[i].number() );
    }

    return atomnums;
}*/

///////////
/////////// Implementation of Amber
///////////

static const RegisterMetaType<Amber> r_amber(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Amber &amber)
{
    //empty class so nothing to stream
    writeHeader(ds, r_amber, 1);

    ds << amber.coul_14scl << amber.lj_14scl;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Amber &amber)
{
    //empty class so nothing to stream

    VersionID v = readHeader(ds, r_amber);

    if (v == 1)
    {
        ds >> amber.coul_14scl >> amber.lj_14scl;
    }
    else if (v == 0)
    {
        amber.coul_14scl = AMBER14COUL;
        amber.lj_14scl = AMBER14LJ;
    }
    else
        throw version_error( v, "1,0", r_amber, CODELOC );

    return ds;
}

/** Constructor */
Amber::Amber() : coul_14scl(AMBER14COUL), lj_14scl(AMBER14LJ)
{}

/** Copy constructor */
Amber::Amber(const Amber &other)
      : coul_14scl(other.coul_14scl), lj_14scl(other.lj_14scl)
{}

/** Destructor */
Amber::~Amber()
{}

/** Copy assignment operator */
Amber& Amber::operator=(const Amber &other)
{
    coul_14scl = other.coul_14scl;
    lj_14scl = other.lj_14scl;
    return *this;
}

/** Comparison operator */
bool Amber::operator==(const Amber &other) const
{
    return coul_14scl == other.coul_14scl and lj_14scl == other.lj_14scl;
}

/** Comparison operator */
bool Amber::operator!=(const Amber &other) const
{
    return not operator==(other);
}

/** Reads the contents of a topfile and associated crdfile and returns a molecule group */
tuple<MoleculeGroup,SpacePtr> Amber::readCrdTop(const QString &crdfile,
        const QString &topfile, QString flag_cutting) const
{
    /**
        See http://ambermd.org/formats.html

        FORMAT(20a4)  (ITITL(i), i=1,20)
    ITITL  : title

    FORMAT(12i6)  NATOM,  NTYPES, NBONH,  MBONA,  NTHETH, MTHETA,
                NPHIH,  MPHIA,  NHPARM, NPARM,  NEXT,   NRES,
                NBONA,  NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA,
                NATYP,  NPHB,   IFPERT, NBPER,  NGPER,  NDPER,
                MBPER,  MGPER,  MDPER,  IFBOX,  NMXRS,  IFCAP
    NATOM  : total number of atoms
    NTYPES : total number of distinct atom types
    NBONH  : number of bonds containing hydrogen
    MBONA  : number of bonds not containing hydrogen
    NTHETH : number of angles containing hydrogen
    MTHETA : number of angles not containing hydrogen
    NPHIH  : number of dihedrals containing hydrogen
    MPHIA  : number of dihedrals not containing hydrogen
    NHPARM : currently not used
    NPARM  : currently not used
    NEXT   : number of excluded atoms
    NRES   : number of residues
    NBONA  : MBONA + number of constraint bonds
    NTHETA : MTHETA + number of constraint angles
    NPHIA  : MPHIA + number of constraint dihedrals
    NUMBND : number of unique bond types
    NUMANG : number of unique angle types
    NPTRA  : number of unique dihedral types
    NATYP  : number of atom types in parameter file, see SOLTY below
    NPHB   : number of distinct 10-12 hydrogen bond pair types
    IFPERT : set to 1 if perturbation info is to be read in
    NBPER  : number of bonds to be perturbed
    NGPER  : number of angles to be perturbed
    NDPER  : number of dihedrals to be perturbed
    MBPER  : number of bonds with atoms completely in perturbed group
    MGPER  : number of angles with atoms completely in perturbed group
    MDPER  : number of dihedrals with atoms completely in perturbed groups
    IFBOX  : set to 1 if standard periodic box, 2 when truncated octahedral
    NMXRS  : number of atoms in the largest residue
    IFCAP  : set to 1 if the CAP option from edit was specified

    FORMAT(20a4)  (IGRAPH(i), i=1,NATOM)
    IGRAPH : the user atoms names

    FORMAT(5E16.8)  (CHRG(i), i=1,NATOM)
    CHRG   : the atom charges.  (Divide by 18.2223 to convert to charge
             in units of the electron charge)

    FORMAT(5E16.8)  (AMASS(i), i=1,NATOM)
    AMASS  : the atom masses

    FORMAT(12I6)  (IAC(i), i=1,NATOM)
    IAC    : index for the atom types involved in Lennard Jones (6-12)
             interactions.  See ICO below.

    FORMAT(12I6)  (NUMEX(i), i=1,NATOM)
    NUMEX  : total number of excluded atoms for atom "i".  See
             NATEX below.

    FORMAT(12I6)  (ICO(i), i=1,NTYPES*NTYPES)
    ICO    : provides the index to the nonbon parameter
             arrays CN1, CN2 and ASOL, BSOL.  All possible 6-12
             or 10-12 atoms type interactions are represented.
             NOTE: A particular atom type can have either a 10-12
             or a 6-12 interaction, but not both.  The index is
             calculated as follows:
               index = ICO(NTYPES*(IAC(i)-1)+IAC(j))
             If index is positive, this is an index into the
             6-12 parameter arrays (CN1 and CN2) otherwise it
             is an index into the 10-12 parameter arrays (ASOL
             and BSOL).

    FORMAT(20A4)  (LABRES(i), i=1,NRES)
    LABRES : the residue labels

    FORMAT(12I6)  (IPRES(i), i=1,NRES)
    IPRES  : atoms in each residue are listed for atom "i" in
             IPRES(i) to IPRES(i+1)-1

    FORMAT(5E16.8)  (RK(i), i=1,NUMBND)
    RK     : force constant for the bonds of each type, kcal/mol

    FORMAT(5E16.8)  (REQ(i), i=1,NUMBND)
    REQ    : the equilibrium bond length for the bonds of each type, angstroms

    FORMAT(5E16.8)  (TK(i), i=1,NUMANG)
    TK     : force constant for the angles of each type, kcal/mol A**2

    FORMAT(5E16.8)  (TEQ(i), i=1,NUMANG)
    TEQ    : the equilibrium angle for the angles of each type, radians

    FORMAT(5E16.8)  (PK(i), i=1,NPTRA)
    PK     : force constant for the dihedrals of each type, kcal/mol

    FORMAT(5E16.8)  (PN(i), i=1,NPTRA)
    PN     : periodicity of the dihedral of a given type

    FORMAT(5E16.8)  (PHASE(i), i=1,NPTRA)
    PHASE  : phase of the dihedral of a given type, radians

    FORMAT(5E16.8)  (SOLTY(i), i=1,NATYP)
    SOLTY  : currently unused (reserved for future use)

    FORMAT(5E16.8)  (CN1(i), i=1,NTYPES*(NTYPES+1)/2)
    LENNARD_JONES_ACOEF    : Lennard Jones r**12 terms for all possible
             atom type interactions, indexed by ICO and IAC; for atom i
             and j where i < j, the index into this array is as follows
             (assuming the value of ICO(index) is positive):
             CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j))).

    FORMAT(5E16.8)  (CN2(i), i=1,NTYPES*(NTYPES+1)/2)
    LENNARD_JONES_BCOEF    : Lennard Jones r**6 terms for all possible
             atom type interactions.  Indexed like CN1 above.

    NOTE: the atom numbers in the following arrays that describe bonds, angles, and dihedrals are
    coordinate array indexes for runtime speed. The true atom number equals the absolute value of
    the number divided by three, plus one. In the case of the dihedrals, if the fourth atom is negative,
    this implies that the dihedral is an improper. If the third atom is negative, this implies that the
    end group interations are to be ignored. End group interactions are ignored, for example, in
    dihedrals of various ring systems (to prevent double counting of 1-4 interactions) and in
    multiterm dihedrals.

    FORMAT(12I6)  (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
    IBH    : atom involved in bond "i", bond contains hydrogen
    JBH    : atom involved in bond "i", bond contains hydrogen
    ICBH   : index into parameter arrays RK and REQ

    FORMAT(12I6)  (IB(i),JB(i),ICB(i), i=1,NBONA)
    IB     : atom involved in bond "i", bond does not contain hydrogen
    JB     : atom involved in bond "i", bond does not contain hydrogen
    ICB    : index into parameter arrays RK and REQ

    FORMAT(12I6)  (ITH(i),JTH(i),KTH(i),ICTH(i), i=1,NTHETH)
    ITH    : atom involved in angle "i", angle contains hydrogen
    JTH    : atom involved in angle "i", angle contains hydrogen
    KTH    : atom involved in angle "i", angle contains hydrogen
    ICTH   : index into parameter arrays TK and TEQ for angle
             ITH(i)-JTH(i)-KTH(i)

    FORMAT(12I6)  (IT(i),JT(i),KT(i),ICT(i), i=1,NTHETA)
    IT     : atom involved in angle "i", angle does not contain hydrogen
    JT     : atom involved in angle "i", angle does not contain hydrogen
    KT     : atom involved in angle "i", angle does not contain hydrogen
    ICT    : index into parameter arrays TK and TEQ for angle
             IT(i)-JT(i)-KT(i)

    FORMAT(12I6)  (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i), i=1,NPHIH)
    IPH    : atom involved in dihedral "i", dihedral contains hydrogen
    JPH    : atom involved in dihedral "i", dihedral contains hydrogen
    KPH    : atom involved in dihedral "i", dihedral contains hydrogen
    LPH    : atom involved in dihedral "i", dihedral contains hydrogen
    ICPH   : index into parameter arrays PK, PN, and PHASE for
             dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i)

    FORMAT(12I6)  (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)
    IP     : atom involved in dihedral "i", dihedral does not contain hydrogen
    JP     : atom involved in dihedral "i", dihedral does not contain hydrogen
    KP     : atom involved in dihedral "i", dihedral does not contain hydrogen
    LP     : atom involved in dihedral "i", dihedral does not contain hydrogen
    ICP    : index into parameter arrays PK, PN, and PHASE for
             dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i).  Note, if the
             periodicity is negative, this implies the following entry
             in the PK, PN, and PHASE arrays is another term in a
             multitermed dihedral.

    FORMAT(12I6)  (NATEX(i), i=1,NEXT)
    NATEX  : the excluded atom list.  To get the excluded list for atom
             "i" you need to traverse the NUMEX list, adding up all
             the previous NUMEX values, since NUMEX(i) holds the number
             of excluded atoms for atom "i", not the index into the
             NATEX list.  Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
             excluded atoms are NATEX(IEXCL) to NATEX(IEXCL+NUMEX(i)).

    FORMAT(5E16.8)  (ASOL(i), i=1,NPHB)
    ASOL   : the value for the r**12 term for hydrogen bonds of all
             possible types.  Index into these arrays is equivalent
             to the CN1 and CN2 arrays, however the index is negative.
             For example, for atoms i and j, with i < j, the index is
             -ICO(NTYPES*(IAC(i)-1+IAC(j)).

    FORMAT(5E16.8)  (BSOL(i), i=1,NPHB)
    BSOL   : the value for the r**10 term for hydrogen bonds of all
             possible types.  Indexed like ASOL.

    FORMAT(5E16.8)  (HBCUT(i), i=1,NPHB)
    HBCUT  : no longer in use

    FORMAT(20A4)  (ISYMBL(i), i=1,NATOM)
    ISYMBL : the AMBER atom types for each atom

    FORMAT(20A4)  (ITREE(i), i=1,NATOM)
    ITREE  : the list of tree joining information, classified into five
             types.  M -- main chain, S -- side chain, B -- branch point,
             3 -- branch into three chains, E -- end of the chain

    FORMAT(12I6)  (JOIN(i), i=1,NATOM)
    JOIN   : tree joining information, potentially used in ancient
             analysis programs.  Currently unused in sander or gibbs.

    FORMAT(12I6)  (IROTAT(i), i = 1, NATOM)
    IROTAT : apparently the last atom that would move if atom i was
             rotated, however the meaning has been lost over time.
             Currently unused in sander or gibbs.

    The following are only present if IFBOX .gt. 0

    FORMAT(12I6)  IPTRES, NSPM, NSPSOL
    IPTRES : final residue that is considered part of the solute,
             reset in sander and gibbs
    NSPM   : total number of molecules
    NSPSOL : the first solvent "molecule"

    FORMAT(12I6)  (NSP(i), i=1,NSPM)
    NSP    : the total number of atoms in each molecule,
             necessary to correctly perform the pressure
             scaling.

    FORMAT(5E16.8)  BETA, BOX(1), BOX(2), BOX(3)
    BETA   : periodic box, angle between the XY and YZ planes in
             degrees.
    BOX    : the periodic box lengths in the X, Y, and Z directions

    The following are only present if IFCAP .gt. 0

    FORMAT(12I6)  NATCAP
    NATCAP : last atom before the start of the cap of waters
             placed by edit

    FORMAT(5E16.8)  CUTCAP, XCAP, YCAP, ZCAP
    CUTCAP : the distance from the center of the cap to the outside
    XCAP   : X coordinate for the center of the cap
    YCAP   : Y coordinate for the center of the cap
    ZCAP   : Z coordinate for the center of the cap

    The following is only present if IFPERT .gt. 0
    Note that the initial state, or equivalently the prep/link/edit state, is represented by lambda=1 and the
    perturbed state, or final state specified in parm, is the lambda=0 state.

    FORMAT(12I6)  (IBPER(i), JBPER(i), i=1,NBPER)
    IBPER  : atoms involved in perturbed bonds
    JBPER  : atoms involved in perturbed bonds

    FORMAT(12I6)  (ICBPER(i), i=1,2*NBPER)
    ICBPER : pointer into the bond parameter arrays RK and REQ for the
             perturbed bonds.  ICBPER(i) represents lambda=1 and
             ICBPER(i+NBPER) represents lambda=0.

    FORMAT(12I6)  (ITPER(i), JTPER(i), KTPER(i), i=1,NGPER)
    IPTER  : atoms involved in perturbed angles
    JTPER  : atoms involved in perturbed angles
    KTPER  : atoms involved in perturbed angles

    FORMAT(12I6)  (ICTPER(i), i=1,2*NGPER)
    ICTPER : pointer into the angle parameter arrays TK and TEQ for
             the perturbed angles.  ICTPER(i) represents lambda=0 and
             ICTPER(i+NGPER) represents lambda=1.

    FORMAT(12I6)  (IPPER(i), JPPER(i), KPPER(i), LPPER(i), i=1,NDPER)
    IPPER  : atoms involved in perturbed dihedrals
    JPPER  : atoms involved in perturbed dihedrals
    KPPER  : atoms involved in perturbed dihedrals
    LPPER  : atoms involved in pertrubed dihedrals

    FORMAT(12I6)  (ICPPER(i), i=1,2*NDPER)
    ICPPER : pointer into the dihedral parameter arrays PK, PN and
             PHASE for the perturbed dihedrals.  ICPPER(i) represents
             lambda=1 and ICPPER(i+NGPER) represents lambda=0.

    FORMAT(20A4)  (LABRES(i), i=1,NRES)
    LABRES : residue names at lambda=0

    FORMAT(20A4)  (IGRPER(i), i=1,NATOM)
    IGRPER : atomic names at lambda=0

    FORMAT(20A4)  (ISMPER(i), i=1,NATOM)
    ISMPER : atomic symbols at lambda=0

    FORMAT(5E16.8)  (ALMPER(i), i=1,NATOM)
    ALMPER : unused currently in gibbs

    FORMAT(12I6)  (IAPER(i), i=1,NATOM)
    IAPER  : IAPER(i) = 1 if the atom is being perturbed

    FORMAT(12I6)  (IACPER(i), i=1,NATOM)
    IACPER : index for the atom types involved in Lennard Jones
             interactions at lambda=0.  Similar to IAC above.
             See ICO above.

    FORMAT(5E16.8)  (CGPER(i), i=1,NATOM)
    CGPER  : atomic charges at lambda=0

    The following is only present if IPOL .eq. 1

    FORMAT(5E18.8) (ATPOL(i), i=1,NATOM)
    ATPOL  : atomic polarizabilities

    The following is only present if IPOL .eq. 1 .and. IFPERT .eq. 1

    FORMAT(5E18.8) (ATPOL1(i), i=1,NATOM)
    ATPOL1 : atomic polarizabilities at lambda = 1 (above is at lambda = 0)
    */

    // Read the contents of the top file
    QFile top_f(topfile);

    if ( not (top_f.exists() and top_f.open(QIODevice::ReadOnly) ) )
    {
        throw SireError::file_error(top_f, CODELOC);
    }

    QTextStream ts(&top_f);

    // TOP file format generated by sleap in Amber-tools 1.4
    // see amber11/AmberTools/src/gleap/mortsrc/ambfmt/prmtop.cpp

    // The following holds the data read from the top file
    int currentFlag = 0;

    FortranFormat currentFormat;

    //QString currentVersion = " ";
    QString title = " ";
    QList<int> pointers;
    QStringList atom_name;
    QList<double> charge;
    QList<double> mass;
    QList<int> atom_type_index;
    QList<int> element;
    QList<int> num_excluded_atoms;
    QList<int> nb_parm_index;
    QStringList res_label;
    QList<int> res_pointer;
    QList<double> bond_force_constant;
    QList<double> bond_equil_value;
    QList<double> ang_force_constant;
    QList<double> ang_equil_value;
    QList<double> dih_force_constant;
    QList<double> dih_periodicity;
    QList<double> dih_phase;
    QList<double> solty;
    QList<double> lj_a_coeff;
    QList<double> lj_b_coeff;
    QList<int> bond_inc_h;
    QList<int> bonds_exc_h;
    QList<int> angs_inc_h;
    QList<int> angs_exc_h;
    QList<int> dihs_inc_h;
    QList<int> dihs_exc_h;
    QList<int> exc_atom_list;
    QList<double> hbond_a_coeff;
    QList<double> hbond_b_coeff;
    QList<double> hbcut;
    QStringList amber_type;
    QStringList tree_chain_class;
    QList<int> join_array;
    QList<int> irotat;
    QList<int> svn_pointers;
    QList<int> atoms_per_mol;
    QList<double> box_dims;
    QStringList radius_set;
    QList<double> radii;
    QList<double> screen;
    QList<double> sceefactor;
    QList<double> scnbfactor;

    QString line = ts.readLine();

    // The first line containts the version.
    // Make sure we know how to read this version
    //processVersionLine(words, currentVersion);
    while (not line.isNull() )
    {
        line = ts.readLine();

        if (line.startsWith("%"))
        {
            // We are reading meta data, can be VERSION, FLAG, FORMAT
            QStringList words = line.split(" ", Qt::SkipEmptyParts);

            if (line.startsWith("%FLAG"))
                processFlagLine(words, currentFlag);
            else if (line.startsWith("%FORMAT"))
                processFormatLine(words, currentFormat);
            else if (line.startsWith("%COMMENT"))
                continue;
            else
            {
                qDebug() << "ERROR" << line;
                throw SireError::program_bug( QObject::tr(
                                                  "Does not know what to do with a '%1' statement")
                                              .arg(words[0]), CODELOC );
            }
        }
        else
        {
            // if not we are reading data. Read it according to the
            // latest FLAG/FORMAT values

            //qDebug() << "FORMAT IS " << currentFormat.type << " FLAG IS " << currentFlag;
            switch ( currentFlag )
            {
                case TITLE:
                    processTitleLine(line, currentFormat, title);
                    break;
                case POINTERS:
                    processIntegerLine(line, currentFormat, pointers);
                    break;
                case ATOM_NAME:
                    processStringLine(line, currentFormat, atom_name);
                    break;
                case CHARGE:
                    processDoubleLine(line, currentFormat, charge);
                    break;
                case MASS:
                    processDoubleLine(line, currentFormat, mass);
                    break;
                case ATOM_TYPE_INDEX:
                    processIntegerLine(line, currentFormat, atom_type_index);
                    break;
                case NUMBER_EXCLUDED_ATOMS:
                    processIntegerLine(line, currentFormat, num_excluded_atoms);
                    break;
                case NONBONDED_PARM_INDEX:
                    processIntegerLine(line, currentFormat, nb_parm_index);
                    break;
                case RESIDUE_LABEL:
                    processStringLine(line, currentFormat, res_label);
                    break;
                case RESIDUE_POINTER:
                    processIntegerLine(line, currentFormat, res_pointer);
                    break;
                case BOND_FORCE_CONSTANT:
                    processDoubleLine(line, currentFormat, bond_force_constant);
                    break;
                case BOND_EQUIL_VALUE:
                    processDoubleLine(line, currentFormat, bond_equil_value);
                    break;
                case ANGLE_FORCE_CONSTANT:
                    processDoubleLine(line, currentFormat, ang_force_constant);
                    break;
                case ANGLE_EQUIL_VALUE:
                    processDoubleLine(line, currentFormat, ang_equil_value);
                    break;
                case DIHEDRAL_FORCE_CONSTANT:
                    processDoubleLine(line, currentFormat, dih_force_constant);
                    break;
                case DIHEDRAL_PERIODICITY:
                    processDoubleLine(line, currentFormat, dih_periodicity);
                    break;
                case DIHEDRAL_PHASE:
                    processDoubleLine(line, currentFormat, dih_phase);
                    break;
                case SOLTY:
                    processDoubleLine(line, currentFormat, solty);
                    break;
                case LENNARD_JONES_ACOEF:
                    processDoubleLine(line, currentFormat, lj_a_coeff);
                    break;
                case LENNARD_JONES_BCOEF:
                    processDoubleLine(line, currentFormat, lj_b_coeff);
                    break;
                case BONDS_INC_HYDROGEN:
                    processIntegerLine(line, currentFormat, bond_inc_h);
                    break;
                case BONDS_WITHOUT_HYDROGEN:
                    processIntegerLine(line, currentFormat, bonds_exc_h);
                    break;
                case ANGLES_INC_HYDROGEN:
                    processIntegerLine(line, currentFormat, angs_inc_h);
                    break;
                case ANGLES_WITHOUT_HYDROGEN:
                    processIntegerLine(line, currentFormat, angs_exc_h);
                    break;
                case DIHEDRALS_INC_HYDROGEN:
                    processIntegerLine(line, currentFormat, dihs_inc_h);
                    break;
                case DIHEDRALS_WITHOUT_HYDROGEN:
                    processIntegerLine(line, currentFormat, dihs_exc_h);
                    break;
                case EXCLUDED_ATOMS_LIST:
                    processIntegerLine(line, currentFormat, exc_atom_list);
                    break;
                case HBOND_ACOEF:
                    processDoubleLine(line, currentFormat, hbond_a_coeff);
                    break;
                case HBOND_BCOEF:
                    processDoubleLine(line, currentFormat, hbond_b_coeff);
                    break;
                case HBCUT:
                    processDoubleLine(line, currentFormat, hbcut);
                    break;
                case AMBER_ATOM_TYPE:
                    processStringLine(line, currentFormat, amber_type);
                    break;
                case TREE_CHAIN_CLASSIFICATION:
                    processStringLine(line, currentFormat, tree_chain_class);
                    break;
                case JOIN_ARRAY:
                    processIntegerLine(line, currentFormat, join_array);
                    break;
                case IROTAT:
                    processIntegerLine(line, currentFormat, irotat);
                    break;
                case SOLVENT_POINTERS:
                    processIntegerLine(line, currentFormat, svn_pointers);
                    break;
                case ATOMS_PER_MOLECULE:
                    processIntegerLine(line, currentFormat, atoms_per_mol);
                    break;
                case BOX_DIMENSIONS:
                    processDoubleLine(line, currentFormat, box_dims);
                    break;
                case RADIUS_SET:
                    processStringLine(line, currentFormat, radius_set);
                    break;
                case RADII:
                    processDoubleLine(line, currentFormat, radii);
                    break;
                case SCREEN:
                    processDoubleLine(line, currentFormat, screen);
                    break;
                case ATOMIC_NUMBER:
                    processIntegerLine(line, currentFormat, element);
                    break;
	        case SCEE_SCALE_FACTOR:
		    processDoubleLine(line, currentFormat, sceefactor);
		    break;
	        case SCNB_SCALE_FACTOR:
		    processDoubleLine(line, currentFormat, scnbfactor);
		    break;
                case UNKNOWN:
                    break;
                default:
                {
                    qDebug() << "PROCESSING UNKNOWN LINE";
                    qDebug() << currentFlag;
                    qDebug() << line;
                    throw SireError::program_bug( QObject::tr(
                                                      "Serious problem with the value of the variable "
                                                      "currentFlag, '%1'").arg(currentFlag),
                                                  CODELOC );
                    break;
                }
            }
        }
    }

    int cutting;

    if(flag_cutting == "perresidue")
        cutting = PERRESIDUE;
    else if (flag_cutting == "peratom")
        cutting = PERATOM;
    else
        throw SireError::program_bug(QObject::tr(
    "The Cutting method has not been correctly specified. Possible choises: perresidue, peratom"),
                                     CODELOC);

    //DEBUG STUFF
    //qDebug() << " POINTERS " << pointers << " ATOMNAME "
    //         << atom_name.size() << " CHARGE " << charge.size()
    //         << " MASS " << mass.size() << " ATOM_TYPE_INDEX "
    //         << atom_type_index.size() << " SCREEN " << screen.size();

    //qDebug() << " Finished reading the top file";

    // Now read the contents of the crd file to get the coordinates
    QFile crd_f(crdfile);

    if ( not (crd_f.exists() and crd_f.open(QIODevice::ReadOnly) ) )
    {
        throw SireError::file_error(crd_f, CODELOC);
    }

    QTextStream ts2(&crd_f);

    // the first line contains the title
    line = ts2.readLine();

    // the second line contains the number of atoms, read as 1I6
    line = ts2.readLine();

    //QList<int> temp;
    //FortranFormat crd_int_format = FortranFormat(1, "I", 6, 0);
    //processIntegerLine(line, crd_int_format,temp);

    // Split by space to get the first number. I unfortunately discovered that
    // leap and sander do not produce exactly the same crd files...
    QStringList temp = line.split(" ", Qt::SkipEmptyParts);
    bool ok;
    int crd_atoms = temp[0].toInt(&ok);

    //qDebug() << " THERE ARE " << crdAtoms << " ATOMS IN THE CRD FILE";
    // Check that this number of atoms is compatible with what is in the top file
    if (pointers[NATOM] != crd_atoms)
        throw SireError::incompatible_error( QObject::tr(
                "The number of atoms in the crd file (%1) does not equal the number "
                "of atoms in the top file (%2)!")
                                             .arg(crd_atoms).arg(pointers[NATOM]), CODELOC );

    QList<double> crd_coords;
    FortranFormat crd_double_format = FortranFormat(6,"E",12,7);

    // Must read crdAtoms / 2 lines, but make sure to round up !
    int crd_coord_lines = ( double (crd_atoms) / 2.0 ) + 0.5 ;

    for (int i=0; i < crd_coord_lines ; i++ )
    {
        line = ts2.readLine();
        processDoubleLine(line, crd_double_format, crd_coords);
    }

    //qDebug() << " THE COORDS ARE " << crd_coords;

    // And now the box dimensions. Skip to the last line of the file, because
    // the crd file could have contained velocities, which are not used for the moment

    while (not ts2.atEnd())
    {
        line = ts2.readLine();
    }

    //qDebug() << "THE LAST LINE IS " << line;

    QList<double> crd_box;

    if ( pointers[IFBOX] != 0 )
    {
        processDoubleLine(line, crd_double_format, crd_box);
        //qDebug() << " THE BOX DIMENSIONS ARE " << crd_box;
    }

    //qDebug() << " Finished reading the crd file";

    // Now create the atoms and molecules etc..

    PropertyName coords_property = PropertyName("coordinates");
    PropertyName charge_property = PropertyName("charge");
    PropertyName element_property = PropertyName("element");
    PropertyName mass_property = PropertyName("mass");
    PropertyName lj_property = PropertyName("LJ");
    PropertyName ambertype_property = PropertyName("ambertype");

    PropertyName connectivity_property = PropertyName("connectivity");
    PropertyName bond_property = PropertyName("bond");
    PropertyName angle_property = PropertyName("angle");
    PropertyName dihedral_property = PropertyName("dihedral");
    PropertyName improper_property = PropertyName("improper");
    PropertyName nb_property = PropertyName("intrascale");

    PropertyName amberparameters_property = PropertyName("amberparameters");

    MoleculeGroup molecules( QString("%1:%2").arg(crdfile,topfile) );

    int molnum = 1;
    int resnum = 1;

    int total_molecules;

    if ( pointers[IFBOX] != 0 )
    {
        // When loading a top file setup with a periodic box,
        // the number of molecules and atoms per molecule
        // has been specified, which makes our life easier
        total_molecules = svn_pointers[NSPM];
    }
    else
    {
        // Otherwise we need to figure out the number of molecules
        // using the information about bonds
        //qDebug() << "Information about the position of atoms in molecules is "
        //            "not present in the top file and must be determined. "
        //            "This may take a while...";
        calcNumberMolecules(total_molecules, atoms_per_mol,
                            bond_inc_h, bonds_exc_h,
                            pointers[NATOM], pointers[NBONH], pointers[MBONA]);
    }

    int idx_inc_h = 0;
    int idx_exc_h = 0;
    int idx_bnd_inc_h = 0;
    int idx_bnd_exc_h = 0;
    int idx_ang_inc_h = 0;
    int idx_ang_exc_h = 0;
    int idx_dih_inc_h = 0;
    int idx_dih_exc_h = 0;

    for (int i=0; i < total_molecules ; i++)
    {
        //qDebug() << " Parameterizing molecule " << i;
        /** First pass, use StructureEditors to build the layout of the molecule*/
        MolStructureEditor molstructeditor;
        int atoms_in_mol = 0;

        while (atoms_in_mol < atoms_per_mol[i])
        {
            int start_atom = res_pointer[resnum - 1];

            // Be careful not to overflow
            int end_atom;
            //qDebug() << " atomStart " << atomStart << " pointers[NRES] "
            //         << pointers[NRES] ;

            if ( resnum < ( pointers[NRES] ) )
                end_atom = res_pointer[resnum - 1 + 1 ] - 1;
            else
                end_atom = pointers[NATOM] ;

            //qDebug() << "Residue " << resnum << " start "
            //         << start_atom << " end " << end_atom;

            // create an empty residue. Use RESIDUE_LABEL for the name
            ResStructureEditor resstructeditor = molstructeditor.add( ResNum(resnum) );
            resstructeditor.rename( ResName( res_label[resnum - 1]) );

            for (int j=start_atom; j <= end_atom; ++j)
            {
                AtomStructureEditor atomstructeditor = molstructeditor.add( AtomNum(j) );
                atomstructeditor.rename( AtomName(atom_name[j -1]) );
                atomstructeditor.reparent( ResNum(resnum) );
            }

            atoms_in_mol += ( end_atom - start_atom ) + 1 ;
            ++resnum;
        }

        // Create cut groups using a per residue or per atom scheme

        if(cutting == PERRESIDUE)
        {

            ResidueCutting residue_cutfunc = ResidueCutting();

            molstructeditor = residue_cutfunc(molstructeditor);
        }

        else if(cutting == PERATOM)
        {

            AtomCutting atom_cutfunc = AtomCutting();

            molstructeditor = atom_cutfunc(molstructeditor);

        }

        Molecule molecule = molstructeditor.commit();

        MolEditor editmol = molecule.edit();

        ConnectivityEditor connectivity = Connectivity(editmol.data()).edit();

        TwoAtomFunctions bondfuncs(editmol);
        ThreeAtomFunctions anglefuncs(editmol);
        FourAtomFunctions dihedralfuncs(editmol);
        FourAtomFunctions improperfuncs(editmol);

        AmberParameters amberparams(editmol);

        CLJNBPairs nbpairs;
        QHash<AtomNum, QList<AtomNum> > atoms14;
        // JM 07/14 Store info about variable scale factors
        QHash<AtomNum, QHash<AtomNum, double> > atoms14sclee;
        QHash<AtomNum, QHash<AtomNum, double> > atoms14sclnb;

        int natoms = editmol.nAtoms();

        for (int i=0; i < natoms ; ++i)
        {
            //qDebug() << " Parameterizing atom " << i;
            // Now that the structure of the molecule has been built, we assign the
            // following atom properties: coordinates, charge, mass, lj , amber_atom_type
            // and element (if element is available)
            AtomEditor editatom = editmol.atom(AtomIdx(i));

            setAtomParameters( editatom, editmol, crd_coords, coords_property,
                               element, element_property,
                               charge, charge_property,
                               mass, mass_property, atom_type_index,
                               nb_parm_index, lj_a_coeff,
                               lj_b_coeff, lj_property,
                               amber_type, ambertype_property, pointers);
        }

        // Find all the bonds that involve an atom in this molecule.
        // Because bonds are indexed in two arrays in the top file
        // (those without hydrogen atoms and those with hydrogen atoms)
        // we have to look in both.

        // The following routines are slow because for each atom in a molecule
        // we look at all arrays to find entries this is inefficient for solvent
        // molecules. One way to improve efficiency could be to construct all the molecules
        // first and maintain an array of atom pointers. Then in a single pass
        // create connectivity, bonds, angles, dihedrals in molecules of matching atoms
        if (natoms > 1)
        {
            //qDebug() << " Setting up connectivity ";
            idx_inc_h = setConnectivity(editmol, pointers[NBONH],
                                        bond_inc_h,
                                        connectivity, connectivity_property, idx_inc_h);

            idx_exc_h = setConnectivity(editmol, pointers[MBONA],
                                        bonds_exc_h,
                                        connectivity, connectivity_property, idx_exc_h);

            // Next all the forcefield terms
            //qDebug() << " Setting up bonds ";
            idx_bnd_inc_h = setBonds(editmol, pointers[NBONH],
                                     bond_inc_h,
                                     bond_force_constant, bond_equil_value,
                                     bondfuncs, bond_property,
                                     amberparams, amberparameters_property,
                                     idx_bnd_inc_h);

            idx_bnd_exc_h = setBonds(editmol, pointers[MBONA],
                                     bonds_exc_h,
                                     bond_force_constant, bond_equil_value,
                                     bondfuncs, bond_property,
                                     amberparams, amberparameters_property,
                                     idx_bnd_exc_h);
        }

        if (natoms > 2)
        {
            //qDebug() << " Setting up angles ";
            idx_ang_inc_h = setAngles(editmol, pointers[NTHETH],
                                      angs_inc_h,
                                      ang_force_constant, ang_equil_value,
                                      anglefuncs, angle_property,
                                      amberparams, amberparameters_property,
                                      idx_ang_inc_h);

            idx_ang_exc_h = setAngles(editmol, pointers[MTHETA],
                                      angs_exc_h,
                                      ang_force_constant, ang_equil_value,
                                      anglefuncs, angle_property,
                                      amberparams, amberparameters_property,
                                      idx_ang_exc_h);
        }

        if (natoms >3)
        {
            //qDebug() << " Setting up dihedrals ";
            idx_dih_inc_h = setDihedrals(editmol, pointers[NPHIH],
                                         dihs_inc_h,
                                         dih_force_constant, dih_periodicity, dih_phase,
                                         sceefactor, scnbfactor,
                                         dihedralfuncs, dihedral_property,
                                         improperfuncs, improper_property,
                                         atoms14, atoms14sclee, atoms14sclnb,
                                         coul_14scl, lj_14scl,
                                         amberparams, amberparameters_property,
                                         idx_dih_inc_h);

            idx_dih_exc_h = setDihedrals(editmol, pointers[MPHIA],
                                         dihs_exc_h,
                                         dih_force_constant, dih_periodicity, dih_phase,
                                         sceefactor, scnbfactor,
                                         dihedralfuncs, dihedral_property,
                                         improperfuncs, improper_property,
                                         atoms14, atoms14sclee, atoms14sclnb,
                                         coul_14scl, lj_14scl,
                                         amberparams, amberparameters_property,
                                         idx_dih_exc_h);
        }

        setNonBondedPairs(editmol, pointers[NEXT],
                          num_excluded_atoms, exc_atom_list,
                          nbpairs, nb_property,
                          atoms14, atoms14sclee, atoms14sclnb);

        molecule = editmol.commit();

        molecules.add(molecule);
        ++molnum;
    }

    // Now the box information
    SpacePtr spce;

    if ( pointers[IFBOX] == 1)
    {
        /** Rectangular box, dimensions read from the crd file */
        Vector dimensions( crd_box[0], crd_box[1], crd_box[2] );

        //qDebug() << "We have a periodic box of dimensions"
        //         << crdBox[0] << crdBox[1] << crdBox[2];

        // PeriodicBox.
        if ( crd_box[3] == 90.0 && crd_box[4] == 90.0 && crd_box[5] == 90.0 )
        {
            spce = PeriodicBox( dimensions );
        }
        // TriclinicBox.
        else
        {
            spce = TriclinicBox(dimensions.x(), dimensions.y(), dimensions.z(),
                                crd_box[3]*degrees, crd_box[4]*degrees, crd_box[5]*degrees);

        }
        //spce = PeriodicBox( Vector ( crdBox[0], crdBox[1], crdBox[2] ) ).asA<Space>() ;
        //qDebug() << " periodic box " << spce.toString() ;
    }
    else if ( pointers[IFBOX] == 2 )
    {
        /** Truncated Octahedral box*/
        throw SireError::incompatible_error( QObject::tr(
                "Sire does not yet support a truncated octahedral box"),
                                             CODELOC );
    }
    else
    {
        /** Default is a non periodic system */
        spce = Cartesian();
    }

    return tuple<MoleculeGroup, SpacePtr>(molecules, spce);
}

/** Write the coordinates of the molecules in the passed MoleculeGroup in the
    passed Space to an Amber7,
    format coordinate/restart file. The passed property map is used to find
    the required properties */
void Amber::writeCrd(const MoleculeGroup &mols, const Space &space, const QString &filename,
                     const PropertyMap &map) const
{
    QFile file(filename);

    if (not file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        throw SireError::file_error(file, CODELOC);
    }

    QTextStream fout(&file);

    fout << "Sire generated restart\n";

    //count up the number of atoms in the molecules
    int total_nats = 0;

    for (auto mol : mols)
    {
        total_nats += mol.molecule().nAtoms();
    }

    //write out the number of atoms
    fout << total_nats <<  "  0.0000000E+00\n";

    //now write out the coordinates as six columns of F12.7
    int nmols = mols.nMolecules();

    const PropertyName coords_property = map["coordinates"];

    int nwritten = 0;
    int nwritten_atoms = total_nats;

    fout.setFieldWidth(12);
    fout.setRealNumberPrecision(7);
    fout.setRealNumberNotation(QTextStream::FixedNotation);

    for (int i=0; i<nmols; ++i)
    {
        Molecule mol = mols[MolIdx(i)].molecule();

        int nats = mol.nAtoms();

        for (int j=0; j<nats; ++j)
        {
            Atom atom = mol.atom(AtomIdx(j));
            Vector coords = atom.property<Vector>(coords_property);

            fout << coords.x() << coords.y() << coords.z();

            nwritten += 3;
            nwritten_atoms -= 1;

            if (nwritten >= 6)
            {
                fout.setFieldWidth(0);
                fout << "\n";
                fout.setFieldWidth(12);
                nwritten = 0;
            }
        }
    }

    if (nwritten < 6)
    {
        fout.setFieldWidth(0);
        fout << "\n";
        fout.setFieldWidth(12);
    }

    //now write out the space
    if (space.isA<PeriodicBox>())
    {
        Vector boxsize = space.asA<PeriodicBox>().dimensions();

        fout << boxsize.x() << boxsize.y() << boxsize.z()
             << 90.0 << 90.0 << 90.0;

        fout.setFieldWidth(0);
        fout << "\n";
        fout.setFieldWidth(12);
    }
    else if (space.isA<TriclinicBox>())
    {
        Vector v0 = space.asA<TriclinicBox>().vector0();
        Vector v1 = space.asA<TriclinicBox>().vector1();
        Vector v2 = space.asA<TriclinicBox>().vector2();

        auto x = v0.magnitude();
        auto y = v1.magnitude();
        auto z = v2.magnitude();

        double alpha = space.asA<TriclinicBox>().alpha();
        double beta  = space.asA<TriclinicBox>().beta();
        double gamma = space.asA<TriclinicBox>().gamma();

        fout << x << y << z << alpha << beta << gamma;

        fout.setFieldWidth(0);
        fout << "\n";
        fout.setFieldWidth(12);
    }

    //ok, all done
    file.close();

    if (nwritten_atoms > 0)
    {
        qDebug() << "WARNING: Number of written atoms is" << nwritten_atoms
                 << "less than the number in the Molecules (" << total_nats << ")";
    }
    else if (total_nats < 0)
    {
        qDebug() << "WARNING: Number of written atoms is" << (-nwritten_atoms)
                 << "more than the number in the Molecules (" << total_nats << ")";
    }
}

const char* Amber::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Amber>() );
}

void Amber::set14Factors(double coul_14, double lj_14)
{
    coul_14scl = coul_14;
    lj_14scl = lj_14;
}

double Amber::coulomb14Factor() const
{
    return coul_14scl;
}

double Amber::lj14Factor() const
{
    return lj_14scl;
}

const char* Amber::what() const
{
    return Amber::typeName();
}


