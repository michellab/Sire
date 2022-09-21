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

#include "zmatrixmaker.h"

#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/moleditor.h"
#include "SireMol/reseditor.h"
#include "SireMol/selector.hpp"
#include "SireMol/connectivity.h"
#include "SireMol/moleculedata.h"
#include "SireMol/bondid.h"

#include "SireMove/zmatrix.h"

#include "SireUnits/convert.h"
#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"
#include "SireMol/errors.h"

using namespace SireIO;
using namespace SireMol;
using namespace SireStream;
using namespace SireMove;
using namespace SireUnits;

//
// Implementation of ZmatrixLine
//

ZmatrixLineTemplate::ZmatrixLineTemplate(const QString &atom, const QString &bond, const QString &angle,
                                         const QString &dihedral, const double &bondDelta, const double &angleDelta,
                                         const double &dihedralDelta)
    : atom(atom), bond(bond), angle(angle), dihedral(dihedral),
      bondDelta(bondDelta), angleDelta(angleDelta), dihedralDelta(dihedralDelta)
{
}

bool ZmatrixLineTemplate::has(const QString &atom)
{
    if (this->atom == atom)
        return true;
    else
        return false;
}

const QString ZmatrixLineTemplate::getAtom()
{
    return this->atom;
}

const QString ZmatrixLineTemplate::getBond()
{
    return this->bond;
}

const QString ZmatrixLineTemplate::getAngle()
{
    return this->angle;
}

const QString ZmatrixLineTemplate::getDihedral()
{
    return this->dihedral;
}

const double ZmatrixLineTemplate::getBondDelta()
{
    return this->bondDelta;
}

const double ZmatrixLineTemplate::getAngleDelta()
{
    return this->angleDelta;
}

const double ZmatrixLineTemplate::getDihedralDelta()
{
    return this->dihedralDelta;
}

void ZmatrixLineTemplate::setBondDelta(double bondDelta)
{
    this->bondDelta = bondDelta;
}

void ZmatrixLineTemplate::setAngleDelta(double angleDelta)
{
    this->angleDelta = angleDelta;
}

void ZmatrixLineTemplate::setDihedralDelta(double dihedralDelta)
{
    this->dihedralDelta = dihedralDelta;
}

QString ZmatrixLineTemplate::toString()
{
    return QObject::tr("ZmatrixLineTemplate: %1 %2 %3 %4 %5 %6 %7")
        .arg(this->getAtom())
        .arg(this->getBond())
        .arg(this->getAngle())
        .arg(this->getDihedral())
        .arg(this->getBondDelta())
        .arg(this->getAngleDelta())
        .arg(this->getDihedralDelta());
}

//
// Implementation of ZmatrixTemplate
//

ZmatrixTemplate::ZmatrixTemplate(const QString &name)
    : name(name)
{
}

QList<ZmatrixLineTemplate> ZmatrixTemplate::getZmatrix()
{
    QList<ZmatrixLineTemplate> zmatrix;
    foreach (ZmatrixLineTemplate zmatline, this->zmatrix)
        zmatrix.append(zmatline);
    return zmatrix;
}

ZmatrixLineTemplate ZmatrixTemplate::getZmatrixLineTemplate(const QString &atom)
{
    if (not this->zmatrix.contains(atom))
        throw SireError::invalid_key(atom, CODELOC);
    else
        return this->zmatrix.value(atom);
}

const QString ZmatrixTemplate::getName()
{
    return this->name;
}

void ZmatrixTemplate::addZmatrixLineTemplate(const QString &atom, const ZmatrixLineTemplate &zmatline)
{
    if (this->zmatrix.contains(atom))
        throw SireError::invalid_key(atom, CODELOC);
    else
        this->zmatrix[atom] = zmatline;
}

void ZmatrixTemplate::setBondDelta(const QString &atom, const QString &bond, double bondDelta)
{
    if (not this->zmatrix.contains(atom))
        throw SireError::invalid_key(atom, CODELOC);
    // Also check that the bond atom is correct
    else if (bond != zmatrix[atom].getBond())
        throw SireError::invalid_key(bond, CODELOC);
    else
        this->zmatrix[atom].setBondDelta(bondDelta);
}

void ZmatrixTemplate::setAngleDelta(const QString &atom, const QString &bond, const QString &angle, double angleDelta)
{
    if (not this->zmatrix.contains(atom))
        throw SireError::invalid_key(atom, CODELOC);
    else if (bond != zmatrix[atom].getBond())
        throw SireError::invalid_key(bond, CODELOC);
    else if (angle != zmatrix[atom].getAngle())
        throw SireError::invalid_key(angle, CODELOC);
    else
        this->zmatrix[atom].setAngleDelta(angleDelta);
}

void ZmatrixTemplate::setDihedralDelta(const QString &atom, const QString &bond,
                                       const QString &angle, const QString &dihedral, double dihedralDelta)
{
    if (not this->zmatrix.contains(atom))
        throw SireError::invalid_key(atom, CODELOC);
    else if (bond != zmatrix[atom].getBond())
        throw SireError::invalid_key(bond, CODELOC);
    else if (angle != zmatrix[atom].getAngle())
        throw SireError::invalid_key(angle, CODELOC);
    else if (dihedral != zmatrix[atom].getDihedral())
        throw SireError::invalid_key(dihedral, CODELOC);
    else
        this->zmatrix[atom].setDihedralDelta(dihedralDelta);
}

//
// Implementation of ZmatrixResidue
//

ZmatrixResidue::ZmatrixResidue(const QString &name)
    : ZmatrixTemplate(name)
{
}

void ZmatrixResidue::setRotation(const double rot)
{
    this->rotate = rot;
}

void ZmatrixResidue::setTranslation(const double trans)
{
    this->translate = trans;
}

const double ZmatrixResidue::getRotation()
{
    return this->rotate;
}

const double ZmatrixResidue::getTranslation()
{
    return this->translate;
}

void ZmatrixResidue::addChain(const QString &name, const ZmatrixTemplate &chain)
{
    if (this->backbone.contains(name))
        throw SireError::invalid_key(name, CODELOC);
    else
        this->backbone[name] = chain;
}

const ZmatrixTemplate ZmatrixResidue::getChain(const QString &name)
{
    if (not this->backbone.contains(name))
        throw SireError::invalid_key(name, CODELOC);
    else
        return this->backbone[name];
}

QList<ZmatrixTemplate> ZmatrixResidue::getChains()
{
    QList<ZmatrixTemplate> chains;
    foreach (ZmatrixTemplate chain, this->backbone)
        chains.append(chain);
    return chains;
}

QStringList ZmatrixResidue::getBBatoms()
{
    return this->bbatoms;
}

void ZmatrixResidue::addBBatom(const QString &bbAtom)
{
    this->bbatoms.append(bbAtom);
}

//
// Helper functions to parse a templates input file
//

static int processVersionLine(QString &line)
{
    QStringList words = line.split(" ", Qt::SkipEmptyParts);
    bool ok;
    int version = words[1].toInt(&ok);
    if (not ok)
        throw SireError::program_bug(QObject::tr(
                                         "Unexpected error while trying to read the version line of the zmatrix template file"),
                                     CODELOC);
    return version;
}

//
// Implementation of ZmatrixMaker
//
static const RegisterMetaType<ZmatrixMaker> r_zmatrixmaker(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ZmatrixMaker &zmatrixmaker)
{
    writeHeader(ds, r_zmatrixmaker, 1);

    SharedDataStream sds(ds);

    //    sds << zmatrixmaker.residues;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ZmatrixMaker &zmatrixmaker)
{
    VersionID v = readHeader(ds, r_zmatrixmaker);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        //        sds >> zmatrixmaker.residues;
    }
    else
        throw version_error(v, "1", r_zmatrixmaker, CODELOC);

    return ds;
}

/** Default constructor */
ZmatrixMaker::ZmatrixMaker()
{
}

//** Destructor */
ZmatrixMaker::~ZmatrixMaker()
{
}

const char *ZmatrixMaker::typeName()
{
    return QMetaType::typeName(qMetaTypeId<ZmatrixMaker>());
}

/** Load the zmatrix templates specified in &templatefile*/
void ZmatrixMaker::loadTemplates(const QString &templatefile)
{
    QFile template_f(templatefile);

    if (not(template_f.exists() and template_f.open(QIODevice::ReadOnly)))
    {
        throw SireError::file_error(template_f, CODELOC);
    }

    QTextStream ts(&template_f);

    QString line = ts.readLine();
    /** The first line contains the version*/
    int version = ::processVersionLine(line);
    if (version != 1)
        throw SireError::process_error(QObject::tr(
                                           "Invalid version of the template, got '%1' but only support '1' ")
                                           .arg(version),
                                       CODELOC);

    /** Holds the dictionnary of chain and residues read from the templates*/
    QHash<QString, ZmatrixTemplate> chains;
    QHash<QString, ZmatrixResidue> residues;

    enum
    {
        CHAIN = 0,
        RESIDUE = 1
    };
    QString current = " "; // the chain/residue currently read
    int currentType = CHAIN;

    /** Now read rest of the file */
    while (not line.isNull())
    {
        line = ts.readLine();
        QStringList words = line.split(" ", Qt::SkipEmptyParts);
        // qDebug() << line;

        if (line.startsWith("chain"))
        {
            ZmatrixTemplate chain = ZmatrixTemplate(words[1]);
            current = chain.getName();
            chains[current] = chain;
            currentType = CHAIN;
        }
        else if (line.startsWith("residue"))
        {
            ZmatrixResidue residue = ZmatrixResidue(words[1]);
            current = residue.getName();
            residues[current] = residue;
            currentType = RESIDUE;
        }
        else if (line.startsWith("bbatom"))
        {
            if (currentType != RESIDUE)
                throw SireError::io_error(QObject::tr(
                                              "There is a problem with the input file %1, "
                                              " a bbatom line should only be defined for a residue template")
                                              .arg(templatefile),
                                          CODELOC);
            residues[current].addBBatom(words[1]);
        }
        else if (line.startsWith("zmatrix"))
        {
            ZmatrixLineTemplate zmatline = ZmatrixLineTemplate(words[1], words[2], words[3], words[4]);
            if (currentType == CHAIN)
                chains[current].addZmatrixLineTemplate(words[1], zmatline);
            else if (currentType == RESIDUE)
                residues[current].addZmatrixLineTemplate(words[1], zmatline);
        }
        else if (line.startsWith("bond"))
        {
            if (currentType == CHAIN)
                chains[current].setBondDelta(words[1], words[2], words[4].toDouble());
            else if (currentType == RESIDUE)
                residues[current].setBondDelta(words[1], words[2], words[4].toDouble());
        }
        else if (line.startsWith("angle"))
        {
            if (currentType == CHAIN)
                chains[current].setAngleDelta(words[1], words[2], words[3], words[5].toDouble());
            else if (currentType == RESIDUE)
                residues[current].setAngleDelta(words[1], words[2], words[3], words[5].toDouble());
        }
        else if (line.startsWith("dihedral"))
        {
            if (currentType == CHAIN)
                chains[current].setDihedralDelta(words[1], words[2], words[3], words[4], words[6].toDouble());
            else if (currentType == RESIDUE)
                residues[current].setDihedralDelta(words[1], words[2], words[3], words[4], words[6].toDouble());
        }
        else if (line.startsWith("rigidbody"))
        {
            if (currentType != RESIDUE)
                throw SireError::io_error(QObject::tr(
                                              " There is a problem with the input file %1, "
                                              " a rigidbody line should only be defined for a residue template")
                                              .arg(templatefile),
                                          CODELOC);
            residues[current].setRotation(words[2].toDouble());
            residues[current].setTranslation(words[4].toDouble());
        }
        else if (line.startsWith("backbone"))
        {
            if (currentType != RESIDUE)
                throw SireError::io_error(QObject::tr(
                                              " There is a problem with the input file %1, "
                                              " a backbone line should only be defined for a residue template")
                                              .arg(templatefile),
                                          CODELOC);
            /** The appropriate chain can be defined for a variety of environments */
            residues[current].addChain(words[1], chains[words[2]]);
            residues[current].addChain(words[3], chains[words[4]]);
            residues[current].addChain(words[5], chains[words[6]]);
            residues[current].addChain(words[7], chains[words[8]]);
        }
    }

    /** Store the new residues */
    // foreach (ZmatrixTemplate chain, chains)
    //   {
    //     qDebug() << " CHAIN " << chain.getName();
    //     QList<ZmatrixLineTemplate> zmatrix = chain.getZmatrix();
    //     foreach ( ZmatrixLineTemplate zmatline, zmatrix)
    // 	qDebug() << zmatline.toString();
    //   }
    foreach (ZmatrixResidue residue, residues)
    {
        // qDebug() << " RESIDUE " << residue.getName();
        // qDebug() << " Rotate " << residue.getRotation();
        // qDebug() << " Translate " << residue.getTranslation();
        // qDebug() << " CHAINS " ;
        // QList<ZmatrixTemplate> chains = residue.getChains();
        // foreach (ZmatrixTemplate chain, chains)
        //	qDebug() << chain.getName();
        // QList<ZmatrixLineTemplate> zmatrix = residue.getZmatrix();
        // foreach ( ZmatrixLineTemplate zmatline, zmatrix)
        //	qDebug() << zmatline.toString();
        QString resname = residue.getName();
        if (not this->residues.contains(resname))
            this->residues[resname] = residue;
    }
}

Molecule ZmatrixMaker::applyTemplates(Molecule &molecule)
{
    PropertyName zmatrix_property = PropertyName("z-matrix");

    // Does it already have a zmatrix_property?
    if (molecule.hasProperty(zmatrix_property))
        return molecule;

    const Connectivity &connectivity = molecule.data()
                                           .property("connectivity")
                                           .asA<Connectivity>();

    MolEditor editmol = molecule.edit();

    ZMatrix zmatrix(editmol);

    int nres = editmol.nResidues();

    for (ResIdx i(0); i < nres; ++i)
    {
        Residue residue = editmol.residue(i);

        /** Look up a residue among the templates*/
        QString resname = residue.name().value();

        if (not this->residues.contains(resname))
            throw SireError::invalid_key(resname, CODELOC);

        ZmatrixResidue restemplate = this->residues[resname];

        /** Is this residue 'first', 'middle' or 'last' or 'single' ?*/
        bool firstbondedwithother = false;
        bool lastbondedwithother = false;

        // nterm has first false and last true
        // cterm has first true and last false
        // middle has first and last true
        // single hast first and last false
        QStringList bbatoms = restemplate.getBBatoms();

        Atom lastatom = residue.select(AtomName(bbatoms.last()));

        QList<BondID> bonds = connectivity.getBonds(lastatom.index());

        foreach (BondID bond, bonds)
        {
            int n_in_residue = 0;

            if (residue.contains(bond.atom0()))
                n_in_residue += 1;

            if (residue.contains(bond.atom1()))
                n_in_residue += 1;

            if (n_in_residue == 1)
            {
                lastbondedwithother = true;
                break;
            }
        }

        Atom firstatom = residue.select(AtomName(bbatoms.first()));
        bonds = connectivity.getBonds(firstatom.index());

        foreach (BondID bond, bonds)
        {
            int n_in_residue = 0;

            if (residue.contains(bond.atom0()))
                n_in_residue += 1;

            if (residue.contains(bond.atom1()))
                n_in_residue += 1;

            if (n_in_residue == 1)
            {
                firstbondedwithother = true;
                break;
            }
        }

        QString position;

        if (firstbondedwithother and lastbondedwithother)
            position = "middle";
        else if (firstbondedwithother)
            position = "last";
        else if (lastbondedwithother)
            position = "first";
        else
            position = "single";

        //qDebug() << " The position of this residue is " << position ;

        /** Need to decide where to store info about rigid body
            translation and rotations of this residue*/

        Selector<Atom> atoms = residue.atoms();

        for (int i = 0; i < atoms.count(); ++i)
        {
            Atom atom = atoms(i);

            /** Backbone atoms do not have a zmatrix line */
            if (bbatoms.contains(atom.name().value()))
                continue;

            /** Try to find matching atom in the residue zmatrix */
            ZmatrixLineTemplate linetemplate;

            try
            {
                linetemplate = restemplate.getZmatrixLineTemplate(atom.name().value());
            }
            catch (SireError::invalid_key)
            {
                /** If this fails, look also in the matching backbone zmatrix*/
                ZmatrixTemplate chain = restemplate.getChain(position);
                linetemplate = chain.getZmatrixLineTemplate(atom.name().value());
            }

            Atom bond = residue.select(AtomName(linetemplate.getBond()));
            Atom angle = residue.select(AtomName(linetemplate.getAngle()));
            Atom dihedral = residue.select(AtomName(linetemplate.getDihedral()));

            double bondDelta = linetemplate.getBondDelta();
            double angleDelta = linetemplate.getAngleDelta();
            double dihedralDelta = linetemplate.getDihedralDelta();

            ZMatrixLine zmatrixline = ZMatrixLine(atom.index(), bond.index(),
                                                  angle.index(), dihedral.index());

            zmatrixline.setBondDelta(bondDelta * angstrom);
            zmatrixline.setAngleDelta(angleDelta * degrees);
            zmatrixline.setDihedralDelta(dihedralDelta * degrees);

            zmatrix.add(zmatrixline);
        }
    }

    editmol.setProperty(zmatrix_property, zmatrix);

    return editmol.commit();
}
