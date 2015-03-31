/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "flexibilitylibrary.h"

#include "SireMol/atom.h"
#include "SireMol/atomeditor.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/selector.hpp"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"

#include "SireUnits/convert.h"
#include "SireUnits/units.h"
#include "SireUnits/dimensions.h"

#include "SireMove/internalmove.h"
#include "SireMove/flexibility.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"
#include "SireMol/errors.h"

using namespace SireIO;
using namespace SireMol;
using namespace SireMove;
using namespace SireStream;
using namespace SireUnits;

//
// Implementation of FlexibilityTemplate
//

static RegisterMetaType<FlexibilityTemplate> r_flextemplate(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds,
                                      const FlexibilityTemplate &flextemplate)
{
    writeHeader(ds, r_flextemplate, 1);
    
    SharedDataStream sds(ds);
    
    sds << flextemplate.name << flextemplate.translation 
        << flextemplate.rotation << flextemplate.maxbondvar
	<< flextemplate.maxanglevar << flextemplate.maxdihedralvar
        << flextemplate.bonds << flextemplate.angles
        << flextemplate.dihedrals;
        
    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, FlexibilityTemplate &flextemplate)
{
    VersionID v = readHeader(ds, r_flextemplate);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> flextemplate.name >> flextemplate.translation
            >> flextemplate.rotation >> flextemplate.maxbondvar
	    >> flextemplate.maxanglevar >> flextemplate.maxdihedralvar
            >> flextemplate.bonds >> flextemplate.angles
            >> flextemplate.dihedrals;
    }
    else
        throw version_error(v, "1", r_flextemplate, CODELOC);
        
    return ds;
}

/** Constructor */
FlexibilityTemplate::FlexibilityTemplate()
                    : translation(0), rotation(0), maxbondvar(0), maxanglevar(0), maxdihedralvar(0)
{}

FlexibilityTemplate::FlexibilityTemplate(const QString &name)
  : name(name), translation(0), rotation(0), maxbondvar(0), maxanglevar(0), maxdihedralvar(0)
{}

/** Copy constructor */
FlexibilityTemplate::FlexibilityTemplate(const FlexibilityTemplate &other)
                    : name(other.name), translation(other.translation),
                      rotation(other.rotation), maxbondvar(other.maxbondvar),
		      maxanglevar(other.maxanglevar), maxdihedralvar(other.maxdihedralvar),
                      bonds(other.bonds), angles(other.angles), 
                      dihedrals(other.dihedrals)
{}

/** Destructor */
FlexibilityTemplate::~FlexibilityTemplate()
{}

const char* FlexibilityTemplate::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FlexibilityTemplate>() );
}

/** Copy assignment operator */
FlexibilityTemplate& FlexibilityTemplate::operator=(const FlexibilityTemplate &other)
{
    if (this != &other)
    {
        name = other.name;
        translation = other.translation;
        rotation = other.rotation;
        maxbondvar = other.maxbondvar;
	maxanglevar = other.maxanglevar;
	maxdihedralvar = other.maxdihedralvar;
        bonds = other.bonds;
        angles = other.angles;
        dihedrals = other.dihedrals;
    }
    
    return *this;
}

/** Comparison operator */
bool FlexibilityTemplate::operator==(const FlexibilityTemplate &other) const
{
    return this == &other or
           (name == other.name and translation == other.translation and
            rotation == other.rotation and maxbondvar == other.maxbondvar and
	    maxanglevar == other.maxanglevar and maxdihedralvar == other.maxdihedralvar and
            bonds == other.bonds and angles == other.angles and
            dihedrals == other.dihedrals);
}

/** Comparison operator */
bool FlexibilityTemplate::operator!=(const FlexibilityTemplate &other) const
{
    return not FlexibilityTemplate::operator==(other);
}

const QString FlexibilityTemplate::getName()
{
    return this->name;
}

void FlexibilityTemplate::setRotation(const Angle &rotation)
{
    this->rotation = rotation;
}

void FlexibilityTemplate::setTranslation(const Length &translation)
{
    this->translation = translation;
}

void FlexibilityTemplate::setMaximumBondVar(int maxvar)
{
    this->maxbondvar = maxvar;
}

void FlexibilityTemplate::setMaximumAngleVar(int maxvar)
{
    this->maxanglevar = maxvar;
}

void FlexibilityTemplate::setMaximumDihedralVar(int maxvar)
{
    this->maxdihedralvar = maxvar;
}


Angle FlexibilityTemplate::getRotation() const
{
    return this->rotation;
}

Length FlexibilityTemplate::getTranslation() const
{
    return this->translation;
}

int FlexibilityTemplate::getMaximumBondVar() const
{
    return this->maxbondvar;
}

int FlexibilityTemplate::getMaximumAngleVar() const
{
    return this->maxanglevar;
}

int FlexibilityTemplate::getMaximumDihedralVar() const
{
    return this->maxdihedralvar;
}


void FlexibilityTemplate::setBondDelta(const BondID &bond, const Length &delta)
{
    bonds.insert(bond, delta);
}

void FlexibilityTemplate::setAngleDelta(const AngleID &angle, const Angle &delta)
{
    angles.insert(angle, delta);
}

void FlexibilityTemplate::setDihedralDelta(const DihedralID &dihedral, const Angle &delta)
{
    dihedrals.insert(dihedral,delta);
}

Length FlexibilityTemplate::getDelta(const BondID &bondid) const
{
    return bonds.value(bondid, Length(0));
}

Angle FlexibilityTemplate::getDelta(const AngleID &angleid) const
{
    return angles.value(angleid, Angle(0));
}

Angle FlexibilityTemplate::getDelta(const DihedralID &dihedralid) const
{
    return dihedrals.value(dihedralid, Angle(0));
}

const QHash<BondID,Length>& FlexibilityTemplate::getBondDeltas() const
{
    return this->bonds;
}

const QHash<AngleID,Angle>& FlexibilityTemplate::getAngleDeltas() const
{
    return this->angles;
}

const QHash<DihedralID,Angle>& FlexibilityTemplate::getDihedralDeltas() const
{
    return this->dihedrals;
}

//
// Helper functions to parse a templates input file
//

static int processVersionLine( QString& line)
{
    QStringList words = line.split(" ", QString::SkipEmptyParts);
    bool ok;
    int version = words[1].toInt(&ok);
    
    if (not ok)
        throw SireError::program_bug( QObject::tr(
                    "Unexpected error while trying to read the version line "
                    "of the zmatrix template file"), CODELOC);

    return version;
}

//
// Implementation of FlexibilityLibrary
//
static const RegisterMetaType<FlexibilityLibrary> r_flexibilitymaker;

/** Serialise to a binary datastream */
QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, 
                                      const FlexibilityLibrary &flexibilitymaker)
{
    writeHeader(ds, r_flexibilitymaker, 1);
    
    SharedDataStream sds(ds);

    sds << flexibilitymaker.templates;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, 
                                      FlexibilityLibrary &flexibilitymaker)
{
    VersionID v = readHeader(ds, r_flexibilitymaker);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> flexibilitymaker.templates;
    }
    else
        throw version_error( v, "1", r_flexibilitymaker, CODELOC );
        
    return ds;
}

/** Default constructor */
FlexibilityLibrary::FlexibilityLibrary()
{}

/** Construct, loading the library from the passed file */
FlexibilityLibrary::FlexibilityLibrary(const QString &file)
{
    this->loadTemplates(file);
}

/** Copy constructor */
FlexibilityLibrary::FlexibilityLibrary(const FlexibilityLibrary &other)
                 : templates(other.templates)
{}

//** Destructor */
FlexibilityLibrary::~FlexibilityLibrary()
{}

/** Copy assignment operator */
FlexibilityLibrary& FlexibilityLibrary::operator=(const FlexibilityLibrary &other)
{
    templates = other.templates;
    return *this;
}

/** Comparison operator */
bool FlexibilityLibrary::operator==(const FlexibilityLibrary &other) const
{
    return templates == other.templates;
}

/** Comparison operator */
bool FlexibilityLibrary::operator!=(const FlexibilityLibrary &other) const
{
    return not FlexibilityLibrary::operator==(other);
}

const char* FlexibilityLibrary::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FlexibilityLibrary>() );
}

/** Add the templates in 'other' into this library */
FlexibilityLibrary& FlexibilityLibrary::operator+=(const FlexibilityLibrary &other)
{
    if (templates.isEmpty())
    {
        templates = other.templates;
    }
    else
    {
        for (QHash<QString,FlexibilityTemplate>::const_iterator 
                                            it = other.templates.constBegin();
             it != other.templates.constEnd();
             ++it)
        {
            templates.insert( it.key(), it.value() );
        }
    }
    
    return *this;
}

/** Return the combination of this library with 'other' */
FlexibilityLibrary FlexibilityLibrary::operator+(const FlexibilityLibrary &other) const
{
    FlexibilityLibrary ret(*this);
    ret += other;
    return ret;
}

/** Add the templates in 'other' into this library */
void FlexibilityLibrary::add(const FlexibilityLibrary &other)
{
    this->operator+=(other);
}

/** Return the template for the specified 'key'

    \throw SireError::invalid_key
*/
const FlexibilityTemplate& FlexibilityLibrary::getTemplate(const QString &key)
{
    QHash<QString,FlexibilityTemplate>::const_iterator it = templates.constFind(key);
    
    if (it == templates.constEnd())
        throw SireError::invalid_key( QObject::tr(
                "Cannot find the template with key \"%1\". Available templates "
                "are %2.")
                    .arg(key, Sire::toString(templates.keys())), CODELOC );

    return it.value();
}

/** Set the template associated with the passed key */
void FlexibilityLibrary::setTemplate(const QString &key, 
                                     const FlexibilityTemplate &tmplate)
{
    templates.insert(key, tmplate);
}

void FlexibilityLibrary::loadTemplates(const QString &templatefile)
{
    QFile template_f(templatefile);
  
    if ( not (template_f.exists() and template_f.open(QIODevice::ReadOnly) ) )
    {
        throw SireError::file_error(template_f, CODELOC);
    }
  
    QTextStream ts(&template_f);
  
    QString line = ts.readLine();

    // The first line contains the version
    int version = ::processVersionLine(line);

    if (version != 1)
        throw SireError::process_error( QObject::tr(
                    "Invalid version of the template, got '%1' but only support '1'")
                        .arg(version), CODELOC);
  
    QString current = " "; // the template currently being read

    QHash <QString,FlexibilityTemplate> new_templates; 

    /** Now read rest of the file */
    while ( not line.isNull() )
    {
        line = ts.readLine();
        QStringList words = line.split(" ", QString::SkipEmptyParts);
        //     qDebug() << line;

        if ( line.startsWith("molecule") )
        {
            // create a new flexibilitytemplate
            FlexibilityTemplate flextemplate = FlexibilityTemplate( words[1] );
            current = flextemplate.getName();
            new_templates[current] = flextemplate;
        }
        else if ( line.startsWith("rigidbody") )
        {
            new_templates[current].setRotation( words[2].toDouble() * degrees );
            new_templates[current].setTranslation( words[4].toDouble() * angstroms );
        }
        else if ( line.startsWith("maximumbondvariables") )
        {
            new_templates[current].setMaximumBondVar( words[1].toInt() );
        }
        else if ( line.startsWith("maximumanglevariables") )
        {
            new_templates[current].setMaximumAngleVar( words[1].toInt() );
        }
        else if ( line.startsWith("maximumdihedralvariables") )
        {
            new_templates[current].setMaximumDihedralVar( words[1].toInt() );
        }	
        else if ( line.startsWith("bond") )
        {
            new_templates[current].setBondDelta( BondID(AtomName(words[1]),
                                                        AtomName(words[2])),
                                                 words[4].toDouble() * angstroms );
        }
        else if ( line.startsWith("angle") )
        {
            new_templates[current].setAngleDelta( AngleID(AtomName(words[1]),
                                                          AtomName(words[2]),
                                                          AtomName(words[3])),
                                                  words[5].toDouble() * degrees );
        }
        else if ( line.startsWith("dihedral") )
        {
            new_templates[current].setDihedralDelta( DihedralID(AtomName(words[1]),
                                                                AtomName(words[2]),
                                                                AtomName(words[3]),
                                                                AtomName(words[4])),
                                                     words[6].toDouble() * degrees );
        }
    }

    foreach (FlexibilityTemplate templ, new_templates)
    {
        QString templname = templ.getName();

        if ( not templates.contains(templname) )
            templates.insert(templname, templ);
    }
}

/** Generate the Flexibility property for the atoms in the passed molecule view */
Flexibility FlexibilityLibrary::getFlexibility(const MoleculeView &molecule) const
{
    AtomSelection selected_atoms = molecule.selection();

    Flexibility flexibility = Flexibility(molecule);
  
    MolName moleculename = molecule.data().name();

    if ( not this->templates.contains(moleculename) )
        throw SireError::invalid_key(QObject::tr("There is no flexibility template for the "
               "molecule with name \"%1\" - available templates are %2.")
                    .arg(moleculename, Sire::toString(templates.keys())), CODELOC);

    FlexibilityTemplate templ = this->templates[moleculename];

    flexibility.setRotation(templ.getRotation());
    flexibility.setTranslation(templ.getTranslation());
    flexibility.setMaximumBondVar(templ.getMaximumBondVar());
    flexibility.setMaximumAngleVar(templ.getMaximumAngleVar());
    flexibility.setMaximumDihedralVar(templ.getMaximumDihedralVar());    

    for (QHash<BondID,Length>::const_iterator it = templ.getBondDeltas().constBegin();
         it != templ.getBondDeltas().constEnd();
         ++it)
    {
        const BondID &bond = it.key();
    
        if ( selected_atoms.selectedAll() or
             (selected_atoms.selected(bond.atom0()) and
              selected_atoms.selected(bond.atom1())) )
        {
            flexibility.add(bond, it.value());
        }
    }

    for (QHash<AngleID,Angle>::const_iterator it = templ.getAngleDeltas().constBegin();
         it != templ.getAngleDeltas().constEnd();
         ++it)
    {
        const AngleID &angle = it.key();
        
        if ( selected_atoms.selectedAll() or
             (selected_atoms.selected(angle.atom0()) and
              selected_atoms.selected(angle.atom1()) and
              selected_atoms.selected(angle.atom2())) )
        {
            flexibility.add(angle, it.value());
        }
    }

    for (QHash<DihedralID,Angle>::const_iterator 
                                    it = templ.getDihedralDeltas().constBegin();
         it != templ.getDihedralDeltas().constEnd();
         ++it)
    {
        const DihedralID &dihedral = it.key();
        
        if ( selected_atoms.selectedAll() or
             (selected_atoms.selected(dihedral.atom0()) and
              selected_atoms.selected(dihedral.atom1()) and
              selected_atoms.selected(dihedral.atom2()) and
              selected_atoms.selected(dihedral.atom3())) )
        {
            flexibility.add(dihedral, it.value());
        }
    }

    return flexibility;
}
