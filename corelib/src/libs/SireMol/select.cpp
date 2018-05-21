/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Christopher Woods
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

#include "SireMol/select.h"
#include "SireMol/parser.h"

#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;
using namespace SireMol;

///////////
/////////// Implementation of SelectEngine
///////////

/** Constructor */
SireMol::parser::SelectEngine::SelectEngine()
{}

/** Destructor */
SireMol::parser::SelectEngine::~SelectEngine()
{}

SelectResult SireMol::parser::SelectEngine::operator()(const MolGroupsBase &molgroups,
                                                       const PropertyMap &map) const
{
    return this->select(molgroups, map);
}

SelectResult SireMol::parser::SelectEngine::operator()(const MoleculeGroup &molgroup,
                                                       const PropertyMap &map) const
{
    return this->select(molgroup, map);
}

SelectResult SireMol::parser::SelectEngine::operator()(const Molecules &molecules,
                                                       const PropertyMap &map) const
{
    return this->select(molecules, map);
}

SelectResult SireMol::parser::SelectEngine::operator()(const MoleculeView &molecule,
                                                       const PropertyMap &map) const
{
    return this->select(molecule, map);
}

SelectResult SireMol::parser::SelectEngine::select(const MoleculeView &molecule,
                                                   const PropertyMap &map) const
{
    ViewsOfMol views = this->selectFromMolecule(molecule, map);
    
    if (views.isEmpty())
        return SelectResult();
    else
    {
        QList<ViewsOfMol> result;
        result.append(views);
        return SelectResult(views);
    }
}

SelectResult SireMol::parser::SelectEngine::select(const Molecules &molecules,
                                                   const PropertyMap &map) const
{
    //loop through in molecule number order so that we are consistent
    //every time we run this search
    const auto molnums = molecules.molNums();
    qSort(molnums);
    
    QList<ViewsOfMol> result;
    
    for (const auto molnum : molnums)
    {
        auto views = this->selectFromMolecule( molecules[molnum], map );
        
        if (not views.isEmpty())
            result.append(views);
    }
    
    return SelectResult(views);
}

SelectResult SireMol::parser::SelectEngine::select(const MoleculeGroup &molgroup,
                                                   const PropertyMap &map) const
{
    QList<ViewsOfMol> result;
    
    for (int i=0; i<molgroup.nMolecules(); ++i)
    {
        auto views = this->selectFromMolecule( molgroup.moleculeAt(i), map );
        
        if (not views.isEmpty())
            result.append(views);
    }
    
    return SelectResult(views);
}

SelectResult SireMol::parser::SelectEngine::select(const MolGroupsBase &molgroups,
                                                   const PropertyMap &map) const
{
    return this->select(molgroups.molecules(), map);
}

///////////
/////////// Implementation of Select
///////////

static const RegisterMetaType<Select> r_select;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const Select &select)
{
    writeHeader(ds, r_select, 1);
    
    SharedDataStream sds(ds);
    
    sds << select.search_string;
    
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, Select &select)
{
    VersionID v = readHeader(ds, r_select);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        QString search_string;
        sds >> search_string;

        select = Select(search_string);
    }
    else
        throw version_error(v, "1", r_select, CODELOC);
    
    return ds;
}

/** Construct an empty selection (will select nothing) */
Select::Select() : ConcreteProperty<Select,Property>()
{}

/** Construct a selection based on the passed string */
Select::Select(const QString &str) : ConcreteProperty<Select,Property>()
{
    e = SireMol::parser::parse(str);
    search_string = str;
}

/** Copy constructor */
Select::Select(const Select &other)
       : ConcreteProperty<Select,Property>(other),
         search_string(other.search_string), e(other.e)
{}

/** Destructor */
Select::~Select()
{}

Select& Select::operator=(const Select &other)
{
    search_string = other.search_string;
    e = other.e;
    Property::operator=(other);
    return *this;
}

bool Select::operator==(const Select &other) const
{
    return search_string == other.search_string and Property::operator==(other);
}

bool Select::operator!=(const Select &other) const
{
    return not Select::operator==(other);
}

Select* Select::clone() const
{
    return new Select(*this);
}

const char* Select::what() const
{
    return Select::typeName();
}

const char* Select::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Select>() );
}

SelectResult Select::operator()(const MolGroupsBase &molgroups,
                                const PropertyMap &map) const
{
    if (e.get() != 0)
        return e->operator()(molgroups,map);
    else
        return SelectResult();
}

SelectResult Select::operator()(const MoleculeGroup &molgroup,
                                const PropertyMap &map) const
{
    if (e.get() != 0)
        return e->operator()(molgroups,map);
    else
        return SelectResult();
}

SelectResult Select::operator()(const Molecules &molecules,
                                const PropertyMap &map) const
{
    if (e.get() != 0)
        return e->operator()(molgroups,map);
    else
        return SelectResult();
}

SelectResult Select::operator()(const MoleculeView &molecule,
                                const PropertyMap &map) const
{
    if (e.get() != 0)
        return e->operator()(molgroups,map);
    else
        return SelectResult();
}

void Select::setToken(const QString &token, const QString &selection)
{
    SireMol::parser::set_token(token, selection);
}

void Select::resetTokens()
{
    SireMol::parser::reset_tokens();
}

QString Select::toString() const
{
    if (search_string.isEmpty())
    {
        return QObject::tr("Select::null");
    }
    else
        return QObject::tr("Select( %1 )").arg(search_string);
}

///////////
/////////// Implementation of SelectResult
///////////

static const RegisterMetaType<SelectResult> r_result;

QDataStream...

