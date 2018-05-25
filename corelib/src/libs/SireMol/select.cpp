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

#include "SireStream/datastream.h"
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

/** Return the pointer to self */
SireMol::parser::SelectEnginePtr SireMol::parser::SelectEngine::self()
{
    return selfptr.lock();
}

SelectResult SireMol::parser::SelectEngine::operator()(const SelectResult &result,
                                                       const PropertyMap &map) const
{
    SelectResult r = this->select(result, map);
    
    if (hasParent())
        return r;
    else
        return this->expand(r);
}

SelectResult SireMol::parser::SelectEngine::operator()(const MolGroupsBase &molgroups,
                                                       const PropertyMap &map) const
{
    return this->operator()( SelectResult(molgroups), map );
}

SelectResult SireMol::parser::SelectEngine::operator()(const MoleculeGroup &molgroup,
                                                       const PropertyMap &map) const
{
    return this->operator()( SelectResult(molgroup), map);
}

SelectResult SireMol::parser::SelectEngine::operator()(const Molecules &molecules,
                                                       const PropertyMap &map) const
{
    return this->operator()( SelectResult(molecules), map);
}

SelectResult SireMol::parser::SelectEngine::operator()(const MoleculeView &molecule,
                                                       const PropertyMap &map) const
{
    return this->operator()( SelectResult(molecule), map);
}

/** Internal function used to set the parent of this engine */
void SireMol::parser::SelectEngine::setParent(SireMol::parser::SelectEnginePtr ptr) const
{
    const_cast<SireMol::parser::SelectEngine*>(this)->parent = ptr;
}

/** Internal function used to make a shared pointer out of the passed pointer */
SireMol::parser::SelectEnginePtr SireMol::parser::SelectEngine::makePtr(SelectEngine *ptr)
{
    SelectEnginePtr p(ptr);
    ptr->selfptr = p;
    return p;
}

/** Return a simplified version of this engine (e.g. remove double-nots etc.) */
SireMol::parser::SelectEnginePtr SireMol::parser::SelectEngine::simplify()
{
    return selfptr.lock();
}

/** Return whether or not this selection depends on coordinates of atoms */
bool SireMol::parser::SelectEngine::usesCoordinates() const
{
    return false;
}

/** Return whether or not this engine has a parent */
bool SireMol::parser::SelectEngine::hasParent() const
{
    return parent.lock().get() != 0;
}

/** Expand the passed SelectResult based on the selection type of this SelectEngine */
SelectResult SireMol::parser::SelectEngine::expand(const SelectResult &results) const
{
    const auto objtyp = this->objectType();
    
    if (objtyp == SelectEngine::COMPLEX or objtyp == SelectEngine::MOLECULE)
    {
        //we don't need to do anything
        return results;
    }

    QList<ViewsOfMol> expanded;

    if (objtyp == SelectEngine::ATOM)
    {
        for (auto result : results)
        {
            expanded.append( ViewsOfMol(result.atoms()) );
        }
    }
    else if (objtyp == SelectEngine::CUTGROUP)
    {
        for (auto result : results)
        {
            expanded.append( ViewsOfMol(result.cutGroups()) );
        }
    }
    else if (objtyp == SelectEngine::RESIDUE)
    {
        for (auto result : results)
        {
            expanded.append( ViewsOfMol(result.residues()) );
        }
    }
    else if (objtyp == SelectEngine::CHAIN)
    {
        for (auto result : results)
        {
            expanded.append( ViewsOfMol(result.chains()) );
        }
    }
    else if (objtyp == SelectEngine::SEGMENT)
    {
        for (auto result : results)
        {
            expanded.append( ViewsOfMol(result.segments()) );
        }
    }
    else
    {
        return results;
    }

    return SelectResult(expanded);
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
        return e->operator()(molgroup,map);
    else
        return SelectResult();
}

SelectResult Select::operator()(const Molecules &molecules,
                                const PropertyMap &map) const
{
    if (e.get() != 0)
        return e->operator()(molecules,map);
    else
        return SelectResult();
}

SelectResult Select::operator()(const MoleculeView &molecule,
                                const PropertyMap &map) const
{
    if (e.get() != 0)
        return e->operator()(molecule,map);
    else
        return SelectResult();
}

SelectResult Select::operator()(const SelectResult &result,
                                const PropertyMap &map) const
{
    if (e.get() != 0)
        return e->operator()(result,map);
    else
        return SelectResult();
}

/** Set a user token that will be substituted for the passed selection, e.g.

    setToken("protein", "molecules with resname /ala/i")
    
    would allow you to use "protein" to refer to any molecules that contain
    residues called /ala/i
    
    Note that the token is set globally for all searches
*/
void Select::setToken(const QString &token, const QString &selection)
{
    SireMol::parser::set_token(token, selection);
}

/** Clear all user-set tokens */
void Select::resetTokens()
{
    SireMol::parser::reset_tokens();
}

QString Select::objectType() const
{
    if (e.get() == 0)
        return QObject::tr("nothing");
    
    switch(e->objectType())
    {
    case SireMol::parser::SelectEngine::COMPLEX:
        return QObject::tr("complex view");
    case SireMol::parser::SelectEngine::ATOM:
        return QObject::tr("atoms");
    case SireMol::parser::SelectEngine::CUTGROUP:
        return QObject::tr("cutgroups");
    case SireMol::parser::SelectEngine::RESIDUE:
        return QObject::tr("residues");
    case SireMol::parser::SelectEngine::CHAIN:
        return QObject::tr("chains");
    case SireMol::parser::SelectEngine::SEGMENT:
        return QObject::tr("segments");
    case SireMol::parser::SelectEngine::MOLECULE:
        return QObject::tr("molecules");
    }
}

QString Select::toString() const
{
    if (search_string.isEmpty())
    {
        return QObject::tr("Select::null");
    }
    else
        return QObject::tr("Select( %1, result=%2 )").arg(search_string).arg(objectType());
}

///////////
/////////// Implementation of SelectResult
///////////

static const RegisterMetaType<SelectResult> r_result;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const SelectResult &result)
{
    writeHeader(ds, r_result, 1);
    
    SharedDataStream sds(ds);
    
    sds << result.molviews << static_cast<const Property&>(result);
    
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, SelectResult &result)
{
    VersionID v = readHeader(ds, r_result);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> result.molviews >> static_cast<Property&>(result);
    }
    else
        throw version_error(v, "1", r_result, CODELOC);
    
    return ds;
}

/** Constructor */
SelectResult::SelectResult() : ConcreteProperty<SelectResult,Property>()
{}

/** Construct from the passed molecules */
SelectResult::SelectResult(const Molecules &molecules)
             : ConcreteProperty<SelectResult,Property>()
{
    auto molnums = molecules.molNums().toList();
    qSort(molnums);
    
    for (const auto molnum : molnums)
    {
        const auto &view = molecules[molnum];
        
        if (not view.isEmpty())
            molviews.append(view);
    }
}

/** Construct from the passed molecules */
SelectResult::SelectResult(const MolGroupsBase &molgroups)
             : ConcreteProperty<SelectResult,Property>()
{
    this->operator=( SelectResult(molgroups.molecules()) );
}

/** Construct from the passed molecules */
SelectResult::SelectResult(const MoleculeGroup &molgroup)
             : ConcreteProperty<SelectResult,Property>()
{
    for (int i=0; i<molgroup.nMolecules(); ++i)
    {
        const auto &view = molgroup.moleculeAt(i);
        
        if (not view.isEmpty())
            molviews.append(view);
    }
}

/** Construct from the passed molecules */
SelectResult::SelectResult(const MoleculeView &molview)
             : ConcreteProperty<SelectResult,Property>()
{
    if (not molview.isEmpty())
    {
        if (molview.isA<ViewsOfMol>())
            molviews.append( molview.asA<ViewsOfMol>() );
        else
            molviews.append( ViewsOfMol(molview) );
    }
}

/** Construct from the passed molecules */
SelectResult::SelectResult(const QList<ViewsOfMol> views)
             : ConcreteProperty<SelectResult,Property>()
{
    molviews = views;
    
    bool remove_empty = false;
    
    for (const auto &view : views)
    {
        if (view.isEmpty())
        {
            remove_empty = true;
            break;
        }
    }
    
    if (remove_empty)
    {
        QMutableListIterator<ViewsOfMol> it(molviews);
        
        while (it.hasNext())
        {
            const auto &view = it.next();
            
            if (view.isEmpty())
                it.remove();
        }
    }
}

/** Copy constructor */
SelectResult::SelectResult(const SelectResult &other)
             : ConcreteProperty<SelectResult,Property>(other),
               molviews(other.molviews)
{}

/** Destructor */
SelectResult::~SelectResult()
{}

/** Copy assignment operator */
SelectResult& SelectResult::operator=(const SelectResult &other)
{
    molviews = other.molviews;
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool SelectResult::operator==(const SelectResult &other) const
{
    return molviews == other.molviews and Property::operator==(other);
}

/** Comparison operator */
bool SelectResult::operator!=(const SelectResult &other) const
{
    return not SelectResult::operator==(other);
}

const char* SelectResult::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SelectResult>() );
}

const char* SelectResult::what() const
{
    return SelectResult::typeName();
}

SelectResult* SelectResult::clone() const
{
    return new SelectResult(*this);
}

/** Return whether or not this is empty */
bool SelectResult::isEmpty() const
{
    return molviews.isEmpty();
}

/** Return the number of views in this result */
int SelectResult::count() const
{
    int total = 0;
    
    for (const auto &view : molviews)
    {
        total += view.nViews();
    }
    
    return total;
}

/** Return the number of views in this result */
int SelectResult::size() const
{
    return this->count();
}

/** Return all of the views in this result, grouped by molecule */
QList<ViewsOfMol> SelectResult::views() const
{
    return molviews;
}

/** Return the ith view in the result. This is automatically converted to
    the right molecule view type */
MolViewPtr SelectResult::operator[](int i) const
{
    i = Index(i).map( this->count() );
    
    int nmol = 0;
    
    for (const auto &view : molviews)
    {
        nmol += 1;
    
        if (i >= view.nViews())
        {
            i -= view.nViews();
        }
        else
        {
            auto mol = view.valueAt(i);
            
            if (mol.selectedAll())
            {
                return MolViewPtr( mol.molecule() );
            }
            else if (mol.nAtoms() == 1)
            {
                return mol.atom();
            }
            else if (mol.nResidues() == 1)
            {
                const auto res = mol.residue();
                
                if (mol.selection().selectedAll(res.number()))
                    return MolViewPtr(res);
            }
            else if (mol.nChains() == 1)
            {
                const auto chain = mol.chain();
                
                if (mol.selection().selectedAll(chain.name()))
                    return MolViewPtr(chain);
            }
            else if (mol.nSegments() == 1)
            {
                const auto seg = mol.segment();
                
                if (mol.selection().selectedAll(seg.name()))
                    return MolViewPtr(seg);
            }
            else if (mol.nCutGroups() == 1)
            {
                const auto cg = mol.cutGroup();
                
                if (mol.selection().selectedAll(cg.name()))
                    return MolViewPtr(cg);
            }
            
            //this is a mixed view
            return MolViewPtr(mol);
        }
    }
    
    throw SireError::program_bug( QObject::tr(
                "Error as we could not find the view..."), CODELOC );
    
    return MolViewPtr();
}

QString SelectResult::toString() const
{
    QStringList lines;

    const int nviews = this->count();

    if (nviews == 0)
        return QObject::tr("SelectResult::empty");

    else if (nviews <= 5)
    {
        for (int i=0; i<nviews; ++i)
        {
            lines.append( QString("%1 : %2").arg(i).arg(this->operator[](i).read().toString()) );
        }
    }
    else
    {
        for (int i=0; i<3; ++i)
        {
            lines.append( QString("%1 : %2").arg(i).arg(this->operator[](i).read().toString()) );
        }
        
        lines.append( "..." );

        for (int i=-2; i<0; ++i)
        {
            int idx = Index(i).map(nviews);
        
            lines.append( QString("%1 : %2").arg(idx)
                                            .arg(this->operator[](idx).read().toString()) );
        }
    }
    
    return QObject::tr("SelectResult{ count() == %1,\n  %2\n}")
                .arg(nviews).arg(lines.join(",\n  "));
}

/** Return a object that can be used to move all of the views in this result */
SelectResultMover SelectResult::move() const
{
    return SelectResultMover(*this);
}

SelectResult::const_iterator SelectResult::begin() const
{
    return molviews.begin();
}

SelectResult::const_iterator SelectResult::end() const
{
    return molviews.end();
}

SelectResult::const_iterator SelectResult::constBegin() const
{
    return begin();
}

SelectResult::const_iterator SelectResult::constEnd() const
{
    return end();
}
///////////
/////////// Implementation of SelectResult
///////////

/** Constructor */
SelectResultMover::SelectResultMover()
                  : ConcreteProperty<SelectResultMover,Property>()
{}

/** Construct from the passed SelectResult */
SelectResultMover::SelectResultMover(const SelectResult &other)
                  : ConcreteProperty<SelectResultMover,Property>()
{
    for (const auto view : other.views())
    {
        molviews.append( view.move() );
    }
}

/** Copy constructor */
SelectResultMover::SelectResultMover(const SelectResultMover &other)
                  : ConcreteProperty<SelectResultMover,Property>(), molviews(other.molviews)
{}

/** Destructor */
SelectResultMover::~SelectResultMover()
{}

/** Copy assignment operator */
SelectResultMover& SelectResultMover::operator=(const SelectResultMover &other)
{
    molviews = other.molviews;
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool SelectResultMover::operator==(const SelectResultMover &other) const
{
    return false;
}

/** Comparison operator */
bool SelectResultMover::operator!=(const SelectResultMover &other) const
{
    return not SelectResultMover::operator==(other);
}

const char* SelectResultMover::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SelectResultMover>() );
}

const char* SelectResultMover::what() const
{
    return SelectResultMover::typeName();
}

SelectResultMover* SelectResultMover::clone() const
{
    return new SelectResultMover(*this);
}

QString SelectResultMover::toString() const
{
    return this->what();
}

/** Translate all of the views by 'delta' */
SelectResultMover& SelectResultMover::translate(const Vector &delta)
{
    for (auto &view : molviews)
    {
        view.translate(delta);
    }
    
    return *this;
}

/** Commit all of the moves */
SelectResult SelectResultMover::commit() const
{
    QList<ViewsOfMol> views;
    
    for (const auto view : molviews)
    {
        views.append( view.commit() );
    }
    
    return SelectResult(views);
}
