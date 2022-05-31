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
#include "SireMol/molecules.h"
#include "SireMol/moleculegroup.h"

#include "SireMol/core.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;
using namespace SireMol;

////////
//////// implementation of parse_error
////////

const char* parse_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<parse_error>() );
}

static const RegisterMetaType<parse_error> r_parse;

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

/** Return if any of the parts in 'molecule' match this engine */
bool SireMol::parser::SelectEngine::matches(const MoleculeView &molecule,
                                            const PropertyMap &map) const
{
    auto r = this->operator()(molecule, map);

    return not r.isEmpty();
}

/** Return if all of the parts in 'molecule' match this engine */
bool SireMol::parser::SelectEngine::matchesAll(const MoleculeView &molecule,
                                               const PropertyMap &map) const
{
    auto r = this->operator()(molecule, map);

    if (r.isEmpty())
        return false;

    if (molecule.nAtoms() > 1)
    {
        return r.contains(molecule);
    }
    else
    {
        return true;
    }
}

SelectResult SireMol::parser::SelectEngine::operator()(const SelectResult &result,
                                                       const PropertyMap &map) const
{
    if (hasParent())
    {
        return this->select(result, map);
    }
    else
    {
        // need to pass in the context so that, e.g. a MolIdx search can
        // know which molecule was referred to
        PropertyMap m(map);
        m.set("_context", result);

        auto r = this->select(result, m);

        return this->expand(r);
    }
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

QString SireMol::parser::SelectEngine::toString() const
{
    return QObject::tr("%1").arg(typeid(*this).name());
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

/** Expand the passed molecule based on the selection type of this SelectEngine */
MolViewPtr SireMol::parser::SelectEngine::expandMol(const MoleculeView &mol) const
{
    switch(this->objectType())
    {
    case SelectEngine::ATOM:
        return mol.atoms();
    case SelectEngine::CUTGROUP:
        return mol.cutGroups();
    case SelectEngine::RESIDUE:
        return mol.residues();
    case SelectEngine::CHAIN:
        return mol.chains();
    case SelectEngine::SEGMENT:
        return mol.segments();
    case SelectEngine::MOLECULE:
        return mol.molecule();
    default:
        return mol;
    }
}

/** Expand the passed SelectResult based on the selection type of this SelectEngine */
SelectResult SireMol::parser::SelectEngine::expand(const SelectResult &results) const
{
    const auto objtyp = this->objectType();

    if (objtyp == SelectEngine::COMPLEX | objtyp == SelectEngine::BOND)
    {
        //we don't need to do anything (or can't do anything for bonds!)
        return results;
    }

    QList<MolViewPtr> expanded;

    if (objtyp == SelectEngine::ATOM)
    {
        for (auto result : results)
        {
            expanded.append( result->atoms() );
        }
    }
    else if (objtyp == SelectEngine::CUTGROUP)
    {
        for (auto result : results)
        {
            expanded.append( result->cutGroups() );
        }
    }
    else if (objtyp == SelectEngine::RESIDUE)
    {
        for (auto result : results)
        {
            expanded.append( result->residues() );
        }
    }
    else if (objtyp == SelectEngine::CHAIN)
    {
        for (auto result : results)
        {
            expanded.append( result->chains() );
        }
    }
    else if (objtyp == SelectEngine::SEGMENT)
    {
        for (auto result : results)
        {
            expanded.append( result->segments() );
        }
    }
    else if (objtyp == SelectEngine::MOLECULE)
    {
        for (auto result : results)
        {
            expanded.append( result->molecule() );
        }
    }
    else
    {
        qDebug() << "UNRECOGNISED TYPE" << objtyp;
        return results;
    }

    return SelectResult(expanded);
}

///////////
/////////// Implementation of Select
///////////

static const RegisterMetaType<Select> r_select;

QDataStream &operator<<(QDataStream &ds, const Select &select)
{
    writeHeader(ds, r_select, 1);

    SharedDataStream sds(ds);

    sds << select.search_string;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Select &select)
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
    case SireMol::parser::SelectEngine::BOND:
        return QObject::tr("bonds");
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

    return QObject::tr("nothing");
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

QDataStream &operator<<(QDataStream &ds, const SelectResult &result)
{
    writeHeader(ds, r_result, 1);

    SharedDataStream sds(ds);

    sds << result.molviews << static_cast<const Property&>(result);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SelectResult &result)
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
    auto molnums = molecules.molNums().values();
    std::sort(molnums.begin(), molnums.end());

    for (const auto &molnum : molnums)
    {
        const auto &view = molecules[molnum];

        if (not view.isEmpty())
            molviews.append(view);
    }
}

SelectResult::SelectResult(const QList<Molecule> &views)
             : ConcreteProperty<SelectResult,Property>()
{
    molviews.reserve(views.count());

    for (const auto &view : views)
    {
        molviews.append(view);
    }
}

SelectResult::SelectResult(const QList< Selector<Atom> > &views)
             : ConcreteProperty<SelectResult,Property>()
{
    molviews.reserve(views.count());

    for (const auto &view : views)
    {
        molviews.append(view);
    }
}

SelectResult::SelectResult(const QList< Selector<Residue> > &views)
             : ConcreteProperty<SelectResult,Property>()
{
    molviews.reserve(views.count());

    for (const auto &view : views)
    {
        molviews.append(view);
    }
}

SelectResult::SelectResult(const QList< Selector<Chain> > &views)
             : ConcreteProperty<SelectResult,Property>()
{
    molviews.reserve(views.count());

    for (const auto &view : views)
    {
        molviews.append(view);
    }
}

SelectResult::SelectResult(const QList< Selector<Segment> > &views)
             : ConcreteProperty<SelectResult,Property>()
{
    molviews.reserve(views.count());

    for (const auto &view : views)
    {
        molviews.append(view);
    }
}

SelectResult::SelectResult(const QList< Selector<CutGroup> > &views)
             : ConcreteProperty<SelectResult,Property>()
{
    molviews.reserve(views.count());

    for (const auto &view : views)
    {
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
        molviews.append(molview);
    }
}

/** Construct from the passed molecules */
SelectResult::SelectResult(const QList<ViewsOfMol> &views)
             : ConcreteProperty<SelectResult,Property>()
{
    molviews.reserve(views.count());

    for (const auto &view : views)
    {
        if (not view.isEmpty())
        {
            molviews.append(view);
        }
    }
}

/** Construct from the passed molecules */
SelectResult::SelectResult(const QList<MolViewPtr> &views)
             : ConcreteProperty<SelectResult,Property>()
{
    molviews.reserve(views.count());

    for (const auto &view : views)
    {
        if (not (view.isNull() or view->isEmpty()))
        {
            molviews.append(view);
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
        total += view->nViews();
    }

    return total;
}

/** Return the number of views in this result */
int SelectResult::size() const
{
    return this->count();
}

/** Return whether or not this set contains views of the molecule with
    number 'molnum' */
bool SelectResult::contains(MolNum molnum) const
{
    for (const auto &molview : molviews)
    {
        if (molview->data().number() == molnum)
            return true;
    }

    return false;
}

/** Return whether or not this set contains all of the atoms in the
    passed molecule */
bool SelectResult::contains(const MoleculeView &mol) const
{
    AtomSelection s;

    QList<AtomSelection> selections;

    for (const auto &molview : molviews)
    {
        if (molview->data().number() == mol.data().number())
        {
            if (s.isNull())
            {
                s = mol.selection();
            }

            auto selection = molview->selection();

            if (selection.contains(s))
            {
                return true;
            }

            // maybe we will see another view of this molecule?
            selections.append(selection);
        }
    }

    if (selections.count() < 2)
        // nope
        return false;

    auto selection = selections.takeFirst();
    selection.unite(selections);

    return selection.contains(s);
}

/** Return all of the views in this result, grouped by molecule */
QList<ViewsOfMol> SelectResult::views() const
{
    QList<ViewsOfMol> v;

    for (const auto &view : molviews)
    {
        v.append(ViewsOfMol(*view));
    }

    return v;
}

/** Return all of the views of the molecule with number 'molnum'. This
    returns an empty set of views if the molecule is not in this set */
ViewsOfMol SelectResult::views(MolNum molnum) const
{
    ViewsOfMol v;

    for (const auto &molview : molviews)
    {
        if (molview->data().number() == molnum)
        {
            for (int i=0; i<molview->nViews(); ++i)
            {
                if (v.isEmpty())
                {
                    v = molview->at(i);
                }
                else
                {
                    v += ViewsOfMol(molview->at(i));
                }
            }
        }
    }

    return v;
}

/** Return the numbers of all molecules whose views are in this set,
    in the order they appear in this set */
QList<MolNum> SelectResult::molNums() const
{
    QList<MolNum> molnums;

    for (const auto &molview : molviews)
    {
        molnums.append(molview->data().number());
    }

    return molnums;
}

/** Return the ith MolViewPtr in the underlying list */
MolViewPtr SelectResult::listAt(int i) const
{
    i = Index(i).map(molviews.count());
    return MolViewPtr(molviews.at(i)->clone());
}

/** Return the number of items in the list */
int SelectResult::listCount() const
{
    return molviews.count();
}

/** Return the ith view in the result. This is automatically converted to
    the right molecule view type */
MolViewPtr SelectResult::operator[](int i) const
{
    i = Index(i).map( this->count() );

    for (const auto &view : molviews)
    {
        if (i >= view->nViews())
        {
            i -= view->nViews();
        }
        else
        {
            return view->at(i);
        }
    }

    throw SireError::program_bug( QObject::tr(
                "Error as we could not find the view..."), CODELOC );

    return MolViewPtr();
}

/** Return any views of the molecule with number 'molnum'. This returns
    an empty molecule if there are no matching molecules */
MolViewPtr SelectResult::operator[](MolNum molnum) const
{
    return this->views(molnum);
}

/** Return the results as a list of MolViewPtrs */
QList<MolViewPtr> SelectResult::toList() const
{
    return molviews;
}

QString SelectResult::toString() const
{
    QStringList lines;

    const int nviews = this->listCount();

    if (nviews == 0)
        return QObject::tr("SelectResult::empty");

    else if (nviews <= 10)
    {
        for (int i=0; i<nviews; ++i)
        {
            lines.append( QString("%1 : %2").arg(i).arg(this->listAt(i).read().toString()) );
        }
    }
    else
    {
        for (int i=0; i<5; ++i)
        {
            lines.append( QString("%1 : %2").arg(i).arg(this->listAt(i).read().toString()) );
        }

        lines.append( "..." );

        for (int i=nviews-5; i<nviews; ++i)
        {
            lines.append( QString("%1 : %2").arg(i).arg(this->listAt(i).read().toString()) );
        }
    }

    return QObject::tr("SelectResult( size=%1,\n%2\n)")
                .arg(nviews).arg(lines.join("\n"));
}

/** Return a object that can be used to move all of the views in this result */
SelectResultMover SelectResult::move() const
{
    return SelectResultMover(*this);
}

/** Return this result as a new molecule group */
MoleculeGroup SelectResult::toGroup() const
{
    return MoleculeGroup(*this);
}

/** Return this result as a new molecule group called 'name' */
MoleculeGroup SelectResult::toGroup(const QString &name) const
{
    return MoleculeGroup(name, *this);
}

/** Return this result as a set of Molecules */
Molecules SelectResult::toMolecules() const
{
    return Molecules(*this);
}

/** Return the result of searching this result with 'search_term' */
SelectResult SelectResult::search(const QString &search_term) const
{
    return Select(search_term)(*this);
}

/** Return a copy of this result with all views joined into single views */
SelectResult SelectResult::join() const
{
    QList<MolViewPtr> result;
    result.reserve(molviews.count());

    for (const auto &mol : molviews)
    {
        result.append( PartialMolecule(*mol).toUnit() );
    }

    return SelectResult(result);
}

/** Return a copy of this result with all views split into individual atoms */
SelectResult SelectResult::atoms() const
{
    QList<MolViewPtr> result;
    result.reserve(molviews.count());

    for (const auto &mol : molviews)
    {
        result.append(mol->atoms());
    }

    return SelectResult(result);
}

/** Return a copy of this result with all views split into individual cutgroups */
SelectResult SelectResult::cutGroups() const
{
    QList<MolViewPtr> result;
    result.reserve(molviews.count());

    for (const auto &mol : molviews)
    {
        result.append(mol->cutGroups());
    }

    return SelectResult(result);
}

/** Return a copy of this result with all views split into individual residues */
SelectResult SelectResult::residues() const
{
    QList<MolViewPtr> result;
    result.reserve(molviews.count());

    for (const auto &mol : molviews)
    {
        result.append(mol->residues());
    }

    return SelectResult(result);
}

/** Return a copy of this result with all views split into individual chains */
SelectResult SelectResult::chains() const
{
    QList<MolViewPtr> result;
    result.reserve(molviews.count());

    for (const auto &mol : molviews)
    {
        result.append(mol->chains());
    }

    return SelectResult(result);
}

/** Return a copy of this result with all views split into individual segments */
SelectResult SelectResult::segments() const
{
    QList<MolViewPtr> result;
    result.reserve(molviews.count());

    for (const auto &mol : molviews)
    {
        result.append(mol->segments());
    }

    return SelectResult(result);
}

/** Return a copy of this result with all views split into individual molecules */
SelectResult SelectResult::molecules() const
{
    QList<MolViewPtr> result;
    result.reserve(molviews.count());

    for (const auto &mol : molviews)
    {
        result.append(mol->molecule());
    }

    return SelectResult(result);
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

/** Return the highest common type (e.g. SireMol::Atom, SireMol::Residue etc)
 *  that suits all of the views in this result
 */
QString SelectResult::getCommonType() const
{
    if (molviews.isEmpty())
        return QString();

    // we need to start from the largest type and work downwards
    // (as some atoms may be whole residues)

    QList<AtomSelection> selections;

    for (const auto &molview : molviews)
    {
        const auto typ = molview->what();

        if (typ == Molecule::typeName() or
            typ == Atom::typeName() or
            typ == Residue::typeName() or
            typ == Chain::typeName() or
            typ == Segment::typeName() or
            typ == CutGroup::typeName())
        {
            selections.append(molview->selection());
        }
        else
        {
            //this is a composite container
            for (int i=0; i<molview->nViews(); ++i)
            {
                selections.append(molview->at(i)->selection());
            }
        }
    }

    bool is_molecule = true;

    for (const auto &s : selections)
    {
        if (not s.isMolecule())
        {
            is_molecule = false;
            break;
        }
    }

    if (is_molecule)
        return Molecule::typeName();

    bool is_segment = true;

    for (const auto &s : selections)
    {
        if (not s.isSegment())
        {
            is_segment = false;
            break;
        }
    }

    if (is_segment)
        return Segment::typeName();

    bool is_chain = true;

    for (const auto &s : selections)
    {
        if (not s.isChain())
        {
            is_chain = false;
            break;
        }
    }

    if (is_chain)
        return Chain::typeName();

    bool is_residue = true;

    for (const auto &s : selections)
    {
        if (not s.isResidue())
        {
            is_residue = false;
            break;
        }
    }

    if (is_residue)
        return Residue::typeName();

    bool is_atom = true;

    for (const auto &s : selections)
    {
        if (not s.isAtom())
        {
            is_atom = false;
            break;
        }
    }

    if (is_atom)
        return Atom::typeName();

    bool is_cutgroup = true;

    for (const auto &s : selections)
    {
        if (not s.isCutGroup())
        {
            is_cutgroup = false;
            break;
        }
    }

    if (is_cutgroup)
        return CutGroup::typeName();

    if (is_cutgroup)
        return CutGroup::typeName();

    // ok, go with the default
    return PartialMolecule::typeName();
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
    for (const auto &view : other.views())
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

    for (const auto &view : molviews)
    {
        views.append( view.commit() );
    }

    return SelectResult(views);
}
