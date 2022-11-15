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

#include "idengine.h"
#include "approx_equal.h"

#include "SireBase/parallel.h"
#include "SireBase/booleanproperty.h"

#include "SireMol/atomelements.h"
#include "SireMol/core.h"

#include "SireMM/selectorbond.h"

#include "SireVol/space.h"
#include "SireVol/cartesian.h"

#include "SireSearch/helper_funcs.h"

#include "tostring.h"

#include <QRegExp>

using namespace SireSearch;
using namespace SireMol;
using namespace SireMM;
using namespace SireBase;
using namespace parser_idengine;

SelectEngine::ObjType _to_obj_type(AST::IDObject obj)
{
    switch(obj)
    {
    case AST::ATOM:
        return SelectEngine::ATOM;
    case AST::CUTGROUP:
        return SelectEngine::CUTGROUP;
    case AST::RESIDUE:
        return SelectEngine::RESIDUE;
    case AST::CHAIN:
        return SelectEngine::CHAIN;
    case AST::SEGMENT:
        return SelectEngine::SEGMENT;
    case AST::MOLECULE:
        return SelectEngine::MOLECULE;
    case AST::BOND:
        return SelectEngine::BOND;
    default:
        return SelectEngine::COMPLEX;
    }
}

MolViewPtr _expand(const MoleculeView &mol, SelectEngine::ObjType typ,
                   const PropertyMap &map)
{
    switch(typ)
    {
    case SelectEngine::ATOM:
        if (mol.isA< Selector<Atom> >())
            return mol;
        else
            return mol.atoms();
    case SelectEngine::CUTGROUP:
        if (mol.isA< Selector<CutGroup> >())
            return mol;
        else
            return mol.cutGroups();
    case SelectEngine::RESIDUE:
        if (mol.isA< Selector<Residue> >())
            return mol;
        else
            return mol.residues();
    case SelectEngine::CHAIN:
        if (mol.isA< Selector<Chain> >())
            return mol;
        else
            return mol.chains();
    case SelectEngine::SEGMENT:
        if (mol.isA< Selector<Segment> >())
            return mol;
        else
            return mol.segments();
    case SelectEngine::BOND:
        if (mol.isA<SelectorBond>())
            return mol;
        else
            return SelectorBond(mol, map);
    case SelectEngine::MOLECULE:
        return mol.molecule();
    default:
        return mol.atoms();
    }
}

MolViewPtr _invert(const MoleculeView &mol, SelectEngine::ObjType typ,
                   const PropertyMap &map)
{
    switch(typ)
    {
    case SelectEngine::ATOM:
        if (mol.isA< Selector<Atom> >())
            return mol.asA< Selector<Atom> >().invert();
        else
            return mol.atoms().invert();
    case SelectEngine::CUTGROUP:
        if (mol.isA< Selector<CutGroup> >())
            return mol.asA< Selector<CutGroup> >().invert();
        else
            return mol.cutGroups().invert();
    case SelectEngine::RESIDUE:
        if (mol.isA< Selector<Residue> >())
            return mol.asA< Selector<Residue> >().invert();
        else
            return mol.residues().invert();
    case SelectEngine::CHAIN:
        if (mol.isA< Selector<Chain> >())
            return mol.asA< Selector<Chain> >().invert();
        else
            return mol.chains().invert();
    case SelectEngine::SEGMENT:
        if (mol.isA< Selector<Segment> >())
            return mol.asA< Selector<Segment> >().invert();
        else
            return mol.segments().invert();
    case SelectEngine::BOND:
        if (mol.isA<SelectorBond>())
            return mol.asA<SelectorBond>().invert(map);
        else
            return SelectorBond(mol, map).invert(map);
    default:
        return mol.atoms().invert();
    }
}

MolViewPtr _invert_and_intersect(const MoleculeView &mol, SelectEngine::ObjType typ,
                                 const MoleculeView &view, const PropertyMap &map)
{
    auto s = _expand(view, typ, map);

    switch(typ)
    {
    case SelectEngine::ATOM:
        if (mol.isA< Selector<Atom> >())
            return mol.asA< Selector<Atom> >().invert().intersection(s->asA< Selector<Atom> >());
        else
            return mol.atoms().invert().intersection(s->asA< Selector<Atom> >());
    case SelectEngine::CUTGROUP:
        if (mol.isA< Selector<Atom> >())
            return mol.asA< Selector<Atom> >().invert().intersection(s->asA< Selector<Atom> >());
        else
            return mol.atoms().invert().intersection(s->asA< Selector<Atom> >());
    case SelectEngine::RESIDUE:
        if (mol.isA< Selector<Atom> >())
            return mol.asA< Selector<Atom> >().invert().intersection(s->asA< Selector<Atom> >());
        else
            return mol.atoms().invert().intersection(s->asA< Selector<Atom> >());
    case SelectEngine::CHAIN:
        if (mol.isA< Selector<Atom> >())
            return mol.asA< Selector<Atom> >().invert().intersection(s->asA< Selector<Atom> >());
        else
            return mol.atoms().invert().intersection(s->asA< Selector<Atom> >());
    case SelectEngine::SEGMENT:
        if (mol.isA< Selector<Atom> >())
            return mol.asA< Selector<Atom> >().invert().intersection(s->asA< Selector<Atom> >());
        else
            return mol.atoms().invert().intersection(s->asA< Selector<Atom> >());
    case SelectEngine::BOND:
        if (mol.isA<SelectorBond>())
            return mol.asA<SelectorBond>().invert(map).intersection(s->asA<SelectorBond>());
        else
            return SelectorBond(mol, map).invert(map).intersection(s->asA<SelectorBond>());
    default:
        return mol.atoms().invert().intersection(s->asA< Selector<Atom> >());
    }
}

MolViewPtr _intersection(const MoleculeView &mol0, const MoleculeView &mol1,
                         SelectEngine::ObjType obj, const PropertyMap &map)
{
    auto s0 = _expand(mol0, obj, map);
    auto s1 = _expand(mol1, obj, map);

    switch(obj)
    {
    case SelectEngine::ATOM:
        return s0->asA< Selector<Atom> >().intersection(s1->asA< Selector<Atom> >());
    case SelectEngine::CUTGROUP:
        return s0->asA< Selector<CutGroup> >().intersection(s1->asA< Selector<CutGroup> >());
    case SelectEngine::RESIDUE:
        return s0->asA< Selector<Residue> >().intersection(s1->asA< Selector<Residue> >());
    case SelectEngine::CHAIN:
        return s0->asA< Selector<Chain> >().intersection(s1->asA< Selector<Chain> >());
    case SelectEngine::SEGMENT:
        return s0->asA< Selector<Segment> >().intersection(s1->asA< Selector<Segment> >());
    case SelectEngine::BOND:
        return s0->asA< SelectorBond >().intersection(s1->asA< SelectorBond >());
    default:
        return s0->molecule();
    }
}

MolViewPtr _unite(const MoleculeView &mol0, const MoleculeView &mol1,
                  SelectEngine::ObjType obj, const PropertyMap &map)
{
    auto s0 = _expand(mol0, obj, map);
    auto s1 = _expand(mol1, obj, map);

    switch(obj)
    {
    case SelectEngine::ATOM:
        return s0->asA< Selector<Atom> >().add(s1->asA< Selector<Atom> >());
    case SelectEngine::CUTGROUP:
        return s0->asA< Selector<CutGroup> >().add(s1->asA< Selector<CutGroup> >());
    case SelectEngine::RESIDUE:
        return s0->asA< Selector<Residue> >().add(s1->asA< Selector<Residue> >());
    case SelectEngine::CHAIN:
        return s0->asA< Selector<Chain> >().add(s1->asA< Selector<Chain> >());
    case SelectEngine::SEGMENT:
        return s0->asA< Selector<Segment> >().add(s1->asA< Selector<Segment> >());
    case SelectEngine::BOND:
        return s0->asA< SelectorBond >().add(s1->asA< SelectorBond >());
    default:
        return s0->molecule();
    }
}

////////
//////// Implementation of the IDNameEngine
////////

IDNameEngine::IDNameEngine() : SelectEngine()
{}

// backport of Qt5 wildcardToRegularExpression
QString _wildcardToRegularExpression(const QString &pattern)
{
    const int wclen = pattern.length();
    QString rx;
    rx.reserve(wclen + wclen / 16);
    int i = 0;
    const QChar *wc = pattern.unicode();

#ifdef Q_OS_WIN
    const QLatin1Char nativePathSeparator('\\');
    const QLatin1String starEscape("[^/\\\\]*");
    const QLatin1String questionMarkEscape("[^/\\\\]");
#else
    const QLatin1Char nativePathSeparator('/');
    const QLatin1String starEscape("[^/]*");
    const QLatin1String questionMarkEscape("[^/]");
#endif

    while (i < wclen) {
        const QChar c = wc[i++];
        switch (c.unicode()) {
        case '*':
            rx += starEscape;
            break;
        case '?':
            rx += questionMarkEscape;
            break;
        case '\\':
#ifdef Q_OS_WIN
        case '/':
            rx += QLatin1String("[/\\\\]");
            break;
#endif
        case '$':
        case '(':
        case ')':
        case '+':
        case '.':
        case '^':
        case '{':
        case '|':
        case '}':
            rx += QLatin1Char('\\');
            rx += c;
            break;
        case '[':
            rx += c;
            // Support for the [!abc] or [!a-c] syntax
            if (i < wclen) {
                if (wc[i] == QLatin1Char('!')) {
                    rx += QLatin1Char('^');
                    ++i;
                }

                if (i < wclen && wc[i] == QLatin1Char(']'))
                    rx += wc[i++];

                while (i < wclen && wc[i] != QLatin1Char(']')) {
                    // The '/' appearing in a character class invalidates the
                    // regular expression parsing. It also concerns '\\' on
                    // Windows OS types.
                    if (wc[i] == QLatin1Char('/') || wc[i] == nativePathSeparator)
                        return rx;
                    if (wc[i] == QLatin1Char('\\'))
                        rx += QLatin1Char('\\');
                    rx += wc[i++];
                }
            }
            break;
        default:
            rx += c;
            break;
        }
    }

    return QRegularExpression::anchoredPattern(rx);
}

QString IDNameEngine::toString() const
{
    if (names.count() > 0)
       return QObject::tr("IDNameEngine(%1 from %2)")
                    .arg(idobject_to_string(obj))
                    .arg(names.join(", "));
    else if (regexps.count() > 0)
        return QObject::tr("IDNameEngine(%1 from %2)")
                    .arg(idobject_to_string(obj))
                    .arg("regexps");
    else
        return QObject::tr("IDNameEngine(%1 from *)")
                    .arg(idobject_to_string(obj));
}

SelectEnginePtr IDNameEngine::construct(IDObject o, NameValues vals)
{
    IDNameEngine *ptr = new IDNameEngine();
    ptr->obj = o;

    try
    {
        for (const auto &val : vals)
        {
            if (val.value.which() == 0)
            {
                RegExpValue v = boost::get<RegExpValue>(val.value);
                QString r = QString::fromStdString(v.value);

                QRegularExpression regexp;

                if (v.is_case_sensitive)
                    regexp = QRegularExpression(_wildcardToRegularExpression(r));
                else
                    regexp = QRegularExpression(_wildcardToRegularExpression(r),
                                                QRegularExpression::CaseInsensitiveOption);

                if (not regexp.isValid())
                {
                    throw SireMol::parse_error( QObject::tr("Failed to interpret the "
                      "regular expression '%1' (escaped version is '%2'). Error is '%3'")
                        .arg( QString::fromStdString(v.value) )
                        .arg(r)
                        .arg(regexp.errorString()), CODELOC );
                }

                //optimise (JIT-compile) the regular expression now as it will
                //be used a lot
                #if QT_VERSION >= QT_VERSION_CHECK(5, 4, 0)
                    regexp.optimize();  // introduced in Qt 5.4
                #endif

                ptr->regexps.append(regexp);
            }
            else if (val.value.which() == 1)
                ptr->names.append( QString::fromStdString(boost::get<std::string>(val.value)) );
        }
    }
    catch(...)
    {
        delete ptr;
        throw;
    }

    return makePtr(ptr);
}

IDNameEngine::~IDNameEngine()
{}

/** Internal function used to see if the passed string matches */
bool IDNameEngine::match(const QString &val) const
{
    //try all of the fixed names
    for (const auto &name : names)
    {
        if (name == val)
        {
            //name matches exactly
            return true;
        }
    }

    //now try all of the regexps
    for (const auto &regexp : regexps)
    {
        auto match = regexp.match(val);

        if (match.hasMatch())
        {
            //we have a regexp match. Make sure we have matched the
            //entire string
            return match.captured(0).length() == val.length();
        }
    }

    return false;
}

template<class T>
Selector<T> get_views(const MolViewPtr &mol);

template<>
Selector<Atom> get_views<Atom>(const MolViewPtr &mol)
{
    return mol->atoms();
}

template<>
Selector<Residue> get_views<Residue>(const MolViewPtr &mol)
{
    return mol->residues();
}

template<>
Selector<Chain> get_views<Chain>(const MolViewPtr &mol)
{
    return mol->chains();
}

template<>
Selector<CutGroup> get_views<CutGroup>(const MolViewPtr &mol)
{
    return mol->cutGroups();
}

template<>
Selector<Segment> get_views<Segment>(const MolViewPtr &mol)
{
    return mol->segments();
}

template<class T>
SelectResult IDNameEngine::searchName(const SelectResult &mols,
                                      bool uses_parallel) const
{
    QList< Selector<T> > matches;

    if (this->names.count() == 1 and this->regexps.isEmpty())
    {
        //this is a simple name match
        const auto id = typename T::Name(this->names.at(0));

        for (const auto &mol : mols)
        {
            auto idxs = mol->data().info().mapNoThrow(id);

            if (not idxs.isEmpty())
            {
                if (mol->selectedAll())
                    matches.append(Selector<T>(mol->data(), idxs));
                else
                {
                    // need to filter out all of the non-matching idxs
                    auto s = mol->selection();

                    QList<typename T::Index> selected_idxs;

                    for (const auto &idx : idxs)
                    {
                        if (s.selected(idx))
                            selected_idxs.append(idx);
                    }

                    if (not selected_idxs.isEmpty())
                        matches.append(Selector<T>(mol->data(), selected_idxs));
                }
            }
        }

        return SelectResult(matches);
    }

    for (const auto &mol : mols)
    {
        QList<qint64> idxs;

        auto views = get_views<T>(mol);

        if (views.nViews() == 0)
            continue;

        const auto n = views.names();

        for (qint64 i=0; i<n.count(); ++i)
        {
            if (match(n[i]))
            {
                idxs.append(i);
            }
        }

        if (idxs.count() == views.nViews())
        {
            if (idxs.count() == 1)
                matches.append(views(0));
            else
                matches.append(views);
        }
        else if (idxs.count() == 1)
        {
            matches.append(views(idxs.at(0)));
        }
        else if (not idxs.isEmpty())
        {
            matches.append(views(idxs));
        }
    }

    return SelectResult(matches);
}

#include "tostring.h"

SelectResult IDNameEngine::searchMolName(const SelectResult &mols,
                                         bool uses_parallel) const
{
    QList<Molecule> matches;

    for (const auto &mol : mols)
    {
        const auto molname = mol->data().name().value();

        if (match(molname))
        {
            matches.append(mol->molecule());
        }
    };

    return SelectResult(matches);
}

SelectResult IDNameEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    bool use_parallel = true;

    if (map["parallel"].hasValue())
    {
        use_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    switch(obj)
    {
    case AST::ATOM:
        return searchName<Atom>(mols, use_parallel);
    case AST::CUTGROUP:
        return searchName<CutGroup>(mols, use_parallel);
    case AST::RESIDUE:
        return searchName<Residue>(mols, use_parallel);
    case AST::CHAIN:
        return searchName<Chain>(mols, use_parallel);
    case AST::SEGMENT:
        return searchName<Segment>(mols, use_parallel);
    case AST::MOLECULE:
        return searchMolName(mols, use_parallel);
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDNameEngine::objectType() const
{
    return _to_obj_type(obj);
}

////////
//////// Implementation of the IDNumberEngine
////////

/** Internal function used to see if the passed integer matches any of the
    range values */
bool IDNumberEngine::match(const int idx) const
{
    for (const auto &val : vals)
    {
        if (val.which() == 0)
        {
            auto v = boost::get<RangeValue>(val);

            auto slice = v.toSlice();

            for (auto it = slice.begin(std::numeric_limits<int>::max());
                 not it.atEnd(); it.next())
            {
                if (it.value() == idx)
                    return true;
            }
        }
        else
        {
            auto v = boost::get<CompareValue>(val);

            switch(v.compare)
            {
            case ID_CMP_LT:
                if (idx < v.value)
                    return true;
                break;
            case ID_CMP_LE:
                if (idx <= v.value)
                    return true;
                break;
            case ID_CMP_AE:
                if (approx_equal(idx, v.value))
                    return true;
                break;
            case ID_CMP_EQ:
                if (idx == v.value)
                    return true;
                break;
            case ID_CMP_NE:
                if (idx != v.value)
                    return true;
                break;
            case ID_CMP_GE:
                if (idx >= v.value)
                    return true;
                break;
            case ID_CMP_GT:
                if (idx > v.value)
                    return true;
                break;
            default:
                return false;
            }
        }
    }

    return false;
}

bool _is_single_value(const RangeValues &vals)
{
    if (vals.size() == 1)
    {
        if (vals[0].which() == 0)
        {
            auto v = boost::get<RangeValue>(vals[0]);

            auto slice = v.toSlice();

            auto it = slice.begin(std::numeric_limits<int>::max());
            it.next();

            if (it.atEnd())
                return true;
        }
    }

    return false;
}

int _to_single_value(const RangeValues &vals)
{
    if (vals.size() == 1)
    {
        if (vals[0].which() == 0)
        {
            auto v = boost::get<RangeValue>(vals[0]);

            auto slice = v.toSlice();

            auto it = slice.begin(std::numeric_limits<int>::max());

            int first = it.value();
            it.next();

            if (it.atEnd())
                return first;
        }
    }

    return 0;
}

template<class T>
SelectResult IDNumberEngine::searchNum(const SelectResult &mols,
                                       bool uses_parallel) const
{
    QList< Selector<T> > matches;

    if (_is_single_value(this->vals))
    {
        //this is a simple name match
        const auto id = typename T::Number(_to_single_value(this->vals));

        for (const auto &mol : mols)
        {
            auto idxs = mol->data().info().mapNoThrow(id);

            if (not idxs.isEmpty())
            {
                if (mol->selectedAll())
                    matches.append(Selector<T>(mol->data(), idxs));
                else
                {
                    // need to filter out all of the non-matching idxs
                    auto s = mol->selection();

                    QList<typename T::Index> selected_idxs;

                    for (const auto &idx : idxs)
                    {
                        if (s.selected(idx))
                            selected_idxs.append(idx);
                    }

                    if (not selected_idxs.isEmpty())
                        matches.append(Selector<T>(mol->data(), selected_idxs));
                }
            }
        }

        return SelectResult(matches);
    }

    for (const auto &mol : mols)
    {
        QList<qint64> idxs;

        auto views = get_views<T>(mol);

        if (views.nViews() == 0)
            continue;

        const auto numbers = views.numbers();

        for (qint64 i=0; i<numbers.count(); ++i)
        {
            const auto &viewnum = numbers[i];

            if (this->match(viewnum))
            {
                idxs.append(i);
            }
        }

        if (idxs.count() == views.nViews())
        {
            if (idxs.count() == 1)
                matches.append(views(0));
            else
                matches.append(views);
        }
        else if (idxs.count() == 1)
        {
            matches.append(views(idxs.at(0)));
        }
        else if (not idxs.isEmpty())
        {
            matches.append(views(idxs));
        }
    }

    return SelectResult(matches);
}

SelectResult IDNumberEngine::searchMolNum(const SelectResult &mols,
                                          bool uses_parallel) const
{
    QList<Molecule> matches;

    for (const auto &mol : mols)
    {
        const auto molnum = mol->data().number().value();

        if (this->match(molnum))
            matches.append(mol->molecule());
    };

    return SelectResult(matches);
}

IDNumberEngine::IDNumberEngine()
{}

SelectEnginePtr IDNumberEngine::construct( IDObject o, RangeValues v )
{
    IDNumberEngine *ptr = new IDNumberEngine();
    ptr->obj = o;
    ptr->vals = v;

    return makePtr(ptr);
}

IDNumberEngine::~IDNumberEngine()
{}

SelectResult IDNumberEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    bool use_parallel = true;

    if (map["parallel"].hasValue())
    {
        use_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    switch(obj)
    {
    case AST::ATOM:
        return searchNum<Atom>(mols, use_parallel);
    case AST::RESIDUE:
        return searchNum<Residue>(mols, use_parallel);
    case AST::CHAIN:
        return searchNum<Chain>(mols, use_parallel);
    case AST::CUTGROUP:
        return searchNum<CutGroup>(mols, use_parallel);
    case AST::SEGMENT:
        return searchNum<Segment>(mols, use_parallel);
    case AST::MOLECULE:
        return searchMolNum(mols, use_parallel);
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDNumberEngine::objectType() const
{
    return _to_obj_type(obj);
}

////////
//////// Implementation of the IDIndexEngine
////////

static int map(int idx, int n)
{
    if (idx >= 0)
        return qMin(idx, n-1);
    else
        return n - qMin(abs(idx), n);
}

/** Internal function used to see if the passed integer matches any of the
    range values */
bool IDIndexEngine::match(int idx, const int count) const
{
    idx = map(idx, count);

    for (const auto &val : vals)
    {
        if (val.which() == 0)
        {
            auto v = boost::get<RangeValue>(val);

            auto slice = v.toSlice();

            for (auto it = slice.begin(count, true); not it.atEnd(); it.next())
            {
                if (idx == it.value())
                    return true;
            }
        }
        else
        {
            auto v = boost::get<CompareValue>(val);

            int value = map(v.value,count);

            if (value < count and value >= 0)
            {
                switch(v.compare)
                {
                case ID_CMP_LT:
                    if (idx < value)
                        return true;
                    break;
                case ID_CMP_LE:
                    if (idx <= value)
                        return true;
                    break;
                case ID_CMP_EQ:
                    if (idx == value)
                        return true;
                    break;
                case ID_CMP_AE:
                    if (approx_equal(idx, value))
                        return true;
                    break;
                case ID_CMP_NE:
                    if (idx != value)
                        return true;
                    break;
                case ID_CMP_GE:
                    if (idx >= value)
                        return true;
                    break;
                case ID_CMP_GT:
                    if (idx > value)
                        return true;
                    break;
                default:
                    return false;
                }
            }
        }
    }

    return false;
}

template<class T>
SelectResult IDIndexEngine::searchIdx(const SelectResult &mols,
                                      bool uses_parallel) const
{
    QList< Selector<T> > matches;

    if (_is_single_value(this->vals))
    {
        //this is a simple name match
        const auto id = typename T::Index(_to_single_value(this->vals));

        for (const auto &mol : mols)
        {
            auto idxs = mol->data().info().mapNoThrow(id);

            if (not idxs.isEmpty())
            {
                if (mol->selectedAll())
                    matches.append(Selector<T>(mol->data(), idxs));
                else
                {
                    // need to filter out all of the non-matching idxs
                    auto s = mol->selection();

                    QList<typename T::Index> selected_idxs;

                    for (const auto &idx : idxs)
                    {
                        if (s.selected(idx))
                            selected_idxs.append(idx);
                    }

                    if (not selected_idxs.isEmpty())
                        matches.append(Selector<T>(mol->data(), selected_idxs));
                }
            }
        }

        return SelectResult(matches);
    }

    for (const auto &mol : mols)
    {
        QList<qint64> idxs;

        auto views = get_views<T>(mol);

        const int count = views.nViews();

        if (views.nViews() == 0)
            continue;

        for (qint64 i=0; i<count; ++i)
        {
            if (this->match(i, count))
            {
                idxs.append(i);
            }
        }

        if (idxs.count() == views.nViews())
        {
            if (idxs.count() == 1)
                matches.append(views(0));
            else
                matches.append(views);
        }
        else if (idxs.count() == 1)
        {
            matches.append(views(idxs.at(0)));
        }
        else if (not idxs.isEmpty())
        {
            matches.append(views(idxs));
        }
    }

    return SelectResult(matches);
}

SelectResult IDIndexEngine::searchMolIdx(const SelectResult &mols,
                                         const SelectResult &context,
                                         bool uses_parallel) const
{
    QList<Molecule> matches;

    int idx = 0;
    int count = context.listCount();

    for (const auto &mol : context)
    {
        if (this->match(idx, count))
        {
            // we have found the molecule in the context - need to see how
            // much of this molecule remains in 'mols'
            if (&mols == &context)
            {
                // parent search
                matches.append(mol->molecule());
            }
            else
            {
                // child search - need to see if this molecule is still
                // in the child
                const auto molnum = mol->data().number();

                for (const auto &m : mols)
                {
                    if (m->data().number() == molnum)
                    {
                        // assume we are talking about this molecule
                        matches.append(m->molecule());
                        break;
                    }
                }
            }
        }

        idx += 1;
    };

    return SelectResult(matches);
}

IDIndexEngine::IDIndexEngine()
{}

SelectEnginePtr IDIndexEngine::construct( IDObject o, RangeValues v )
{
    IDIndexEngine *ptr = new IDIndexEngine();
    ptr->obj = o;
    ptr->vals = v;

    return makePtr(ptr);
}

IDIndexEngine::~IDIndexEngine()
{}

SelectResult IDIndexEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;

    bool use_parallel = true;

    if (map["parallel"].hasValue())
    {
        use_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    switch(obj)
    {
    case AST::ATOM:
        return searchIdx<Atom>(mols, use_parallel);
    case AST::CUTGROUP:
        return searchIdx<CutGroup>(mols, use_parallel);
    case AST::RESIDUE:
        return searchIdx<Residue>(mols, use_parallel);
    case AST::CHAIN:
        return searchIdx<Chain>(mols, use_parallel);
    case AST::SEGMENT:
        return searchIdx<Segment>(mols, use_parallel);
    case AST::MOLECULE:
    {
        if (this->hasParent())
        {
            auto context = map["_context"].value().asA<SelectResult>();
            return searchMolIdx(mols, context, use_parallel);
        }
        else
        {
            return searchMolIdx(mols, mols, use_parallel);
        }
    }
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDIndexEngine::objectType() const
{
    return _to_obj_type(obj);
}

////////
//////// Implementation of the IDAndEngine
////////

IDAndEngine::IDAndEngine()
{}

SelectEnginePtr IDAndEngine::construct(SelectEnginePtr p0, SelectEnginePtr p1)
{
    IDAndEngine *ptr = new IDAndEngine();

    auto p = makePtr(ptr);

    if (p0)
        p0->setParent(p);

    if (p1)
        p1->setParent(p);

    ptr->part0 = p0;
    ptr->part1 = p1;

    return p;
}

IDAndEngine::~IDAndEngine()
{}

SelectResult IDAndEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    if (part0.get() == 0 and part1.get() == 0)
        return SelectResult();
    else if (part0.get() == 0)
        return part1->operator()(mols, map);
    else if (part1.get() == 0)
        return part0->operator()(mols, map);

    // need to get the object type - this is the smallest of the two parts
    auto obj = qMin(part0->objectType(), part1->objectType());

    // perform the search from left to right, so search the left hand side first...
    QList<MolViewPtr> left = part0->operator()(mols, map).toList();

    QHash<MolNum, int> molnum_to_idx;
    QList<MolNum> ordered_molnums;
    molnum_to_idx.reserve(left.count());

    for (int i=0; i<left.count(); ++i)
    {
        auto molnum = left[i]->data().number();

        if (molnum_to_idx.contains(molnum))
        {
            auto idx = molnum_to_idx[molnum];
            left[i] = _intersection(*(left[idx]), *(left[i]), obj, map);
        }
        else
        {
            molnum_to_idx[molnum] = i;
            ordered_molnums.append(molnum);
        }
    }

    // now perform the right hand side on the result of the left hand side
    QSet<MolNum> seen;

    for (const auto &mol : part1->operator()(mols, map))
    {
        auto molnum = mol->data().number();

        if (molnum_to_idx.contains(molnum))
        {
            auto idx = molnum_to_idx[molnum];

            if (not left[idx]->isEmpty())
                left[idx] = _intersection(*(left[idx]), *mol, obj, map);

            if (not (seen.contains(molnum) or left[idx]->isEmpty()))
            {
                seen.insert(molnum);
            }
        }
    }

    QList<MolViewPtr> result;

    for (const auto &molnum : ordered_molnums)
    {
        if (seen.contains(molnum))
        {
            auto idx = molnum_to_idx[molnum];

            if (not left[idx]->isEmpty())
                result.append(left[idx]);
        }
    }

    return SelectResult(result);
}

SelectEnginePtr IDAndEngine::simplify()
{
    if (part0.get())
        part0 = part0->simplify();

    if (part1.get())
        part1 = part1->simplify();

    if (part0.get() == 0)
        return part1;
    else if (part1.get() == 0)
        return part0;
    else
        return selfptr.lock();
}

SelectEngine::ObjType IDAndEngine::objectType() const
{
    if (part0.get() and part1.get())
    {
        auto o0 = part0->objectType();
        auto o1 = part1->objectType();

        if (o0 == o1)
            return o0;
        else
        {
            //the object type is always the smallest, e.g. atom and residue == atom
            return qMin(o0,o1);
        }
    }
    else if (part0.get())
    {
        return part0->objectType();
    }
    else if (part1.get())
    {
        return part1->objectType();
    }
    else
        return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDOrEngine
////////

IDOrEngine::IDOrEngine()
{}

SelectEnginePtr IDOrEngine::construct(SelectEnginePtr part0, SelectEnginePtr part1)
{
    IDOrEngine *ptr = new IDOrEngine();

    auto p = makePtr(ptr);

    if (part0)
        part0->setParent(p);

    if (part1)
        part1->setParent(p);

    ptr->parts.append(part0);
    ptr->parts.append(part1);

    return p;
}

SelectEnginePtr IDOrEngine::construct(QList<SelectEnginePtr> parts)
{
    IDOrEngine *ptr = new IDOrEngine();

    auto p = makePtr(ptr);

    for (auto &part : parts)
    {
        if (part)
            part->setParent(p);
    }

    ptr->parts = parts;

    return p;
}

IDOrEngine::~IDOrEngine()
{}

SelectResult IDOrEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    if (parts.isEmpty())
        return SelectResult();
    else if (parts.count() == 1)
        return parts[0]->operator()(mols, map);

    QList<MolViewPtr> result;

    QHash<MolNum, int> molnum_to_idx;

    ObjType obj = parts[0]->objectType();

    for (const auto &mol : parts[0]->operator()(mols, map))
    {
        auto molnum = mol->data().number();

        if (molnum_to_idx.contains(molnum))
        {
            int idx = molnum_to_idx[molnum];
            result[idx] = _unite(*(result[idx]), *mol, obj, map);
        }
        else
        {
            molnum_to_idx[molnum] = result.count();
            result.append(mol);
        }
    }

    for (int i=1; i<parts.count(); ++i)
    {
        obj = qMin(obj, parts[i]->objectType());

        for (const auto &mol : parts[i]->operator()(mols, map))
        {
            auto molnum = mol->data().number();

            if (molnum_to_idx.contains(molnum))
            {
                int idx = molnum_to_idx[molnum];
                result[idx] = _unite(*(result[idx]), *mol, obj, map);
            }
            else
            {
                molnum_to_idx[molnum] = result.count();
                result.append(mol);
            }
        }
    }

    return SelectResult(result);
}

SelectEnginePtr IDOrEngine::simplify()
{
    for (auto &part : parts)
    {
        part = part->simplify();
    }

    return selfptr.lock();
}

SelectEngine::ObjType IDOrEngine::objectType() const
{
    bool set = false;
    auto o = SelectEngine::COMPLEX;

    for (auto part : parts)
    {
        if (part.get())
        {
            if (not set)
            {
                o = part->objectType();
                set = true;
            }
            else if (o != part->objectType())
            {
                //the object type is always the smallest, e.g. atom or residue == atom
                o = qMin(o, part->objectType());
            }
        }
    }

    return o;
}

////////
//////// Implementation of the IDNotEngine
////////

IDNotEngine::IDNotEngine()
{}

SelectEnginePtr IDNotEngine::construct(SelectEnginePtr part)
{
    IDNotEngine *ptr = new IDNotEngine();
    auto p = makePtr(ptr);

    if (part)
        part->setParent(p);

    ptr->part = part;

    return p;
}

IDNotEngine::~IDNotEngine()
{}

SelectResult IDNotEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    if (not part.get())
        return SelectResult();

    QList<MolViewPtr> result;

    // first, make the selection
    auto selected = part->operator()(mols, map).toList();

    // now index the selection by molecule number
    QHash<MolNum, int> molnum_to_idx;
    molnum_to_idx.reserve(selected.count());

    for (int i=0; i<selected.count(); ++i)
    {
        const auto molnum = selected[i]->data().number();

        if (molnum_to_idx.contains(molnum))
        {
            selected[molnum_to_idx[molnum]] =
                    MolViewPtr(
                        PartialMolecule(selected[i]->data(),
                                        selected[molnum_to_idx[molnum]]->selection() +
                                        selected[i]->selection()) );
        }
        else
        {
            molnum_to_idx.insert(molnum, i);
        }
    }

    auto typ = part->objectType();

    // now go through the views and see if they have made it into the
    // selection
    for (const auto &mol : mols)
    {
        const auto molnum = mol->data().number();

        if (molnum_to_idx.contains(molnum))
        {
            auto s = selected[molnum_to_idx[molnum]];

            if (not s->selectedAll())
            {
                if (mol->selectedAll())
                {
                    result.append(_invert(*s, typ, map));
                }
                else
                {
                    result.append(_invert_and_intersect(*s, typ, *mol, map));
                }
            }
        }
        else
        {
            // this wasn't selected, so need to add the whole view
            result.append(_expand(*mol, typ, map));
        }
    }

    return SelectResult(result);
}

SelectEnginePtr IDNotEngine::simplify()
{
    if (part.get())
        part = part->simplify();

    return selfptr.lock();
}

SelectEngine::ObjType IDNotEngine::objectType() const
{
    if (part)
        return part->objectType();
    else
        return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDJoinEngine
////////

IDJoinEngine::IDJoinEngine()
{}

SelectEnginePtr IDJoinEngine::construct(SelectEnginePtr part)
{
    IDJoinEngine *ptr = new IDJoinEngine();
    auto p = makePtr(ptr);

    if (part)
        part->setParent(p);

    ptr->part = part;

    return p;
}

IDJoinEngine::~IDJoinEngine()
{}

SelectResult IDJoinEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    if (not part.get())
        return SelectResult();

    //first, select the parts...
    auto selected = part->operator()(mols, map);

    QList<MolViewPtr> result;
    result.reserve(selected.count());

    for (const auto &mol : selected)
    {
        result.append( PartialMolecule(*mol).toUnit() );
    }

    return SelectResult(result);
}

SelectEnginePtr IDJoinEngine::simplify()
{
    if (part.get())
        part = part->simplify();

    return selfptr.lock();
}

SelectEngine::ObjType IDJoinEngine::objectType() const
{
    return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDSubScriptEngine
////////

IDSubScriptEngine::IDSubScriptEngine()
{}

SelectEnginePtr IDSubScriptEngine::construct(SelectEnginePtr part, const RangeValue &val)
{
    IDSubScriptEngine *ptr = new IDSubScriptEngine();
    auto p = makePtr(ptr);

    if (part)
        part->setParent(p);

    ptr->part = part;
    ptr->val = val;

    return p;
}

IDSubScriptEngine::~IDSubScriptEngine()
{}

SelectResult IDSubScriptEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    if (not part.get())
        return SelectResult();

    //first, select the parts...
    QList<MolViewPtr> all;
    auto obj = part->objectType();

    for (auto &mol : part->operator()(mols, map).toList())
    {
        all += _expand(*mol, obj, map)->toList();
    }

    const int nviews = all.count();

    auto addView = [](const MoleculeView &view, QList<MolViewPtr> &result,
                      SelectEngine::ObjType obj,
                      const PropertyMap &map)
    {
        const int molnum = view.data().number();

        for (int i=0; i<result.count(); ++i)
        {
            if (result[i]->data().number() == molnum)
            {
                result[i] = _unite(*(result[i]), view, obj, map);
                return;
            }
        }

        result.append(view);
    };

    QList<MolViewPtr> result;

    //now get the range of views to return
    auto slice = val.toSlice();

    for (auto it = slice.begin(nviews, true); not it.atEnd(); it.next())
    {
        addView( *(all[it.value()]), result, obj, map );
    }

    return SelectResult(result);
}

SelectEnginePtr IDSubScriptEngine::simplify()
{
    if (part.get())
        part = part->simplify();

    return selfptr.lock();
}

SelectEngine::ObjType IDSubScriptEngine::objectType() const
{
    if (part)
        return part->objectType();
    else
        return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDMass
////////

QString IDMass::toString() const
{
    return QObject::tr("%1")
            .arg((value * units).toString());
}

SelectEnginePtr IDMass::toEngine() const
{
    IDCmpMass m;
    m.compare = ID_CMP_AE;
    m.value = *this;

    return m.toEngine();
}

QString IDCmpMass::toString() const
{
    return QObject::tr("mass %1 %2")
            .arg(idcomparison_to_string(compare))
            .arg((value.value * value.units).toString());
}

SelectEnginePtr IDCmpMass::toEngine() const
{
    return IDMassEngine::construct(ATOM, compare, value.value * value.units);
}

QString IDObjMass::toString() const
{
    return QObject::tr("%1 mass %2")
            .arg(idobject_to_string(name))
            .arg((value.value * value.units).toString());
}

QString IDObjCmpMass::toString() const
{
    return QObject::tr("%1 mass %2 %3")
            .arg(idobject_to_string(name))
            .arg(idcomparison_to_string(compare))
            .arg((value.value * value.units).toString());
}

SelectEnginePtr IDObjMass::toEngine() const
{
    IDObjCmpMass m;
    m.name = name;
    m.compare = ID_CMP_AE;
    m.value = value;

    return m.toEngine();
}

SelectEnginePtr IDObjCmpMass::toEngine() const
{
    return IDMassEngine::construct(name, compare, value.value * value.units);
}

////////
//////// Implementation of the IDCharge
////////

QString IDCharge::toString() const
{
    return QObject::tr("IDCharge( %1 )")
            .arg((value * units).toString());
}

SelectEnginePtr IDCharge::toEngine() const
{
    IDCmpCharge c;
    c.compare = ID_CMP_AE;
    c.value = *this;

    return c.toEngine();
}

QString IDCmpCharge::toString() const
{
    return QObject::tr("charge %1 %2")
            .arg(idcomparison_to_string(compare))
            .arg((value.value * value.units).toString());
}

SelectEnginePtr IDCmpCharge::toEngine() const
{
    return IDChargeEngine::construct(ATOM, compare, value.value * value.units);
}

QString IDObjCharge::toString() const
{
    return QObject::tr("%1 charge %2")
            .arg(idobject_to_string(name))
            .arg((value.value * value.units).toString());
}

QString IDObjCmpCharge::toString() const
{
    return QObject::tr("%1 charge %2 %3")
            .arg(idobject_to_string(name))
            .arg(idcomparison_to_string(compare))
            .arg((value.value * value.units).toString());
}

SelectEnginePtr IDObjCharge::toEngine() const
{
    IDObjCmpCharge c;
    c.name = name;
    c.compare = ID_CMP_AE;
    c.value = value;

    return c.toEngine();
}

SelectEnginePtr IDObjCmpCharge::toEngine() const
{
    return IDChargeEngine::construct(name, compare, value.value * value.units);
}

////////
//////// Implementation of the IDMassEngine
////////

IDMassEngine::IDMassEngine()
{}

SelectEngine::ObjType IDMassEngine::objectType() const
{
    return _to_obj_type(obj);
}

SelectEnginePtr IDMassEngine::construct(IDObject obj, IDComparison compare,
                                        SireUnits::Dimension::MolarMass value)
{
    IDMassEngine *ptr = new IDMassEngine();
    auto p = makePtr(ptr);

    ptr->obj = obj;
    ptr->compare = compare;
    ptr->value = value;

    return p;
}

IDMassEngine::~IDMassEngine()
{}

#include <functional>

std::function<bool (double, double)> _get_compare(IDComparison compare)
{
    switch(compare)
    {
    case ID_CMP_LT:
        return std::less<double>();
    case ID_CMP_LE:
        return std::less_equal<double>();
    case ID_CMP_EQ:
        return std::equal_to<double>();
    case ID_CMP_AE:
        return SireSearch::approx_equal;
    case ID_CMP_NE:
        return std::not_equal_to<double>();
    case ID_CMP_GE:
        return std::greater_equal<double>();
    case ID_CMP_GT:
        return std::greater<double>();
    default:
        return std::equal_to<double>();
    }
}

template<class T>
SelectResult IDMassEngine::select_t(const SelectResult &mols,
                                    const PropertyMap &map) const
{
    if (mols.count() == 0)
        return SelectResult();

    QList<MolViewPtr> ret;

    auto compare_func = _get_compare(compare);

    for (const auto &mol : mols)
    {
        auto views = get_views<T>(mol);

        if (views.isEmpty())
            continue;

        QList<qint64> idxs;

        for (qint64 i=0; i<views.nViews(); ++i)
        {
            if (compare_func(views(i).evaluate().mass(map).value(), value.value()))
            {
                idxs.append(i);
            }
        }

        if (idxs.count() > 0)
        {
            if (idxs.count() == views.nViews())
            {
                ret.append(views);
            }
            else
            {
                ret.append(views(idxs));
            }
        }
    }

    return SelectResult(ret);
}

SelectResult IDMassEngine::select_bonds(const SelectResult &mols,
                                        const PropertyMap &map) const
{
    if (mols.count() == 0)
        return SelectResult();

    QList<MolViewPtr> ret;

    auto compare_func = _get_compare(compare);

    for (const auto &mol : mols)
    {
        auto bonds = SelectorBond(*mol, map);

        if (bonds.isEmpty())
            continue;

        QList<qint64> idxs;

        for (qint64 i=0; i<bonds.nViews(); ++i)
        {
            if (compare_func(bonds(i).evaluate().mass(map).value(), value.value()))
            {
                idxs.append(i);
            }
        }

        if (idxs.count() > 0)
        {
            if (idxs.count() == bonds.nViews())
            {
                ret.append(bonds);
            }
            else
            {
                ret.append(bonds(idxs));
            }
        }
    }

    return SelectResult(ret);
}

SelectResult IDMassEngine::select_mols(const SelectResult &mols,
                                       const PropertyMap &map) const
{
    if (mols.count() == 0)
        return SelectResult();

    QList<MolViewPtr> ret;

    auto compare_func = _get_compare(compare);

    for (const auto &mol : mols)
    {
        auto molecule = mol->molecule();

        if (compare_func(molecule.evaluate().mass(map).value(), value.value()))
        {
            ret.append(molecule);
        }
    }

    return SelectResult(ret);
}

SelectResult IDMassEngine::select(const SelectResult &mols,
                                  const PropertyMap &map) const
{
    switch (obj)
    {
    case AST::ATOM:
        return this->select_t<Atom>(mols, map);
    case AST::RESIDUE:
        return this->select_t<Residue>(mols, map);
    case AST::CHAIN:
        return this->select_t<Chain>(mols, map);
    case AST::SEGMENT:
        return this->select_t<Segment>(mols, map);
    case AST::CUTGROUP:
        return this->select_t<CutGroup>(mols, map);
    case AST::MOLECULE:
        return this->select_mols(mols, map);
    case AST::BOND:
        return this->select_bonds(mols, map);
    default:
        throw SireError::invalid_key(QObject::tr("Unsupported search object!"));
    }

    return SelectResult();
}

////////
//////// Implementation of the IDChargeEngine
////////

IDChargeEngine::IDChargeEngine()
{}

SelectEngine::ObjType IDChargeEngine::objectType() const
{
    return _to_obj_type(obj);
}

SelectEnginePtr IDChargeEngine::construct(IDObject obj, IDComparison compare,
                                          SireUnits::Dimension::Charge value)
{
    IDChargeEngine *ptr = new IDChargeEngine();
    auto p = makePtr(ptr);

    ptr->obj = obj;
    ptr->compare = compare;
    ptr->value = value;

    return p;
}

IDChargeEngine::~IDChargeEngine()
{}

template<class T>
SelectResult IDChargeEngine::select_t(const SelectResult &mols,
                                      const PropertyMap &map) const
{
    if (mols.count() == 0)
        return SelectResult();

    QList<MolViewPtr> ret;

    auto compare_func = _get_compare(compare);

    for (const auto &mol : mols)
    {
        auto views = get_views<T>(mol);

        if (views.isEmpty())
            continue;

        QList<qint64> idxs;

        for (qint64 i=0; i<views.nViews(); ++i)
        {
            if (compare_func(views(i).evaluate().charge(map).value(), value.value()))
            {
                idxs.append(i);
            }
        }

        if (idxs.count() > 0)
        {
            if (idxs.count() == views.nViews())
            {
                ret.append(views);
            }
            else
            {
                ret.append(views(idxs));
            }
        }
    }

    return SelectResult(ret);
}

SelectResult IDChargeEngine::select_bonds(const SelectResult &mols,
                                          const PropertyMap &map) const
{
    if (mols.count() == 0)
        return SelectResult();

    QList<MolViewPtr> ret;

    auto compare_func = _get_compare(compare);

    for (const auto &mol : mols)
    {
        auto bonds = SelectorBond(*mol, map);

        if (bonds.isEmpty())
            continue;

        QList<qint64> idxs;

        for (qint64 i=0; i<bonds.nViews(); ++i)
        {
            if (compare_func(bonds(i).evaluate().charge(map).value(), value.value()))
            {
                idxs.append(i);
            }
        }

        if (idxs.count() > 0)
        {
            if (idxs.count() == bonds.nViews())
            {
                ret.append(bonds);
            }
            else
            {
                ret.append(bonds(idxs));
            }
        }
    }

    return SelectResult(ret);
}

SelectResult IDChargeEngine::select_mols(const SelectResult &mols,
                                         const PropertyMap &map) const
{
    if (mols.count() == 0)
        return SelectResult();

    QList<MolViewPtr> ret;

    auto compare_func = _get_compare(compare);

    for (const auto &mol : mols)
    {
        auto molecule = mol->molecule();

        if (compare_func(molecule.evaluate().charge(map).value(), value.value()))
        {
            ret.append(molecule);
        }
    }

    return SelectResult(ret);
}

SelectResult IDChargeEngine::select(const SelectResult &mols,
                                    const PropertyMap &map) const
{
    switch (obj)
    {
    case AST::ATOM:
        return this->select_t<Atom>(mols, map);
    case AST::RESIDUE:
        return this->select_t<Residue>(mols, map);
    case AST::CHAIN:
        return this->select_t<Chain>(mols, map);
    case AST::SEGMENT:
        return this->select_t<Segment>(mols, map);
    case AST::CUTGROUP:
        return this->select_t<CutGroup>(mols, map);
    case AST::MOLECULE:
        return this->select_mols(mols, map);
    case AST::BOND:
        return this->select_bonds(mols, map);
    default:
        throw SireError::invalid_key(QObject::tr("Unsupported search object!"));
    }

    return SelectResult();
}

////////
//////// Implementation of the IDPropertyEngine
////////

IDPropertyEngine::IDPropertyEngine()
{}

SelectEnginePtr IDPropertyEngine::construct(const IDObject &name,
                                            const QString &property,
                                            const IDComparison &compare,
                                            const QString &value)
{
    IDPropertyEngine *ptr = new IDPropertyEngine();
    auto p = makePtr(ptr);

    ptr->name = name;
    ptr->property = property;
    ptr->compare = compare;
    ptr->value = value;

    return p;
}

IDPropertyEngine::~IDPropertyEngine()
{}

/** Return whether or not the passed string is True.
 *
 *  If 'pytrue_only' is true, then this only
 *  tests against the Python 'True' value as anything else risks
 *  causing false matches against numeric values if we actually
 *  want to find values that are equal to 1 (versus equal to 4)
 *
 *  Otherwise this compares against any case comparison of
 *  'true', plus the integer not zero (0) or float not zero (0.0)
*/
inline bool _is_true(const QString &value, bool pytrue_only=true)
{
    if (pytrue_only)
        return value == "True";
    else
        return value.toLower() == "true" or (value != "0" and value != "0.0");
}

/** Return whether or not the passed string is False
 *
 *  If 'pyfalse_only' is true, then this only
 *  tests against the Python 'False' value as anything else risks
 *  causing false matches against numeric values.
 *
 *  Otherwise this compares against any case comparison of
 *  'false', plus the integer zero (0) or float zero (0.0)
*/
inline bool _is_false(const QString &value, bool pyfalse_only=true)
{
    if (pyfalse_only)
        return value == "False";
    else
        return value.toLower() == "false" or value == "0" or value == "0.0";
}

bool _compare_equal(const QString &left,
                    const IDComparison &compare,
                    const QString &right)
{
    // if they are strings, then this should be ok
    if (left == right)
        return true;

    if (compare == ID_CMP_EQ and _is_true(right))
    {
        // we are just asking if 'left' is anything other than false
        if (not _is_false(left, false))
        {
            return true;
        }
    }

    // could they both be numbers?
    {
        bool left_ok;
        double left_num = left.toDouble(&left_ok);

        bool right_ok;
        double right_num = right.toDouble(&right_ok);

        if (left_ok and right_ok)
        {
            auto compare_func = _get_compare(compare);
            return compare_func(left_num, right_num);
        }
    }

    // could they both be booleans?
    try
    {
        BooleanProperty left_bool(left);
        BooleanProperty right_bool(right);

        return left_bool.value() == right_bool.value();
    }
    catch(...)
    {
        // not bools - this is ok.
    }

    if (compare == ID_CMP_AE)
    {
        // allow simplified, case-insensitive comparison
        return left.toLower().simplified() == right.toLower().simplified();
    }

    return false;
}

bool _compare(const QString &left,
              const IDComparison &compare,
              const QString &right)
{
    if (compare == ID_CMP_NE)
    {
        return not _compare(left, ID_CMP_EQ, right);
    }
    else if (compare == ID_CMP_EQ or compare == ID_CMP_AE)
    {
        return _compare_equal(left, compare, right);
    }
    else
    {
        // all other comparisons must be numeric. Check if these are both
        // numbers
        bool ok;

        double left_num = left.toDouble(&ok);

        if (!ok)
            return false;

        double right_num = right.toDouble(&ok);

        if (!ok)
            return false;

        auto compare_func = _get_compare(compare);

        return compare_func(left_num, right_num);
    }
}

template<class T>
MolViewPtr _select_property_(const MolViewPtr &mol,
                             const PropertyName &property,
                             const IDComparison &compare,
                             const QString &value)
{
    auto views = get_views<T>(mol);

    QList<qint64> idxs;

    for (int i=0; i<views.count(); ++i)
    {
        const auto &view = views(i);

        if (view.hasProperty(property))
        {
            try
            {
                auto p = view.propertyAsProperty(property)->asAString();

                if (_compare(p, compare, value))
                {
                    idxs.append(i);
                }
            }
            catch(...)
            {
                if (compare == ID_CMP_EQ and _is_true(value))
                {
                    idxs.append(i);
                }
            }
        }
    }

    if (idxs.count() == views.count())
        //everything has been selected
        return views;
    else if (not idxs.isEmpty())
        return views(idxs);
    else
        return Selector<T>();
}

MolViewPtr _select_property_bond(const SelectorBond &bonds,
                                 const PropertyName &property,
                                 const IDComparison &compare,
                                 const QString &value)
{
    QList<qint64> idxs;

    for (int i=0; i<bonds.count(); ++i)
    {
        auto bond = bonds(i);

        if (bond.hasProperty(property))
        {
            try
            {
                if (_compare(bond.property(property).asAString(), compare, value))
                {
                    idxs.append(i);
                }
            }
            catch(...)
            {
                if (compare == ID_CMP_EQ and _is_true(value))
                {
                    idxs.append(i);
                }
            }
        }
    }

    if (idxs.count() == bonds.count())
        return bonds;
    else if (not idxs.isEmpty())
        return bonds(idxs);
    else
        return SelectorBond();
}

MolViewPtr _select_property(const MolViewPtr &mol,
                            const IDObject &name,
                            const PropertyName &property,
                            const IDComparison &compare,
                            const QString &value,
                            const PropertyMap &map)
{
    switch(name)
    {
    case AST::ATOM:
        return _select_property_<Atom>(mol, property, compare, value);
    case AST::RESIDUE:
        return _select_property_<Residue>(mol, property, compare, value);
    case AST::CHAIN:
        return _select_property_<Chain>(mol, property, compare, value);
    case AST::SEGMENT:
        return _select_property_<Segment>(mol, property, compare, value);
    case AST::CUTGROUP:
        return _select_property_<CutGroup>(mol, property, compare, value);
    case AST::BOND:
    {
        if (mol->isA<SelectorBond>())
            return _select_property_bond(mol->asA<SelectorBond>(), property, compare, value);
        else
            return _select_property_bond(SelectorBond(*mol, map), property, compare, value);
    }
    default:
        qDebug() << "UNRECOGNISED" << name;
        return MolViewPtr();
    }
}

SelectResult IDPropertyEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    QList<MolViewPtr> ret;

    const auto p = map[this->property];

    if (name == AST::MOLECULE)
    {
        for (const auto &mol : mols)
        {
            auto m = mol->molecule();

            if (m.hasProperty(p))
            {
                try
                {
                    if (_compare(m.property(p).asAString(), this->compare, this->value))
                    {
                        ret.append(m);
                    }
                }
                catch(...)
                {
                    // this isn't a property that can be converted to a string.
                    if (this->compare == ID_CMP_EQ and _is_true(value))
                    {
                        // we are only asking if this property exists
                        ret.append(m);
                    }
                }
            }
        }
    }
    else
    {
        for (const auto &mol : mols)
        {
            try
            {
                auto selected = _select_property(mol, this->name, p, this->compare, this->value, map);

                if (not selected->isEmpty())
                    ret.append(selected);
            }
            catch(...)
            {
                // this isn't a property that can be converted to a string
            }
        }
    }

    return SelectResult(ret);
}

SelectEnginePtr IDPropertyEngine::simplify()
{
    return selfptr.lock();
}

SelectEngine::ObjType IDPropertyEngine::objectType() const
{
    return _to_obj_type(name);
}

////////
//////// Implementation of the IDBondEngine
////////

IDBondEngine::IDBondEngine()
{}

SelectEnginePtr IDBondEngine::construct( IDBondToken from_token,
                                         SelectEnginePtr from_value,
                                         IDBondToken to_token,
                                         SelectEnginePtr to_value )
{
    IDBondEngine *ptr = new IDBondEngine();
    auto p = makePtr(ptr);

    if (to_token == ID_BOND_TO)
    {
        if (from_token != ID_BOND_FROM)
        {
            throw SireMol::parse_error(QObject::tr(
                "Invalid syntax: Should be 'bonds from X to Y', not "
                "'bonds %1 X %2 Y'")
                    .arg(AST::idbondtoken_to_string(from_token))
                    .arg(AST::idbondtoken_to_string(to_token)), CODELOC);
        }
    }
    else if (to_token != ID_BOND_UNKNOWN)
    {
        throw SireMol::parse_error(QObject::tr(
            "Invalid 'to' token, %1. Should only be 'to'.")
                .arg(AST::idbondtoken_to_string(to_token)),
                    CODELOC);
    }

    if (from_token == ID_BOND_FROM and to_token != ID_BOND_TO)
    {
        throw SireMol::parse_error(QObject::tr(
            "Invalid syntax: Should be 'bonds from X to Y', not "
            "'bonds %1 X %2 Y'")
                .arg(AST::idbondtoken_to_string(from_token))
                .arg(AST::idbondtoken_to_string(to_token)), CODELOC);
    }

    if (from_value)
        from_value->setParent(p);
    else
        throw SireMol::parse_error(QObject::tr(
            "Missing first (from) group to search..."), CODELOC);

    if (to_value)
    {
        if (to_token != ID_BOND_TO)
            throw SireMol::parse_error(QObject::tr(
                "Should not have a 'to' group if not running a 'to' match!"),
                    CODELOC);

        ptr->to_token = ID_BOND_TO;
        ptr->to_value = to_value;
        to_value->setParent(p);
    }
    else if (to_token != ID_BOND_UNKNOWN)
    {
        throw SireMol::parse_error(QObject::tr(
            "Missing 'to' group when running a 'to' search!"), CODELOC);
    }

    ptr->from_token = from_token;
    ptr->from_value = from_value;

    return p;
}

IDBondEngine::~IDBondEngine()
{}

#include "SireMM/selectorbond.h"
using namespace SireMM;

SelectResult IDBondEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    if (not from_value.get() or from_token == ID_BOND_UNKNOWN)
        return SelectResult();

    SelectResult result = from_value->operator()(mols, map);

    if (result.count() == 0)
        return SelectResult();

    SelectResult to_result;

    if (from_token == ID_BOND_FROM and to_token == ID_BOND_TO)
    {
        to_result = to_value->operator()(mols, map);

        if (to_result.count() == 0)
            return SelectResult();
    }

    QList<MolViewPtr> ret;

    for (const auto &mol : result)
    {
        if (from_token == ID_BOND_TO)
        {
            // only bonds to this object, not wholly contained in this object
            auto from = mol->atoms();
            auto to = from.invert();
            auto result = SelectorBond(from, to, map);

            ret.append(result);
        }
        else if (from_token == ID_BOND_FROM and to_token == ID_BOND_TO)
        {
            for (const auto &to_mol : to_result)
            {
                if (mol->data().number() == to_mol->data().number())
                    ret.append(SelectorBond(mol->atoms(), to_mol->atoms(), map));
            }
        }
    }

    return SelectResult(ret);
}

SelectEnginePtr IDBondEngine::simplify()
{
    if (from_value.get())
        from_value = from_value->simplify();

    if (to_value.get())
        to_value = to_value->simplify();

    return selfptr.lock();
}

SelectEngine::ObjType IDBondEngine::objectType() const
{
    return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDWithEngine
////////

IDWithEngine::IDWithEngine()
{}

SelectEnginePtr IDWithEngine::construct( SelectEnginePtr part0,
                                         IDToken token,
                                         SelectEnginePtr part1)
{
    IDWithEngine *ptr = new IDWithEngine();
    auto p = makePtr(ptr);

    if (part0.get() == 0)
        part0 = IDAllEngine::construct();

    if (part1.get() == 0)
        part1 = IDAllEngine::construct();

    part0->setParent(p);
    part1->setParent(p);

    ptr->part0 = part0;
    ptr->token = token;
    ptr->part1 = part1;

    return p;
}

IDWithEngine::~IDWithEngine()
{}

SelectResult IDWithEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    QList<MolViewPtr> result;

    auto first = part0;
    auto second = part1;

    if (token == AST::ID_IN)
    {
        // an "in" search is an inverted "with" search
        first = part1;
        second = part0;
    }

    // we need to expand bonds manually
    const bool is_bond = part0->objectType() == SelectEngine::BOND;

    // need to expand the result so that we have the right type of unit
    // when we try to match the second part
    for (const auto &mol : first->expand(first->operator()(mols, map)))
    {
        const auto units = mol->toList();

        QList<qint64> matches;
        matches.reserve(units.count());

        for (int i=0; i<units.count(); ++i)
        {
            if (second->matches(*(units[i]), map))
            {
                matches.append(i);
            }
        }

        if (not matches.isEmpty())
        {
            if (matches.count() == units.count())
            {
                if (is_bond)
                {
                    auto bonds = SelectorBond(mol, map);
                    if (not bonds.isEmpty())
                        result.append(bonds);
                }
                else
                    result.append(mol);
            }
            else
            {
                //rejoin the matches into the appropriate Selector
                if (is_bond)
                {
                    auto bonds = SelectorBond(mol->operator[](matches), map);
                    if (not bonds.isEmpty())
                        result.append(bonds);
                }
                else
                    result.append(mol->operator[](matches));
            }
        }
    }

    return SelectResult(result);
}

SelectEnginePtr IDWithEngine::simplify()
{
    if (part0.get())
        part0 = part0->simplify();

    if (part1.get())
        part1 = part1->simplify();

    return selfptr.lock();
}

SelectEngine::ObjType IDWithEngine::objectType() const
{
    if (part0.get())
        return part0->objectType();
    else if (part1.get())
        return part1->objectType();
    else
        return SelectEngine::ATOM;
}

////////
//////// Implementation of the IDElementEngine
////////

IDElementEngine::IDElementEngine()
{}

SelectEnginePtr IDElementEngine::construct(const std::vector<SireMol::Element> &vals)
{
    IDElementEngine *ptr = new IDElementEngine();
    auto p = makePtr(ptr);

    for (const auto &value : vals)
    {
        ptr->elements.insert(value);
    }

    return p;
}

IDElementEngine::~IDElementEngine()
{}

SelectResult IDElementEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    QList<MolViewPtr> result;

    const auto element_property = map["element"];

    for (const auto &mol : mols)
    {
        const auto atoms = mol->atoms();

        QList<qint64> matches;
        matches.reserve(atoms.count());

        try
        {
            for (int i=0; i<atoms.count(); ++i)
            {
                const Atom &atom = atoms(i);

                const Element &element = atom.property<Element>(element_property);

                if (this->elements.contains(element))
                {
                    matches.append(i);
                }
            }
        }
        catch(...)
        {
            // no matching element property - this would be the same
            // for all atoms
        }

        if (not matches.isEmpty())
        {
            if (matches.count() == atoms.count())
            {
                // all atoms matched
                result.append(atoms);
            }
            else
            {
                result.append(atoms(matches));
            }
        }
    }

    return SelectResult(result);
}

SelectEngine::ObjType IDElementEngine::objectType() const
{
    return SelectEngine::ATOM;
}

////////
//////// Implementation of the IDDistanceEngine
////////

IDDistanceEngine::IDDistanceEngine()
{}

SelectEnginePtr IDDistanceEngine::construct(IDObject obj, IDCoordType typ,
                                            SireUnits::Dimension::Length distance,
                                            SelectEnginePtr part)
{
    IDDistanceEngine *ptr = new IDDistanceEngine();
    auto p = makePtr(ptr);

    if (part)
        part->setParent(p);

    ptr->part = part;
    ptr->obj = obj;
    ptr->typ = typ;
    ptr->distance = distance.value();

    return p;
}

SelectEnginePtr IDDistanceEngine::construct(IDObject obj,
                                            SireUnits::Dimension::Length distance,
                                            SelectEnginePtr part)
{
    return IDDistanceEngine::construct(obj, ID_COORD_CLOSEST, distance, part);
}

IDDistanceEngine::~IDDistanceEngine()
{}

SelectEnginePtr IDDistanceEngine::simplify()
{
    if (part.get() != 0)
        part->simplify();

    return selfptr.lock();
}

QVector<Vector> _get_coords(const Selector<Atom> &atoms,
                            const PropertyName &coords_property,
                            const Space &space)
{
    QVector<Vector> coords;

    coords.reserve(atoms.count());

    for (int i=0; i<atoms.count(); ++i)
    {
        const Atom &atom = atoms(i);

        coords.append(space.getMinimumImage(
                            atom.property<Vector>(coords_property),
                            Vector(0)));
    }

    return coords;
}

Vector _get_point(const AABox &box, const IDCoordType &typ)
{
    Vector center = box.center();

    switch(typ)
    {
        case ID_COORD_CENTER:
        case ID_COORD_CENTER_X:
        case ID_COORD_CENTER_Y:
        case ID_COORD_CENTER_Z:
        case ID_COORD_X:
        case ID_COORD_Y:
        case ID_COORD_Z:
            return center;
        case ID_COORD_MAX:
            return box.maxCoords();
        case ID_COORD_MIN:
            return box.minCoords();
        case ID_COORD_MAX_X:
            center.setX(box.maxCoords().x());
            return center;
        case ID_COORD_MAX_Y:
            center.setY(box.maxCoords().y());
            return center;
        case ID_COORD_MAX_Z:
            center.setZ(box.maxCoords().z());
            return center;
        case ID_COORD_MIN_X:
            center.setX(box.minCoords().x());
            return center;
        case ID_COORD_MIN_Y:
            center.setY(box.minCoords().y());
            return center;
        case ID_COORD_MIN_Z:
            center.setZ(box.minCoords().z());
            return center;
        default:
            qDebug() << "UNRECOGNISED ID TYP" << typ;
            return center;
    }
}

bool _is_within(const QVector<Vector> &coords,
                const CoordGroup &ref_group,
                double distance,
                const IDCoordType &typ, const Space &space)
{
    CoordGroup group(coords);

    if (typ == ID_COORD_CLOSEST)
    {
        return space.minimumDistance(group, ref_group) < distance;
    }

    Vector point = _get_point(group.aaBox(), typ);

    return space.minimumDistance(CoordGroup(1, point), ref_group) < distance;
}

SelectResult IDDistanceEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    //first, get the objects against where the distance is calculated
    if (part.get() == 0)
        return SelectResult();

    QList<MolViewPtr> ret;

    const auto refmols = part->operator()(mols, map);

    if (refmols.isEmpty())
        //nothing against which to compare
        return SelectResult();

    const auto coords_property = map["coordinates"];

    SireVol::SpacePtr space = SireVol::Cartesian();

    if (map["space"].hasValue())
    {
        space = map["space"].value().asA<SireVol::Space>();
    }

    // merge all atoms in 'refmols' into a single array of coordinates
    QVector<Vector> ref_coords;

    for (const auto &mol : refmols)
    {
        ref_coords += _get_coords(mol->atoms(), coords_property, *space);
    }

    CoordGroup ref_group(ref_coords);

    for (const auto &mol : mols)
    {
        // expand this molecule into the views that are requested
        const auto expanded = _expand(*mol, this->objectType(), map);

        QList<qint64> idxs;
        idxs.reserve(expanded->count());

        for (int i=0; i<expanded->count(); ++i)
        {
            const auto view = expanded->operator[](i);

            QVector<Vector> coords = _get_coords(view->atoms(),
                                                 coords_property, *space);

            if (_is_within(coords, ref_group, distance, typ, *space))
            {
                idxs.append(i);
            }
        }

        if (idxs.count() == expanded->count())
        {
            ret.append(expanded);
        }
        else if (not idxs.isEmpty())
        {
            ret.append(expanded->operator[](idxs));
        }
    }

    return SelectResult(ret);
}

SelectEngine::ObjType IDDistanceEngine::objectType() const
{
    return _to_obj_type(obj);
}

////////
//////// Implementation of the IDDistanceVectorEngine
////////

IDDistanceVectorEngine::IDDistanceVectorEngine()
{}

SelectEnginePtr IDDistanceVectorEngine::construct(IDObject obj, IDCoordType typ,
                                                  SireUnits::Dimension::Length distance,
                                                  VectorValue position)
{
    IDDistanceVectorEngine *ptr = new IDDistanceVectorEngine();
    auto p = makePtr(ptr);

    ptr->position = position;
    ptr->obj = obj;
    ptr->typ = typ;
    ptr->distance = distance.value();

    return p;
}

SelectEnginePtr IDDistanceVectorEngine::construct(IDObject obj,
                                                  SireUnits::Dimension::Length distance,
                                                  VectorValue position)
{
    return IDDistanceVectorEngine::construct(obj, ID_COORD_CLOSEST, distance, position);
}

IDDistanceVectorEngine::~IDDistanceVectorEngine()
{}

SelectEnginePtr IDDistanceVectorEngine::simplify()
{
    return selfptr.lock();
}

SelectResult IDDistanceVectorEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    QList<MolViewPtr> ret;

    const auto coords_property = map["coordinates"];

    SireVol::SpacePtr space = SireVol::Cartesian();

    if (map["space"].hasValue())
    {
        space = map["space"].value().asA<SireVol::Space>();
    }

    // turn the reference point into a single-point CoordGroup
    Vector point((position.x.value * position.x.unit).to(SireUnits::angstrom),
                 (position.y.value * position.y.unit).to(SireUnits::angstrom),
                 (position.z.value * position.z.unit).to(SireUnits::angstrom));

    QVector<Vector> ref_coords(1, point);
    CoordGroup ref_group(ref_coords);

    for (const auto &mol : mols)
    {
        // expand this molecule into the views that are requested
        const auto expanded = _expand(*mol, this->objectType(), map);

        QList<qint64> idxs;
        idxs.reserve(expanded->count());

        for (int i=0; i<expanded->count(); ++i)
        {
            const auto view = expanded->operator[](i);

            QVector<Vector> coords = _get_coords(view->atoms(),
                                                 coords_property, *space);

            if (_is_within(coords, ref_group, distance, typ, *space))
            {
                idxs.append(i);
            }
        }

        if (idxs.count() == expanded->count())
        {
            ret.append(expanded);
        }
        else if (not idxs.isEmpty())
        {
            ret.append(expanded->operator[](idxs));
        }
    }

    return SelectResult(ret);
}

SelectEngine::ObjType IDDistanceVectorEngine::objectType() const
{
    return _to_obj_type(obj);
}

////////
//////// Implementation of the IDAllEngine
////////

IDAllEngine::IDAllEngine()
{}

SelectEnginePtr IDAllEngine::construct(IDObject obj)
{
    IDAllEngine *ptr = new IDAllEngine();
    auto p = makePtr(ptr);

    ptr->obj = obj;

    return p;
}

IDAllEngine::~IDAllEngine()
{}

SelectResult _get_bonds(const SelectResult &mols,
                        const PropertyMap &map)
{
    QList<MolViewPtr> result;

    for (const auto &mol : mols)
    {
        SelectorBond b(*mol, map);

        if (not b.isEmpty())
            result.append(b);
    }

    return result;
}

SelectResult IDAllEngine::select(const SelectResult &mols,
                                 const PropertyMap &map) const
{
    switch (obj)
    {
    case AST::ATOM:
        return mols.atoms();
    case AST::RESIDUE:
        return mols.residues();
    case AST::CHAIN:
        return mols.chains();
    case AST::SEGMENT:
        return mols.segments();
    case AST::CUTGROUP:
        return mols.cutGroups();
    case AST::MOLECULE:
        return mols.molecules();
    case AST::BOND:
        return _get_bonds(mols, map);
    default:
        return mols.molecules();
    }
}

SelectEngine::ObjType IDAllEngine::objectType() const
{
    return _to_obj_type(obj);
}

////////
//////// Implementation of the IDWaterEngine
////////

IDWaterEngine::IDWaterEngine()
{}

SelectEnginePtr IDWaterEngine::construct()
{
    IDWaterEngine *ptr = new IDWaterEngine();
    auto p = makePtr(ptr);

    return p;
}

IDWaterEngine::~IDWaterEngine()
{}

SelectResult IDWaterEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    QList<MolViewPtr> result;

    for (const auto &molview : mols)
    {
        // Counters for the number of hydrogens, oxygens, and protons in the molecule.
        int num_hydrogen = 0;
        int num_oxygen = 0;
        int num_protons = 0;

        // Convert to a molecule.
        auto molecule = molview->molecule();

        // Skip if there is no element property.
        if (not molecule.hasProperty(map["element"]))
            continue;

        // Extract the element property.
        const auto &elements = molecule.property(map["element"]).asA<AtomElements>();

        // Whether this a water molecule.
        bool is_water = true;

        // Loop over all cut-groups associated with the elements.
        for (int i=0; i<elements.nCutGroups(); ++i)
        {
            // Create the cut-group index.
            CGIdx cg(i);

            // Extract the data for this cut-group.
            auto data = elements.constData(cg);

            // Loop over all atoms in this cut-group.
            for (int j=0; j<elements.nAtoms(cg); ++j)
            {
                // Get the element.
                const auto element = data[j];

                // Update the number of protons.
                num_protons += element.nProtons();

                // Hydrogen.
                if (element.nProtons() == 1)
                    num_hydrogen++;
                // Oxygen.
                else if (element.nProtons() == 8)
                    num_oxygen++;

                // Not a water molecule, abort!
                if (num_oxygen > 1 or
                    num_hydrogen > 2 or
                    num_protons > 10)
                {
                    is_water = false;
                    break;
                }
            }

            // Break out of inner loop.
            if (not is_water)
                break;
        }

        // If this is a water molecule, append the result.
        if (is_water and num_oxygen == 1 and num_hydrogen == 2 and num_protons == 10)
            result.append(molecule);
    }

    return SelectResult(result);
}

SelectEngine::ObjType IDWaterEngine::objectType() const
{
    return SelectEngine::MOLECULE;
}


////////
//////// Implementation of the IDProteinEngine
////////

QString IDProtein::toString() const
{
    return QString("protein");
}

SelectEnginePtr IDProtein::toEngine() const
{
    return IDProteinEngine::construct();
}

IDProteinEngine::IDProteinEngine()
{}

SelectEnginePtr IDProteinEngine::construct()
{
    IDProteinEngine *ptr = new IDProteinEngine();
    auto p = makePtr(ptr);

    return p;
}

IDProteinEngine::~IDProteinEngine()
{}


SelectResult IDProteinEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    QList<MolViewPtr> result;

    auto min_res = get_min_protein_residues();
    auto resnames = get_protein_residue_names();

    for (const auto &mol : mols)
    {
        auto molinfo = mol->data().info();

        if (molinfo.nResidues() < min_res)
            continue;

        int nres = 0;
        bool is_protein = false;

        for (int i=0; i<molinfo.nResidues(); ++i)
        {
            auto name = molinfo.name(ResIdx(i)).value().toLower();

            if (resnames.contains(name))
            {
                nres += 1;

                if (nres >= min_res)
                {
                    is_protein = true;
                    break;
                }
            }
        }

        if (is_protein)
            result.append(mol->molecule());
    }

    return SelectResult(result);
}

SelectEngine::ObjType IDProteinEngine::objectType() const
{
    return SelectEngine::MOLECULE;
}

////////
//////// Implementation of the IDPerturbableEngine
////////

IDPerturbableEngine::IDPerturbableEngine()
{}

SelectEnginePtr IDPerturbableEngine::construct()
{
    IDPerturbableEngine *ptr = new IDPerturbableEngine();
    auto p = makePtr(ptr);

    return p;
}

IDPerturbableEngine::~IDPerturbableEngine()
{}

SelectResult IDPerturbableEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    QList<MolViewPtr> result;

    for (const auto &molview : mols)
    {
        // Convert to a molecule.
        auto molecule = molview->molecule();

        // Check whether this molecule is flagged as being perturbable.
        if (molecule.hasProperty(map["is_perturbable"]))
            result.append(molecule);
    }

    return SelectResult(result);
}

SelectEngine::ObjType IDPerturbableEngine::objectType() const
{
    return SelectEngine::MOLECULE;
}

////////
//////// Implementation of the IDCountEngine
////////

IDCountEngine::IDCountEngine()
{}

SelectEnginePtr IDCount::toEngine() const
{
    return IDCountEngine::construct(object.toEngine(), compare, value);
}

QString IDCount::toString() const
{
    return QObject::tr("count( %1 ) %2 %3")
            .arg(object.toString())
            .arg(idcomparison_to_string(compare))
            .arg(value);
}

SelectEnginePtr IDCountEngine::construct(SelectEnginePtr object,
                                         IDComparison compare, int value)
{
    IDCountEngine *ptr = new IDCountEngine();
    auto p = makePtr(ptr);

    if (object)
        object->setParent(p);

    ptr->object = object;
    ptr->compare = compare;
    ptr->value = value;

    return p;
}

IDCountEngine::~IDCountEngine()
{}

std::function<int (const MoleculeView&)> _get_count(SelectEngine::ObjType obj)
{
    switch(obj)
    {
    case SelectEngine::ATOM:
        return [](const MoleculeView &view){ return view.nAtoms(); };
    case SelectEngine::RESIDUE:
        return [](const MoleculeView &view){ return view.nResidues(); };
    case SelectEngine::CHAIN:
        return [](const MoleculeView &view){ return view.nChains(); };
    case SelectEngine::SEGMENT:
        return [](const MoleculeView &view){ return view.nSegments(); };
    case SelectEngine::CUTGROUP:
        return [](const MoleculeView &view){ return view.nCutGroups(); };
    case SelectEngine::BOND:
        return [](const MoleculeView &view){ return SelectorBond(view).count(); };
    default:
        return [](const MoleculeView&){ return 0; };
    }
}

SelectResult IDCountEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    QList<MolViewPtr> result;

    // first select our expression from mols
    if (not object)
        return result;

    auto compare_func = _get_compare(compare);
    auto count_func = _get_count(object->objectType());

    auto items = object->operator()(mols, map);

    if (items.isEmpty())
        return SelectResult();

    for (const auto &molview : items.toList())
    {
        if (compare_func(count_func(*molview), value))
        {
            result.append(molview);
        }
    }

    return SelectResult(result);
}

SelectEngine::ObjType IDCountEngine::objectType() const
{
    if (object)
        return object->objectType();
    else
        return SelectEngine::COMPLEX;
}

SelectEnginePtr IDCountEngine::simplify()
{
    if (object.get())
        object = object->simplify();

    return selfptr.lock();
}
