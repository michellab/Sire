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

#include "SireBase/parallel.h"
#include "SireBase/booleanproperty.h"

#include "SireMol/atomelements.h"

#include "SireVol/space.h"
#include "SireVol/cartesian.h"

#include "tostring.h"

#include <QRegExp>

using namespace SireSearch;
using namespace SireMol;
using namespace SireBase;
using namespace parser_idengine;

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

/** Function used to find the matching parts in mols using 'searchFunction' */
template<class T>
SelectResult searchForMatch(const SelectResult &mols, const T &searchFunction, bool uses_parallel)
{
    //now loop through all of the molecules and find the matching atoms
    SelectResult::Container result;

    if (uses_parallel)
    {
        const auto molviews = mols.views();
        QVector<ViewsOfMol> tmpresult(molviews.count());

        tbb::parallel_for( tbb::blocked_range<int>(0,molviews.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                tmpresult[i] = searchFunction(molviews[i]);
            }
        });

        for (const auto &r : tmpresult)
        {
            if (not r.isEmpty())
                result.append(r);
        }
    }
    else
    {
        for (const auto &mol : mols.views())
        {
            auto match = searchFunction(mol);

            if (not match.isEmpty())
                result.append(match);
        }
    }

    return result;
}

/** Internal function used by nearly all of the ID selections to return a ViewsOfMol
    with the parts selected from 'mol' that match 'selectFunc'. The passed 'countFunc'
    and 'getSelectedFunc' are used to get the number of items from the molecule and
    see what is selected, while 'uses_parallel' sets whether or not parallel selection
    should be performed */
template<class IDX, class OBJ, class T, class U, class V>
ViewsOfMol genericSelect(const ViewsOfMol &mol, const T &selectFunc,
                         const U &countFunc, const V &getSelectedFunc,
                         bool uses_parallel)
{
    const auto molinfo = mol.data().info();

    QList<IDX> selected;

    if (mol.selectedAll())
    {
        const int count = countFunc(molinfo);

        if (uses_parallel)
        {
            QMutex mutex;

            tbb::parallel_for( tbb::blocked_range<int>(0,count),
                               [&](const tbb::blocked_range<int> &r)
            {
                QList<IDX> thread_selected;

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const IDX idx(i);

                    if (selectFunc(molinfo, idx))
                        thread_selected.append(idx);
                }

                QMutexLocker lkr(&mutex);
                selected += thread_selected;
            });
        }
        else
        {
            for (int i=0; i<count; ++i)
            {
                const IDX idx(i);

                if (selectFunc(molinfo, idx))
                    selected.append(idx);
            }
        }
    }
    else
    {
        const auto viewed = getSelectedFunc(mol);

        if (uses_parallel)
        {
            QMutex mutex;

            tbb::parallel_for( tbb::blocked_range<int>(0,viewed.count()),
                               [&](const tbb::blocked_range<int> &r)
            {
                QList<IDX> thread_selected;

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto idx = viewed[i];

                    if (selectFunc(molinfo, idx))
                        thread_selected.append(idx);
                }

                QMutexLocker lkr(&mutex);
                selected += thread_selected;
            });
        }
        else
        {
            for (const auto &idx : viewed)
            {
                if (selectFunc(molinfo, idx))
                    selected.append(idx);
            }
        }
    }

    if (selected.isEmpty())
        //nothing matched
        return ViewsOfMol();
    else if (selected.count() == 1)
        //only a single object matched
        return ViewsOfMol( OBJ(mol.data(),selected[0]) );
    else if (selected.count() == countFunc(molinfo))
        //the entire molecule matched
        return ViewsOfMol( mol.molecule() );
    else
    {
        //a subset of the molecule matches
        return ViewsOfMol( mol.data(),
                           Selector<OBJ>(mol.data(),selected).selection() );
    }
}

static auto countAtoms = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nAtoms();
};

static auto getSelectedAtoms = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedAtoms();
};

static auto countCutGroups = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nCutGroups();
};

static auto getSelectedCutGroups = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedCutGroups();
};

static auto countResidues = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nResidues();
};

static auto getSelectedResidues = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedResidues();
};

static auto countChains = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nChains();
};

static auto getSelectedChains = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedChains();
};

static auto countSegments = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nSegments();
};

static auto getSelectedSegments = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedSegments();
};

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

/** Function used to find all of the atoms that match by name */
SelectResult IDNameEngine::selectAtoms(const SelectResult &mols, bool uses_parallel) const
{
    auto matchAtom = [&](const MoleculeInfoData &molinfo, const AtomIdx &idx)
    {
        return this->match(molinfo.name(idx).value());
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<AtomIdx,Atom>(mol, matchAtom, countAtoms,
                                             getSelectedAtoms, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectCutGroups(const SelectResult &mols, bool uses_parallel) const
{
    auto matchCutGroup = [&](const MoleculeInfoData &molinfo, const CGIdx &idx)
    {
        return this->match(molinfo.name(idx).value());
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<CGIdx,CutGroup>(mol, matchCutGroup, countCutGroups,
                                               getSelectedCutGroups, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectResidues(const SelectResult &mols, bool uses_parallel) const
{
    auto matchResidue = [&](const MoleculeInfoData &molinfo, const ResIdx &idx)
    {
        return this->match(molinfo.name(idx).value());
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<ResIdx,Residue>(mol, matchResidue, countResidues,
                                               getSelectedResidues, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectChains(const SelectResult &mols, bool uses_parallel) const
{
    auto matchChain = [&](const MoleculeInfoData &molinfo, const ChainIdx &idx)
    {
        return this->match(molinfo.name(idx).value());
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<ChainIdx,Chain>(mol, matchChain, countChains,
                                               getSelectedChains, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectSegments(const SelectResult &mols, bool uses_parallel) const
{
    auto matchSegment = [&](const MoleculeInfoData &molinfo, const SegIdx &idx)
    {
        return this->match(molinfo.name(idx).value());
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<SegIdx,Segment>(mol, matchSegment, countSegments,
                                               getSelectedSegments, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectMolecules(const SelectResult &mols, bool uses_parallel) const
{
    // function that finds all of the molecules that have been selected
    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        const auto molname = mol.data().name().value();

        //try all of the fixed names
        for (const auto &name : names)
        {
            if (name == molname)
            {
                //name matches exactly
                return ViewsOfMol(mol.molecule());
            }
        }

        //now try all of the regexps
        for (const auto &regexp : regexps)
        {
            auto match = regexp.match(molname);

            if (match.hasMatch())
            {
                //we have a regexp match :-)
                return ViewsOfMol(mol.molecule());
            }
        }

        //no match
        return ViewsOfMol();
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;

    bool uses_parallel = true;

    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    switch(obj)
    {
    case AST::ATOM:
        return selectAtoms(mols, uses_parallel);
    case AST::CUTGROUP:
        return selectCutGroups(mols, uses_parallel);
    case AST::RESIDUE:
        return selectResidues(mols, uses_parallel);
    case AST::CHAIN:
        return selectChains(mols, uses_parallel);
    case AST::SEGMENT:
        return selectSegments(mols, uses_parallel);
    case AST::MOLECULE:
        return selectMolecules(mols, uses_parallel);
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDNameEngine::objectType() const
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
    default:
        return SelectEngine::COMPLEX;
    }
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

            if (v.step == 0)
                v.step = 1;

            if (v.step < 0)
                v.step *= -1;

            if (v.start <= v.end)
            {
                for (int i=v.start; i<=v.end; i+=v.step)
                {
                    if (i == idx)
                        return true;
                }
            }
            else
            {
                for (int i=v.start; i>=v.end; i-=v.step)
                {
                    if (i == idx)
                        return true;
                }
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

/** Function used to find all of the atoms that match by number */
SelectResult IDNumberEngine::selectAtoms(const SelectResult &mols, bool uses_parallel) const
{
    auto matchAtom = [&](const MoleculeInfoData &molinfo, const AtomIdx &idx)
    {
        const auto number = molinfo.number(idx).value();
        return match(number);
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<AtomIdx,Atom>(mol, matchAtom, countAtoms,
                                             getSelectedAtoms, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

/** Function used to find all of the residues that match by number */
SelectResult IDNumberEngine::selectResidues(const SelectResult &mols, bool uses_parallel) const
{
    auto matchResidue = [&](const MoleculeInfoData &molinfo, const ResIdx &idx)
    {
        const auto number = molinfo.number(idx).value();
        return match(number);
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<ResIdx,Residue>(mol, matchResidue, countResidues,
                                               getSelectedResidues, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

/** Function used to find all of the molecules that match by number */
SelectResult IDNumberEngine::selectMolecules(const SelectResult &mols, bool uses_parallel) const
{
    // function that finds all of the molecules that have been selected
    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        const auto molnum = mol.data().number().value();

        if (match(molnum))
            return ViewsOfMol(mol.molecule());
        else
            return ViewsOfMol();
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNumberEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;

    bool uses_parallel = true;

    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    switch(obj)
    {
    case AST::ATOM:
        return selectAtoms(mols, uses_parallel);
    case AST::RESIDUE:
        return selectResidues(mols, uses_parallel);
    case AST::MOLECULE:
        return selectMolecules(mols, uses_parallel);
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDNumberEngine::objectType() const
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
    default:
        return SelectEngine::COMPLEX;
    }
}

////////
//////// Implementation of the IDIndexEngine
////////

static int map(int idx, int n)
{
    if (idx >= 0)
        return idx;
    else if (idx < 0 and idx >= -n)
        return n + idx;
    else
        return -1;
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

            int start = map(v.start,count);
            int end = map(v.end,count);

            //only loop if the range is valid
            if (start < count and end < count and start >= 0 and end >= 0)
            {
                if (start <= end)
                {
                    for (int i=start; i<=end; i+=v.step)
                    {
                        if (i == idx)
                            return true;
                    }
                }
                else
                {
                    for (int i=start; i>=end; i-=v.step)
                    {
                        if (i == idx)
                            return true;
                    }
                }
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

SelectResult IDIndexEngine::selectAtoms(const SelectResult &mols, bool uses_parallel) const
{
    auto matchAtom = [&](const MoleculeInfoData &molinfo, const AtomIdx &idx)
    {
        return match(idx.value(), countAtoms(molinfo));
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<AtomIdx,Atom>(mol, matchAtom, countAtoms,
                                             getSelectedAtoms, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDIndexEngine::selectCutGroups(const SelectResult &mols, bool uses_parallel) const
{
    auto matchCutGroup = [&](const MoleculeInfoData &molinfo, const CGIdx &idx)
    {
        return match(idx.value(), countCutGroups(molinfo));
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<CGIdx,CutGroup>(mol, matchCutGroup, countCutGroups,
                                               getSelectedCutGroups, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDIndexEngine::selectResidues(const SelectResult &mols, bool uses_parallel) const
{
    auto matchResidue = [&](const MoleculeInfoData &molinfo, const ResIdx &idx)
    {
        return match(idx.value(), countResidues(molinfo));
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<ResIdx,Residue>(mol, matchResidue, countResidues,
                                               getSelectedResidues, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDIndexEngine::selectChains(const SelectResult &mols, bool uses_parallel) const
{
    auto matchChain = [&](const MoleculeInfoData &molinfo, const ChainIdx &idx)
    {
        return match(idx.value(), countChains(molinfo));
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<ChainIdx,Chain>(mol, matchChain, countChains,
                                               getSelectedChains, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDIndexEngine::selectSegments(const SelectResult &mols, bool uses_parallel) const
{
    auto matchSegment = [&](const MoleculeInfoData &molinfo, const SegIdx &idx)
    {
        return match(idx.value(), countSegments(molinfo));
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<SegIdx,Segment>(mol, matchSegment, countSegments,
                                               getSelectedSegments, uses_parallel);
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDIndexEngine::selectMolecules(const SelectResult &mols, bool uses_parallel) const
{
    //have to use the index order from 'mols'
    SelectResult::Container result;

    const auto m = mols.views();

    for (int i=0; i<m.count(); ++i)
    {
        if (match(i,m.count()))
            result.append( ViewsOfMol(m[i].molecule()) );
    }

    return result;
}

SelectResult IDIndexEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;

    bool uses_parallel = true;

    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    switch(obj)
    {
    case AST::ATOM:
        return selectAtoms(mols, uses_parallel);
    case AST::CUTGROUP:
        return selectCutGroups(mols, uses_parallel);
    case AST::RESIDUE:
        return selectResidues(mols, uses_parallel);
    case AST::CHAIN:
        return selectChains(mols, uses_parallel);
    case AST::SEGMENT:
        return selectSegments(mols, uses_parallel);
    case AST::MOLECULE:
        return selectMolecules(mols, uses_parallel);
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDIndexEngine::objectType() const
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
    default:
        return SelectEngine::COMPLEX;
    }
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

    SelectResult::Container result;

    bool uses_parallel = true;

    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    SelectResult result0, result1;

    //first get the results of both sub-matches
    if (uses_parallel)
    {
        tbb::parallel_invoke( [&](){ result0 = part0->operator()(mols, map); },
                              [&](){ result1 = part1->operator()(mols, map); } );
    }
    else
    {
        result0 = part0->operator()(mols, map);
        result1 = part1->operator()(mols, map);
    }

    //match up the results
    for (const auto &mol0 : result0.views())
    {
        const auto molnum = mol0.data().number();

        //is this molecule in the other match?
        for (const auto &mol1 : result1.views())
        {
            if (mol1.data().number() == molnum)
            {
                //yes - unite these two views
                auto selected_atoms = mol0.selection();
                selected_atoms = selected_atoms.intersect(mol1.selection());

                if (not selected_atoms.isEmpty())
                    result.append( ViewsOfMol(mol0.data(),selected_atoms) );
            }
        }
    }

    return result;
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

/** Simple internal function that expands a mol into the specified type of views */
static ViewsOfMol expand(const ViewsOfMol &mol, SelectEngine::ObjType obj)
{
    switch(obj)
    {
    case SelectEngine::ATOM:
        return ViewsOfMol(mol.atoms());
    case SelectEngine::CUTGROUP:
        return ViewsOfMol(mol.cutGroups());
    case SelectEngine::RESIDUE:
        return ViewsOfMol(mol.residues());
    case SelectEngine::CHAIN:
        return ViewsOfMol(mol.chains());
    case SelectEngine::SEGMENT:
        return ViewsOfMol(mol.segments());
    case SelectEngine::MOLECULE:
        return ViewsOfMol(mol.molecule());
    default:
        return mol;
    }
}

SelectResult IDOrEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    if (parts.isEmpty())
        return SelectResult();
    else if (parts.count() == 1)
        return parts[0]->operator()(mols, map);

    bool uses_parallel = true;

    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    QVector<SelectResult> results(parts.count());

    //first get the results of the sub-matches
    if (uses_parallel)
    {
        tbb::parallel_for( tbb::blocked_range<int>(0,parts.count(),1),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                auto r = parts.at(i)->operator()(mols, map);
                results[i] = r;
            }
        });
    }
    else
    {
        for (int i=0; i<parts.count(); ++i)
        {
            results[i] = parts.at(i)->operator()(mols, map);
        }
    }

    //see if all parts use the same object type
    bool same_type = true;

    for (int i=1; i<parts.count(); ++i)
    {
        if (parts.at(i)->objectType() != parts.at(0)->objectType())
        {
            same_type = false;
            break;
        }
    }

    //combine all of the results
    QMap<MolNum,ViewsOfMol> result;

    for (int i=0; i<results.count(); ++i)
    {
        for (const auto &mol : results[i].views())
        {
            const auto molnum = mol.data().number();

            auto it = result.find(molnum);

            if (it == result.end())
            {
                if (same_type)
                    result.insert(molnum, mol);
                else
                    result.insert(molnum, ::expand(mol,parts.at(i)->objectType()));
            }
            else
            {
                if (same_type)
                    it.value() = ViewsOfMol((it.value() + mol).join());
                else
                    it.value().add( ::expand(mol,parts.at(i)->objectType()).selections() );
            }
        }
    }

    return SelectResult( result.values() );
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
                return SelectEngine::COMPLEX;
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

    bool uses_parallel = true;

    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    //first, select the parts...
    auto selected = part->operator()(mols, map);

    //function that returns the parts of 'original' that are not in 'selected'
    auto invert = [&](const ViewsOfMol &original, const ViewsOfMol &selected)
    {
        if (selected.selectedAll())
            return ViewsOfMol();
        else if (selected.isEmpty())
            return original;

        auto inverted = original.selection() - selected.selection();

        return ViewsOfMol(original.data(),inverted);
    };

    //now go through and find from 'mols' what is not in 'selected'
    QList<ViewsOfMol> result;

    if (uses_parallel)
    {
        QVector<ViewsOfMol> resultmols(mols.count());
        const auto molviews = mols.views();

        tbb::parallel_for( tbb::blocked_range<int>(0,molviews.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                const auto mol = molviews.at(i);

                //find the molecule in the set of selected
                const auto selectmol = selected.views(mol.data().number());

                resultmols[i] = invert(mol, selectmol);
            }
        });

        for (const auto &mol : resultmols)
        {
            if (not mol.isEmpty())
                result.append(mol);
        }
    }
    else
    {
        for (const auto &mol : mols.views())
        {
            auto inverted = invert(mol, selected.views(mol.data().number()));

            if (not inverted.isEmpty())
                result.append(inverted);
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

    //do we need to do anything?
    bool needs_joining = false;

    for (const auto &molview : selected.views())
    {
        if (molview.nViews() > 1)
        {
            needs_joining = true;
            break;
        }
    }

    if (not needs_joining)
        return selected;

    QList<ViewsOfMol> result;

    for (const auto &molview : selected.views())
    {
        if (molview.nViews() > 1)
            result.append( molview.join() );
        else
            result.append( molview );
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
    auto selected = this->expand( part->operator()(mols, map) );

    //how many views are there?
    const int nviews = selected.count();

    //now get the range of views to return
    const int start = ::map(val.start,nviews);
    const int end = ::map(val.end,nviews);
    const int step = val.step;

    auto addView = [](const MoleculeView &view, QList<ViewsOfMol> &result)
    {
        const int molnum = view.data().number();

        for (auto mol : result)
        {
            if (mol.data().number() == molnum)
            {
                mol += view;
                return;
            }
        }

        result.append( ViewsOfMol(view) );
    };

    QList<ViewsOfMol> result;

    //only loop if the range is valid
    if (start < nviews and end < nviews and start >= 0 and end >= 0)
    {
        if (start <= end)
        {
            for (int i=start; i<=end; i+=step)
            {
                addView( selected[i], result );
            }
        }
        else
        {
            for (int i=start; i>=end; i-=step)
            {
                addView( selected[i], result );
            }
        }
    }

    //re-join the views if this is a simple view type
    if (this->objectType() != SelectEngine::COMPLEX)
    {
        for (auto &mol : result)
        {
            mol = mol.join();
        }
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

    if (from_value)
        from_value->setParent(p);

    if (to_value)
        to_value->setParent(p);

    ptr->from_token = from_token;
    ptr->from_value = from_value;
    ptr->to_token = to_token;
    ptr->to_value = to_value;

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

    auto result = from_value->operator()(mols, map);

    if (result.count() == 0)
        return SelectResult();

    if (result.count() > 1)
        throw SireMol::parse_error(QObject::tr(
            "Can only extract bonds from matches that return a single molecule."),
                CODELOC);

    auto mol = result[0];

    if (from_token == ID_BOND_WITHIN)
    {
        /*if (to_value.get() or to_token != ID_BOND_UNKNOWN)
            throw SireMol::parse_error(QObject::tr(
                "'bonds in value' cannot be followed by any other terms."),
                    CODELOC);*/
        auto result = SelectorBond(*mol, map);

        return SelectResult(result);
    }

    return SelectResult();
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

SelectEnginePtr IDWithEngine::construct( IDObject obj, IDToken token, SelectEnginePtr part)
{
    IDWithEngine *ptr = new IDWithEngine();
    auto p = makePtr(ptr);

    if (part)
        part->setParent(p);

    ptr->part = part;
    ptr->obj = obj;
    ptr->token = token;

    return p;
}

IDWithEngine::~IDWithEngine()
{}

SelectResult IDWithEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    if (not part.get())
        return SelectResult();

    //first, select the parts...
    auto selected = part->operator()(mols, map);

    const auto objtype = this->objectType();

    if (objtype == SelectEngine::COMPLEX)
        return selected;

    QList<ViewsOfMol> result;

    if (objtype == SelectEngine::ATOM or token == AST::ID_WITH)
    {
        //we just need to expand and then re-join the parts
        selected = this->expand(selected);

        for (const auto &mol : selected.views())
        {
            result.append( mol.join() );
        }
    }
    else if (token == AST::ID_IN)
    {
        // 'in' will not include partially selected parts, so only completely
        // selected molecules, chains, residues etc.
        for (auto mol : selected.views())
        {
            //get the original molecule
            auto selected = mols.views(mol.data().number()).selection();

            //is the whole molecule selected?
            if (selected.selectedAll())
            {
                //yes - by definition it will contain all parts
                result.append( mol.join() );
            }
            else
            {
                auto in_selected = selected;
                in_selected = in_selected.selectNone();

                if (objtype == SelectEngine::CUTGROUP)
                {
                    for (auto cgidx : mol.selection().selectedCutGroups())
                    {
                        if (selected.selectedAll(cgidx))
                            in_selected = in_selected.select(cgidx);
                    }
                }
                else if (objtype == SelectEngine::RESIDUE)
                {
                    for (auto residx : mol.selection().selectedResidues())
                    {
                        if (selected.selectedAll(residx))
                            in_selected = in_selected.select(residx);
                    }
                }
                else if (objtype == SelectEngine::CHAIN)
                {
                    for (auto chainidx : mol.selection().selectedChains())
                    {
                        if (selected.selectedAll(chainidx))
                            in_selected = in_selected.select(chainidx);
                    }
                }
                else if (objtype == SelectEngine::SEGMENT)
                {
                    for (auto segidx : mol.selection().selectedSegments())
                    {
                        if (selected.selectedAll(segidx))
                            in_selected = in_selected.select(segidx);
                    }
                }

                result.append( ViewsOfMol(mol.data(),in_selected) );
            }
        }
    }
    else
    {
        throw SireError::program_bug( QObject::tr(
                "Invalid 'with' token? %1").arg(AST::idtoken_to_string(token)), CODELOC );
    }

    return SelectResult(result);
}

SelectEnginePtr IDWithEngine::simplify()
{
    if (part.get())
        part = part->simplify();

    return selfptr.lock();
}

SelectEngine::ObjType IDWithEngine::objectType() const
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
    default:
        return SelectEngine::COMPLEX;
    }
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
    const auto element_property = map["element"];

    bool uses_parallel = true;

    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    QList<ViewsOfMol> result;

    //function that extracts all of the atoms with matching elements
    auto getAtomsWithElement = [&](const ViewsOfMol &mol)
    {
        const auto &moldata = mol.data();
        const auto &molinfo = moldata.info();

        if (not moldata.hasProperty(element_property))
            return ViewsOfMol();

        const auto &prop = moldata.property(element_property);

        if (not prop.isA<AtomElements>())
            return ViewsOfMol();

        const auto &elms = prop.asA<AtomElements>();

        auto selected = mol.selection();

        //deselect any atoms that are not the right element
        for (const auto &atomidx : selected.selectedAtoms())
        {
            const auto element = elms[ molinfo.cgAtomIdx(atomidx) ];

            if (not elements.contains(element))
            {
                selected = selected.deselect(atomidx);
            }
        }

        if (not selected.selectedNone())
        {
            return ViewsOfMol(moldata,selected);
        }
        else
            return ViewsOfMol();
    };

    if (uses_parallel)
    {
        const auto molviews = mols.views();
        QVector<ViewsOfMol> matches(molviews.count());

        tbb::parallel_for( tbb::blocked_range<int>(0,molviews.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                auto match = getAtomsWithElement(molviews.at(i));

                if (not match.isEmpty())
                    matches[i] = match;
            }
        });

        for (const auto &match : matches)
        {
            if (not match.isEmpty())
                result.append(match);
        }
    }
    else
    {
        for (const auto &mol : mols.views())
        {
            auto match = getAtomsWithElement(mol);

            if (not match.isEmpty())
                result.append(match);
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

SelectResult IDDistanceEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    //first, get the objects against where the distance is calculated
    if (part.get() == 0)
        return SelectResult();

    const auto refmols = part->operator()(mols, map);

    if (refmols.isEmpty())
        //nothing against which to compare
        return SelectResult();

    const auto coords_property = map["coordinates"];

    bool uses_parallel = true;

    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    SireVol::SpacePtr space = SireVol::Cartesian();

    if (map["space"].hasValue())
    {
        space = map["space"].value().asA<SireVol::Space>();
    }

    //function that gets the reference coordinates of the passed view
    auto getCoords = [&](const AtomSelection &selection, const MoleculeInfoData &molinfo,
                         const AtomCoords &coords)
    {
        Vector minval(0), maxval(0), center(0);

        auto atoms = selection.selectedAtoms();

        if (atoms.isEmpty())
            return Vector(0);

        minval = coords[ molinfo.cgAtomIdx(atoms.at(0)) ];
        maxval = minval;
        center = minval;

        //now find the maximum and minimum of all other coordinates
        for (int i=1; i<atoms.count(); ++i)
        {
            const auto c = coords[ molinfo.cgAtomIdx(atoms.at(i)) ];
            minval.setMin(c);
            maxval.setMax(c);
        }

        center = minval + 0.5*(maxval-minval);

        switch (typ)
        {
        case ID_COORD_MAX:
            return maxval;
        case ID_COORD_MIN:
            return minval;
        case ID_COORD_MAX_X:
            center.setX(maxval.x());
            return center;
        case ID_COORD_MAX_Y:
            center.setY(maxval.y());
            return center;
        case ID_COORD_MAX_Z:
            center.setZ(maxval.z());
            return center;
        case ID_COORD_MIN_X:
            center.setX(minval.x());
            return center;
        case ID_COORD_MIN_Y:
            center.setY(minval.y());
            return center;
        case ID_COORD_MIN_Z:
            center.setZ(minval.z());
            return center;
        default:
            return center;
        }
    };

    const double dist2 = distance*distance;

    //function that tests whether the passed point is within the
    //distance of the atoms in the passed molecule - this could
    //be made significantly more efficient if we use the rapid distance
    //calculation abilities of Space and CoordGroup better...
    auto isWithin = [&](const Vector &point, const ViewsOfMol &mol)
    {
        const auto &moldata = mol.data();

        if (not moldata.hasProperty(coords_property))
            return false;

        const auto &prop = moldata.property(coords_property);

        if (not prop.isA<AtomCoords>())
            return false;

        const auto &coords = prop.asA<AtomCoords>();
        const auto &molinfo = moldata.info();

        if (mol.selectedAll())
        {
            //loop over all coordinates
            for (int i=0; i<molinfo.nAtoms(); ++i)
            {
                if (space.read().calcDist2(point, coords[molinfo.cgAtomIdx(AtomIdx(i))]) < dist2)
                {
                    return true;
                }
            }
        }
        else
        {
            for (const auto &atom : mol.selection().selectedAtoms())
            {
                if (space.read().calcDist2(point, coords[molinfo.cgAtomIdx(atom)]) < dist2)
                {
                    return true;
                }
            }
        }

        return false;
    };

    //function that gets the views within the specified distance
    auto getWithin = [&](const ViewsOfMol &searchmol)
    {
        const auto &moldata = searchmol.data();
        const auto &molinfo = moldata.info();

        //get the atoms that are to be selected
        auto selected_atoms = searchmol.selection();
        selected_atoms = selected_atoms.selectNone();

        //now get the coordinates of the atoms
        if (not moldata.hasProperty(coords_property))
            return ViewsOfMol();

        const auto &prop = moldata.property(coords_property);

        if (not prop.isA<AtomCoords>())
            return ViewsOfMol();

        const auto &coords = prop.asA<AtomCoords>();

        //we need to loop over the expanded molecule
        ViewsOfMol mol = this->expandMol(searchmol);

        //loop over each view in 'mol'
        for (int i=0; i<mol.nViews(); ++i)
        {
            bool within_distance = false;

            if (typ == ID_COORD_CLOSEST)
            {
                //need to compare all atoms in this view...
                for (const auto &atom : mol.viewAt(i).selectedAtoms())
                {
                    const auto point = coords[ molinfo.cgAtomIdx(atom) ];

                    //now loop over all atoms in selected molecules and see if they
                    //are within the distance
                    for (const auto &refmol : refmols)
                    {
                        bool is_within = isWithin(point, refmol);

                        if (is_within)
                        {
                            within_distance = true;
                            break;
                        }
                    }

                    if (within_distance)
                        break;
                }
            }
            else
            {
                //get the comparison coordinates for this view
                const auto point = getCoords(mol.viewAt(i), molinfo, coords);

                //now loop over all atoms in selected molecules and see if they
                //are within the distance
                for (const auto &refmol : refmols)
                {
                    bool is_within = isWithin(point, refmol);

                    if (is_within)
                    {
                        within_distance = true;
                        break;
                    }
                }
            }

            if (within_distance)
            {
                selected_atoms = selected_atoms.select(mol.viewAt(i));

                if (selected_atoms.selectedAll())
                    break;
            }
        }

        if (selected_atoms.isEmpty())
            return ViewsOfMol();
        else
            return ViewsOfMol(moldata, selected_atoms);
    };

    QList<ViewsOfMol> result;

    if (uses_parallel)
    {
        const auto molviews = mols.views();
        QVector<ViewsOfMol> matched(molviews.count());

        tbb::parallel_for( tbb::blocked_range<int>(0,molviews.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                auto match = getWithin(molviews.at(i));

                if (not match.isEmpty())
                    matched[i] = match;
            }
        });

        for (const auto &match : matched)
        {
            if (not match.isEmpty())
                result.append(match);
        }
    }
    else
    {
        for (const auto &mol : mols.views())
        {
            auto match = getWithin(mol);

            if (not match.isEmpty())
                result.append(match);
        }
    }

    return SelectResult(result);
}

SelectEngine::ObjType IDDistanceEngine::objectType() const
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
    default:
        return SelectEngine::COMPLEX;
    }
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
    // Extract the x,y,z components of the position (implicitly converted to Angstrom).
    double x = position.x.value * position.x.unit;
    double y = position.y.value * position.y.unit;
    double z = position.z.value * position.z.unit;

    // Create a reference point using the vector components.
    Vector ref_point(x, y, z);

    const auto coords_property = map["coordinates"];

    bool uses_parallel = true;

    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    SireVol::SpacePtr space = SireVol::Cartesian();

    if (map["space"].hasValue())
    {
        space = map["space"].value().asA<SireVol::Space>();
    }

    //function that gets the reference coordinates of the passed view
    auto getCoords = [&](const AtomSelection &selection, const MoleculeInfoData &molinfo,
                         const AtomCoords &coords)
    {
        Vector minval(0), maxval(0), center(0);

        auto atoms = selection.selectedAtoms();

        if (atoms.isEmpty())
            return Vector(0);

        minval = coords[ molinfo.cgAtomIdx(atoms.at(0)) ];
        maxval = minval;
        center = minval;

        //now find the maximum and minimum of all other coordinates
        for (int i=1; i<atoms.count(); ++i)
        {
            const auto c = coords[ molinfo.cgAtomIdx(atoms.at(i)) ];
            minval.setMin(c);
            maxval.setMax(c);
        }

        center = minval + 0.5*(maxval-minval);

        switch (typ)
        {
        case ID_COORD_MAX:
            return maxval;
        case ID_COORD_MIN:
            return minval;
        case ID_COORD_MAX_X:
            center.setX(maxval.x());
            return center;
        case ID_COORD_MAX_Y:
            center.setY(maxval.y());
            return center;
        case ID_COORD_MAX_Z:
            center.setZ(maxval.z());
            return center;
        case ID_COORD_MIN_X:
            center.setX(minval.x());
            return center;
        case ID_COORD_MIN_Y:
            center.setY(minval.y());
            return center;
        case ID_COORD_MIN_Z:
            center.setZ(minval.z());
            return center;
        default:
            return center;
        }
    };

    const double dist2 = distance*distance;

    //function that gets the views within the specified distance
    auto getWithin = [&](const ViewsOfMol &searchmol)
    {
        const auto &moldata = searchmol.data();
        const auto &molinfo = moldata.info();

        //get the atoms that are to be selected
        auto selected_atoms = searchmol.selection();
        selected_atoms = selected_atoms.selectNone();

        //now get the coordinates of the atoms
        if (not moldata.hasProperty(coords_property))
            return ViewsOfMol();

        const auto &prop = moldata.property(coords_property);

        if (not prop.isA<AtomCoords>())
            return ViewsOfMol();

        const auto &coords = prop.asA<AtomCoords>();

        //we need to loop over the expanded molecule
        ViewsOfMol mol = this->expandMol(searchmol);

        //loop over each view in 'mol'
        for (int i=0; i<mol.nViews(); ++i)
        {
            bool within_distance = false;

            if (typ == ID_COORD_CLOSEST)
            {
                //need to compare all atoms in this view...
                for (const auto &atom : mol.viewAt(i).selectedAtoms())
                {
                    const auto point = coords[ molinfo.cgAtomIdx(atom) ];

                    if (space.read().calcDist2(point, ref_point) < dist2)
                    {
                        within_distance = true;
                        break;
                    }
                }
            }
            else
            {
                //get the comparison coordinates for this view
                const auto point = getCoords(mol.viewAt(i), molinfo, coords);

                if (space.read().calcDist2(point, ref_point) < dist2)
                {
                    within_distance = true;
                    break;
                }
            }

            if (within_distance)
            {
                selected_atoms = selected_atoms.select(mol.viewAt(i));

                if (selected_atoms.selectedAll())
                    break;
            }
        }

        if (selected_atoms.isEmpty())
            return ViewsOfMol();
        else
            return ViewsOfMol(moldata, selected_atoms);
    };

    QList<ViewsOfMol> result;

    if (uses_parallel)
    {
        const auto molviews = mols.views();
        QVector<ViewsOfMol> matched(molviews.count());

        tbb::parallel_for( tbb::blocked_range<int>(0,molviews.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                auto match = getWithin(molviews.at(i));

                if (not match.isEmpty())
                    matched[i] = match;
            }
        });

        for (const auto &match : matched)
        {
            if (not match.isEmpty())
                result.append(match);
        }
    }
    else
    {
        for (const auto &mol : mols.views())
        {
            auto match = getWithin(mol);

            if (not match.isEmpty())
                result.append(match);
        }
    }

    return SelectResult(result);
}

SelectEngine::ObjType IDDistanceVectorEngine::objectType() const
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
    default:
        return SelectEngine::COMPLEX;
    }
}

////////
//////// Implementation of the IDAllEngine
////////

IDAllEngine::IDAllEngine()
{}

SelectEnginePtr IDAllEngine::construct()
{
    IDAllEngine *ptr = new IDAllEngine();
    auto p = makePtr(ptr);

    return p;
}

IDAllEngine::~IDAllEngine()
{}

SelectResult IDAllEngine::select(const SelectResult &mols, const PropertyMap&) const
{
    return mols;
}

SelectEngine::ObjType IDAllEngine::objectType() const
{
    return SelectEngine::COMPLEX;
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
    QList<ViewsOfMol> result;

    for (const auto &molview : mols.views())
    {
        // Counters for the number of hydrogens, oxygens, and protons in the molecule.
        int num_hydrogen = 0;
        int num_oxygen = 0;
        int num_protons = 0;

        // Convert to a molecule.
        auto molecule = molview.molecule();

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
                    num_protons > 10 or
                    num_protons == 0)
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
        if (is_water)
            result.append(molecule);
    }

    return SelectResult(result);
}

SelectEngine::ObjType IDWaterEngine::objectType() const
{
    return SelectEngine::COMPLEX;
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
    QList<ViewsOfMol> result;

    for (const auto &molview : mols.views())
    {
        // Convert to a molecule.
        auto molecule = molview.molecule();

        // Check whether this molecule is flagged as being perturbable.
        if (molecule.hasProperty(map["is_perturbable"]))
            result.append(molecule);
    }

    return SelectResult(result);
}

SelectEngine::ObjType IDPerturbableEngine::objectType() const
{
    return SelectEngine::COMPLEX;
}
