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

#include "viewsofmol.h"
#include "partialmolecule.h"

#include "editor.hpp"
#include "mover.hpp"
#include "selector.hpp"

#include "molecule.h"
#include "segment.h"
#include "chain.h"
#include "residue.h"
#include "cutgroup.h"
#include "atom.h"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <QDebug>

using namespace SireMol;
using namespace SireStream;

RegisterMetaType<ViewsOfMol> r_molviews;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                       const ViewsOfMol &molviews)
{
    writeHeader(ds, r_molviews, 1);

    SharedDataStream sds(ds);
    
    sds << molviews.selected_atoms << molviews.views
        << static_cast<const MoleculeView&>(molviews);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, 
                                       ViewsOfMol &molviews)
{
    VersionID v = readHeader(ds, r_molviews);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> molviews.selected_atoms >> molviews.views
            >> static_cast<MoleculeView&>(molviews);
    }
    else
        throw version_error(v, "1", r_molviews, CODELOC);

    return ds;
}

/** Null constructor */
ViewsOfMol::ViewsOfMol() : ConcreteProperty<ViewsOfMol,MoleculeView>()
{}

/** Construct an empty view of the molecule whose data
    is in 'moldata' */
ViewsOfMol::ViewsOfMol(const MoleculeData &moldata)
           : ConcreteProperty<ViewsOfMol,MoleculeView>(moldata), 
             selected_atoms(moldata)
{}

/** Construct the view of the passed molecule */
ViewsOfMol::ViewsOfMol(const MoleculeData &moldata,
                       const AtomSelection &molview)
           : ConcreteProperty<ViewsOfMol,MoleculeView>(moldata), 
             selected_atoms(molview)
{
    selected_atoms.assertCompatibleWith(moldata);
}

/** Construct the views of the passed molecule */
ViewsOfMol::ViewsOfMol(const MoleculeData &moldata,
                       const QList<AtomSelection> &molviews)
           : ConcreteProperty<ViewsOfMol,MoleculeView>(moldata)
{
    if (molviews.isEmpty())
    {
        selected_atoms = AtomSelection(moldata);
    }
    else if (molviews.count() == 1)
    {
        selected_atoms = molviews.first();
        selected_atoms.assertCompatibleWith(moldata);
    }
    else
    {
        views = molviews;
        views.at(0).assertCompatibleWith(moldata);
        
        selected_atoms = views.at(0);
        selected_atoms.unite(views);
    }
}

/** Construct just a single view of a molecule */
ViewsOfMol::ViewsOfMol(const MoleculeView &view)
           : ConcreteProperty<ViewsOfMol,MoleculeView>(view), 
             selected_atoms(view.selection())
{}

/** Set this equal to the multiple views held in 'selection' */
template<class T>
void ViewsOfMol::setEqualTo(const Selector<T> &selection)
{
    MoleculeView::operator=(selection);
    selected_atoms = AtomSelection(this->data());
    views.clear();

    int nviews = selection.count();
    
    if (nviews == 1)
    { 
        selected_atoms = selection(0).selection();
    }
    else if (nviews > 1)
    {
        for (int i=0; i<nviews; ++i)
        {
            views.append( selection(i).selection() );
        }
        
        selected_atoms = views.at(0);
        selected_atoms.unite(views);
    }
}

/** Construct from the set of atoms */
ViewsOfMol::ViewsOfMol(const Selector<Atom> &atoms)
           : ConcreteProperty<ViewsOfMol,MoleculeView>()
{
    this->setEqualTo(atoms);
}

/** Construct from the set of CutGroups */
ViewsOfMol::ViewsOfMol(const Selector<CutGroup> &cgroups)
           : ConcreteProperty<ViewsOfMol,MoleculeView>()
{
    this->setEqualTo(cgroups);
}

/** Construct from the set of residues */
ViewsOfMol::ViewsOfMol(const Selector<Residue> &residues)
           : ConcreteProperty<ViewsOfMol,MoleculeView>()
{
    this->setEqualTo(residues);
}

/** Construct from the set of chains */
ViewsOfMol::ViewsOfMol(const Selector<Chain> &chains)
           : ConcreteProperty<ViewsOfMol,MoleculeView>()
{
    this->setEqualTo(chains);
}

/** Construct from the set of segments */
ViewsOfMol::ViewsOfMol(const Selector<Segment> &segments)
           : ConcreteProperty<ViewsOfMol,MoleculeView>()
{
    this->setEqualTo(segments);
}

/** Copy constructor */
ViewsOfMol::ViewsOfMol(const ViewsOfMol &other)
           : ConcreteProperty<ViewsOfMol,MoleculeView>(other), 
             selected_atoms(other.selected_atoms),
             views(other.views)
{}

/** Destructor */
ViewsOfMol::~ViewsOfMol()
{}

/** Copy assignment operator */
ViewsOfMol& ViewsOfMol::operator=(const ViewsOfMol &other)
{
    MoleculeView::operator=(other);
    selected_atoms = other.selected_atoms;
    views = other.views;
    
    return *this;
}

/** Copy assignment operator */
ViewsOfMol& ViewsOfMol::operator=(const MoleculeView &view)
{
    MoleculeView::operator=(view);
    selected_atoms = view.selection();
    views.clear();
    
    return *this;
}

/** Copy assignment operator from a set of atoms */
ViewsOfMol& ViewsOfMol::operator=(const Selector<Atom> &atoms)
{
    this->setEqualTo(atoms);
    return *this;
}

/** Copy assignment operator from a set of CutGroups */
ViewsOfMol& ViewsOfMol::operator=(const Selector<CutGroup> &cgroups)
{
    this->setEqualTo(cgroups);
    return *this;
}

/** Copy assignment operator from a set of residues */
ViewsOfMol& ViewsOfMol::operator=(const Selector<Residue> &residues)
{
    this->setEqualTo(residues);
    return *this;
}

/** Copy assignment operator from a set of chains */
ViewsOfMol& ViewsOfMol::operator=(const Selector<Chain> &chains)
{
    this->setEqualTo(chains);
    return *this;
}

/** Copy assignment operator from a set of segments */
ViewsOfMol& ViewsOfMol::operator=(const Selector<Segment> &segments)
{
    this->setEqualTo(segments);
    return *this;
}

/** Comparison operator */
bool ViewsOfMol::operator==(const ViewsOfMol &other) const
{
    return MoleculeView::operator==(other) and
           selected_atoms == other.selected_atoms and
           views == other.views;
}

/** Comparison operator */
bool ViewsOfMol::operator!=(const ViewsOfMol &other) const
{
    return MoleculeView::operator!=(other) or
           selected_atoms != other.selected_atoms or
           views != other.views;
}

/** Return a string representation of these views */
QString ViewsOfMol::toString() const
{
    return QObject::tr( "ViewsOfMol( %1 : %2 : nViews() == %3 )" )
                .arg( this->name() )
                .arg( this->number() )
                .arg( this->nViews() );
}

/** Return whether or not this is empty */
bool ViewsOfMol::isEmpty() const
{
    return selected_atoms.selectedNone();
}

/** Return whether or not this contains all atoms in the molecule */
bool ViewsOfMol::selectedAll() const
{
    return selected_atoms.selectedAll();
}

/** Return the number of views in this set */
int ViewsOfMol::nViews() const
{
    if (views.isEmpty())
    {
        if (selected_atoms.selectedNone())
            return 0;
        else
            return 1;
    }
    else
        return views.count();
}

/** Return the ith view in this set

    \throw SireError::invalid_index
*/
MolViewPtr ViewsOfMol::operator[](int i) const
{
    i = Index(i).map( this->nViews() );
    
    if ( i == 0 and views.isEmpty() )
        return PartialMolecule(*d, selected_atoms).toUnit();
    else
        return PartialMolecule(*d, views.at(i)).toUnit();
}

/** Return the name of the molecule being viewed */
const MolName& ViewsOfMol::name() const
{
    return d->name();
}

/** Return the number of the molecule being viewed */
MolNum ViewsOfMol::number() const
{
    return d->number();
}

/** Return the version of the molecule being viewed */
quint64 ViewsOfMol::version() const
{
    return d->version();
}

/** Return the version of the property at key 'key'

    \throw SireBase::missing_property
*/
quint64 ViewsOfMol::version(const PropertyName &key) const
{
    return d->version(key);
}

/** Synonym for ViewsOfMol::selection(i)

    \throw SireError::invalid_index
*/
const AtomSelection& ViewsOfMol::viewAt(int i) const
{
    return this->selection(i);
}

PartialMolecule ViewsOfMol::valueAt(int i) const
{
    return PartialMolecule(*this, this->viewAt(i));
}

/** Add the view 'view' to this set - this adds the 
    view even if it already exists in this set

    \throw SireError::incompatible_error
*/
void ViewsOfMol::add(const AtomSelection &view)
{
    int nviews = this->nViews();
    
    if (nviews == 0)
        this->operator=(ViewsOfMol(*d, view));
    else if (nviews == 1)
    {
        view.assertCompatibleWith(*d);
    
        views.append(selected_atoms);
        views.append(view);
        
        if (not selected_atoms.selectedAll())
            selected_atoms = selected_atoms.unite(view);
    }
    else
    {
        view.assertCompatibleWith(*d);
        
        views.append(view);
        
        if (not selected_atoms.selectedAll())
            selected_atoms = selected_atoms.unite(view);
    }
}

/** Return the views where 'views' have been added to 
    the set - this duplicates any views that already
    exist
    
    \throw SireError::incompatible_error
*/
void ViewsOfMol::add(const QList<AtomSelection> &views)
{
    if (views.isEmpty())
        return;
    else if (views.count() == 1)
    {
        this->add(views.first());
        return;
    }

    ViewsOfMol new_views(*this);
    
    foreach (const AtomSelection &view, views)
    {
        new_views.add(view);
    }
    
    this->operator=(new_views);
}

/** Add the view 'view' to this set, only if it doesn't
    already exist - this returns whether the view has
    been added
    
    \throw SireError::incompatible_error
*/
bool ViewsOfMol::addIfUnique(const AtomSelection &view)
{
    if (this->contains(view))
        return false;
    else
    {
        this->add(view);
        return true;
    }
}

/** Add the views in 'views' that don't already exist
    in this set - this returns the views that were added
    (or an empty list if nothing was added)
    
    \throw SireError::incompatible_error
*/
QList<AtomSelection> ViewsOfMol::addIfUnique(const QList<AtomSelection> &views)
{
    if (views.isEmpty())
        return views;
    else if (views.count() == 1)
    {
        if (this->addIfUnique(views.first()))
            return views;
        else
            return QList<AtomSelection>();
    }
    
    QList<AtomSelection> added_views;
    ViewsOfMol new_views(*this);
    
    foreach (const AtomSelection &view, views)
    {
        if (new_views.addIfUnique(view))
            added_views.append(view);
    }
    
    this->operator=(new_views);
    
    return added_views;
}

/** Synonym for ViewsOfMol::addIfUnique(view)
      
    \throw SireError::incompatible_error
*/
bool ViewsOfMol::unite(const AtomSelection &view)
{
    return this->addIfUnique(view);
}

/** Synonym for ViewsOfMol:addIfUnique(views) 
      
    \throw SireError::incompatible_error
*/
QList<AtomSelection> ViewsOfMol::unite(const QList<AtomSelection> &views)
{
    return this->addIfUnique(views);
}

/** Remove the ith view from this set - this returns
    the view that was removed

    \throw SireError::invalid_index
*/
AtomSelection ViewsOfMol::removeAt(int i)
{
    i = Index(i).map( this->nViews() );

    AtomSelection removed_view;

    if (views.isEmpty())
    {
        //there is only one view in this set - and it is
        //the total selection
        removed_view = selected_atoms;
        this->removeAll();
    }
    else
    {
        removed_view = views.takeAt(i);
        
        if (views.count() < 2)
        {
            selected_atoms = views.first();
            views.clear();
        }
        else
        {
            selected_atoms = views.at(0);
            selected_atoms.unite(views);
        }
    }

    return removed_view;
}

/** Remove the view 'view' from this set, if 
    any copies exist. This only removes the first
    copy of this view from this set, if multiple
    copies of this view exist. This returns whether 
    any copies were removed from this set.
    
    \throw SireError::incompatible_error
*/
bool ViewsOfMol::remove(const AtomSelection &view)
{
    if (not selected_atoms.contains(view))
        return false;
        
    if (views.count() == 0)
    {
        if (selected_atoms == view)
        {
            this->removeAll();
            return true;
        }
        else
            return false;
    }
    else
    {
        int i = views.indexOf(view);
        
        if (i >= 0)
        {
            views.removeAt(i);
            
            if (views.count() < 2)
            {
                selected_atoms = views.first();
                views.clear();
            }
            else
            {
                selected_atoms = views.at(0);
                selected_atoms.unite(views);
            }
            
            return true;
        }
        else
            return false;
    }
}

/** Remove the views in 'views' from this set. This removes only
    the first copy of the view if multiple exist in this set.
    This returns the views that were successfully removed. */
QList<AtomSelection> ViewsOfMol::remove(const QList<AtomSelection> &views)
{
    if (views.isEmpty())
        return views;
    else if (views.count() == 1)
    {
        if (this->remove(views.first()))
            return views;
        else
            return QList<AtomSelection>();
    }

    ViewsOfMol new_views(*this);
    QList<AtomSelection> removed_views;
    
    foreach (const AtomSelection &view, views)
    {
        if (new_views.remove(view))
            removed_views.append(view);
    }
    
    this->operator=(new_views);
    return removed_views;
}

/** Remove all copies of 'view' from this set, if any
    copies of this view are contained in this set. This
    return whether any views were removed

    \throw SireError::incompatible_error
*/
bool ViewsOfMol::removeAll(const AtomSelection &view)
{
    if (not selected_atoms.contains(view))
        return false;
    
    if (views.isEmpty())
    {
        if (selected_atoms == view)
        {
            this->removeAll();
            return true;
        }
        else
            return false;
    }
    else
    {
        if (views.removeAll(view) > 0)
        {
            if (views.isEmpty())
                this->removeAll();
            
            else if (views.count() == 1)
            {
                selected_atoms = views.first();
                views.clear();
            }
            else
            {
                selected_atoms = views.first();
                selected_atoms.unite(views);
            }
            
            return true;
        }
        else
            return false;
    }
}

/** Remove all copies of all of the views in 'views'. This removes
    all copies of any duplicated views in this set, and returns
    a list of all of the views that were successfully removed. */
QList<AtomSelection> ViewsOfMol::removeAll(const QList<AtomSelection> &views)
{
    if (views.isEmpty())
        return views;
    else if (views.count() == 1)
    {
        if (this->removeAll(views.first()))
            return views;
        else
            return QList<AtomSelection>();
    }
    
    ViewsOfMol new_views;
    QList<AtomSelection> removed_views;
    
    foreach (const AtomSelection &view, views)
    {
        if (new_views.removeAll(view))
            removed_views.append(view);
    }
    
    this->operator=(new_views);
    
    return removed_views;
}

/** Remove all duplicate views from this set - this returns
    copies of all of the removed views (multiple times
    if a view is contained more than twice in this set) */
QList<AtomSelection> ViewsOfMol::removeDuplicates()
{
    if (views.isEmpty())
        return views;
        
    QList<AtomSelection> removed_views;

    int nviews = views.count();
    
    //compare all pairs of views
    int i = 0;
    
    while (i < nviews-1)
    {
        //need to take a pointer as removing views would
        //screw up the reference
        const AtomSelection *view0 = &(views.at(i));
        
        int j = i+1;
        
        while (j < nviews)
        {
            const AtomSelection &view1 = views.at(j);
            
            if (*view0 == view1)
            {
                //this is a duplicate view - remove the later view
                removed_views.append(view1);
                views.removeAt(j);
                --nviews;
                
                //need to retake the pointer to view0 as removing the
                //item may have moved view0 in memory
                view0 = &(views.at(i));
            }
            else
                ++j;
        }
        
        ++i;
    }
    
    return removed_views;
}

/** Remove all views from this set */
void ViewsOfMol::removeAll()
{
    views.clear();
    selected_atoms = selected_atoms.selectNone();
}

/** Adds the other views onto this set.

    \throw SireError::incompatible_error
*/
ViewsOfMol ViewsOfMol::operator+(const ViewsOfMol &other) const
{
    ViewsOfMol ret(*this);
    ret.add(other);
    return ret;
}

/** Subtract the other views from this set

    \throw SireError::incompatible_error
*/
ViewsOfMol ViewsOfMol::operator-(const ViewsOfMol &other) const
{
    ViewsOfMol ret(*this);
    ret.remove(other);
    return ret;
}

/** Synonym for ViewsOfMol::add(views), except that this
    set is returned
    
    \throw SireError::incompatible_error
*/
ViewsOfMol& ViewsOfMol::operator+=(const ViewsOfMol &other)
{
    this->add(other);
    return *this;
}

/** Synonym for ViewsOfMol::remove(views), except that this
    set is returned 
    
    \throw SireError::incompatible_error
*/
ViewsOfMol& ViewsOfMol::operator-=(const ViewsOfMol &other)
{
    this->remove(other);
    return *this;
}

/** Return the molecule that is the union of all of the
    views in this set */
PartialMolecule ViewsOfMol::join() const
{
    return PartialMolecule(*d, selected_atoms);
}

/** Return the molecule that is the union of all of the
    views in this set */
PartialMolecule ViewsOfMol::all() const
{
    return this->join();
}

/** Return the molecule that contains these views */
Molecule ViewsOfMol::molecule() const
{
    return Molecule(*d);
}

/** Return the mover that can move all of the atoms 
    in all of the views */
Mover<ViewsOfMol> ViewsOfMol::move() const
{
    return Mover<ViewsOfMol>(*this);
}

/** Return all of the atoms selected across all of the 
    views in this set */
AtomSelection ViewsOfMol::selection() const
{
    return selected_atoms;
}

/** Return the atoms selected in the ith view in this set
    
    \throw SireError:invalid_index
*/
const AtomSelection& ViewsOfMol::selection(int i) const
{
    i = Index(i).map( this->nViews() );
    
    if (views.isEmpty())
        return selected_atoms;
    else
        return views.at(i);
}

/** Return all of the selections in this set of views */
QList<AtomSelection> ViewsOfMol::selections() const
{
    if (views.isEmpty())
    {
        if (selected_atoms.isEmpty())
            return QList<AtomSelection>();
        else
        {
            QList<AtomSelection> allviews;
            allviews.append(selected_atoms);
            return allviews;
        }
    }
    else
        return views;
}

/** Return the mover that can move all of the atoms
    in the ith view in this set 
    
    \throw SireError::invalid_index
*/
Mover<ViewsOfMol> ViewsOfMol::move(int i) const
{
    return Mover<ViewsOfMol>(*this, selection(i));
}

/** Return an evaluator that can evaluate properties
    over all of the atoms in the views */
Evaluator ViewsOfMol::evaluate() const
{
    return Evaluator(*this);
}

/** Return an evaluator that can evaluate properties
    over the atoms in the ith view of this set 
    
    \throw SireError::invalid_index
*/
Evaluator ViewsOfMol::evaluate(int i) const
{
    return Evaluator(*this, selection(i));
}

/** Return whether or not any of the views contains
    the atom at index 'atomidx' 
    
    \throw SireError::invalid_index
*/
bool ViewsOfMol::contains(AtomIdx atomidx) const
{
    return selected_atoms.contains(atomidx);
}

/** Return whether or not the views between them contain
    all of the atoms identified by 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
bool ViewsOfMol::contains(const AtomID &atomid) const
{
    return selected_atoms.contains(atomid);
}

/** Return whether or not the views between them contain
    all of the atoms in the selection 'selection'
    
    \throw SireError::incompatible_error
*/
bool ViewsOfMol::contains(const AtomSelection &selection) const
{
    return selected_atoms.contains(selection);
}

/** Return whether or not the views between them contain
    all of the atoms in all of the selections in 'selections'
    
    \throw SireError::incompatible_error
*/
bool ViewsOfMol::contains(const QList<AtomSelection> &selections) const
{
    foreach (const AtomSelection &selection, selections)
    {
        if (not selected_atoms.contains(selection))
            return false;
    }
    
    return true;
}

/** Return whether or not any of the views contains any
    of the atoms identified by 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
bool ViewsOfMol::intersects(const AtomID &atomid) const
{
    return selected_atoms.intersects(atomid);
}

/** Return whether or not any of the views contains any
    of the atoms in 'selection'
    
    \throw SireError::incompatible_error
*/
bool ViewsOfMol::intersects(const AtomSelection &selection) const
{
    return selected_atoms.intersects(selection);
}

/** Return the index of the view 'selection' in this set, 
    searching forward from the index position 'from'. This
    returns -1 if this view is not present in this set
    
    \throw SireError::incompatible_error
*/
int ViewsOfMol::indexOf(const AtomSelection &selection, int from) const
{
    if (selection.isEmpty() or this->isEmpty())
        return -1;

    if (from != 0)
        from = Index(from).map(this->nViews());

    selected_atoms.assertCompatibleWith(selection);
    
    if (views.isEmpty())
    {
        return (selection == selected_atoms) - 1;
    }
    else
    {
        return views.indexOf(selection, from);
    }
}

/** Assert that none of the views contain the same atoms
    (i.e. there is no overlap between each view)
    
    \throw SireMol::duplicate_atom
*/
void ViewsOfMol::assertNoOverlap() const
{
    int nviews = this->nViews();

    if (nviews < 2)
        return;
        
    //compare each pair of views...
    for (int i=0; i<nviews-1; ++i)
    {
        const AtomSelection &view0 = views.at(i);
    
        for (int j=i+1; j<nviews; ++j)
        {
            const AtomSelection &view1 = views.at(j);
            
            if (view0.intersects(view1))
                //there are some common atoms...!
                throw SireMol::duplicate_atom( QObject::tr(
                    "The views (%1 and %2) of molecule %3 (%4) contain "
                    "overlapping atoms.")
                        .arg(i).arg(j)
                        .arg(this->name()).arg(this->number()),
                            CODELOC );
        }
    }
}

const char* ViewsOfMol::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ViewsOfMol>() );
}

namespace SireMol
{
    ////////
    //////// Explicitly instantiate templates
    ////////

    template class Mover<ViewsOfMol>;
}

ViewsOfMol* ViewsOfMol::clone() const
{
    return new ViewsOfMol(*this);
}
