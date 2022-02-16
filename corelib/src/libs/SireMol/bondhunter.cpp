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

#include <QMutex>

#include "bondhunter.h"

#include "connectivity.h"
#include "moleculeview.h"
#include "atomselection.h"

#include "atomelements.h"
#include "atomcoords.h"
#include "atom.h"
#include "molecule.h"
#include "moleculedata.h"
#include "moleculeinfodata.h"

#include "mover.hpp"
#include "selector.hpp"

#include "SireVol/coordgroup.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireVol;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implementation of BondHunter
/////////

static const RegisterMetaType<BondHunter> r_hunter( MAGIC_ONLY,
                                                        "SireMol::BondHunter" );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const BondHunter &hunter)
{
    writeHeader(ds, r_hunter, 0);

    ds << static_cast<const Property&>(hunter);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, BondHunter &hunter)
{
    VersionID v = readHeader(ds, r_hunter);

    if (v == 0)
    {
        ds >> static_cast<Property&>(hunter);
    }
    else
        throw version_error(v, "0", r_hunter, CODELOC);

    return ds;
}

/** Constructor */
BondHunter::BondHunter() : Property()
{}

/** Copy constructor */
BondHunter::BondHunter(const BondHunter &other)
           : Property(other)
{}

/** Destructor */
BondHunter::~BondHunter()
{}

/////////
///////// Implementation of NullBondHunter
/////////

static const RegisterMetaType<NullBondHunter> r_nullhunter;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const NullBondHunter &hunter)
{
    writeHeader(ds, r_nullhunter, 1);

    ds << static_cast<const BondHunter&>(hunter);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, NullBondHunter &hunter)
{
    VersionID v = readHeader(ds, r_nullhunter);

    if (v == 1)
    {
        ds >> static_cast<BondHunter&>(hunter);
    }
    else
        throw version_error(v, "1", r_nullhunter, CODELOC);

    return ds;
}

/** Constructor */
NullBondHunter::NullBondHunter() : ConcreteProperty<NullBondHunter,BondHunter>()
{}

/** Copy constructor */
NullBondHunter::NullBondHunter(const NullBondHunter &other)
               : ConcreteProperty<NullBondHunter,BondHunter>(other)
{}

/** Destructor */
NullBondHunter::~NullBondHunter()
{}

const char* NullBondHunter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullBondHunter>() );
}

const char* NullBondHunter::what() const
{
    return NullBondHunter::typeName();
}

/** Return an empty connectivity object that is ready to be populated
    by the passed molecule */
Connectivity NullBondHunter::operator()(const MoleculeView &molview,
                                        const PropertyMap &map) const
{
    return Connectivity(molview.data());
}

const NullBondHunter& BondHunter::null()
{
    return *(create_shared_null<NullBondHunter>());
}

/////////
///////// Implementation of CovalentBondHunter
/////////

PropertyName CovalentBondHunterParameters::coords_param("coordinates");
PropertyName CovalentBondHunterParameters::elements_param("element");

CovalentBondHunterParameters CovalentBondHunter::params;

static const RegisterMetaType<CovalentBondHunter> r_covalenthunter;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CovalentBondHunter &hunter)
{
    writeHeader(ds, r_covalenthunter, 1);

    ds << hunter.tol << static_cast<const BondHunter&>(hunter);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CovalentBondHunter &hunter)
{
    VersionID v = readHeader(ds, r_covalenthunter);

    if (v == 1)
    {
        ds >> hunter.tol >> static_cast<BondHunter&>(hunter);
    }
    else
        throw version_error(v, "1", r_covalenthunter, CODELOC);

    return ds;
}

/** Construct a CovalentBondHunter with a specified tolerance */
CovalentBondHunter::CovalentBondHunter(double tolerance)
              : ConcreteProperty<CovalentBondHunter,BondHunter>(),
                tol(tolerance)
{}

/** Copy constructor */
CovalentBondHunter::CovalentBondHunter(const CovalentBondHunter &other)
              : ConcreteProperty<CovalentBondHunter,BondHunter>(other),
                tol(other.tol)
{}

/** Destructor */
CovalentBondHunter::~CovalentBondHunter()
{}

/** Copy assignment operator */
CovalentBondHunter& CovalentBondHunter::operator=(const CovalentBondHunter &other)
{
    tol = other.tol;
    BondHunter::operator=(other);

    return *this;
}

/** Comparison operator */
bool CovalentBondHunter::operator==(const CovalentBondHunter &other) const
{
    return tol == other.tol;
}

/** Comparison operator */
bool CovalentBondHunter::operator!=(const CovalentBondHunter &other) const
{
    return tol != other.tol;
}

void addAllIntraBonds(ConnectivityEditor &connectivity, CGIdx cgidx,
                      const CoordGroup &coords,
                      const AtomElements::Array &elements,
                      double tolerance)
{
    int nats = coords.count();
    BOOST_ASSERT(elements.count() == nats );

    const Vector *coords_array = coords.constData();
    const Element *elements_array = elements.constData();

    for (Index i(0); i<nats-1; ++i)
    {
        const Vector &v0 = coords_array[i];
        double r0 = elements_array[i].bondOrderRadius();

        for (Index j(i+1); j<nats; ++j)
        {
            const Vector &v1 = coords_array[j];
            double r1 = elements_array[j].bondOrderRadius();

            if ( SireMaths::pow_2( tolerance*(r0+r1) ) > Vector::distance2(v0,v1) )
            {
                connectivity.connect( CGAtomIdx(cgidx,i), CGAtomIdx(cgidx,j) );
            }
        }
    }
}

void addSomeIntraBonds(ConnectivityEditor &connectivity, CGIdx cgidx,
                       const CoordGroup &coords, const AtomElements::Array &elements,
                       const QSet<Index> &selected_atoms, double tolerance)
{
    BOOST_ASSERT( coords.count() == elements.count() );

    const Vector *coords_array = coords.constData();
    const Element *elements_array = elements.constData();

    for (QSet<Index>::const_iterator it0 = selected_atoms.constBegin();
         it0 != selected_atoms.constEnd();
         ++it0)
    {
        const Vector &v0 = coords_array[*it0];
        double r0 = elements_array[*it0].bondOrderRadius();

        QSet<Index>::const_iterator it1 = it0;
        for (++it1;
             it1 != selected_atoms.constEnd();
             ++it1)
        {
            const Vector &v1 = coords_array[*it1];
            double r1 = elements_array[*it1].bondOrderRadius();

            if ( SireMaths::pow_2(tolerance*(r0+r1)) > Vector::distance2(v0,v1) )
            {
                connectivity.connect( CGAtomIdx(cgidx,*it0), CGAtomIdx(cgidx,*it1) );
            }
        }
    }
}

//the largest covalent radius is 2.4 A, so if the CoordGroups are more
//than 5 A apart then they are not bonded!
const double max_radius2 = 25;

void addAllInterBonds(ConnectivityEditor &connectivity,
                      CGIdx cgidx0, CGIdx cgidx1,
                      const CoordGroup &coords0, const AtomElements::Array &elements0,
                      const CoordGroup &coords1, const AtomElements::Array &elements1,
                      double tolerance)
{
    if ( Vector::distance2(coords0.aaBox().center(), coords1.aaBox().center())
             > max_radius2 )
    {
        //the CoordGroups are too far apart to contain any bonded atoms
        return;
    }

    int nats0 = coords0.count();
    int nats1 = coords1.count();

    BOOST_ASSERT( nats0 == elements0.count() );
    BOOST_ASSERT( nats1 == elements1.count() );

    const Vector *coords0_array = coords0.constData();
    const Vector *coords1_array = coords1.constData();

    const Element *elements0_array = elements0.constData();
    const Element *elements1_array = elements1.constData();

    for (Index i(0); i<nats0; ++i)
    {
        const Vector &v0 = coords0_array[i];
        double r0 = elements0_array[i].bondOrderRadius();

        for (Index j(0); j<nats1; ++j)
        {
            const Vector &v1 = coords1_array[j];
            double r1 = elements1_array[j].bondOrderRadius();

            if ( SireMaths::pow_2( tolerance * (r0+r1) ) > Vector::distance2(v0,v1) )
            {
                connectivity.connect( CGAtomIdx(cgidx0,i), CGAtomIdx(cgidx1,j) );
            }
        }
    }
}

void addSomeInterBonds(ConnectivityEditor &connectivity, CGIdx cgidx0, CGIdx cgidx1,
                       const CoordGroup &coords0, const AtomElements::Array &elements0,
                       const CoordGroup &coords1, const AtomElements::Array &elements1,
                       const QSet<Index> &selected_atoms1, double tolerance)
{
    int nats0 = coords0.count();
    BOOST_ASSERT( elements0.count() == nats0 );
    BOOST_ASSERT( coords1.count() == elements1.count() );

    const Vector *coords0_array = coords0.constData();
    const Element *elements0_array = elements0.constData();

    const Vector *coords1_array = coords1.constData();
    const Element *elements1_array = elements1.constData();

    for (QSet<Index>::const_iterator it = selected_atoms1.constBegin();
         it != selected_atoms1.constEnd();
         ++it)
    {
        const Vector &v1 = coords1_array[*it];
        double r1 = elements1_array[*it].bondOrderRadius();

        for (Index i(0); i<nats0; ++i)
        {
            const Vector &v0 = coords0_array[i];
            double r0 = elements0_array[i].bondOrderRadius();

            if ( SireMaths::pow_2(tolerance*(r0+r1)) > Vector::distance2(v0,v1) )
            {
                connectivity.connect( CGAtomIdx(cgidx0,i), CGAtomIdx(cgidx1,*it) );
            }
        }
    }
}

void addSomeInterBonds(ConnectivityEditor &connectivity, CGIdx cgidx0, CGIdx cgidx1,
                       const CoordGroup &coords0, const AtomElements::Array &elements0,
                       const CoordGroup &coords1, const AtomElements::Array &elements1,
                       const QSet<Index> &selected_atoms0,
                       const QSet<Index> &selected_atoms1, double tolerance)
{
    BOOST_ASSERT( coords0.count() == elements0.count() );
    BOOST_ASSERT( coords1.count() == elements1.count() );

    const Vector *coords0_array = coords0.constData();
    const Element *elements0_array = elements0.constData();

    const Vector *coords1_array = coords1.constData();
    const Element *elements1_array = elements1.constData();

    for (QSet<Index>::const_iterator it0 = selected_atoms0.constBegin();
         it0 != selected_atoms0.constEnd();
         ++it0)
    {
        const Vector &v0 = coords0_array[*it0];
        double r0 = elements0_array[*it0].bondOrderRadius();

        for (QSet<Index>::const_iterator it1 = selected_atoms1.constBegin();
             it1 != selected_atoms1.constEnd();
             ++it1)
        {
            const Vector &v1 = coords1_array[*it1];
            double r1 = elements1_array[*it1].bondOrderRadius();

            if ( SireMaths::pow_2(tolerance*(r0+r1)) > Vector::distance2(v0,v1) )
            {
                connectivity.connect( CGAtomIdx(cgidx0,*it0), CGAtomIdx(cgidx1,*it1) );
            }
        }
    }
}

/** Return the connectivity for the atoms viewed in the view 'molview'.
    This uses the 'coordinates' property to find the coordinates of the
    atoms, and the 'element' property find the elements, from which the
    VDW radius is taken

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Connectivity CovalentBondHunter::operator()(const MoleculeView &molview,
                                            const PropertyMap &map) const
{
    // Get the names of the "coordinates" and "element" properties, with the
    // passed map taking precedence. Don't use the CovalentBondHunterParameters
    // since it doesn't allow any flexibility in the parameter naming.
    auto coords_prop = map["coordinates"];
    auto  elems_prop = map["element"];

    // Get the required properties
    const Property &coords_property = molview.data().property(coords_prop);

    const AtomCoords &coords = coords_property.asA<AtomCoords>();
    const CoordGroup *coords_array = coords.constData();

    const Property &elements_property = molview.data().property(elems_prop);

    const AtomElements &elements = elements_property.asA<AtomElements>();
    const AtomElements::Array *elements_array = elements.constData();

    ConnectivityEditor connectivity = Connectivity(molview.data()).edit();

    AtomSelection selected_atoms = molview.selection();

    if (selected_atoms.selectedAll())
    {
        int ncgroups = coords.count();

        for (CGIdx i(0); i<ncgroups; ++i)
        {
            //add the intra-cutgroup bonds
            addAllIntraBonds(connectivity, i, coords_array[i], elements_array[i], tol);

            for (CGIdx j(i+1); j<ncgroups; ++j)
            {
                addAllInterBonds(connectivity, i, j,
                                 coords_array[i], elements_array[i],
                                 coords_array[j], elements_array[j], tol);
            }
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        int ncgroups = coords.count();

        for (CGIdx i(0); i<ncgroups; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                addAllIntraBonds(connectivity, i, coords_array[i],
                                 elements_array[i], tol);

                for (CGIdx j(i+1); j<ncgroups; ++j)
                {
                    if (selected_atoms.selectedAll(j))
                    {
                        addAllInterBonds(connectivity, i, j,
                                         coords_array[i], elements_array[i],
                                         coords_array[j], elements_array[j], tol);
                    }
                    else
                    {
                        addSomeInterBonds(connectivity, i, j,
                                          coords_array[i], elements_array[i],
                                          coords_array[j], elements_array[j],
                                          selected_atoms.selectedAtoms(j), tol);
                    }
                }
            }
            else
            {
                QSet<Index> atoms0 = selected_atoms.selectedAtoms(i);

                addSomeIntraBonds(connectivity, i, coords_array[i],
                                  elements_array[i], atoms0, tol);

                for (CGIdx j(i+1); j<ncgroups; ++j)
                {
                    if (selected_atoms.selectedAll(j))
                    {
                        addSomeInterBonds(connectivity, j, i,
                                          coords_array[j], elements_array[j],
                                          coords_array[i], elements_array[i],
                                          atoms0, tol);
                    }
                    else
                    {
                        addSomeInterBonds(connectivity, i, j,
                                          coords_array[i], elements_array[i],
                                          coords_array[j], elements_array[j],
                                          atoms0, selected_atoms.selectedAtoms(j), tol);
                    }
                }
            }
        }
    }
    else
    {
        QList<CGIdx> selected_cgroups = selected_atoms.selectedCutGroups();
        int ncgroups = selected_cgroups.count();

        for (int ii=0; ii<ncgroups; ++ii)
        {
            CGIdx i = selected_cgroups.at(ii);

            if (selected_atoms.selectedAll(i))
            {
                addAllIntraBonds(connectivity, i, coords_array[i],
                                 elements_array[i], tol);

                for (int jj=ii+1; jj<ncgroups; ++jj)
                {
                    CGIdx j = selected_cgroups.at(jj);

                    if (selected_atoms.selectedAll(j))
                    {
                        addAllInterBonds(connectivity, i, j,
                                         coords_array[i], elements_array[i],
                                         coords_array[j], elements_array[j], tol);
                    }
                    else
                    {
                        addSomeInterBonds(connectivity, i, j,
                                          coords_array[i], elements_array[i],
                                          coords_array[j], elements_array[j],
                                          selected_atoms.selectedAtoms(j), tol);
                    }
                }
            }
            else
            {
                QSet<Index> atoms0 = selected_atoms.selectedAtoms(i);

                addSomeIntraBonds(connectivity, i, coords_array[i],
                                  elements_array[i], atoms0, tol);

                for (int jj=ii+1; jj<ncgroups; ++jj)
                {
                    CGIdx j = selected_cgroups.at(jj);

                    if (selected_atoms.selectedAll(j))
                    {
                        addSomeInterBonds(connectivity, j, i,
                                          coords_array[j], elements_array[j],
                                          coords_array[i], elements_array[i],
                                          atoms0, tol);
                    }
                    else
                    {
                        addSomeInterBonds(connectivity, i, j,
                                          coords_array[i], elements_array[i],
                                          coords_array[j], elements_array[j],
                                          atoms0, selected_atoms.selectedAtoms(j), tol);
                    }
                }
            }
        }
    }

    return connectivity;
}

const char* CovalentBondHunter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CovalentBondHunter>() );
}

/////////
///////// Implementation of ChemicalBondHunter
/////////

static const RegisterMetaType<ChemicalBondHunter> r_chemicalhunter;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ChemicalBondHunter &hunter)
{
    writeHeader(ds, r_chemicalhunter, 1);

    ds << static_cast<const CovalentBondHunter&>(hunter);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ChemicalBondHunter &hunter)
{
    VersionID v = readHeader(ds, r_chemicalhunter);

    if (v == 1)
    {
        ds >> static_cast<ChemicalBondHunter&>(hunter);
    }
    else
        throw version_error(v, "1", r_covalenthunter, CODELOC);

    return ds;
}

/** Construct with the default tolerance */
ChemicalBondHunter::ChemicalBondHunter()
                   : ConcreteProperty<ChemicalBondHunter,CovalentBondHunter>()
{}

/** Construct with the specified tolerance */
ChemicalBondHunter::ChemicalBondHunter(double tolerance)
                   : ConcreteProperty<ChemicalBondHunter,CovalentBondHunter>(tolerance)
{}

/** Copy constructor */
ChemicalBondHunter::ChemicalBondHunter(const ChemicalBondHunter &other)
                   : ConcreteProperty<ChemicalBondHunter,CovalentBondHunter>(other)
{}

/** Destructor */
ChemicalBondHunter::~ChemicalBondHunter()
{}

/** Return the connectivity for the atoms viewed in the view 'molview'.
    This uses the 'coordinates' property to find the coordinates of the
    atoms, and the 'element' property find the elements, from which the
    VDW radius is taken, and also to ensure that no atom has too many bonds.

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Connectivity ChemicalBondHunter::operator()(const MoleculeView &molview,
                                            const PropertyMap &map) const
{
    //get the connectivity using the covalent radii
    Connectivity connectivity = CovalentBondHunter::operator()(molview, map);

    const Property &elements_property = molview.data()
                                               .property( parameters().element() );

    const AtomElements &elements = elements_property.asA<AtomElements>();
    const AtomElements::Array *elements_array = elements.constData();

    const MoleculeInfoData &molinfo = molview.data().info();

    //now check that no atom has too many bonds...
    QHash<AtomIdx,int> overbonded_atoms;

    AtomSelection selected_atoms = molview.selection();

    if (selected_atoms.selectedAllCutGroups())
    {
        int ngroups = elements.count();

        for (CGIdx i(0); i<ngroups; ++i)
        {
            const AtomElements::Array &group_elements = elements_array[i];
            const Element *group_elements_array = group_elements.constData();

            if (selected_atoms.selectedAll(i))
            {
                int nats = group_elements.count();

                for (Index j(0); j<nats; ++j)
                {
                    AtomIdx atomidx = molinfo.atomIdx(CGAtomIdx(i,j));

                    if ( connectivity.nConnections(atomidx) >
                         int(group_elements_array[j].maxBonds()) )
                    {
                        overbonded_atoms.insert(atomidx,
                                                group_elements_array[j].maxBonds());
                    }
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    AtomIdx atomidx = molinfo.atomIdx(CGAtomIdx(i,j));

                    if ( connectivity.nConnections(atomidx) >
                         int(group_elements_array[j].maxBonds()) )
                    {
                        overbonded_atoms.insert(atomidx,
                                                group_elements_array[j].maxBonds());
                    }
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const AtomElements::Array &group_elements = elements_array[i];
            const Element *group_elements_array = group_elements.constData();

            if (selected_atoms.selectedAll(i))
            {
                int nats = group_elements.count();

                for (Index j(0); j<nats; ++j)
                {
                    AtomIdx atomidx = molinfo.atomIdx(CGAtomIdx(i,j));

                    if ( connectivity.nConnections(atomidx) >
                         int(group_elements_array[j].maxBonds()) )
                    {
                        overbonded_atoms.insert(atomidx,
                                                group_elements_array[j].maxBonds());
                    }
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    AtomIdx atomidx = molinfo.atomIdx(CGAtomIdx(i,j));

                    if ( connectivity.nConnections(atomidx) >
                         int(group_elements_array[j].maxBonds()) )
                    {
                        overbonded_atoms.insert(atomidx,
                                                group_elements_array[j].maxBonds());
                    }
                }
            }
        }
    }

    //overbonded atoms now contains the indicies of all of atoms with
    //too many bonds, and also the maximum number of bonds this atom
    //should have!

    if (overbonded_atoms.isEmpty())
        return connectivity;

    ConnectivityEditor editor = connectivity.edit();

    for (QHash<AtomIdx,int>::const_iterator it = overbonded_atoms.constBegin();
         it != overbonded_atoms.constEnd();
         ++it)
    {
        AtomIdx atom = it.key();
        int maxbonds = it.value();

        //do we still have too many bonds?
        int nbonds = editor.nConnections(atom);

        if (nbonds <= maxbonds)
            continue;

        //get the coordinates of this atom
        Vector v0 = Molecule(molview.data()).atom(atom)
                                            .property<Vector>(parameters().coordinates());

        //calculate all of the distances to the bonded atoms
        QSet<AtomIdx> connected_atoms = editor.connectionsTo(it.key());

        QMap<double, AtomIdx> distances2;

        foreach (AtomIdx bonded_atom, connected_atoms)
        {
            Vector v1 = Molecule(molview.data())
                                          .atom(bonded_atom)
                                          .property<Vector>(parameters().coordinates());

            distances2.insert( Vector::distance2(v0,v1), bonded_atom );
        }

        //remove the n longest bonds
        QMap<double,AtomIdx>::const_iterator it2 = distances2.constEnd();

        for (int i=nbonds; i>maxbonds; --i)
        {
            --it2;
            editor.disconnect( atom, it2.value() );
        }
    }

    return editor;
}

const char* ChemicalBondHunter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ChemicalBondHunter>() );
}
