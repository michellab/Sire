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

#include "mover.h"
#include "atomcoords.h"
#include "connectivity.h"
#include "weightfunction.h"

#include "atommatcher.h"
#include "atommatchers.h"
#include "bondid.h"
#include "angleid.h"
#include "dihedralid.h"
#include "improperid.h"

#include "tostring.h"

#include "SireMaths/align.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/axisset.h"
#include "SireMaths/vectorproperty.h"
#include "SireMaths/rotate.h"

#include "SireVol/coordgroup.h"
#include "SireVol/space.h"

#include "SireUnits/units.h"

#include "SireMol/errors.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireMaths;
using namespace SireVol;
using namespace SireUnits;

namespace SireMol
{
    /** Return the AxisSet needed to move 'view1' so that it is aligned against
        the atoms in 'view0' */
    Transform getAlignment(const MoleculeView &view0, const PropertyMap &map0,
                                          const MoleculeView &view1, const PropertyMap &map1,
                                          const AtomMatcher &matcher, bool fit)
    {
        const AtomCoords &coords0 = view0.data().property( map0["coordinates"] )
                                                .asA<AtomCoords>();
        
        const AtomCoords &coords1 = view1.data().property( map1["coordinates"] )
                                                .asA<AtomCoords>();
        
        QHash<AtomIdx,AtomIdx> map = matcher.match(view0, map0, view1, map1);
        
        if (map.isEmpty())
        {
            //there are no matching atoms - we can't do anything
            qDebug() << "No matching atoms!";
            return Transform();
        }
        
        QVector<Vector> p(map.count()), q(map.count());
        
        int n = 0;
        
        for (QHash<AtomIdx,AtomIdx>::const_iterator it = map.constBegin();
             it != map.constEnd();
             ++it)
        {
            p[n] = coords0.at( view0.data().info().cgAtomIdx(it.key()) );
            q[n] = coords1.at( view1.data().info().cgAtomIdx(it.value()) );
            
            /*qDebug() << n << view0.data().info().name(it.key()).toString() << "=="
                          << view1.data().info().name(it.value()).toString();
            qDebug() << n << p[n].toString() << q[n].toString();*/
            
            n += 1;
        }
        
        Transform t = SireMaths::getAlignment(p, q, fit);
        
        return t;
    }

    /** Return the AxisSet needed to move 'view1' so that it is aligned against
        the atoms in 'view0' */
    Transform getAlignment(const MoleculeView &view0,
                                          const MoleculeView &view1, bool fit)
    {
        return getAlignment(view0, PropertyMap(), view1, PropertyMap(),
                            AtomIdxMatcher(), fit);
    }

    /** Return the AxisSet needed to move 'view1' so that it is aligned against
        the atoms in 'view0' */
    Transform getAlignment(const MoleculeView &view0, const MoleculeView &view1,
                                          const PropertyMap &map, bool fit)
    {
        return getAlignment(view0, map, view1, map, AtomIdxMatcher(), fit);
    }

    /** Return the AxisSet needed to move 'view1' so that it is aligned against
        the atoms in 'view0' */
    Transform getAlignment(const MoleculeView &view0, const PropertyMap &map0,
                                          const MoleculeView &view1, const PropertyMap &map1,
                                          bool fit)
    {
        return getAlignment(view0, map0, view1, map1, AtomIdxMatcher(), fit);
    }

    /** Return the AxisSet needed to move 'view1' so that it is aligned against
        the atoms in 'view0' */
    Transform getAlignment(const MoleculeView &view0, const MoleculeView &view1,
                                          const AtomMatcher &matcher, bool fit)
    {
        return getAlignment(view0, PropertyMap(), view1, PropertyMap(), matcher, fit);
    }

    /** Return the AxisSet needed to move 'view1' so that it is aligned against
        the atoms in 'view0' */
    Transform getAlignment(const MoleculeView &view0, const MoleculeView &view1,
                                          const AtomMatcher &matcher,
                                          const PropertyMap &map, bool fit)
    {
        return getAlignment(view0, map, view1, map, matcher, fit);
    }

} // end of namespace SireMol

/** Constructor */
MoverBase::MoverBase()
{}

/** Construct a MoverBase that moves the specified atoms */
MoverBase::MoverBase(const AtomSelection &selected_atoms)
          : movable_atoms(selected_atoms)
{}

/** Copy constructor */
MoverBase::MoverBase(const MoverBase &other)
          : movable_atoms(other.movable_atoms)
{}

/** Destructor */
MoverBase::~MoverBase()
{}

/** Copy assignment operator */
MoverBase& MoverBase::operator=(const MoverBase &other)
{
    movable_atoms = other.movable_atoms;
    return *this;
}

/** Set the movable atoms of this Mover */
void MoverBase::setMovableAtoms(const AtomSelection &selection)
{
    movable_atoms = selection;
}

/** Translate the selected atoms from 'coords' by 'delta'.
    This function assumes that 'selected_atoms' is compatible
    with 'coords' */
void MoverBase::translate(AtomCoords &coords,
                          const AtomSelection &selected_atoms,
                          const Vector &delta)
{
    if (delta.isZero() or selected_atoms.selectedNone())
        return;

    int ncg = coords.count();

    if (selected_atoms.selectedAll())
    {
        //we are moving everything!
        coords.translate(delta);
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        //we are moving all of the CutGroups
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.translate(i, delta);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);

                CoordGroupEditor editor = coords.constData()[i].edit();

                foreach (const Index atom, atoms_to_move)
                {
                    editor.translate(atom, delta);
                }

                coords.set(i, editor.commit());
            }
        }
    }
    else
    {
        //we are moving only some CutGroups
        QList<CGIdx> cg_to_move = selected_atoms.selectedCutGroups();

        foreach (CGIdx i, cg_to_move)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.translate(i, delta);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);

                CoordGroupEditor editor = coords.constData()[i].edit();

                foreach (const Index &atom, atoms_to_move)
                {
                    editor.translate(atom, delta);
                }

                coords.set(i, editor.commit());
            }
        }
    }
}

/** Rotate the coordinates (in 'coords') of the specified selected
    atoms using the rotation matrix 'rotmat' about the point 'point'.
    This function assumes that coords and selected_atoms are compatible */
void MoverBase::rotate(AtomCoords &coords,
                       const AtomSelection &selected_atoms,
                       const Matrix &rotmat,
                       const Vector &point)
{
    if (selected_atoms.selectedNone())
        return;

    int ncg = coords.count();

    if (selected_atoms.selectedAll())
    {
        //we are rotating everything
        coords.rotate(rotmat,point);
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        //we are rotating every CutGroup
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.rotate(i, rotmat,point);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);

                CoordGroupEditor editor = coords.constData()[i].edit();

                foreach (const Index &atom, atoms_to_move)
                {
                    editor.rotate(atom, rotmat, point);
                }

                coords.set(i, editor.commit());
            }
        }
    }
    else
    {
        QList<CGIdx> cg_to_move = selected_atoms.selectedCutGroups();

        foreach (CGIdx i, cg_to_move)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.rotate(i, rotmat, point);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);
                CoordGroupEditor editor = coords.constData()[i].edit();

                foreach (const Index &atom, atoms_to_move)
                {
                    editor.rotate(atom, rotmat, point);
                }

                coords.set(i, editor.commit());
            }
        }
    }
}

/** Rotate the coordinates (in 'coords') of the specified selected
    atoms using the quaternion 'quat' about the point 'point'.
    This function assumes that coords and selected_atoms are compatible */
void MoverBase::rotate(AtomCoords &coords,
                       const AtomSelection &selected_atoms,
                       const Quaternion &quat, const Vector &point)
{
    MoverBase::rotate(coords, selected_atoms, quat.toMatrix(), point);
}

/** Transform the coordinates (in 'coords') of the specified selected
    atoms using the transformation 't'.
    This function assumes that coords and selected_atoms are compatible */
void MoverBase::transform(AtomCoords &coords,
                          const AtomSelection &selected_atoms,
                          const Transform &t)
{
    if (selected_atoms.selectedNone() or t.isZero())
        return;

    int ncg = coords.count();

    if (selected_atoms.selectedAll())
    {
        //we are rotating everything
        coords.transform(t);
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        //we are rotating every CutGroup
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.transform(i, t);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);

                CoordGroupEditor editor = coords.constData()[i].edit();

                foreach (const Index &atom, atoms_to_move)
                {
                    editor.transform(atom, t);
                }

                coords.set(i, editor.commit());
            }
        }
    }
    else
    {
        QList<CGIdx> cg_to_move = selected_atoms.selectedCutGroups();

        foreach (CGIdx i, cg_to_move)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.transform(i, t);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);
                CoordGroupEditor editor = coords.constData()[i].edit();

                foreach (const Index &atom, atoms_to_move)
                {
                    editor.transform(atom, t);
                }

                coords.set(i, editor.commit());
            }
        }
    }
}

/** Map the selected atoms from 'coords' into the axis set in 'axes'

    This function assumes that coords and selected_atoms are compatible!
*/
void MoverBase::mapInto(AtomCoords &coords,
                        const AtomSelection &selected_atoms,
                        const AxisSet &axes)
{
    if (selected_atoms.selectedNone())
        return;

    int ncg = coords.count();

    if (selected_atoms.selectedAll())
    {
        coords.mapInto(axes);
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.mapInto(i, axes);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);

                CoordGroupEditor editor = coords.constData()[i].edit();

                foreach (const Index &atom, atoms_to_move)
                {
                    editor.mapInto(atom, axes);
                }

                coords.set(i, editor.commit());
            }
        }
    }
    else
    {
        QList<CGIdx> cg_to_move = selected_atoms.selectedCutGroups();

        foreach (CGIdx i, cg_to_move)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.mapInto(i, axes);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);

                CoordGroupEditor editor = coords.constData()[i].edit();

                foreach (const Index &atom, atoms_to_move)
                {
                    editor.mapInto(atom, axes);
                }

                coords.set(i, editor.commit());
            }
        }
    }
}

/** Map the selected atoms from 'coords' from the frame 'from_frame'
    to the frame 'to_frame'

    This function assumes that coords and selected_atoms are compatible!
*/
void MoverBase::changeFrame(AtomCoords &coords,
                            const AtomSelection &selected_atoms,
                            const AxisSet &from_frame,
                            const AxisSet &to_frame)
{
    if (selected_atoms.selectedNone())
        return;

    int ncg = coords.count();

    if (selected_atoms.selectedAll())
    {
        //we are moving everything
        coords.changeFrame(from_frame, to_frame);
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.changeFrame(i, from_frame, to_frame);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);

                CoordGroupEditor editor = coords.constData()[i];

                foreach (const Index &atom, atoms_to_move)
                {
                    editor.changeFrame(atom, from_frame, to_frame);
                }

                coords.set(i, editor.commit());
            }
        }
    }
    else
    {
        QList<CGIdx> cg_to_move = selected_atoms.selectedCutGroups();

        foreach (CGIdx i, cg_to_move)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords.changeFrame(i, from_frame, to_frame);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);

                CoordGroupEditor editor = coords.constData()[i].edit();

                foreach (const Index &atom, atoms_to_move)
                {
                    editor.changeFrame(atom, from_frame, to_frame);
                }

                coords.set(i, editor.commit());
            }
        }
    }
}

/** Translate the selected atoms in the molecule whose data is in 'moldata'
    by 'delta', using 'coord_property' to get the coordinates to
    be translated. This function assumes that selected_atoms
    and moldata are compatible.

    \throw SireBase::missing_property
*/
void MoverBase::translate(MoleculeData &moldata,
                          const AtomSelection &selected_atoms,
                          const Vector &delta,
                          const PropertyMap &map)
{
    if (delta.isZero())
        return;

    //which property contains the coordinates?
    PropertyName coord_property = map["coordinates"];

    //get the current coordinates
    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //translate the coordinates of the selected atoms
    MoverBase::translate(coords, selected_atoms, delta);

    //set the new property
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);

    //if we have translated all atoms, then update the center point
    //of the molecule, if one has been set
    if (selected_atoms.selectedAll())
    {
        PropertyName center_property = map["center"];
        if (center_property.hasSource() and moldata.hasProperty(center_property))
        {
            Vector center = moldata.property(center_property).asA<VectorProperty>();
            moldata.setProperty(center_property.source(), VectorProperty( center + delta ));
        }
    }
}

/** Rotate the selected atoms in the molecule whose data
    is in 'moldata' using the rotation matrix 'rotmat'
    around the point 'point', using 'coord_property'
    to find the property containing the coordinates
    to be rotated.

    This function assumes that moldata and selected_atoms
    are compatible.

    \throw SireBase::missing_property
*/
void MoverBase::rotate(MoleculeData &moldata,
                       const AtomSelection &selected_atoms,
                       const Matrix &rotmat,
                       const Vector &point,
                       const PropertyMap &map)
{
    //get the name of the property that contains the coordinates
    PropertyName coord_property = map["coordinates"];

    //get the coordinates to be rotated
    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //rotate the coordinates
    MoverBase::rotate(coords, selected_atoms, rotmat, point);

    //set the new property
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);

    //if we have rotated all atoms, then update the center point
    //of the molecule, if one has been set
    if (selected_atoms.selectedAll())
    {
        PropertyName center_property = map["center"];
        if (center_property.hasSource() and moldata.hasProperty(center_property))
        {
            Vector center = moldata.property(center_property).asA<VectorProperty>();
            
            if (center != point)
                moldata.setProperty(center_property.source(),
                                    VectorProperty(SireMaths::rotate(center, rotmat, point)));
        }
    }
}

/** Transform the selected atoms in the molecule whose data
    is in 'moldata' using the transformation in 't', using 'coord_property'
    to find the property containing the coordinates
    to be rotated.

    This function assumes that moldata and selected_atoms
    are compatible.

    \throw SireBase::missing_property
*/
void MoverBase::transform(MoleculeData &moldata,
                          const AtomSelection &selected_atoms,
                          const Transform &t, const PropertyMap &map)
{
    //get the name of the property that contains the coordinates
    PropertyName coord_property = map["coordinates"];

    //get the coordinates to be rotated
    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //transform the coordinates
    MoverBase::transform(coords, selected_atoms, t);

    //set the new property
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);

    //if we have rotated all atoms, then update the center point
    //of the molecule, if one has been set
    if (selected_atoms.selectedAll())
    {
        PropertyName center_property = map["center"];
        if (center_property.hasSource() and moldata.hasProperty(center_property))
        {
            Vector center = moldata.property(center_property).asA<VectorProperty>();
            
            moldata.setProperty(center_property.source(),
                                VectorProperty(t.apply(center)));
        }
    }
}

/** This function maps the selected atoms from their current
    (cartesian) coordinate frame into the coordinate frame
    described by 'axes'. This function assumes that
    moldata and 'selected_atoms' are compatible

    \throw SireBase::missing_property
*/
void MoverBase::mapInto(MoleculeData &moldata,
                        const AtomSelection &selected_atoms,
                        const AxisSet &axes,
                        const PropertyMap &map)
{
    //get the name of the property that holds the coordinates
    PropertyName coord_property = map["coordinates"];

    //get the coordinates to be mapped
    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //map the coordinates
    MoverBase::mapInto(coords, selected_atoms, axes);

    //save the new coordinates
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);

    //if we have mapped all atoms, then update the center point
    //of the molecule, if one has been set
    if (selected_atoms.selectedAll())
    {
        PropertyName center_property = map["center"];
        if (center_property.hasSource() and moldata.hasProperty(center_property))
        {
            Vector center = moldata.property(center_property).asA<VectorProperty>();
            moldata.setProperty(center_property.source(),
                                VectorProperty(axes.fromIdentity(center)));
        }
    }
}

/** This function maps the selected atoms from the frame
    'from_frame' into the coordinate frame
    described by 'to_frame'. This function assumes that
    moldata and 'selected_atoms' are compatible

    \throw SireBase::missing_property
*/
void MoverBase::changeFrame(MoleculeData &moldata,
                            const AtomSelection &selected_atoms,
                            const AxisSet &from_frame,
                            const AxisSet &to_frame,
                            const PropertyMap &map)
{
    //get the name of the property that holds the coordinates
    PropertyName coord_property = map["coordinates"];

    //get the coordinates to be mapped
    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //map the coordinates
    MoverBase::changeFrame(coords, selected_atoms, from_frame, to_frame);

    //save the new coordinates
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);

    //if we changed the frame of all atoms, then update the center point
    //of the molecule, if one has been set
    if (selected_atoms.selectedAll())
    {
        PropertyName center_property = map["center"];
        if (center_property.hasSource() and moldata.hasProperty(center_property))
        {
            Vector center = moldata.property(center_property).asA<VectorProperty>();
            moldata.setProperty(center_property.source(),
                                VectorProperty(from_frame.toFrame(to_frame,center)));
        }
    }
}

/** Map the atoms we are allowed to move into the passed axes,
    finding the coordinates using the passed property map

    \throw SireBase::missing_property
*/
void MoverBase::mapInto(MoleculeData &moldata,
                        const AxisSet &axes,
                        const PropertyMap &map) const
{
    MoverBase::mapInto(moldata, movable_atoms, axes, map);
}

/** Map the atoms we are allowed to move from the coordinate frame
    'from_frame' to the coordinate frame 'to_frame', finding the 
    coordinates using the passed property map
    
    \throw SireBase::missing_property
*/
void MoverBase::changeFrame(MoleculeData &moldata,
                            const AxisSet &from_frame,
                            const AxisSet &to_frame,
                            const PropertyMap &map) const
{
    MoverBase::changeFrame(moldata, movable_atoms, from_frame,
                           to_frame, map);
}

/** Translate atoms we are allowed to move from the molecule whose
    data is in 'moldata' by 'delta', finding the coordinates
    using 'coord_property'

    \throw SireBase::missing_property
*/
void MoverBase::translate(MoleculeData &moldata,
                          const Vector &delta,
                          const PropertyMap &map) const
{
    MoverBase::translate(moldata, movable_atoms, delta, map);
}

/** Rotate the atoms we are allowed to move from the molecule whose
    data is in 'moldata' using the quaternion 'quat' about the
    point 'point', finding the coordinates using 'coord_property'

    \throw SireBase::missing_property
*/
void MoverBase::rotate(MoleculeData &moldata,
                       const Quaternion &quat,
                       const Vector &point,
                       const PropertyMap &map) const
{
    MoverBase::rotate(moldata, movable_atoms, quat.toMatrix(),
                      point, map);
}

/** Rotate the atoms we are allowed to move from the molecule whose
    data is in 'moldata' using the rotation matrix 'rotmat' about the
    point 'point', finding the coordinates using 'coord_property'

    \throw SireBase::missing_property
*/
void MoverBase::rotate(MoleculeData &moldata,
                       const Matrix &rotmat,
                       const Vector &point,
                       const PropertyMap &map) const
{
    MoverBase::rotate(moldata, movable_atoms, rotmat,
                      point, map);
}

/** Transform the atoms we are allowed to move from the molecule whose
    data is in 'moldata' using the transformation 't',
    finding the coordinates using 'coord_property'

    \throw SireBase::missing_property
*/
void MoverBase::transform(MoleculeData &moldata,
                          const Transform &t,
                          const PropertyMap &map) const
{
    MoverBase::transform(moldata, movable_atoms, t, map);
}

/** Apply anchors to the groups - this clears a group if it
    contains an anchor atom */
static void applyAnchors(const AtomSelection &anchors,
                         tuple<AtomSelection,AtomSelection> &groups)
{
    if (groups.get<0>().intersects(anchors))
    {
        groups.get<0>().deselectAll();
    }
    
    if (groups.get<1>().intersects(anchors))
    {
        groups.get<1>().deselectAll();
    }
}

static const AtomSelection& getAnchors(const MoleculeData &moldata,
                                       const PropertyMap &map)
{
    return moldata.property(map["anchors"]).asA<AtomSelection>();
}

/** Change the length of the bond identified by 'bond' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    The bond is labelled;

    atom0--atom1

    The molecule is split into two about this bond, i.e.
    atom0 and everything it is bonded to is in group0, while
    atom1 and everything it is bonded to is in group1.

    The two groups are then translated along the vector atom0->atom1

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const BondID &bond,
                       SireUnits::Dimension::Length delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;

    //get the indicies of the two atoms of the bond
    tuple<AtomIdx,AtomIdx> atomidxs = bond.map(moldata.info());

    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();

    //get the connectivity property that is used to split
    //the molecule into two parts
    const Connectivity &connectivity =
            moldata.property(map["connectivity"]).asA<Connectivity>();

    //split the molecule into the two parts that are
    //going to move - the two groups are only able to
    //contain the atoms that are in 'movable_atoms'
    tuple<AtomSelection,AtomSelection> groups =
                        connectivity.split(atom0, atom1, movable_atoms);

    //see if there are any anchors that must be applied to
    //the section of molecule
    if (map.specified("anchors"))
    {
        applyAnchors( getAnchors(moldata, map), groups );
    }

    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();

    double weight0, weight1;

    if (group0.isEmpty() and group1.isEmpty())
    {
        if (map.specified("anchors"))
        {
            throw SireMol::anchor_error( QObject::tr(
                "Splitting the molecule into two about %1 has resulted "
                "in two groups that are both anchored (anchors = %2)")
                    .arg(bond.toString(), 
                         Sire::toString( getAnchors(moldata, map).selectedAtoms()) ),
                             CODELOC );
        }
    
        throw SireError::program_bug( QObject::tr(
            "Splitting the molecule about %1 has resulted "
            "in two empty groups!").arg(bond.toString()), CODELOC );
    }
    else if (group0.isEmpty())
    {
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        //get the weighting function that is used to weight the two
        //sides of the move
        const WeightFunction &weightfunc =
            moldata.property(map["weight function"],
                             WeightFunction::null()).asA<WeightFunction>();

        weight0 = 1 - weightfunc(moldata, group0, group1, map);
        weight1 = 1 - weightfunc(moldata, group1, group0, map);
    }

    //now get property containing the coordinates of the atoms
    PropertyName coord_property = map["coordinates"];

    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //use these coordinates to calculate the unit vector that
    //points along the bond
    Vector unit_vec = (coords[moldata.info().cgAtomIdx(atom1)] -
                       coords[moldata.info().cgAtomIdx(atom0)]).normalise();

    //scale the vector by 'delta'
    unit_vec *= delta;

    //now translate the groups along this vector by their weighted
    //amount of delta
    if (weight0 != 0)
        MoverBase::translate(coords, group0, -weight0 * unit_vec);

    if (weight1 != 0)
        MoverBase::translate(coords, group1, weight1 * unit_vec);

    //save the new coordinates
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);
}

/** Change the size of the angle identified by 'angle' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    The angle is labelled;

    atom0     atom2
         \    /
         atom1

    The molecule is split by the atom0-..-atom2 bond, i.e. atom0
    is in group0, atom2 is in group1 and atom1 is not in any
    group and is not moved.

    The two groups are rotated around the vector perpendicular to
    atom0->atom1 and atom2->atom1, about the point atom1

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const AngleID &angle,
                       SireUnits::Dimension::Angle delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;

    //get the indicies of the atoms in the angle
    tuple<AtomIdx,AtomIdx,AtomIdx> atomidxs = angle.map(moldata.info());

    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();
    AtomIdx atom2 = atomidxs.get<2>();

    //get the connectivity that is used to split the
    //molecule into two parts
    const Connectivity &connectivity =
            moldata.property(map["connectivity"]).asA<Connectivity>();

    //split the molecule into the two moving parts
    tuple<AtomSelection,AtomSelection> groups =
                          connectivity.split(atom0, atom1, atom2,
                                             movable_atoms);

    //see if there are any anchors that hold part of the
    //molecule stationary
    if (map.specified("anchors"))
        applyAnchors( getAnchors(moldata, map), groups );


    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();

    double weight0, weight1;

    if (group0.isEmpty() and group1.isEmpty())
    {
        if (map.specified("anchors"))
        {
            throw SireMol::anchor_error( QObject::tr(
                "Splitting the molecule into two about %1 has resulted "
                "in two groups that are both anchored (anchors = %2)")
                    .arg(angle.toString(), 
                         Sire::toString( getAnchors(moldata, map).selectedAtoms()) ),
                             CODELOC );
        }
    
        throw SireError::program_bug( QObject::tr(
            "Splitting the molecule about the %1 has resulted "
            "in two empty groups!").arg(angle.toString()), CODELOC );
    }
    else if (group0.isEmpty())
    {
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        //get the weighting function that is used to weight the
        //two sides of the move
        const WeightFunction &weightfunc =
                moldata.property(map["weight function"],
                                 WeightFunction::null()).asA<WeightFunction>();

        weight0 = 1 - weightfunc(moldata, group0, group1, map);
        weight1 = 1 - weightfunc(moldata, group1, group0, map);
    }

    //get the coordinates that are to be changed
    PropertyName coord_property = map["coordinates"];
    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //get the coordinates of the three atoms that comprise the angle
    const Vector &coords0 = coords[moldata.info().cgAtomIdx(atom0)];
    const Vector &coords1 = coords[moldata.info().cgAtomIdx(atom1)];
    const Vector &coords2 = coords[moldata.info().cgAtomIdx(atom2)];

    //get the vector perpendicular to the angle
    // This function contains lots of checks to ensure that a parallel
    // vector is always returned, even for co-linear or co-located atoms
    Vector perp = Vector::cross( coords2-coords0, coords1-coords0 );
    
    //rotate the two groups
    if (weight0 != 0)
        MoverBase::rotate(coords, group0,
                          Quaternion(-weight0*delta, perp), coords1);

    if (weight1 != 0)
        MoverBase::rotate(coords, group1,
                          Quaternion(weight1*delta, perp), coords1);

    //save the new coordinates
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);
}

/** Change the size of the dihedral identified by 'dihedra' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    The dihedral is labelled;

    atom0           atom3
        \           /
        atom1--atom2

    The molecule is split by the atom0-..-atom3 bond (i.e.
    atom0 is in group0, atom3 is in group1, while atom1 and atom2
    are not in any group and are not moved).

    The groups are then rotated about the vector atom1->atom2

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const DihedralID &dihedral,
                       SireUnits::Dimension::Angle delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;

    //get the indicies of the atoms that comprise this dihedral
    tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx> atomidxs =
                                               dihedral.map(moldata.info());

    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();
    AtomIdx atom2 = atomidxs.get<2>();
    AtomIdx atom3 = atomidxs.get<3>();

    //now get the connectivity of the molecule
    const Connectivity &connectivity =
              moldata.property(map["connectivity"]).asA<Connectivity>();

    tuple<AtomSelection,AtomSelection> groups = 
                        connectivity.split(atom0, atom1, 
                                           atom2, atom3, movable_atoms);

    //see if there are any anchors that hold part of the
    //molecule stationary
    if (map.specified("anchors"))
        applyAnchors( getAnchors(moldata, map), groups );

    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();

    double weight0, weight1;

    if (group0.isEmpty() and group1.isEmpty())
    {
        if (map.specified("anchors"))
        {
            throw SireMol::anchor_error( QObject::tr(
                "Splitting the molecule into two about %1 has resulted "
                "in two groups that are both anchored (anchors = %2)")
                    .arg(dihedral.toString(), 
                         Sire::toString( getAnchors(moldata, map).selectedAtoms()) ),
                             CODELOC );
        }
    
        throw SireError::program_bug( QObject::tr(
            "Splitting the molecule about the %1 has resulted "
            "in two empty groups!").arg(dihedral.toString()), CODELOC );
    }
    else if (group0.isEmpty())
    {
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        const WeightFunction &weightfunc =
                      moldata.property(map["weight function"],
                                       WeightFunction::null()).asA<WeightFunction>();

        weight0 = 1 - weightfunc(moldata, group0, group1, map);
        weight1 = 1 - weightfunc(moldata, group1, group0, map);
    }

    //get the coordinates to be moved
    PropertyName coord_property = map["coordinates"];
    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //get the coordinates of the central two atoms of the dihedral
    const Vector &coords1 = coords[moldata.info().cgAtomIdx(atom1)];
    const Vector &coords2 = coords[moldata.info().cgAtomIdx(atom2)];

    //get the vector about which the two parts of the molecule
    //are rotated
    Vector dihvec = coords2 - coords1;

    //now rotate the two parts of the molecule
    if (weight0 != 0)
        MoverBase::rotate(coords, group0,
                          Quaternion(-weight0*delta, dihvec), coords1);

    if (weight1 != 0)
        MoverBase::rotate(coords, group1,
                          Quaternion(weight1*delta, dihvec), coords2);

    // See if there are dihedrals whose rotation should be synchronized to this dihedral
    // (improper torsions)
    // Get the list from the property map
    // For each dihedral, find the atom that must be rotated and with which point/axis/weight
    // Call MoverBase::rotate accordingly

    //save the new coordinates
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);
}

/** Change the size of the dihedral identified by 'bond' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    The dihedral is labelled;

      X--atom0--atom1--X

    The molecule is split by the atom0-atom1 bond, and the two
    groups are rotated about the vector connecting atom0->atom1
    about the point atom0 for group0, and atom1 for group1

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const BondID &bond,
                       SireUnits::Dimension::Angle delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;

    //get the indicies of the atoms that comprise this dihedral
    tuple<AtomIdx,AtomIdx> atomidxs = bond.map(moldata.info());

    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();

    //now get the connectivity of the molecule
    const Connectivity &connectivity =
              moldata.property(map["connectivity"]).asA<Connectivity>();

    tuple<AtomSelection,AtomSelection> groups = 
                            connectivity.split(atom0, atom1, movable_atoms);

    //see if there are any anchors that hold part of the
    //molecule stationary
    if (map.specified("anchors"))
        applyAnchors( getAnchors(moldata, map), groups );

    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();

    double weight0, weight1;

    if (group0.isEmpty() and group1.isEmpty())
    {
        if (map.specified("anchors"))
        {
            throw SireMol::anchor_error( QObject::tr(
                "Splitting the molecule into two about %1 has resulted "
                "in two groups that are both anchored (anchors = %2)")
                    .arg(bond.toString(), 
                         Sire::toString( getAnchors(moldata, map).selectedAtoms()) ),
                             CODELOC );
        }
    
        throw SireError::program_bug( QObject::tr(
            "Splitting the molecule about the %1 has resulted "
            "in two empty groups!").arg(bond.toString()), CODELOC );
    }
    else if (group0.isEmpty())
    {
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        const WeightFunction &weightfunc =
                     moldata.property(map["weight function"],
                                      WeightFunction::null()).asA<WeightFunction>();

        weight0 = 1 - weightfunc(moldata, group0, group1, map);
        weight1 = 1 - weightfunc(moldata, group1, group0, map);
    }

    //get the coordinates to be moved
    PropertyName coord_property = map["coordinates"];
    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //get the coordinates of the central two atoms of the dihedral
    const Vector &coords0 = coords[moldata.info().cgAtomIdx(atom0)];
    const Vector &coords1 = coords[moldata.info().cgAtomIdx(atom1)];

    //get the vector about which the two parts of the molecule
    //are rotated
    Vector dihvec = coords1 - coords0;

    //now rotate the two parts of the molecule
    if (weight0 != 0)
        MoverBase::rotate(coords, group0,
                          Quaternion(-weight0*delta, dihvec), coords0);

    if (weight1 != 0)
        MoverBase::rotate(coords, group1,
                          Quaternion(weight1*delta, dihvec), coords1);

    //save the new coordinates
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);
}

/** Change the size of the improper identified by 'improper' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    An improper is labelled;

                 atom2
                /
    atom0--atom1
                \
                 atom3

    The molecule is split into two along the atom0-atom1 bond
    (i.e. atom0 is in one group, while atom1, atom2 and atom3 are
    in the other). The groups are then rotated around the vector
    from atom2->atom3, about the point 'atom1'

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const ImproperID &improper,
                       SireUnits::Dimension::Angle delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;

    //get the indicies of the atoms that comprise this dihedral
    tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx> atomidxs =
                                               improper.map(moldata.info());

    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();
    AtomIdx atom2 = atomidxs.get<2>();
    AtomIdx atom3 = atomidxs.get<3>();

    //now get the connectivity of the molecule
    const Connectivity &connectivity =
            moldata.property(map["connectivity"]).asA<Connectivity>();

    tuple<AtomSelection,AtomSelection> groups =
                      connectivity.split(atom0, atom1, movable_atoms);

    //see if there are any anchors that hold part of the
    //molecule stationary
    if (map.specified("anchors"))
        applyAnchors( getAnchors(moldata, map), groups );

    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();

    double weight0, weight1;

    if (group0.isEmpty() and group1.isEmpty())
    {
        if (map.specified("anchors"))
        {
            throw SireMol::anchor_error( QObject::tr(
                "Splitting the molecule into two about %1 has resulted "
                "in two groups that are both anchored (anchors = %2)")
                    .arg(improper.toString(), 
                         Sire::toString( getAnchors(moldata, map).selectedAtoms()) ),
                             CODELOC );
        }
    
        throw SireError::program_bug( QObject::tr(
            "Splitting the molecule about the %1 has resulted "
            "in two empty groups!").arg(improper.toString()), CODELOC );
    }
    else if (group0.isEmpty())
    {
        BOOST_ASSERT( not group1.isEmpty() );
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        const WeightFunction &weightfunc =
                   moldata.property(map["weight function"],
                                    WeightFunction::null()).asA<WeightFunction>();

        weight0 = 1 - weightfunc(moldata, group0, group1, map);
        weight1 = 1 - weightfunc(moldata, group1, group0, map);
    }

    //get the coordinates to be moved
    PropertyName coord_property = map["coordinates"];
    AtomCoords coords = moldata.property(coord_property).asA<AtomCoords>();

    //get the coordinates of the last three atoms of the improper
    const Vector &coords1 = coords[moldata.info().cgAtomIdx(atom1)];
    const Vector &coords2 = coords[moldata.info().cgAtomIdx(atom2)];
    const Vector &coords3 = coords[moldata.info().cgAtomIdx(atom3)];

    //get the vector from atom2 to atom3
    Vector impvec = coords3 - coords2;

    //now rotate the two parts of the molecule
    if (weight0 != 0)
        MoverBase::rotate(coords, group0,
                          Quaternion(-weight0*delta, impvec), coords1);

    if (weight1 != 0)
        MoverBase::rotate(coords, group1,
                          Quaternion(weight1*delta, impvec), coords1);

    //save the new coordinates
    if (coord_property.hasSource())
        moldata.setProperty(coord_property.source(), coords);
}

/** Set the length of the bond identified by 'bond' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    The bond is labelled;

    atom0--atom1

    The molecule is split into two about this bond, i.e.
    atom0 and everything it is bonded to is in group0, while
    atom1 and everything it is bonded to is in group1.

    The two groups are then translated along the vector atom0->atom1

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::set(MoleculeData &moldata, const BondID &bond,
                    SireUnits::Dimension::Length value,
                    const PropertyMap &map) const
{
    SireUnits::Dimension::Length current_value( bond.size(moldata,map) );
    this->change(moldata, bond, value - current_value, map);
}

/** Set the size of the angle identified by 'angle' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    The angle is labelled;

    atom0     atom2
         \    /
         atom1

    The molecule is split by the atom0-..-atom2 bond, i.e. atom0
    is in group0, atom2 is in group1 and atom1 is not in any
    group and is not moved.

    The two groups are rotated around the vector perpendicular to
    atom0->atom1 and atom2->atom1, about the point atom1

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::set(MoleculeData &moldata, const AngleID &angle,
                    SireUnits::Dimension::Angle value,
                    const PropertyMap &map) const
{
    SireUnits::Dimension::Angle current_value = angle.size(moldata, map);

    this->change(moldata, angle, value - current_value, map);
}

/** Set the size of the dihedral identified by 'dihedra' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    The dihedral is labelled;

    atom0           atom3
        \           /
        atom1--atom2

    The molecule is split by the atom0-..-atom3 bond (i.e.
    atom0 is in group0, atom3 is in group1, while atom1 and atom2
    are not in any group and are not moved).

    The groups are then rotated about the vector atom1->atom2

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::set(MoleculeData &moldata, const DihedralID &dihedral,
                    SireUnits::Dimension::Angle value,
                    const PropertyMap &map) const
{
    SireUnits::Dimension::Angle current_value = dihedral.size(moldata, map);

    this->change(moldata, dihedral, value - current_value, map);
}

/** Set the size of the dihedral identified by 'dihedral' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    The dihedral is labelled;

    atom0           atom3
        \           /
        atom1--atom2

    The molecule is split by the atom0-..-atom3 bond (i.e.
    atom0 is in group0, atom3 is in group1, while atom1 and atom2
    are not in any group and are not moved).

    The groups are then rotated about the vector atom1->atom2

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::setAll(MoleculeData &moldata, const DihedralID &dihedral,
                       SireUnits::Dimension::Angle value,
                       const PropertyMap &map) const
{
    SireUnits::Dimension::Angle current_value = dihedral.size(moldata, map);

    this->change(moldata, BondID(dihedral.atom1(), dihedral.atom2()),
                 value - current_value, map);
}

/** Set the size of the improper identified by 'improper' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties

    This only moves the movable atoms in this view, and
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.

    An improper is labelled;

                 atom2
                /
    atom0--atom1
                \
                 atom3

    The molecule is split into two along the atom0-atom1 bond
    (i.e. atom0 is in one group, while atom1, atom2 and atom3 are
    in the other). The groups are then rotated around the vector
    from atom2->atom3, about the point 'atom1'

    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::set(MoleculeData &moldata, const ImproperID &improper,
                    SireUnits::Dimension::Angle value,
                    const PropertyMap &map) const
{
    SireUnits::Dimension::Angle current_value = improper.size(moldata, map);

    this->change(moldata, improper, value - current_value, map);
}
