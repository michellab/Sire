/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREMOL_GEOMETRYPERTURBATION_H
#define SIREMOL_GEOMETRYPERTURBATION_H

#include "perturbation.h"
#include "bondid.h"
#include "angleid.h"
#include "dihedralid.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class GeometryPerturbation;
class GeometryPerturbations;
class BondPerturbation;
class AnglePerturbation;
class DihedralPerturbation;
class NullGeometryPerturbation;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::GeometryPerturbation&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::GeometryPerturbation&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::GeometryPerturbations&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::GeometryPerturbations&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::BondPerturbation&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::BondPerturbation&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AnglePerturbation&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AnglePerturbation&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::DihedralPerturbation&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::DihedralPerturbation&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::NullGeometryPerturbation&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::NullGeometryPerturbation&);

namespace SireMol
{

template<class T> class Mover;

/** This is the base class of all geometry perturbations.
    These are perturbations that affect the geometry
    of a molecule 
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT GeometryPerturbation : public Perturbation
{

friend QDataStream& ::operator<<(QDataStream&, const Perturbation&);
friend QDataStream& ::operator>>(QDataStream&, Perturbation&);

friend class GeometryPerturbations;

public:
    GeometryPerturbation(const GeometryPerturbation &other);
    
    ~GeometryPerturbation();
    
    static const char* typeName();
    
    QSet<QString> requiredProperties() const;

    static const NullGeometryPerturbation& null();

protected:
    GeometryPerturbation(const PropertyMap &map = PropertyMap());
    GeometryPerturbation(const SireCAS::Expression &expression,
                         const PropertyMap &map = PropertyMap());
    
    GeometryPerturbation& operator=(const GeometryPerturbation &other);
    
    bool operator==(const GeometryPerturbation &other) const;
    bool operator!=(const GeometryPerturbation &other) const;
    
    void perturbMolecule(MolEditor &molecule, const SireCAS::Values &values) const;
    
    virtual void perturbMolecule(Mover<Molecule> &molecule,
                                 const SireCAS::Values &values) const=0;
};

class SIREMOL_EXPORT NullGeometryPerturbation
        : public SireBase::ConcreteProperty<NullGeometryPerturbation,GeometryPerturbation>
{

friend QDataStream& ::operator<<(QDataStream&, const NullGeometryPerturbation&);
friend QDataStream& ::operator>>(QDataStream&, NullGeometryPerturbation&);

public:
    NullGeometryPerturbation();
    NullGeometryPerturbation(const NullGeometryPerturbation &other);
    ~NullGeometryPerturbation();
    
    NullGeometryPerturbation& operator=(const NullGeometryPerturbation &other);
    
    bool operator==(const NullGeometryPerturbation &other) const;
    bool operator!=(const NullGeometryPerturbation &other) const;
    
    static const char* typeName();

    QSet<Symbol> requiredSymbols() const;
    QSet<QString> requiredProperties() const;
    
    bool wouldChange(const Molecule&, const SireCAS::Values&) const;
    void perturbMolecule(MolEditor&, const SireCAS::Values&) const;
    
protected:
    void perturbMolecule(Mover<Molecule> &molecule, const SireCAS::Values &values) const;
};

typedef SireBase::PropPtr<GeometryPerturbation> GeomPertPtr;

/** This class holds a collection of geometry perturbations */
class SIREMOL_EXPORT GeometryPerturbations
        : public SireBase::ConcreteProperty<GeometryPerturbations,GeometryPerturbation>
{

friend QDataStream& ::operator<<(QDataStream&, const GeometryPerturbations&);
friend QDataStream& ::operator>>(QDataStream&, GeometryPerturbations&);

public:
    GeometryPerturbations();
    
    GeometryPerturbations(const GeometryPerturbation &perturbation);
    GeometryPerturbations(const QList<GeomPertPtr> &perturbations);
                          
    GeometryPerturbations(const GeometryPerturbations &other);
    
    ~GeometryPerturbations();
    
    static const char* typeName();
    
    GeometryPerturbations& operator=(const GeometryPerturbations &other);
    
    bool operator==(const GeometryPerturbations &other) const;
    bool operator!=(const GeometryPerturbations &other) const;
    
    QString toString() const;
    
    QList<GeomPertPtr> perturbations() const;

    PerturbationPtr recreate(const SireCAS::Expression &mapping_function) const;
    PerturbationPtr recreate(const PropertyMap &map) const;
    PerturbationPtr recreate(const SireCAS::Expression &mapping_function,
                             const PropertyMap &map) const;
    
    PerturbationPtr substitute(const SireCAS::Identities &identities) const;
    PerturbationPtr substitute(const SireCAS::Symbol &old_symbol,
                               const SireCAS::Symbol &new_symbol) const;
    
    QList<PerturbationPtr> children() const;
    
    QSet<Symbol> requiredSymbols() const;
    QSet<QString> requiredProperties() const;

    bool wouldChange(const Molecule &molecule, const Values &values) const;
    
protected:
    void perturbMolecule(Mover<Molecule> &molecule, const SireCAS::Values &values) const;

private:
    /** The perturbations to be applied */
    QList<GeomPertPtr> perts;
};

/** This perturbation moves a bond between two lengths.

    This uses the "anchors" property to anchor parts
    of the molecule, the "weight function" property
    to weight the motion of the parts of the molecule,
    and the "coordinates" property to get the coordinates
    to move

    @author Christopher Woods
*/
class SIREMOL_EXPORT BondPerturbation
        : public SireBase::ConcreteProperty<BondPerturbation,GeometryPerturbation>
{

friend QDataStream& ::operator<<(QDataStream&, const BondPerturbation&);
friend QDataStream& ::operator>>(QDataStream&, BondPerturbation&);

public:
    BondPerturbation();
    
    BondPerturbation(const BondID &bond, 
                     const SireUnits::Dimension::Length &start,
                     const SireUnits::Dimension::Length &end,
                     const PropertyMap &map = PropertyMap());
                     
    BondPerturbation(const BondID &bond, 
                     const SireUnits::Dimension::Length &start,
                     const SireUnits::Dimension::Length &end,
                     const SireCAS::Expression &mapping_function,
                     const PropertyMap &map = PropertyMap());
    
    BondPerturbation(const AtomID &atom0, const AtomID &atom1,
                     const SireUnits::Dimension::Length &start,
                     const SireUnits::Dimension::Length &end,
                     const PropertyMap &map = PropertyMap());
                     
    BondPerturbation(const AtomID &atom0, const AtomID &atom1,
                     const SireUnits::Dimension::Length &start,
                     const SireUnits::Dimension::Length &end,
                     const SireCAS::Expression &mapping_function,
                     const PropertyMap &map = PropertyMap());
    
    BondPerturbation(const BondPerturbation &other);
    
    ~BondPerturbation();
    
    static const char* typeName();
    
    BondPerturbation& operator=(const BondPerturbation &other);
    
    bool operator==(const BondPerturbation &other) const;
    bool operator!=(const BondPerturbation &other) const;
    
    QString toString() const;
    
    const BondID& bond() const;
    
    const SireUnits::Dimension::Length& start() const;
    const SireUnits::Dimension::Length& end() const;

    bool wouldChange(const Molecule &molecule, const Values &values) const;

protected:
    void perturbMolecule(Mover<Molecule> &molecule, const Values &values) const;

private:
    /** The ID for the bond */
    BondID bondid;
    
    /** The start and end lengths of the bond */
    SireUnits::Dimension::Length start_size, end_size;
};

/** This perturbation moves an angle between two sizes.

    This uses the "anchors" property to anchor parts
    of the molecule, the "weight function" property
    to weight the motion of the parts of the molecule,
    and the "coordinates" property to get the coordinates
    to move

    @author Christopher Woods
*/
class SIREMOL_EXPORT AnglePerturbation
        : public SireBase::ConcreteProperty<AnglePerturbation,GeometryPerturbation>
{

friend QDataStream& ::operator<<(QDataStream&, const AnglePerturbation&);
friend QDataStream& ::operator>>(QDataStream&, AnglePerturbation&);

public:
    AnglePerturbation();
    
    AnglePerturbation(const AngleID &angle, 
                      const SireUnits::Dimension::Angle &start,
                      const SireUnits::Dimension::Angle &end,
                      const PropertyMap &map = PropertyMap());
                     
    AnglePerturbation(const AngleID &angle, 
                      const SireUnits::Dimension::Angle &start,
                      const SireUnits::Dimension::Angle &end,
                      const SireCAS::Expression &mapping_function,
                      const PropertyMap &map = PropertyMap());
    
    AnglePerturbation(const AtomID &atom0, const AtomID &atom1, 
                      const AtomID &atom2,
                      const SireUnits::Dimension::Angle &start,
                      const SireUnits::Dimension::Angle &end,
                      const PropertyMap &map = PropertyMap());
                     
    AnglePerturbation(const AtomID &atom0, const AtomID &atom1,
                      const AtomID &atom2,
                      const SireUnits::Dimension::Angle &start,
                      const SireUnits::Dimension::Angle &end,
                      const SireCAS::Expression &mapping_function,
                      const PropertyMap &map = PropertyMap());
     
    AnglePerturbation(const AnglePerturbation &other);
    
    ~AnglePerturbation();
    
    static const char* typeName();
    
    AnglePerturbation& operator=(const AnglePerturbation &other);
    
    bool operator==(const AnglePerturbation &other) const;
    bool operator!=(const AnglePerturbation &other) const;
    
    QString toString() const;
    
    const AngleID& angle() const;
    
    const SireUnits::Dimension::Angle& start() const;
    const SireUnits::Dimension::Angle& end() const;

    bool wouldChange(const Molecule &molecule, const Values &values) const;

protected:
    void perturbMolecule(Mover<Molecule> &molecule, const Values &values) const;

private:
    /** The ID for the angle */
    AngleID angleid;
    
    /** The start and end sizes of the angle */
    SireUnits::Dimension::Angle start_size, end_size;
};

/** This perturbation moves a dihedral between two sizes.

    This uses the "anchors" property to anchor parts
    of the molecule, the "weight function" property
    to weight the motion of the parts of the molecule,
    and the "coordinates" property to get the coordinates
    to move

    @author Christopher Woods
*/
class SIREMOL_EXPORT DihedralPerturbation
        : public SireBase::ConcreteProperty<DihedralPerturbation,GeometryPerturbation>
{

friend QDataStream& ::operator<<(QDataStream&, const DihedralPerturbation&);
friend QDataStream& ::operator>>(QDataStream&, DihedralPerturbation&);

public:
    DihedralPerturbation();
    
    DihedralPerturbation(const DihedralID &dihedral, 
                         const SireUnits::Dimension::Angle &start,
                         const SireUnits::Dimension::Angle &end,
                         const PropertyMap &map = PropertyMap());
                     
    DihedralPerturbation(const DihedralID &dihedral, 
                         const SireUnits::Dimension::Angle &start,
                         const SireUnits::Dimension::Angle &end,
                         const SireCAS::Expression &mapping_function,
                         const PropertyMap &map = PropertyMap());
    
    DihedralPerturbation(const AtomID &atom0, const AtomID &atom1, 
                         const AtomID &atom2, const AtomID &atom3,
                         const SireUnits::Dimension::Angle &start,
                         const SireUnits::Dimension::Angle &end,
                         const PropertyMap &map = PropertyMap());
                     
    DihedralPerturbation(const AtomID &atom0, const AtomID &atom1,
                         const AtomID &atom2, const AtomID &atom3,
                         const SireUnits::Dimension::Angle &start,
                         const SireUnits::Dimension::Angle &end,
                         const SireCAS::Expression &mapping_function,
                         const PropertyMap &map = PropertyMap());
     
    DihedralPerturbation(const DihedralPerturbation &other);
    
    ~DihedralPerturbation();
    
    static const char* typeName();
    
    DihedralPerturbation& operator=(const DihedralPerturbation &other);
    
    bool operator==(const DihedralPerturbation &other) const;
    bool operator!=(const DihedralPerturbation &other) const;
    
    QString toString() const;
    
    const DihedralID& dihedral() const;
    
    const SireUnits::Dimension::Angle& start() const;
    const SireUnits::Dimension::Angle& end() const;

    bool wouldChange(const Molecule &molecule, const Values &values) const;

protected:
    void perturbMolecule(Mover<Molecule> &molecule, const Values &values) const;

private:
    /** The ID for the dihedral */
    DihedralID dihedralid;
    
    /** The start and end sizes of the angle */
    SireUnits::Dimension::Angle start_size, end_size;
};

}

Q_DECLARE_METATYPE( SireMol::NullGeometryPerturbation )
Q_DECLARE_METATYPE( SireMol::GeometryPerturbations )
Q_DECLARE_METATYPE( SireMol::BondPerturbation )
Q_DECLARE_METATYPE( SireMol::AnglePerturbation )
Q_DECLARE_METATYPE( SireMol::DihedralPerturbation )

SIRE_EXPOSE_CLASS( SireMol::GeometryPerturbation )
SIRE_EXPOSE_CLASS( SireMol::NullGeometryPerturbation )
SIRE_EXPOSE_CLASS( SireMol::GeometryPerturbations )
SIRE_EXPOSE_CLASS( SireMol::BondPerturbation )
SIRE_EXPOSE_CLASS( SireMol::AnglePerturbation )
SIRE_EXPOSE_CLASS( SireMol::DihedralPerturbation )

SIRE_END_HEADER

#endif
