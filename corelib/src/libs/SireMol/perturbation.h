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

#ifndef SIREMOL_PERTURBATION_H
#define SIREMOL_PERTURBATION_H

#include <QList>

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireCAS/expression.h"
#include "SireCAS/symbol.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Perturbation;
class NullPerturbation;

class Perturbations;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Perturbation&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Perturbation&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::NullPerturbation&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::NullPerturbation&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Perturbations&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Perturbations&);

namespace SireCAS
{
class Identities;
class Values;
}

namespace SireMol
{

class Molecule;
class MolEditor;

using SireBase::PropertyMap;

using SireCAS::Expression;
using SireCAS::Symbol;
using SireCAS::Values;

class SIREMOL_EXPORT PerturbationSymbols
{
public:
    PerturbationSymbols();
    ~PerturbationSymbols();

    const Symbol& lambda() const;
    const Symbol& initial() const;
    const Symbol& final() const;

private:
    Symbol lam, init, finl;
};

typedef SireBase::PropPtr<Perturbation> PerturbationPtr;

/** This is the base class of all perturbation objects. A Perturbation
    is a rule for changing a property of a molecule with respect
    to a driving (reaction) coordinate. Perturbations can be used
    to implement single topology free energy calculations
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT Perturbation : public SireBase::Property
{

friend class Perturbations;   // so can call perturbMolecule directly

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Perturbation&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, Perturbation&);

public:
    Perturbation();
    Perturbation(const Perturbation &other);
    
    virtual ~Perturbation();
    
    static const char* typeName();
    
    virtual Perturbation* clone() const=0;
    
    virtual PerturbationPtr recreate() const;
    virtual PerturbationPtr recreate(const SireCAS::Expression &mapping_function) const;
    virtual PerturbationPtr recreate(const PropertyMap &map) const;
    virtual PerturbationPtr recreate(const SireCAS::Expression &mapping_function,
                                     const PropertyMap &map) const;
    
    virtual PerturbationPtr substitute(const Symbol &old_symbol,
                                       const Symbol &new_symbol) const;
    
    virtual PerturbationPtr substitute(const SireCAS::Identities &identities) const;
    
    virtual QList<PerturbationPtr> children() const;
    
    virtual QSet<Symbol> requiredSymbols() const;
    virtual QSet<QString> requiredProperties() const=0;
    
    static const Expression& defaultFunction();
    static const PerturbationSymbols& symbols();
    
    const Expression& mappingFunction() const;
    
    const PropertyMap& propertyMap() const;
    
    Molecule perturb(const Molecule &molecule, const Values &values) const;

    virtual bool wouldChange(const Molecule &molecule, const Values &values) const=0;

    static const NullPerturbation& null();

protected:
    Perturbation(const PropertyMap &map);
    Perturbation(const Expression &equation, const PropertyMap &map = PropertyMap());

    Perturbation& operator=(const Perturbation &other);
    
    virtual void perturbMolecule(MolEditor &molecule, const Values &values) const=0;
    
    bool operator==(const Perturbation &other) const;
    bool operator!=(const Perturbation &other) const;

private:
    /** The equation used to control the perturbation */
    Expression mapping_eqn;
    
    /** The property map used to find the properties 
        used in this perturbation */
    PropertyMap map;
};

/** This is a null perturbation that does nothing */
class SIREMOL_EXPORT NullPerturbation
        : public SireBase::ConcreteProperty<NullPerturbation,Perturbation>
{
friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const NullPerturbation&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, NullPerturbation&);

public:
    NullPerturbation();
    NullPerturbation(const NullPerturbation &other);
    
    ~NullPerturbation();
    
    static const char* typeName();
    
    NullPerturbation& operator=(const NullPerturbation &other);
    
    bool operator==(const NullPerturbation &other) const;
    bool operator!=(const NullPerturbation &other) const;
    
    QSet<Symbol> requiredSymbols() const;
    QSet<QString> requiredProperties() const;
    
    bool wouldChange(const Molecule &molecule, const Values &values) const;
    
protected:
    void perturbMolecule(MolEditor &molecule, const Values &values) const;
};

/** This class holds a collection of perturbations that can
    be applied to a molecule
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT Perturbations 
        : public SireBase::ConcreteProperty<Perturbations,Perturbation>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Perturbations&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, Perturbations&);

public:
    Perturbations();
    Perturbations(const Perturbation &perturbation);
    Perturbations(const QList<PerturbationPtr> &perturbations);
    
    Perturbations(const Perturbations &other);
    
    ~Perturbations();
    
    Perturbations& operator=(const Perturbations &other);
    
    bool operator==(const Perturbations &other) const;
    bool operator!=(const Perturbations &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    QList<PerturbationPtr> perturbations() const;

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
    void perturbMolecule(MolEditor &molecule, const Values &values) const;

private:
    void makeSane();
    
    /** All of the perturbations used in this object */
    QList<PerturbationPtr> perts;
};

}

Q_DECLARE_METATYPE( SireMol::NullPerturbation )
Q_DECLARE_METATYPE( SireMol::Perturbations )

SIRE_EXPOSE_CLASS( SireMol::Perturbation )
SIRE_EXPOSE_CLASS( SireMol::NullPerturbation )
SIRE_EXPOSE_CLASS( SireMol::Perturbations )
SIRE_EXPOSE_CLASS( SireMol::PerturbationSymbols )

SIRE_EXPOSE_PROPERTY( SireMol::PerturbationPtr, SireMol::Perturbation )

SIRE_END_HEADER

#endif
