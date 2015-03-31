/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREMM_INTERNALPERTURBATION_H
#define SIREMM_INTERNALPERTURBATION_H

#include "SireMol/perturbation.h"
#include "SireMol/atomidx.h"
#include "SireMol/atomidentifier.h"

#include "SireCAS/expression.h"
#include "SireCAS/identities.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class InternalPerturbation;
class TwoAtomPerturbation;
class ThreeAtomPerturbation;
class FourAtomPerturbation;
}

QDataStream& operator<<(QDataStream&, const SireMM::InternalPerturbation&);
QDataStream& operator>>(QDataStream&, SireMM::InternalPerturbation&);

QDataStream& operator<<(QDataStream&, const SireMM::TwoAtomPerturbation&);
QDataStream& operator>>(QDataStream&, SireMM::TwoAtomPerturbation&);

QDataStream& operator<<(QDataStream&, const SireMM::ThreeAtomPerturbation&);
QDataStream& operator>>(QDataStream&, SireMM::ThreeAtomPerturbation&);

QDataStream& operator<<(QDataStream&, const SireMM::FourAtomPerturbation&);
QDataStream& operator>>(QDataStream&, SireMM::FourAtomPerturbation&);

namespace SireMM
{

using SireBase::PropertyMap;

using SireMol::Perturbation;
using SireMol::AtomID;

/** This is the base class of all internal perturbations - these
    are perturbations that change the internal potential of a
    molecule (e.g. the bond, angle and dihedral parameters) 

    Internal perturbation work by applying the mapping function
    to the identities in initialForm() and finalForm() and then
    substituting these identities into baseExpression(), e.g.
    
    initialForm() :=>  k = k_{initial},  r0 = r0_{initial}
    finalForm()   :=>  k = k_{final},    r0 = r0_{final}
    
    baseExpression() :=>  k * (r - r0)**2
    
    mappingFunction() :=>  (1-lam)*initial + lam*final
    
    This will result in the the perturbing function being
    
    perturbFunction() :=> [ (1-lam)*k_{initial} + lam*k_{final} ] *
                            (r - [(1-lam)*r0_{initial} + lam*r0_{final}])**2
                            
    equally, if
    
    initialForm() :=> f = 3 * cos(5 phi)
    finalForm()   :=> f = 5 * cos(8 phi)
    
    baseExpression() :=>  f
    
    mappingFunction() :=> (1-lam)*initial + lam*final
    
    This will result in
    
    perturbFunction() :=> (1-lam)*(3 * cos(5 phi)) + lam * (5 * cos(8 phi))
    
    @author Christopher Woods
*/
class SIREMM_EXPORT InternalPerturbation : public Perturbation
{

friend QDataStream& ::operator<<(QDataStream&, const InternalPerturbation&);
friend QDataStream& ::operator>>(QDataStream&, InternalPerturbation&);

public:
    InternalPerturbation();
    
    InternalPerturbation(const SireCAS::Expression &initial_function,
                         const SireCAS::Expression &final_function,
                         const SireCAS::Expression &mapping_function,
                         const PropertyMap &map);
    
    InternalPerturbation(const SireCAS::Expression &base_expression,
                         const SireCAS::Identities &initial_forms,
                         const SireCAS::Identities &final_forms,
                         const SireCAS::Expression &mapping_function,
                         const PropertyMap &map);
    
    InternalPerturbation(const InternalPerturbation &other);
    
    ~InternalPerturbation();
    
    static const char* typeName();
    
    const SireCAS::Expression& baseExpression() const;
    const SireCAS::Expression& perturbExpression() const;
    
    const SireCAS::Identities& initialForms() const;
    const SireCAS::Identities& finalForms() const;
    
    SireMol::PerturbationPtr recreate(const SireCAS::Expression &expression) const;
    SireMol::PerturbationPtr recreate(const SireCAS::Expression &expression,
                                      const PropertyMap &map) const;
    
    SireMol::PerturbationPtr substitute(const SireCAS::Identities &identities) const;
    
protected:
    InternalPerturbation& operator=(const InternalPerturbation &other);
    
    bool operator==(const InternalPerturbation &other) const;
    bool operator!=(const InternalPerturbation &other) const;

private:
    void buildPerturbExpression();

    /** The base expression - this is combined with the mapping function,
        and identities in 'initial_form' and 'final_form' to give the
        function that is perturbed */
    SireCAS::Expression base_expression;
    
    /** Here is the function that is perturbed */
    SireCAS::Expression perturb_expression;
    
    /** The initial identities */
    SireCAS::Identities initial_forms;
    
    /** The final identities */
    SireCAS::Identities final_forms;
};

/** This class represents a perturbation that maps the two-atom potential
    function using a perturbation function
  
    For example, the perturbation function for a bond could be;
    
    E_{r,lambda} = [ (1-lambda) k_b + lambda k_f ] * 
                      [ ((1-lambda) r0_b + lambda r0_f) - r ]^2
                      
    The perturbation will insert the value of lambda into this 
    expression and set the molecules bond function to the resulting
    expression, e.g at lambda=0
    
    E_{r,0} = k_b * (r0_b - r)^2
    
    and at lambda=1
    
    E_{r,1} = k_f * (r0_f - r)^2
        
    @author Christopher Woods
*/
class SIREMM_EXPORT TwoAtomPerturbation
        : public SireBase::ConcreteProperty<TwoAtomPerturbation,InternalPerturbation>
{

friend QDataStream& ::operator<<(QDataStream&, const TwoAtomPerturbation&);
friend QDataStream& ::operator>>(QDataStream&, TwoAtomPerturbation&);

public:
    TwoAtomPerturbation();

    TwoAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                        const SireCAS::Expression &initial_form,
                        const SireCAS::Expression &final_form,
                        const PropertyMap &map = PropertyMap());
    TwoAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                        const SireCAS::Expression &initial_form,
                        const SireCAS::Expression &final_form,
                        const SireCAS::Expression &mapping_function,
                        const PropertyMap &map = PropertyMap());

    TwoAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                        const SireCAS::Expression &base_expression,
                        const SireCAS::Identities &initial_forms,
                        const SireCAS::Identities &final_forms,
                        const PropertyMap &map = PropertyMap());
    TwoAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                        const SireCAS::Expression &base_expression,
                        const SireCAS::Identities &initial_forms,
                        const SireCAS::Identities &final_forms,
                        const SireCAS::Expression &mapping_function,
                        const PropertyMap &map = PropertyMap());
                        
    TwoAtomPerturbation(const TwoAtomPerturbation &other);
    
    ~TwoAtomPerturbation();
    
    TwoAtomPerturbation& operator=(const TwoAtomPerturbation &other);
    
    bool operator==(const TwoAtomPerturbation &other) const;
    bool operator!=(const TwoAtomPerturbation &other) const;
    
    static const char* typeName();
    
    const AtomID& atom0() const;
    const AtomID& atom1() const;
    
    QString toString() const;

    QSet<QString> requiredProperties() const;
    
    bool wouldChange(const SireMol::Molecule &molecule, 
                     const SireCAS::Values &values) const;

protected:
    void perturbMolecule(SireMol::MolEditor &molecule, 
                         const SireCAS::Values &values) const;
    
private:
    /** The identifiers of the two atoms */
    SireMol::AtomIdentifier atm0, atm1;
};

/** This class represents a perturbation that maps the three-atom potential
    function using a perturbation function
  
    For example, the perturbation function for an angle could be;
    
    E_{theta,lambda} = [ (1-lambda) k_b + lambda k_f ] * 
                       [ ((1-lambda) theta0_b + lambda theta0_f) - theta ]^2
                      
    The perturbation will insert the value of lambda into this 
    expression and set the molecules angle function to the resulting
    expression, e.g at lambda=0
    
    E_{theta,0} = k_b * (theta0_b - theta)^2
    
    and at lambda=1
    
    E_{theta,1} = k_f * (theta0_f - theta)^2
        
    @author Christopher Woods
*/
class SIREMM_EXPORT ThreeAtomPerturbation
        : public SireBase::ConcreteProperty<ThreeAtomPerturbation,InternalPerturbation>
{

friend QDataStream& ::operator<<(QDataStream&, const ThreeAtomPerturbation&);
friend QDataStream& ::operator>>(QDataStream&, ThreeAtomPerturbation&);

public:
    ThreeAtomPerturbation();

    ThreeAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                          const AtomID &atom2,
                          const SireCAS::Expression &initial_form,
                          const SireCAS::Expression &final_form,
                          const PropertyMap &map = PropertyMap());
    ThreeAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                          const AtomID &atom2,
                          const SireCAS::Expression &initial_form,
                          const SireCAS::Expression &final_form,
                          const SireCAS::Expression &mapping_function,
                          const PropertyMap &map = PropertyMap());

    ThreeAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                          const AtomID &atom2,
                          const SireCAS::Expression &base_expression,
                          const SireCAS::Identities &initial_forms,
                          const SireCAS::Identities &final_forms,
                          const PropertyMap &map = PropertyMap());
    ThreeAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                          const AtomID &atom2,
                          const SireCAS::Expression &base_expression,
                          const SireCAS::Identities &initial_forms,
                          const SireCAS::Identities &final_forms,
                          const SireCAS::Expression &mapping_function,
                          const PropertyMap &map = PropertyMap());
                        
    ThreeAtomPerturbation(const ThreeAtomPerturbation &other);
    
    ~ThreeAtomPerturbation();
    
    ThreeAtomPerturbation& operator=(const ThreeAtomPerturbation &other);
    
    bool operator==(const ThreeAtomPerturbation &other) const;
    bool operator!=(const ThreeAtomPerturbation &other) const;
    
    static const char* typeName();
    
    const AtomID& atom0() const;
    const AtomID& atom1() const;
    const AtomID& atom2() const;
    
    QString toString() const;

    QSet<QString> requiredProperties() const;
    
    bool wouldChange(const SireMol::Molecule &molecule, 
                     const SireCAS::Values &values) const;

protected:
    void perturbMolecule(SireMol::MolEditor &molecule, 
                         const SireCAS::Values &values) const;
    
private:
    /** The identifiers of the three atoms */
    SireMol::AtomIdentifier atm0, atm1, atm2;
};

/** This class represents a perturbation that maps the four-atom potential
    function using a perturbation function
  
    For example, the perturbation function for a dihedral could be;
    
    E_{phi,lambda} = (1-lambda)(cos(5 phi)) + lambda (cos(3 phi))
                      
    The perturbation will insert the value of lambda into this 
    expression and set the molecules bond function to the resulting
    expression, e.g at lambda=0
    
    E_{phi,0} = cos(5 phi)
    
    and at lambda=1
    
    E_{phi,1} = cos(3 phi)
        
    @author Christopher Woods
*/
class SIREMM_EXPORT FourAtomPerturbation
        : public SireBase::ConcreteProperty<FourAtomPerturbation,InternalPerturbation>
{

friend QDataStream& ::operator<<(QDataStream&, const FourAtomPerturbation&);
friend QDataStream& ::operator>>(QDataStream&, FourAtomPerturbation&);

public:
    FourAtomPerturbation();

    FourAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                         const AtomID &atom2, const AtomID &atom3,
                         const SireCAS::Expression &initial_form,
                         const SireCAS::Expression &final_form,
                         const PropertyMap &map = PropertyMap());
    FourAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                         const AtomID &atom2, const AtomID &atom3,
                         const SireCAS::Expression &initial_form,
                         const SireCAS::Expression &final_form,
                         const SireCAS::Expression &mapping_function,
                         const PropertyMap &map = PropertyMap());

    FourAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                         const AtomID &atom2, const AtomID &atom3,
                         const SireCAS::Expression &base_expression,
                         const SireCAS::Identities &initial_forms,
                         const SireCAS::Identities &final_forms,
                         const PropertyMap &map = PropertyMap());
    FourAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                         const AtomID &atom2, const AtomID &atom3,
                         const SireCAS::Expression &base_expression,
                         const SireCAS::Identities &initial_forms,
                         const SireCAS::Identities &final_forms,
                         const SireCAS::Expression &mapping_function,
                         const PropertyMap &map = PropertyMap());
                        
    FourAtomPerturbation(const FourAtomPerturbation &other);
    
    ~FourAtomPerturbation();
    
    FourAtomPerturbation& operator=(const FourAtomPerturbation &other);
    
    bool operator==(const FourAtomPerturbation &other) const;
    bool operator!=(const FourAtomPerturbation &other) const;
    
    static const char* typeName();
    
    const AtomID& atom0() const;
    const AtomID& atom1() const;
    const AtomID& atom2() const;
    const AtomID& atom3() const;
    
    QString toString() const;

    QSet<QString> requiredProperties() const;
    
    bool wouldChange(const SireMol::Molecule &molecule, 
                     const SireCAS::Values &values) const;

protected:
    void perturbMolecule(SireMol::MolEditor &molecule, 
                         const SireCAS::Values &values) const;
    
private:
    /** The identifiers of the four atoms */
    SireMol::AtomIdentifier atm0, atm1, atm2, atm3;
};

}

Q_DECLARE_METATYPE( SireMM::TwoAtomPerturbation )
Q_DECLARE_METATYPE( SireMM::ThreeAtomPerturbation )
Q_DECLARE_METATYPE( SireMM::FourAtomPerturbation )

SIRE_EXPOSE_CLASS( SireMM::InternalPerturbation )
SIRE_EXPOSE_CLASS( SireMM::TwoAtomPerturbation )
SIRE_EXPOSE_CLASS( SireMM::ThreeAtomPerturbation )
SIRE_EXPOSE_CLASS( SireMM::FourAtomPerturbation )

SIRE_END_HEADER

#endif
