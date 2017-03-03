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

#ifndef SIREMM_INTERNALPARAMETERS_H
#define SIREMM_INTERNALPARAMETERS_H

#include <QVector>

#include "twoatomfunctions.h"
#include "threeatomfunctions.h"
#include "fouratomfunctions.h"

#include "SireBase/refcountdata.h"
#include "SireBase/shareddatapointer.hpp"

#include "SireMol/cgidx.h"

#include "SireFF/detail/atomiccoords3d.h"

SIRE_BEGIN_HEADER

namespace SireMM
{

class InternalParameters;
class InternalParameters3D;

class GroupInternalParameters;

}

QDataStream& operator<<(QDataStream&, const SireMM::InternalParameters&);
QDataStream& operator>>(QDataStream&, SireMM::InternalParameters&);

QDataStream& operator<<(QDataStream&, const SireMM::InternalParameters3D&);
QDataStream& operator>>(QDataStream&, SireMM::InternalParameters3D&);

QDataStream& operator<<(QDataStream&, const SireMM::GroupInternalParameters&);
QDataStream& operator>>(QDataStream&, SireMM::GroupInternalParameters&);

namespace SireBase
{
class PropertyName;
}

namespace SireMol
{
class AtomSelection;
class PartialMolecule;
}

namespace SireMM
{

using SireVol::CoordGroupArray;

using SireMol::CGIdx;
using SireMol::AtomSelection;
using SireMol::PartialMolecule;

using SireBase::PropertyName;

namespace detail
{

/** This internal class holds the ID of four CutGroups */
class CGIDQuad
{
public:
    CGIDQuad();
    
    CGIDQuad(CGIdx cgidx0);
    CGIDQuad(CGIdx cgidx0, CGIdx cgidx1);
    CGIDQuad(CGIdx cgidx0, CGIdx cgidx1, CGIdx cgidx2);
    CGIDQuad(CGIdx cgidx0, CGIdx cgidx1, CGIdx cgidx2, CGIdx cgidx3);
    
    CGIDQuad(const CGIDQuad &other);
    
    ~CGIDQuad();
    
    CGIDQuad& operator=(const CGIDQuad &other);
    
    bool operator==(const CGIDQuad &other) const;
    bool operator!=(const CGIDQuad &other) const;

    bool isSingleCutGroup() const;
    bool isDoubleCutGroup() const;
    bool isTripleCutGroup() const;
    bool isQuadrupleCutGroup() const;
    
    CGIdx cgIdx0() const;
    CGIdx cgIdx1() const;
    CGIdx cgIdx2() const;
    CGIdx cgIdx3() const;

    uint hash() const;

private:
    CGIdx cgidxs[4];
};

inline uint qHash(const CGIDQuad &idquad)
{
    return idquad.hash();
}

/** Internal class used to hold the non-physical parameters
    (impropers and Urey-Bradley) */
class GroupInternalNonPhysParameters : public SireBase::RefCountData
{
public:
    GroupInternalNonPhysParameters();
    GroupInternalNonPhysParameters(const GroupInternalNonPhysParameters &other);
    
    ~GroupInternalNonPhysParameters();
    
    GroupInternalNonPhysParameters& operator=(
                              const GroupInternalNonPhysParameters &other);

    bool operator==(const GroupInternalNonPhysParameters &other) const;
    bool operator!=(const GroupInternalNonPhysParameters &other) const;
    
    /** The array of improper parameters */
    QVector<FourAtomFunction> improper_params;

    /** The improper forces, acting on theta and phi */
    QVector<FourAtomFunction> improper_theta_forces, improper_phi_forces;

    /** The array of Urey-Bradley parameters and forces */
    QVector<TwoAtomFunction> ub_params, ub_forces;
};

/** Internal class used to hold the cross-term parameters
    (stretch-stretch, stretch-bend etc.) */
class GroupInternalCrossParameters : public SireBase::RefCountData
{
public:
    GroupInternalCrossParameters();
    GroupInternalCrossParameters(const GroupInternalCrossParameters &other);
    
    ~GroupInternalCrossParameters();
    
    GroupInternalCrossParameters& operator=(
                            const GroupInternalCrossParameters &other);

    bool operator==(const GroupInternalCrossParameters &other) const;
    bool operator!=(const GroupInternalCrossParameters &other) const;

    /** The array of stretch-stretch parameters */
    QVector<ThreeAtomFunction> stretch_stretch_params;
    
    /** The stretch-stretch forces */
    QVector<ThreeAtomFunction> stretch_stretch_r01_forces;
    /** The stretch-stretch forces */
    QVector<ThreeAtomFunction> stretch_stretch_r21_forces;

    /** The array of stretch-bend parameters */
    QVector<ThreeAtomFunction> stretch_bend_params;

    /** The stretch-bend forces */
    QVector<ThreeAtomFunction> stretch_bend_theta_forces;
    /** The stretch-bend forces */
    QVector<ThreeAtomFunction> stretch_bend_r01_forces;
    /** The stretch-bend forces */
    QVector<ThreeAtomFunction> stretch_bend_r21_forces;
    
    /** The array of bend-bend parameters */
    QVector<FourAtomFunction> bend_bend_params;

    /** The bend-bend forces */
    QVector<FourAtomFunction> bend_bend_theta012_forces;
    /** The bend-bend forces */
    QVector<FourAtomFunction> bend_bend_theta213_forces;
    /** The bend-bend forces */
    QVector<FourAtomFunction> bend_bend_theta310_forces;

    /** The array of stretch-bend-torsion parameters */
    QVector<FourAtomFunction> stretch_bend_torsion_params;
    
    /** The stretch-bend-torsion forces */
    QVector<FourAtomFunction> stretch_bend_torsion_phi_forces;
    /** The stretch-bend-torsion forces */
    QVector<FourAtomFunction> stretch_bend_torsion_r01_forces;
    /** The stretch-bend-torsion forces */
    QVector<FourAtomFunction> stretch_bend_torsion_r12_forces;
    /** The stretch-bend-torsion forces */
    QVector<FourAtomFunction> stretch_bend_torsion_r32_forces;
    /** The stretch-bend-torsion forces */
    QVector<FourAtomFunction> stretch_bend_torsion_r03_forces;
    /** The stretch-bend-torsion forces */
    QVector<FourAtomFunction> stretch_bend_torsion_theta012_forces;
    /** The stretch-bend-torsion forces */
    QVector<FourAtomFunction> stretch_bend_torsion_theta321_forces;
};

/** Internal class used to hold the private data of GroupInternalParameters */
class GroupInternalParametersData : public SireBase::RefCountData
{
public:
    GroupInternalParametersData();
    GroupInternalParametersData(const CGIDQuad &cgids);
    
    GroupInternalParametersData(const GroupInternalParametersData &other);
    
    ~GroupInternalParametersData();
    
    GroupInternalParametersData& operator=(
                            const GroupInternalParametersData &other);

    bool operator==(const GroupInternalParametersData &other) const;
    bool operator!=(const GroupInternalParametersData &other) const;

    bool hasCrossTerms() const;
    bool hasNonPhysicalParameters() const;

    /** ID of the (up to) four CutGroups */
    CGIDQuad idquad;

    /** The array of bond parameters and forces */
    QVector<TwoAtomFunction> bond_params, bond_forces;

    /** The array of angle parameters and forces */
    QVector<ThreeAtomFunction> angle_params, angle_forces;

    /** The array of dihedral parameters and forces */
    QVector<FourAtomFunction> dihedral_params, dihedral_forces;

    /** Shared pointer to the non-physical terms (impropers
        and Urey-Bradley) - no all groups will have these
        terms, so it is best to hide them to save space */
    SireBase::SharedDataPointer<GroupInternalNonPhysParameters> nonphys_terms;
    
    /** Shared pointer to the cross-term internal parameters
        (double-pointer so as to save memory, as cross-terms
        tend not to be used with very large molecules) */
    SireBase::SharedDataPointer<GroupInternalCrossParameters> cross_terms;
};

}

/** This class holds all of the internal parameters for one group
    combination within a molecule
    
    There are several types of internal parameters, defined
    by the type of internal used to provide the coordinates,
    and the quantities calculated from those coordinates that
    can be used in the function
    
    Bond            : Input Bond - function uses interatomic distance (1-2), r
    Angle           : Input Angle - function uses angle (1-2-3), theta
    Dihedral        : Input Dihedral - function uses torsion (1-2-3-4), phi
    
    Improper        : Input Improper - function uses either torsion angle
                                       (1-3-4-2), phi, or out of plane angle, theta
                                       
    Urey-Bradley    : Input Angle - function uses distance (1-3), r
    
    Stretch-Stretch : Input Angle - function uses distances (1-2), r12, and 
                                    (3-2), r32
                                    
    Stretch-Bend    : Input Angle - function uses distances angle (1-2-3), theta,
                                    and distances (1-2), r12, and (3-2) r32
                                    
    Bend-Bend       : Input Improper - function uses angles (1-2-3), (3-2-4), (4-2-1),
                                        theta123, theta324, theta421 
    
    Stretch-Bend    : Input Dihedral - function uses torsion (1-2-3-4), phi,
      -Torsion                         distances (1-2), (2-3), (3-4), (1-4)
                                       r12, r23, r34, r14 and angles
                                       (1-2-3) and (2-3-4), theta123, theta234
    
    @author Christopher Woods
*/
class SIREMM_EXPORT GroupInternalParameters
{

friend QDataStream& ::operator<<(QDataStream&, const GroupInternalParameters&);
friend QDataStream& ::operator>>(QDataStream&, GroupInternalParameters&);

friend class InternalParameters;  // so can call editing functions

public:
    GroupInternalParameters();
    
    GroupInternalParameters(const GroupInternalParameters &other);
    
    ~GroupInternalParameters();
    
    GroupInternalParameters& operator=(const GroupInternalParameters &other);

    bool operator==(const GroupInternalParameters &other) const;
    bool operator!=(const GroupInternalParameters &other) const;

    bool isEmpty() const;

    bool hasPhysicalParameters() const;
    bool hasNonPhysicalParameters() const;
    bool hasCrossTerms() const;
    
    bool isSingleCutGroup() const;
    bool isDoubleCutGroup() const;
    bool isTripleCutGroup() const;
    bool isQuadrupleCutGroup() const;

    bool isSingleCutGroup(CGIdx cgidx0) const;
    bool isDoubleCutGroup(CGIdx cgidx0, CGIdx cgidx1) const;
    bool isTripleCutGroup(CGIdx cgidx0, CGIdx cgidx1, CGIdx cgidx2) const;
    bool isQuadrupleCutGroup(CGIdx cgidx0, CGIdx cgidx1,
                             CGIdx cgidx2, CGIdx cgidx3) const;

    bool refersTo(CGIdx cgidx) const;
    bool refersTo(const QSet<CGIdx> &cgidxs) const;

    CGIdx cgIdx0() const;
    CGIdx cgIdx1() const;
    CGIdx cgIdx2() const;
    CGIdx cgIdx3() const;

    const QVector<TwoAtomFunction>& bondPotential() const;
    const QVector<TwoAtomFunction>& bondForces() const;
    
    const QVector<ThreeAtomFunction>& anglePotential() const;
    const QVector<ThreeAtomFunction>& angleForces() const;
    
    const QVector<FourAtomFunction>& dihedralPotential() const;
    const QVector<FourAtomFunction>& dihedralForces() const;
    
    const QVector<FourAtomFunction>& improperPotential() const;
    const QVector<FourAtomFunction>& improper_Theta_Forces() const;
    const QVector<FourAtomFunction>& improper_Phi_Forces() const;
    
    const QVector<TwoAtomFunction>& ureyBradleyPotential() const;
    const QVector<TwoAtomFunction>& ureyBradleyForces() const;
    
    const QVector<ThreeAtomFunction>& stretchStretchPotential() const;
    const QVector<ThreeAtomFunction>& stretchStretch_R01_Forces() const;
    const QVector<ThreeAtomFunction>& stretchStretch_R21_Forces() const;
    
    const QVector<ThreeAtomFunction>& stretchBendPotential() const;
    const QVector<ThreeAtomFunction>& stretchBend_Theta_Forces() const;
    const QVector<ThreeAtomFunction>& stretchBend_R01_Forces() const;
    const QVector<ThreeAtomFunction>& stretchBend_R21_Forces() const;
    
    const QVector<FourAtomFunction>& bendBendPotential() const;
    const QVector<FourAtomFunction>& bendBend_Theta012_Forces() const;
    const QVector<FourAtomFunction>& bendBend_Theta213_Forces() const;
    const QVector<FourAtomFunction>& bendBend_Theta310_Forces() const;
    
    const QVector<FourAtomFunction>& stretchBendTorsionPotential() const;
    const QVector<FourAtomFunction>& stretchBendTorsion_Phi_Forces() const;
    const QVector<FourAtomFunction>& stretchBendTorsion_R01_Forces() const;
    const QVector<FourAtomFunction>& stretchBendTorsion_R12_Forces() const;
    const QVector<FourAtomFunction>& stretchBendTorsion_R32_Forces() const;
    const QVector<FourAtomFunction>& stretchBendTorsion_R03_Forces() const;
    const QVector<FourAtomFunction>& stretchBendTorsion_Theta012_Forces() const;
    const QVector<FourAtomFunction>& stretchBendTorsion_Theta321_Forces() const;

protected:
    //these editing functions can be called only by InternalParameters
    GroupInternalParameters(const detail::CGIDQuad &cgids);

    void setBondPotential(const QVector<TwoAtomFunction> &potential,
                          const QVector<TwoAtomFunction> &forces);
                           
    void setAnglePotential(const QVector<ThreeAtomFunction> &potential,
                           const QVector<ThreeAtomFunction> &forces);
    
    void setDihedralPotential(const QVector<FourAtomFunction> &potential,
                              const QVector<FourAtomFunction> &forces);
    
    void setImproperPotential(const QVector<FourAtomFunction> &potential,
                              const QVector<FourAtomFunction> &theta_forces,
                              const QVector<FourAtomFunction> &phi_forces);
                               
    void setUreyBradleyPotential(const QVector<TwoAtomFunction> &potential,
                                 const QVector<TwoAtomFunction> &forces);
    
    void setStretchStretchPotential(const QVector<ThreeAtomFunction> &potential,
                                    const QVector<ThreeAtomFunction> &r01_forces,
                                    const QVector<ThreeAtomFunction> &r21_forces);
                                     
    void setStretchBendPotential(const QVector<ThreeAtomFunction> &potential,
                                 const QVector<ThreeAtomFunction> &theta_forces,
                                 const QVector<ThreeAtomFunction> &r01_forces,
                                 const QVector<ThreeAtomFunction> &r21_forces);
                                  
    void setBendBendPotential(const QVector<FourAtomFunction> &potential,
                              const QVector<FourAtomFunction> &theta012_forces,
                              const QVector<FourAtomFunction> &theta213_forces,
                              const QVector<FourAtomFunction> &theta310_forces);

    void setStretchBendTorsionPotential(const QVector<FourAtomFunction> &potential,
                                     const QVector<FourAtomFunction> &phi_forces,
                                     const QVector<FourAtomFunction> &r01_forces,
                                     const QVector<FourAtomFunction> &r12_forces,
                                     const QVector<FourAtomFunction> &r32_forces,
                                     const QVector<FourAtomFunction> &theta012_forces,
                                     const QVector<FourAtomFunction> &theta321_forces);

private:
    /** Implicitly shared pointer to the data of this class */
    SireBase::SharedDataPointer<detail::GroupInternalParametersData> d;
};

/** This is the base of all of the classes used to hold the symbols
    for the internal parameters */
class SIREMM_EXPORT InternalSymbolsBase
{
public:
    ~InternalSymbolsBase();

    operator QSet<Symbol>() const
    {
        return symbols;
    }

protected:
    InternalSymbolsBase();
    
    QSet<Symbol> symbols;
};

/** This class holds the symbols required for the bond 
    and Urey-Bradley parameters. These are functions
    that act between two atoms, using the distance
    between the atoms (r) as the input variable */
class SIREMM_EXPORT BondSymbols : public InternalSymbolsBase
{
public:
    BondSymbols();
    ~BondSymbols();
    
    const Symbol& r() const;
    
private:
    /** The symbol for the length of a bond */
    Symbol r_;
};

/** This class holds the symbols required for the angle parameters.
    These are functions of the angle 0-1-2 (theta) of three
    atoms (where atom 1 is the central atom of the angle) */
class SIREMM_EXPORT AngleSymbols : public InternalSymbolsBase
{
public:
    AngleSymbols();
    ~AngleSymbols();
    
    const Symbol& theta() const;
    
private:
    /** The symbol for the size of the angle */
    Symbol theta_;
};

/** This class holds the symbols required for the dihedral parameters.
    These are functions of the dihedral (0-1-2-3) (phi) made 
    by four atoms, of atoms 0 and 3 about the bond between
    atoms 1 and 2 (measured clockwise) */
class SIREMM_EXPORT DihedralSymbols : public InternalSymbolsBase
{
public:
    DihedralSymbols();
    ~DihedralSymbols();
    
    const Symbol& phi() const;
    
private:
    /** The symbol for the size of the torsion */
    Symbol phi_;
};

/** This class holds the symbols required for the improper parameters.
    These are functions of the dihedral (0-1-2-3) (phi) and also
    of the improper angle made between the bond 0-1 and the plane
    formed by atoms 1-2-3 (theta) */
class SIREMM_EXPORT ImproperSymbols : public InternalSymbolsBase
{
public:
    ImproperSymbols();
    ~ImproperSymbols();
    
    const Symbol& theta() const;
    const Symbol& phi() const;
    
private:
    /** The symbol for the angle made by the improper bond
        with the plane of the other atoms */
    Symbol theta_;
    
    /** The symbol for the size of the improper torsion */
    Symbol phi_;
};

/** This class holds the symbols required for the stretch-stretch parameters.
    These are functions of the bond distances among three atoms, 0 1 and 2,
    so distances between atoms 0-1 (r01), 2-1 (r21) and 1-2 (r12) */
class SIREMM_EXPORT StretchStretchSymbols : public InternalSymbolsBase
{
public:
    StretchStretchSymbols();
    ~StretchStretchSymbols();
    
    const Symbol& r01() const;
    const Symbol& r21() const;
    const Symbol& r12() const;
    
private:
    /** The symbol for the length of the bond 0->1 */
    Symbol r01_;
    
    /** The symbol for the length of the bond 1<-2 */
    Symbol r21_;
    
    /** The symbol for the length of the bond 1->2 */
    Symbol r12_;
};

/** This class holds the symbols required for the stretch-bend parameters.
    These are a function of the angle of atoms 0-1-2 (theta) and
    the distances between each pair of atoms (r01, r21 and r12) */
class SIREMM_EXPORT StretchBendSymbols : public InternalSymbolsBase
{
public:
    StretchBendSymbols();
    ~StretchBendSymbols();
    
    const Symbol& theta() const;
    const Symbol& r01() const;
    const Symbol& r21() const;
    const Symbol& r12() const;
    
private:
    /** The symbol for the size of the angle */
    Symbol theta_;

    /** The symbol for the length of the bond 0->1 */
    Symbol r01_;
    
    /** The symbol for the length of the bond 1<-2 */
    Symbol r21_;
    
    /** The symbol for the length of the bond 1->2 */
    Symbol r12_;
};

/** This class holds the symbols required for the bend-bend parameters.
    These are functions of the three angles within the four atoms
    0, 1, 2, 3, where atom 1 is in the middle. The angles are
    therefore 0-1-2 (theta012), 2-1-3 (theta213) and 3-1-0 (theta310) */
class SIREMM_EXPORT BendBendSymbols : public InternalSymbolsBase
{
public:
    BendBendSymbols();
    ~BendBendSymbols();
    
    const Symbol& theta012() const;
    const Symbol& theta213() const;
    const Symbol& theta310() const;
    
private:
    /** The symbol for the size of the angle 0->1<-2 */
    Symbol theta012_;
    /** The symbol for the size of the angle 2->1<-3 */
    Symbol theta213_;
    /** The symbol for the size of the angle 3->1<-0 */
    Symbol theta310_;
};

/** This class holds the symbols required for the stretch-bend-torsion parameters
    These are functions over four atoms, 0, 1, 2, 3 and 4, of the dihedral
    formed over the four atoms 0-1-2-3 (phi), the angles 0-1-2 (theta012)
    and 3-2-1 (theta321) and the distances between atoms, 0-1 (r01),
    1-2 (r12), 3-2 (r32) and 0-3 (r03)
*/
class SIREMM_EXPORT StretchBendTorsionSymbols : public InternalSymbolsBase
{
public:
    StretchBendTorsionSymbols();
    ~StretchBendTorsionSymbols();
    
    const Symbol& phi() const;
    
    const Symbol& theta012() const;
    const Symbol& theta321() const;
    
    const Symbol& r01() const;
    const Symbol& r12() const;
    const Symbol& r32() const;
    const Symbol& r03() const;
    
private:
    /** The symbol for the size of the torsion */
    Symbol phi_;

    /** The symbol for the angle 0->1<-2 */
    Symbol theta012_;
    /** The symbol for the angle 3->2<-1 */
    Symbol theta321_;
    
    /** The symbol for the length of the bond 0->1 */
    Symbol r01_;
    /** The symbol for the length of the bond 1->2 */
    Symbol r12_;
    /** The symbol for the length of the bond 3->2 */
    Symbol r32_;
    /** The symbol for the distance between the 0 and 3 atoms */
    Symbol r03_;
};

/** This class holds all of the symbols used by all of the
    internal parameters */
class SIREMM_EXPORT InternalSymbols : public InternalSymbolsBase
{
public:
    InternalSymbols();
    ~InternalSymbols();

    const BondSymbols& bond() const;
    const AngleSymbols& angle() const;
    const DihedralSymbols& dihedral() const;
    
    const ImproperSymbols& improper() const;
    const BondSymbols& ureyBradley() const;
    
    const StretchStretchSymbols& stretchStretch() const;
    const StretchBendSymbols& stretchBend() const;
    const BendBendSymbols& bendBend() const;
    const StretchBendTorsionSymbols& stretchBendTorsion() const;

private:
    BondSymbols bond_;
    AngleSymbols angle_;
    DihedralSymbols dihedral_;
    ImproperSymbols improper_;
    BondSymbols ureybradley_;
    StretchStretchSymbols stretchstretch_;
    StretchBendSymbols stretchbend_;
    BendBendSymbols bendbend_;
    StretchBendTorsionSymbols stretchbendtorsion_;
};

/** This class holds the internal parameters for a molecule
    (bond, angle, dihedral, improper, Urey-Bradley, 
     stretch-stretch, stretch-bend, bend-bend, stretch-bend-torsion)
     
    The parameters are held in groups, that correspond to the 
    set of CutGroups that contain the atoms, e.g.
    
    group 0 contains all of the parameters that act on
    internals that act only within CutGroup 0
    
    group 0,1 contains all of the parameters that act only
    between atoms in groups 0 and 1
    
    group 0,1,2 contains all of the parameters that act only
    between atoms in groups 0, 1 and 2
     
    @author Christopher Woods
*/
class SIREMM_EXPORT InternalParameters
{

friend QDataStream& ::operator<<(QDataStream&, const InternalParameters&);
friend QDataStream& ::operator>>(QDataStream&, InternalParameters&);

public:
    InternalParameters();
    
    InternalParameters(const PartialMolecule &molecule,
                       const PropertyName &bond_params,
                       const PropertyName &angle_params,
                       const PropertyName &dihedral_params,
                       const PropertyName &improper_params,
                       const PropertyName &ub_params,
                       const PropertyName &ss_params,
                       const PropertyName &sb_params,
                       const PropertyName &bb_params,
                       const PropertyName &sbt_params,
                       bool isstrict);
    
    InternalParameters(const InternalParameters &other);
    
    ~InternalParameters();
    
    static const char* typeName();
    
    const char* what() const
    {
        return InternalParameters::typeName();
    }
    
    InternalParameters& operator=(const InternalParameters &other);
    
    bool operator==(const InternalParameters &other) const;
    bool operator!=(const InternalParameters &other) const;
    
    const InternalSymbols& symbols() const;
    
    bool isEmpty() const;
    
    bool hasBondParameters() const;
    bool hasAngleParameters() const;
    bool hasDihedralParameters() const;
    
    bool hasPhysicalParameters() const;
    
    bool hasImproperParameters() const;
    bool hasUreyBradleyParameters() const;
    
    bool hasNonPhysicalParameters() const;
    
    bool hasStretchStretchParameters() const;
    bool hasStretchBendParameters() const;
    bool hasBendBendParameters() const;
    bool hasStretchBendTorsionParameters() const;
    
    bool hasCrossTerms() const;
    
    bool changedAllGroups(const InternalParameters &other) const;

    void addChangedGroups(const InternalParameters &other, 
                          QSet<quint32> &changed_groups) const;
                          
    QSet<quint32> getChangedGroups(const InternalParameters &other) const;
    
    InternalParameters applyMask(const QSet<quint32> &cgidxs) const;

    const QVector<GroupInternalParameters>& groupParameters() const;
    
    QVector<GroupInternalParameters> groupParameters(quint32 cgidx) const;
    QVector<GroupInternalParameters> groupParameters(const QSet<quint32> &cgidxs) const;
    
private:
    qint32 getIndex(CGIdx cgidx0) const;
    qint32 getIndex(CGIdx cgidx0, CGIdx cgidx1) const;
    qint32 getIndex(CGIdx cgidx0, CGIdx cgidx1,
                    CGIdx cgidx2) const;
    qint32 getIndex(CGIdx cgidx0, CGIdx cgidx1,
                    CGIdx cgidx2, CGIdx cgidx3) const;

    qint32 getIndex(const detail::CGIDQuad &cgids) const;
    
    const GroupInternalParameters& getGroup(CGIdx cgidx0) const;
    const GroupInternalParameters& getGroup(CGIdx cgidx0, CGIdx cgidx1) const;
    const GroupInternalParameters& getGroup(CGIdx cgidx0, CGIdx cgidx1,
                                            CGIdx cgidx2) const;
    const GroupInternalParameters& getGroup(CGIdx cgidx0, CGIdx cgidx1,
                                            CGIdx cgidx2, CGIdx cgidx3) const;
    
    qint32 addGroup(const detail::CGIDQuad &cgids);
    
    GroupInternalParameters& getGroup(const detail::CGIDQuad &cgids,
                                      QHash<detail::CGIDQuad,qint32> &cached_groups);
    
    void assertContainsOnly(const QSet<Symbol> &have_symbols,
                            const QSet<Symbol> &test_symbols) const; 
    
    void addBonds(const TwoAtomFunctions &bondparams,
                  QHash<detail::CGIDQuad,qint32> &cached_groups);
    void addAngles(const ThreeAtomFunctions &angleparams,
                   QHash<detail::CGIDQuad,qint32> &cached_groups);
    void addDihedrals(const FourAtomFunctions &dihedralparams,
                      QHash<detail::CGIDQuad,qint32> &cached_groups);

    void addImpropers(const FourAtomFunctions &improperparams,
                      QHash<detail::CGIDQuad,qint32> &cached_groups);
    void addUBs(const TwoAtomFunctions &ubparams,
                QHash<detail::CGIDQuad,qint32> &cached_groups);

    void addSSs(const ThreeAtomFunctions &ssparams,
                QHash<detail::CGIDQuad,qint32> &cached_groups);
    void addSBs(const ThreeAtomFunctions &sbparams,
                QHash<detail::CGIDQuad,qint32> &cached_groups);
    void addBBs(const FourAtomFunctions &bbparams,
                QHash<detail::CGIDQuad,qint32> &cached_groups);
    void addSBTs(const FourAtomFunctions &sbtparams,
                QHash<detail::CGIDQuad,qint32> &cached_groups);
    
    bool containsOnly(const QSet<quint32> &cgidxs) const;
    
    void updateState();
    void reindex();
    
    /** enum giving the state of these parameters */
    enum { EMPTY        = 0x0000,
           HAS_BOND     = 0x0001,
           HAS_ANGLE    = 0x0002,
           HAS_DIHEDRAL = 0x0004,
           HAS_PHYSICAL = 0x000f,
           HAS_IMPROPER = 0x0010,
           HAS_UB       = 0x0020,
           HAS_NONPHYS  = 0x00f0,
           HAS_SS       = 0x0100,
           HAS_SB       = 0x0200,
           HAS_BB       = 0x0400,
           HAS_SBT      = 0x0800,
           HAS_CROSS    = 0x0f00 };
    
    /** The current state of the parameters */
    uint state;
    
    /** All of the groups of internal parameters */
    QVector<GroupInternalParameters> group_params;
    
    /** The indicies of the groups that contains the parameters
        that involve a particular CutGroup, indexed by the CGIdx
        of that CutGroup */
    QHash< CGIdx, QSet<qint32> > groups_by_cgidx;
    
    /** All of the symbols used by the internal functions */
    static InternalSymbols function_symbols;
};

/** This class holds intramolecular bonding parameters for 3D molecules
    (so it also contains the 3D coordinates of the molecule)
    
    @author Christopher Woods
*/
class SIREMM_EXPORT InternalParameters3D
              : public InternalParameters, protected SireFF::detail::AtomicCoords3D
{

friend QDataStream& ::operator<<(QDataStream&, const InternalParameters3D&);
friend QDataStream& ::operator>>(QDataStream&, InternalParameters3D&);

public:
    InternalParameters3D();
    
    InternalParameters3D(const PartialMolecule &molecule,
                         const PropertyName &coords_property,
                         const PropertyName &bond_params,
                         const PropertyName &angle_params,
                         const PropertyName &dihedral_params,
                         const PropertyName &improper_params,
                         const PropertyName &ub_params,
                         const PropertyName &ss_params,
                         const PropertyName &sb_params,
                         const PropertyName &bb_params,
                         const PropertyName &sbt_params,
                         bool isstrict);
    
    InternalParameters3D(const InternalParameters3D &other);                     
    
    ~InternalParameters3D();
    
    static const char* typeName();
    
    const char* what() const
    {
        return InternalParameters3D::typeName();
    }
    
    InternalParameters3D& operator=(const InternalParameters3D &other);
    
    bool operator==(const InternalParameters3D &other) const;
    bool operator!=(const InternalParameters3D &other) const;
    
    const CoordGroupArray& atomicCoordinates() const;
    
    void setAtomicCoordinates(const AtomicCoords3D &coords);
    
    int nCutGroups() const;
    
    bool changedAllGroups(const InternalParameters3D &other) const;

    void addChangedGroups(const InternalParameters3D &other, 
                          QSet<quint32> &changed_groups) const;
                          
    QSet<quint32> getChangedGroups(const InternalParameters3D &other) const;
    
    InternalParameters3D applyMask(const QSet<quint32> &cgidxs) const;

private:
    InternalParameters3D(const AtomicCoords3D &coords,
                         const InternalParameters &params);
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

///////
/////// Inline functions for InternalParameters
///////

/** Return whether or not this is empty (contains no parameters) */
inline bool InternalParameters::isEmpty() const
{
    return state == EMPTY;
}

/** Return whether or not there are any bond parameters */
inline bool InternalParameters::hasBondParameters() const
{
    return state & HAS_BOND;
}

/** Return whether or not there are any angle parameters */
inline bool InternalParameters::hasAngleParameters() const
{
    return state & HAS_ANGLE;
}

/** Return whether or not there are any dihedral parameters */
inline bool InternalParameters::hasDihedralParameters() const
{
    return state & HAS_DIHEDRAL;
}

/** Return whether or not there are any physical (bond, angle 
    or dihedral) parameters */
inline bool InternalParameters::hasPhysicalParameters() const
{
    return state & HAS_PHYSICAL;
}

/** Return whether or not there are any improper parameters */
inline bool InternalParameters::hasImproperParameters() const
{
    return state & HAS_IMPROPER;
}

/** Return whether or not there are any Urey-Bradley parameters */
inline bool InternalParameters::hasUreyBradleyParameters() const
{
    return state & HAS_UB;
}

/** Return whether or not there are any non-physical (improper or
    Urey-Bradley) parameters */
inline bool InternalParameters::hasNonPhysicalParameters() const
{
    return state & HAS_NONPHYS;
}

/** Return whether or not there are any stretch-stretch parameters */
inline bool InternalParameters::hasStretchStretchParameters() const
{
    return state & HAS_SS;
}

/** Return whether or not there are any stretch-bend parameters */
inline bool InternalParameters::hasStretchBendParameters() const
{
    return state & HAS_SB;
}

/** Return whether or not there are any bend-bend parameters */
inline bool InternalParameters::hasBendBendParameters() const
{
    return state & HAS_BB;
}

/** Return whether or not there are any stretch-bend-torsion parameters */
inline bool InternalParameters::hasStretchBendTorsionParameters() const
{
    return state & HAS_SBT;
}

/** Return whether or not there are any cross (stretch-stretch etc.) parameters */
inline bool InternalParameters::hasCrossTerms() const
{
    return state & HAS_CROSS;
}

#endif

}

Q_DECLARE_METATYPE( SireMM::InternalParameters );
Q_DECLARE_METATYPE( SireMM::InternalParameters3D );

SIRE_EXPOSE_CLASS( SireMM::GroupInternalParameters )
SIRE_EXPOSE_CLASS( SireMM::InternalParameters )
SIRE_EXPOSE_CLASS( SireMM::InternalParameters3D )

SIRE_EXPOSE_CLASS( SireMM::InternalSymbolsBase )
SIRE_EXPOSE_CLASS( SireMM::BondSymbols )
SIRE_EXPOSE_CLASS( SireMM::AngleSymbols )
SIRE_EXPOSE_CLASS( SireMM::DihedralSymbols )
SIRE_EXPOSE_CLASS( SireMM::ImproperSymbols )
SIRE_EXPOSE_CLASS( SireMM::StretchStretchSymbols )
SIRE_EXPOSE_CLASS( SireMM::StretchBendSymbols )
SIRE_EXPOSE_CLASS( SireMM::BendBendSymbols )
SIRE_EXPOSE_CLASS( SireMM::StretchBendTorsionSymbols )
SIRE_EXPOSE_CLASS( SireMM::InternalSymbols )

SIRE_END_HEADER

#endif
