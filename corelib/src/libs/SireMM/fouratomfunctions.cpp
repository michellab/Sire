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

#include "sireglobal.h"

#include "fouratomfunctions.h"

#include "SireCAS/symbols.h"

#include "SireMol/moleculeinfodata.h"
#include "SireMol/atomselection.h"
#include "SireMol/atommatcher.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMM::detail;

using namespace SireBase;
using namespace SireCAS;
using namespace SireMol;
using namespace SireStream;

//////
////// Implementation of FourAtomFunction
//////

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const FourAtomFunction &fouratomfunc)
{
    ds << fouratomfunc.atm0 << fouratomfunc.atm1
       << fouratomfunc.atm2 << fouratomfunc.atm3
       << static_cast<const AtomFunction&>(fouratomfunc);
       
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      FourAtomFunction &fouratomfunc)
{
    ds >> fouratomfunc.atm0 >> fouratomfunc.atm1
       >> fouratomfunc.atm2 >> fouratomfunc.atm3
       >> static_cast<AtomFunction&>(fouratomfunc);
       
    return ds;
}

/** Constructor */
FourAtomFunction::FourAtomFunction() : AtomFunction()
{}

/** Construct for the specified pair of atoms with the specified function */
FourAtomFunction::FourAtomFunction(const CGAtomIdx &atom0, const CGAtomIdx &atom1,
                                   const CGAtomIdx &atom2, const CGAtomIdx &atom3,
                                   const SireCAS::Expression &function)
                : AtomFunction(function),
                  atm0(atom0), atm1(atom1), atm2(atom2), atm3(atom3)
{}
  
/** Copy constructor */            
FourAtomFunction::FourAtomFunction(const FourAtomFunction &other)
                : AtomFunction(other),
                  atm0(other.atm0), atm1(other.atm1),
                  atm2(other.atm2), atm3(other.atm3)
{}

/** Destructor */
FourAtomFunction::~FourAtomFunction()
{}

/** Copy assignment operator */
FourAtomFunction& FourAtomFunction::operator=(const FourAtomFunction &other)
{
    AtomFunction::operator=(other);
    atm0 = other.atm0;
    atm1 = other.atm1;
    atm2 = other.atm2;
    atm3 = other.atm3;
    
    return *this;
}

/** Comparison operator */
bool FourAtomFunction::operator==(const FourAtomFunction &other) const
{
    return atm0 == other.atm0 and atm1 == other.atm1 and
           atm2 == other.atm2 and atm3 == other.atm3 and
           AtomFunction::operator==(other);
}

/** Comparison operator */
bool FourAtomFunction::operator!=(const FourAtomFunction &other) const
{
    return atm0 != other.atm0 or atm1 != other.atm1 or
           atm2 != other.atm2 or atm3 != other.atm3 or
           AtomFunction::operator!=(other);
}

/** Return a string representation */
QString FourAtomFunction::toString() const
{
    return QObject::tr("FourAtomFunction( %1 <- %2 - %3 -> %4 : %5 )")
                .arg(atm0.toString(), atm1.toString(),
                     atm2.toString(), atm3.toString())
                .arg(this->function().toString());
}

//////
////// Implementation of detail::IDQuad
//////

QDataStream& operator<<(QDataStream &ds, const IDQuad &idquad)
{
    ds << idquad.atom0 << idquad.atom1
       << idquad.atom2 << idquad.atom3;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds, IDQuad &idquad)
{
    ds >> idquad.atom0 >> idquad.atom1
       >> idquad.atom2 >> idquad.atom3;
       
    return ds;
}

IDQuad::IDQuad(quint32 atm0, quint32 atm1, quint32 atm2, quint32 atm3) 
       : atom0(atm0), atom1(atm1), atom2(atm2), atom3(atm3)
{}

IDQuad::IDQuad(const IDQuad &other)
       : atom0(other.atom0), atom1(other.atom1), 
         atom2(other.atom2), atom3(other.atom3)
{}

IDQuad::~IDQuad()
{}

IDQuad& IDQuad::operator=(const IDQuad &other)
{
    atom0 = other.atom0;
    atom1 = other.atom1;
    atom2 = other.atom2;
    atom3 = other.atom3;
    
    return *this;
}

bool IDQuad::operator==(const IDQuad &other) const
{
    return atom0 == other.atom0 and atom1 == other.atom1 and
           atom2 == other.atom2 and atom3 == other.atom3;
}

bool IDQuad::operator!=(const IDQuad &other) const
{
    return atom0 != other.atom0 or atom1 != other.atom1 or
           atom2 != other.atom2 or atom3 != other.atom3;
}

//////
////// Implementation of FourAtomFunctions
//////

static const RegisterMetaType<FourAtomFunctions> r_fouratomfuncs;

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const FourAtomFunctions &fouratomfuncs)
{
    writeHeader(ds, r_fouratomfuncs, 1);
    
    SharedDataStream sds(ds);
    
    sds << fouratomfuncs.potentials_by_atoms
        << static_cast<const AtomFunctions&>(fouratomfuncs);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      FourAtomFunctions &fouratomfuncs)
{
    VersionID v = readHeader(ds, r_fouratomfuncs);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> fouratomfuncs.potentials_by_atoms
            >> static_cast<AtomFunctions&>(fouratomfuncs);
    }
    else
        throw version_error(v, "1", r_fouratomfuncs, CODELOC);
        
    return ds;
}

/** Constructor */
FourAtomFunctions::FourAtomFunctions()
                 : ConcreteProperty<FourAtomFunctions,AtomFunctions>()
{}

/** Construct the container to hold the set of four-atom functions
    for the molecule whose data is in 'moldata' */
FourAtomFunctions::FourAtomFunctions(const MoleculeData &moldata)
                 : ConcreteProperty<FourAtomFunctions,AtomFunctions>(moldata)
{}

/** Construct the container to hold the set of four-atom functions
    for the molecule whose layout information is in 'molinfo' */
FourAtomFunctions::FourAtomFunctions(const MoleculeInfoData &molinfo)
                 : ConcreteProperty<FourAtomFunctions,AtomFunctions>(molinfo)
{}

/** Copy constructor */
FourAtomFunctions::FourAtomFunctions(const FourAtomFunctions &other)
                 : ConcreteProperty<FourAtomFunctions,AtomFunctions>(other),
                   potentials_by_atoms(other.potentials_by_atoms)
{}

/** Destructor */
FourAtomFunctions::~FourAtomFunctions()
{}

/** Copy assignment operator */
FourAtomFunctions& FourAtomFunctions::operator=(const FourAtomFunctions &other)
{
    AtomFunctions::operator=(other);
    potentials_by_atoms = other.potentials_by_atoms;
    
    return *this;
}

/** Comparison operator */
bool FourAtomFunctions::operator==(const FourAtomFunctions &other) const
{
    return AtomFunctions::operator==(other) and
           potentials_by_atoms == other.potentials_by_atoms;
}

/** Comparison operator */
bool FourAtomFunctions::operator!=(const FourAtomFunctions &other) const
{
    return AtomFunctions::operator!=(other) or
           potentials_by_atoms != other.potentials_by_atoms;
}

/** Return a string representation */
QString FourAtomFunctions::toString() const
{
    return QObject::tr("FourAtomFunctions( nFunctions() == %1 )")
                .arg(potentials_by_atoms.count());
}

/** Set the potential energy function used by atoms 'atom0' to 'atom3'
    to be equal to 'expression' - this replaces any existing expression
    
    If you want to add an improper then use  
    
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
void FourAtomFunctions::set(AtomIdx atom0, AtomIdx atom1, 
                            AtomIdx atom2, AtomIdx atom3,
                            const Expression &expression)
{
    quint32 atm0 = atom0.map( info().nAtoms() );
    quint32 atm1 = atom1.map( info().nAtoms() );
    quint32 atm2 = atom2.map( info().nAtoms() );
    quint32 atm3 = atom3.map( info().nAtoms() );

    if (atm0 == atm1 or atm0 == atm2 or atm0 == atm3 or
        atm1 == atm2 or atm1 == atm3 or atm2 == atm3 )
        throw SireMol::duplicate_atom( QObject::tr( 
            "You cannot add a function that acts between the same atoms! "
            "(%1-%2-%3-%4)")
                .arg(atm0).arg(atm1).arg(atm2).arg(atm3), CODELOC );

    potentials_by_atoms.insert( IDQuad(atm0,atm1,atm2,atm3), expression );
    AtomFunctions::addSymbols(expression.symbols());
}
         
/** Set the potential energy function used by atoms 'atom0' to 'atom3'
    to be equal to 'expression' - this replaces any existing expression 
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
void FourAtomFunctions::set(const AtomID &atom0, const AtomID &atom1,
                            const AtomID &atom2, const AtomID &atom3,
                            const Expression &expression)
{
    this->set( info().atomIdx(atom0), info().atomIdx(atom1), 
               info().atomIdx(atom2), info().atomIdx(atom3), expression );
}

/** Set the potential energy function used for the dihedral identified by 'dihedralid'
    to be equal to 'expression' - this replaces any existing expression

    Note that this will replace *any* equivalent dihedral, so this will replace
    both 1-2-3-4 and 4-3-2-1 with the new expression
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
void FourAtomFunctions::set(const DihedralID &dihedralid, const Expression &expression)
{
    AtomIdx atom0 = info().atomIdx(dihedralid.atom0());
    AtomIdx atom1 = info().atomIdx(dihedralid.atom1());
    AtomIdx atom2 = info().atomIdx(dihedralid.atom2());
    AtomIdx atom3 = info().atomIdx(dihedralid.atom3());
    
    this->clear(atom0, atom1, atom2, atom3);
    this->clear(atom3, atom2, atom1, atom0);
    
    this->set(atom0, atom1, atom2, atom3, expression);
}

/** Set the potential energy function used for the improper identified by 'improperid'
    to be equal to 'expression' - this replaces any existing expression
    
    Note that this replaces any existing improper with this function,
    e.g. this replaces both 1-2-3-4 and 1-2-4-3 with this function
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
void FourAtomFunctions::set(const ImproperID &improperid, const Expression &expression)
{
    AtomIdx atom0 = info().atomIdx(improperid.atom0());
    AtomIdx atom1 = info().atomIdx(improperid.atom1());
    AtomIdx atom2 = info().atomIdx(improperid.atom2());
    AtomIdx atom3 = info().atomIdx(improperid.atom3());
    
    this->clear(atom0, atom1, atom2, atom3);
    this->clear(atom0, atom1, atom3, atom2);
    
    this->set(atom0, atom1, atom2, atom3, expression);
}

/** Check if any of the symbols in 'symbols' need to be removed... */
void FourAtomFunctions::removeSymbols(QSet<Symbol> symbols)
{
    for (QHash<IDQuad,Expression>::const_iterator 
                                    it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd();
         ++it)
    {
        if (symbols.isEmpty())
            return;
            
        symbols.subtract(it.value().symbols());
    }
    
    //the only remaining symbols are ones that no longer exist
    //in this set
    AtomFunctions::removeSymbols(symbols);
}

/** Clear any function that acts between the atoms 'atom0' to 'atom2'

    \throw SireError::invalid_index
*/
void FourAtomFunctions::clear(AtomIdx atom0, AtomIdx atom1, 
                              AtomIdx atom2, AtomIdx atom3)
{
    quint32 atm0 = atom0.map( info().nAtoms() );
    quint32 atm1 = atom1.map( info().nAtoms() );
    quint32 atm2 = atom2.map( info().nAtoms() );
    quint32 atm3 = atom3.map( info().nAtoms() );
        
    FourAtomFunctions::removeSymbols( potentials_by_atoms
                                       .take( IDQuad(atm0,atm1,atm2,atm3) ).symbols() );
}

/** Clear all functions that involve the atom 'atom' 

    \throw SireError::invalid_index
*/
void FourAtomFunctions::clear(AtomIdx atom)
{
    quint32 atm = atom.map(info().nAtoms());
    
    QList<IDQuad> keys = potentials_by_atoms.keys();
    
    foreach (const IDQuad &key, keys)
    {
        if (key.atom0 == atm or key.atom1 == atm or 
            key.atom2 == atm or key.atom3 == atm)
        {
            FourAtomFunctions::removeSymbols( potentials_by_atoms
                                               .take(key).symbols() );
        }
    }
}

/** Clear any function that acts on the atoms identified by 'atom'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void FourAtomFunctions::clear(const AtomID &atom)
{
    QList<AtomIdx> atomidxs = atom.map(info());
    
    foreach (AtomIdx atomidx, atomidxs)
    {
        this->clear(atomidx);
    }
}

/** Clear any function that acts between the atoms 'atom0' to 'atom3'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void FourAtomFunctions::clear(const AtomID &atom0, const AtomID &atom1,
                              const AtomID &atom2, const AtomID &atom3)
{
    QList<AtomIdx> atoms0 = atom0.map(info());
    QList<AtomIdx> atoms1 = atom1.map(info());
    QList<AtomIdx> atoms2 = atom2.map(info());
    QList<AtomIdx> atoms3 = atom3.map(info());
    
    foreach (AtomIdx atm0, atoms0)
    {
        foreach (AtomIdx atm1, atoms1)
        {
            foreach (AtomIdx atm2, atoms2)
            {
                foreach (AtomIdx atm3, atoms3)
                {
                    this->clear(atm0, atm1, atm2, atm3);
                }
            }
        }
    }
}

/** Clear the potential that acts over the dihedral identified by 'dihedralid'

    This clears all matching dihedrals, so 1-2-3-4 and 4-3-2-1

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void FourAtomFunctions::clear(const DihedralID &dihedralid)
{
    this->clear( dihedralid.atom0(), dihedralid.atom1(),
                 dihedralid.atom2(), dihedralid.atom3() );
                 
    this->clear( dihedralid.atom3(), dihedralid.atom2(),
                 dihedralid.atom1(), dihedralid.atom0() );
}

/** Clear the potential that acts over the improper identified by 'improperid'

    This clears all matching impropers, so 1-2-3-4 and 1-2-4-3

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void FourAtomFunctions::clear(const ImproperID &improperid)
{
    this->clear( improperid.atom0(), improperid.atom1(),
                 improperid.atom2(), improperid.atom3() );
                 
    this->clear( improperid.atom0(), improperid.atom1(),
                 improperid.atom3(), improperid.atom2() );
}

/** Completely clear all of the functions from this set */
void FourAtomFunctions::clear()
{
    potentials_by_atoms.clear();
    AtomFunctions::removeSymbols();
}

/** Perform the substitutions contained in 'identities' in all of 
    the expressions in this set. This could be useful if you have
    defined these expressions with respect to a lambda parameter,
    and now want to set that value of lambda */
void FourAtomFunctions::substitute(const Identities &identities)
{
    AtomFunctions::removeSymbols();

    for (QHash<IDQuad,Expression>::iterator it = potentials_by_atoms.begin();
         it != potentials_by_atoms.end();
         ++it)
    {
        it.value() = it.value().substitute(identities);
        AtomFunctions::addSymbols(it.value().symbols());
    }
}

/** Return whether or not this is empty (has no potentials for any internals) */
bool FourAtomFunctions::isEmpty() const
{
    return potentials_by_atoms.isEmpty();
}

/** Return the function acting between the atoms 'atom0' to 'atom3'.
    This returns an empty expression if there is no expression between
    these atoms 
    
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
Expression FourAtomFunctions::potential(AtomIdx atom0, AtomIdx atom1,
                                        AtomIdx atom2, AtomIdx atom3) const
{
    quint32 atm0 = atom0.map( info().nAtoms() );
    quint32 atm1 = atom1.map( info().nAtoms() );
    quint32 atm2 = atom2.map( info().nAtoms() );
    quint32 atm3 = atom3.map( info().nAtoms() );
    
    if (atm0 == atm1 or atm0 == atm2 or atm0 == atm3 or
        atm1 == atm2 or atm1 == atm3 or atm2 == atm3)
        throw SireMol::duplicate_atom( QObject::tr(
            "There is no potential that acts between the same atoms! "
            "(%1-%2-%3-%4)")
                .arg(atm0).arg(atm1).arg(atm2).arg(atm3), CODELOC );
                
    return potentials_by_atoms.value( IDQuad(atm0,atm1,atm2,atm3) );
}

/** Return the function acting between the atoms 'atom0' to 'atom3'.
    This returns an empty expression if there is no expression between
    these atoms 
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression FourAtomFunctions::potential(const AtomID &atom0, 
                                        const AtomID &atom1,
                                        const AtomID &atom2,
                                        const AtomID &atom3) const
{
    return this->potential( info().atomIdx(atom0), info().atomIdx(atom1),
                            info().atomIdx(atom2), info().atomIdx(atom3) );
}

/** Return the function acting on the dihedral identified by 'dihedralid'.
    This returns an empty expression if there is no expression on
    this dihedral
    
    Note that this searches first for 1-2-3-4, then if no function
    is found, it returns the function for 4-3-2-1
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression FourAtomFunctions::potential(const DihedralID &dihedralid) const
{
    AtomIdx atom0 = info().atomIdx( dihedralid.atom0() );
    AtomIdx atom1 = info().atomIdx( dihedralid.atom1() );
    AtomIdx atom2 = info().atomIdx( dihedralid.atom2() );
    AtomIdx atom3 = info().atomIdx( dihedralid.atom3() );
    
    quint32 atm0 = atom0.map( info().nAtoms() );
    quint32 atm1 = atom1.map( info().nAtoms() );
    quint32 atm2 = atom2.map( info().nAtoms() );
    quint32 atm3 = atom3.map( info().nAtoms() );

    if (potentials_by_atoms.contains(IDQuad(atm0,atm1,atm2,atm3)))
    {
        return potentials_by_atoms.value(IDQuad(atm0,atm1,atm2,atm3));
    }
    else
        return potentials_by_atoms.value( IDQuad(atm3,atm2,atm1,atm0) );
}

/** Return the function acting on the improper identified by 'improperid'.
    This returns an empty expression if there is no expression on
    this improper
   
    Note that this searches first for 1-2-3-4, then if no function
    is found, it returns the function for 1-2-4-3

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression FourAtomFunctions::potential(const ImproperID &improperid) const
{
    AtomIdx atom0 = info().atomIdx( improperid.atom0() );
    AtomIdx atom1 = info().atomIdx( improperid.atom1() );
    AtomIdx atom2 = info().atomIdx( improperid.atom2() );
    AtomIdx atom3 = info().atomIdx( improperid.atom3() );
    
    quint32 atm0 = atom0.map( info().nAtoms() );
    quint32 atm1 = atom1.map( info().nAtoms() );
    quint32 atm2 = atom2.map( info().nAtoms() );
    quint32 atm3 = atom3.map( info().nAtoms() );

    if (potentials_by_atoms.contains(IDQuad(atm0,atm1,atm2,atm3)))
    {
        return potentials_by_atoms.value(IDQuad(atm0,atm1,atm2,atm3));
    }
    else
        return potentials_by_atoms.value( IDQuad(atm0,atm1,atm3,atm2) );
}

/** Return the force (derivative of the potential with respect to 'symbol')
    between the atoms 'atom0' to 'atom3'
    
    \throw SireError::invalid_index
*/
Expression FourAtomFunctions::force(AtomIdx atom0, AtomIdx atom1, 
                                    AtomIdx atom2, AtomIdx atom3,
                                    const Symbol &symbol) const
{
    return -(this->potential(atom0,atom1,atom2,atom3).differentiate(symbol));
}
                                 
/** Return the force (derivative of the potential with respect to 'symbol')
    between the atoms 'atom0' to 'atom3'
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression FourAtomFunctions::force(const AtomID &atom0, const AtomID &atom1,
                                    const AtomID &atom2, const AtomID &atom3,
                                    const Symbol &symbol) const
{
    return -(this->potential(atom0,atom1,atom2,atom3).differentiate(symbol));
}

/** Return the force (derivative of the potential with respect to 'symbol')
    on the dihedral identified by 'dihedralid'
    
    Note that this searches 1-2-3-4 first, then searches for 4-3-2-1
    if no function is found
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression FourAtomFunctions::force(const DihedralID &dihedralid, 
                                    const Symbol &symbol) const
{
    return -(this->potential(dihedralid).differentiate(symbol));
}

/** Return the force (derivative of the potential with respect to 'symbol')
    on the improper identified by 'improperid'
    
    Note that this searches 1-2-3-4 first, then searches for 1-2-4-3
    if no function is found

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression FourAtomFunctions::force(const ImproperID &improperid, 
                                    const Symbol &symbol) const
{
    return -(this->potential(improperid).differentiate(symbol));
}

/** Return the potential energy functions acting between the identified
    quads of atoms */
QVector<FourAtomFunction> FourAtomFunctions::potentials() const
{
    QVector<FourAtomFunction> funcs( potentials_by_atoms.count() );
    
    FourAtomFunction *funcs_array = funcs.data();
    
    int i = 0;
    
    for (QHash<IDQuad,Expression>::const_iterator
                                    it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd();
         ++it)
    {
        funcs_array[i] = FourAtomFunction( info().cgAtomIdx( AtomIdx(it.key().atom0) ),
                                           info().cgAtomIdx( AtomIdx(it.key().atom1) ),
                                           info().cgAtomIdx( AtomIdx(it.key().atom2) ),
                                           info().cgAtomIdx( AtomIdx(it.key().atom3) ),
                                           it.value() );
    
        ++i;
    }
    
    return funcs;
}

/** Return the force functions acting between the identified
    quads of atoms, for the given symbol */
QVector<FourAtomFunction> FourAtomFunctions::forces(const Symbol &symbol) const
{
    QVector<FourAtomFunction> forces;
    forces.reserve(potentials_by_atoms.count());
    
    for (QHash<IDQuad,Expression>::const_iterator
                                    it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd();
         ++it)
    {
        Expression force = it.value().differentiate(symbol);
        
        if (not force.isZero())
        {
            forces.append( FourAtomFunction( 
                              info().cgAtomIdx( AtomIdx(it.key().atom0) ),
                              info().cgAtomIdx( AtomIdx(it.key().atom1) ),
                              info().cgAtomIdx( AtomIdx(it.key().atom2) ),
                              info().cgAtomIdx( AtomIdx(it.key().atom3) ),
                              -force ) );
        }
    }
    
    return forces;
}

/** Return the set of functions where only functions that involve the 
    atoms in 'selected_atoms' are included. If 'isstrict' is true, then
    only include functions where all of the atoms are in 'selected_atoms',
    while if 'isstrict' is false, include functions where at least one
    atom is in 'selected_atoms' */
FourAtomFunctions FourAtomFunctions::includeOnly(const AtomSelection &selected_atoms,
                                                 bool isstrict) const
{
    FourAtomFunctions ret(*this);

    QMutableHashIterator<IDQuad,Expression> it(ret.potentials_by_atoms);

    if (isstrict)
    {
        while (it.hasNext())
        {
            it.next();
            
            if (not (selected_atoms.selected(AtomIdx(it.key().atom0)) and
                     selected_atoms.selected(AtomIdx(it.key().atom1)) and
                     selected_atoms.selected(AtomIdx(it.key().atom2)) and
                     selected_atoms.selected(AtomIdx(it.key().atom3)) ) )
            {
                it.remove();
            }
        }
    }
    else
    {
        while (it.hasNext())
        {
            it.next();
            
            if (not (selected_atoms.selected(AtomIdx(it.key().atom0)) or
                     selected_atoms.selected(AtomIdx(it.key().atom1)) or
                     selected_atoms.selected(AtomIdx(it.key().atom2)) or
                     selected_atoms.selected(AtomIdx(it.key().atom3)) ) )
            {
                it.remove();
            }
        }
    }
    
    return ret;
}

/** Return the number of functions in this set */
int FourAtomFunctions::nFunctions() const
{
    return potentials_by_atoms.count();
}

/** Return a copy of this property that has been made to be compatible
    with the molecule layout in 'molinfo' - this uses the atom matching
    functions in 'atommatcher' to match atoms from the current molecule
    to the atoms in the molecule whose layout is in 'molinfo'

    This only copies the FourAtomFunction for pairs of atoms that
    are successfully matched - it does not copy functions for atoms
    that are not matched. Use FouAtomFunctions::nFunctions() to check
    if the number of functions in the returned set is the same as
    the number in this set, if you want to ensure that all of the 
    functions have been copied.
    
    \throw SireError::incompatible_error
*/
PropertyPtr 
FourAtomFunctions::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                           const AtomMatcher &atommatcher) const
{
    if (not atommatcher.changesOrder(this->info(), molinfo))
    {
        //the order of the atoms remains the same - this means that the 
        //AtomIdx indicies are still valid
        FourAtomFunctions ret(molinfo);
        ret.potentials_by_atoms = this->potentials_by_atoms;
        return ret;
    }

    QHash<AtomIdx,AtomIdx> matched_atoms = atommatcher.match(this->info(), molinfo);
    
    FourAtomFunctions ret(molinfo);
    
    for (QHash<IDQuad,Expression>::const_iterator it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd();
         ++it)
    {
        AtomIdx new_atom0 = matched_atoms.value( AtomIdx(it.key().atom0), AtomIdx(-1) );
        AtomIdx new_atom1 = matched_atoms.value( AtomIdx(it.key().atom1), AtomIdx(-1) );
        AtomIdx new_atom2 = matched_atoms.value( AtomIdx(it.key().atom2), AtomIdx(-1) );
        AtomIdx new_atom3 = matched_atoms.value( AtomIdx(it.key().atom3), AtomIdx(-1) );

        if (new_atom0 == -1 or new_atom1 == -1 or new_atom2 == -1 or new_atom3 == -1)
            continue;
        
        ret.set( new_atom0, new_atom1, new_atom2, new_atom3, it.value() );
    }
    
    return ret;
}

const char* FourAtomFunctions::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FourAtomFunctions>() );
}
