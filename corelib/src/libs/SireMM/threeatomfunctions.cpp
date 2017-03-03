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

#include "threeatomfunctions.h"

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
////// Implementation of ThreeAtomFunction
//////

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const ThreeAtomFunction &threeatomfunc)
{
    ds << threeatomfunc.atm0 << threeatomfunc.atm1
       << threeatomfunc.atm2
       << static_cast<const AtomFunction&>(threeatomfunc);
       
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      ThreeAtomFunction &threeatomfunc)
{
    ds >> threeatomfunc.atm0 >> threeatomfunc.atm1
       >> threeatomfunc.atm2
       >> static_cast<AtomFunction&>(threeatomfunc);
       
    return ds;
}

/** Constructor */
ThreeAtomFunction::ThreeAtomFunction() : AtomFunction()
{}

/** Construct for the specified pair of atoms with the specified function */
ThreeAtomFunction::ThreeAtomFunction(const CGAtomIdx &atom0, const CGAtomIdx &atom1,
                                     const CGAtomIdx &atom2,
                                     const SireCAS::Expression &function)
                : AtomFunction(function),
                  atm0(atom0), atm1(atom1), atm2(atom2)
{}
  
/** Copy constructor */            
ThreeAtomFunction::ThreeAtomFunction(const ThreeAtomFunction &other)
                : AtomFunction(other),
                  atm0(other.atm0), atm1(other.atm1),
                  atm2(other.atm2)
{}

/** Destructor */
ThreeAtomFunction::~ThreeAtomFunction()
{}

/** Copy assignment operator */
ThreeAtomFunction& ThreeAtomFunction::operator=(const ThreeAtomFunction &other)
{
    AtomFunction::operator=(other);
    atm0 = other.atm0;
    atm1 = other.atm1;
    atm2 = other.atm2;
    
    return *this;
}

/** Comparison operator */
bool ThreeAtomFunction::operator==(const ThreeAtomFunction &other) const
{
    return atm0 == other.atm0 and atm1 == other.atm1 and
           atm2 == other.atm2 and
           AtomFunction::operator==(other);
}

/** Comparison operator */
bool ThreeAtomFunction::operator!=(const ThreeAtomFunction &other) const
{
    return atm0 != other.atm0 or atm1 != other.atm1 or
           atm2 != other.atm2 or
           AtomFunction::operator!=(other);
}

/** Return a string representation */
QString ThreeAtomFunction::toString() const
{
    return QObject::tr("ThreeAtomFunction( %1 <- %2 -> %3 : %4 )")
                .arg(atm0.toString(), atm1.toString(), atm2.toString(),
                     this->function().toString());
}

//////
////// Implementation of detail::IDTriple
//////

QDataStream& operator<<(QDataStream &ds, const IDTriple &idtriple)
{
    ds << idtriple.atom0 << idtriple.atom1
       << idtriple.atom2;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds, IDTriple &idtriple)
{
    ds >> idtriple.atom0 >> idtriple.atom1
       >> idtriple.atom2;
       
    return ds;
}

IDTriple::IDTriple(quint32 atm0, quint32 atm1, quint32 atm2) 
       : atom0(atm0), atom1(atm1), atom2(atm2)
{
    if (atm0 > atm2)
    {
        qSwap(atom0,atom2);
    }
}

IDTriple::IDTriple(const IDTriple &other)
       : atom0(other.atom0), atom1(other.atom1), atom2(other.atom2)
{}

IDTriple::~IDTriple()
{}

IDTriple& IDTriple::operator=(const IDTriple &other)
{
    atom0 = other.atom0;
    atom1 = other.atom1;
    atom2 = other.atom2;
    
    return *this;
}

bool IDTriple::operator==(const IDTriple &other) const
{
    return atom0 == other.atom0 and atom1 == other.atom1 and
           atom2 == other.atom2;
}

bool IDTriple::operator!=(const IDTriple &other) const
{
    return atom0 != other.atom0 or atom1 != other.atom1 or
           atom2 != other.atom2;
}

//////
////// Implementation of ThreeAtomFunctions
//////

static const RegisterMetaType<ThreeAtomFunctions> r_threeatomfuncs;

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const ThreeAtomFunctions &threeatomfuncs)
{
    writeHeader(ds, r_threeatomfuncs, 1);
    
    SharedDataStream sds(ds);
    
    sds << threeatomfuncs.potentials_by_atoms
        << static_cast<const AtomFunctions&>(threeatomfuncs);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      ThreeAtomFunctions &threeatomfuncs)
{
    VersionID v = readHeader(ds, r_threeatomfuncs);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> threeatomfuncs.potentials_by_atoms
            >> static_cast<AtomFunctions&>(threeatomfuncs);
    }
    else
        throw version_error(v, "1", r_threeatomfuncs, CODELOC);
        
    return ds;
}

/** Constructor */
ThreeAtomFunctions::ThreeAtomFunctions()
                 : ConcreteProperty<ThreeAtomFunctions,AtomFunctions>()
{}

/** Construct the container to hold the set of three-atom functions
    for the molecule whose data is in 'moldata' */
ThreeAtomFunctions::ThreeAtomFunctions(const MoleculeData &moldata)
                 : ConcreteProperty<ThreeAtomFunctions,AtomFunctions>(moldata)
{}

/** Construct the container to hold the set of three-atom functions
    for the molecule whose layout information is in 'molinfo' */
ThreeAtomFunctions::ThreeAtomFunctions(const MoleculeInfoData &molinfo)
                 : ConcreteProperty<ThreeAtomFunctions,AtomFunctions>(molinfo)
{}

/** Copy constructor */
ThreeAtomFunctions::ThreeAtomFunctions(const ThreeAtomFunctions &other)
                 : ConcreteProperty<ThreeAtomFunctions,AtomFunctions>(other),
                   potentials_by_atoms(other.potentials_by_atoms)
{}

/** Destructor */
ThreeAtomFunctions::~ThreeAtomFunctions()
{}

/** Copy assignment operator */
ThreeAtomFunctions& ThreeAtomFunctions::operator=(const ThreeAtomFunctions &other)
{
    AtomFunctions::operator=(other);
    potentials_by_atoms = other.potentials_by_atoms;
    
    return *this;
}

/** Comparison operator */
bool ThreeAtomFunctions::operator==(const ThreeAtomFunctions &other) const
{
    return AtomFunctions::operator==(other) and
           potentials_by_atoms == other.potentials_by_atoms;
}

/** Comparison operator */
bool ThreeAtomFunctions::operator!=(const ThreeAtomFunctions &other) const
{
    return AtomFunctions::operator!=(other) or
           potentials_by_atoms != other.potentials_by_atoms;
}

/** Return a string representation */
QString ThreeAtomFunctions::toString() const
{
    return QObject::tr("ThreeAtomFunctions( nFunctions() == %1 )")
                .arg(potentials_by_atoms.count());
}

/** Set the potential energy function used by atoms 'atom0', 'atom1' and 'atom2'
    to be equal to 'expression' - this replaces any existing expression 
    
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
void ThreeAtomFunctions::set(AtomIdx atom0, AtomIdx atom1, 
                             AtomIdx atom2, const Expression &expression)
{
    quint32 atm0 = atom0.map( info().nAtoms() );
    quint32 atm1 = atom1.map( info().nAtoms() );
    quint32 atm2 = atom2.map( info().nAtoms() );

    if (atm0 == atm1 or atm0 == atm2 or atm1 == atm2)
        throw SireMol::duplicate_atom( QObject::tr( 
            "You cannot add a function that acts between the same atoms! "
            "(%1-%2-%3)")
                .arg(atm0).arg(atm1).arg(atm2), CODELOC );

    potentials_by_atoms.insert( IDTriple(atm0,atm1,atm2), expression );
    AtomFunctions::addSymbols(expression.symbols());
}
         
/** Set the potential energy function used by atoms 'atom0' to 'atom2'
    to be equal to 'expression' - this replaces any existing expression 
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
void ThreeAtomFunctions::set(const AtomID &atom0, const AtomID &atom1,
                             const AtomID &atom2,
                             const Expression &expression)
{
    this->set( info().atomIdx(atom0), info().atomIdx(atom1), 
               info().atomIdx(atom2), expression );
}

/** Set the potential energy function used for the angle identified by 'angleid'
    to be equal to 'expression' - this replaces any existing expression
    
    This replaces both 1-2-3 and 3-2-1
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
void ThreeAtomFunctions::set(const AngleID &angleid, const Expression &expression)
{
    AtomIdx atom0 = info().atomIdx(angleid.atom0());
    AtomIdx atom1 = info().atomIdx(angleid.atom1());
    AtomIdx atom2 = info().atomIdx(angleid.atom2());
    
    this->clear(atom0, atom1, atom2);
    this->clear(atom2, atom1, atom0);
    
    this->set( atom0, atom1, atom2, expression );
}

/** Check if any of the symbols in 'symbols' need to be removed... */
void ThreeAtomFunctions::removeSymbols(QSet<Symbol> symbols)
{
    for (QHash<IDTriple,Expression>::const_iterator 
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
void ThreeAtomFunctions::clear(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2)
{
    quint32 atm0 = atom0.map( info().nAtoms() );
    quint32 atm1 = atom1.map( info().nAtoms() );
    quint32 atm2 = atom2.map( info().nAtoms() );
        
    ThreeAtomFunctions::removeSymbols( potentials_by_atoms
                                        .take( IDTriple(atm0,atm1,atm2) ).symbols() );
}

/** Clear all functions that involve the atom 'atom' 

    \throw SireError::invalid_index
*/
void ThreeAtomFunctions::clear(AtomIdx atom)
{
    quint32 atm = atom.map(info().nAtoms());
    
    QList<IDTriple> keys = potentials_by_atoms.keys();
    
    foreach (const IDTriple &key, keys)
    {
        if (key.atom0 == atm or key.atom1 == atm or key.atom2 == atm)
        {
            ThreeAtomFunctions::removeSymbols( potentials_by_atoms
                                               .take(key).symbols() );
        }
    }
}

/** Clear any function that acts on the atoms identified by 'atom'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void ThreeAtomFunctions::clear(const AtomID &atom)
{
    QList<AtomIdx> atomidxs = atom.map(info());
    
    foreach (AtomIdx atomidx, atomidxs)
    {
        this->clear(atomidx);
    }
}

/** Clear any function that acts between the atoms 'atom0' to 'atom2'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void ThreeAtomFunctions::clear(const AtomID &atom0, const AtomID &atom1,
                               const AtomID &atom2)
{
    QList<AtomIdx> atoms0 = atom0.map(info());
    QList<AtomIdx> atoms1 = atom1.map(info());
    QList<AtomIdx> atoms2 = atom2.map(info());
    
    foreach (AtomIdx atm0, atoms0)
    {
        foreach (AtomIdx atm1, atoms1)
        {
            foreach (AtomIdx atm2, atoms2)
            {
                this->clear(atm0, atm1, atm2);
            }
        }
    }
}

/** Clear the potential that acts over the angle identified by 'angleid'

    This clears both 1-2-3 and 3-2-1

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void ThreeAtomFunctions::clear(const AngleID &angleid)
{
    this->clear( angleid.atom0(), angleid.atom1(),
                 angleid.atom2() );
                 
    this->clear( angleid.atom2(), angleid.atom1(),
                 angleid.atom0() );
}

/** Completely clear all of the functions from this set */
void ThreeAtomFunctions::clear()
{
    potentials_by_atoms.clear();
    AtomFunctions::removeSymbols();
}

/** Perform the substitutions contained in 'identities' in all of 
    the expressions in this set. This could be useful if you have
    defined these expressions with respect to a lambda parameter,
    and now want to set that value of lambda */
void ThreeAtomFunctions::substitute(const Identities &identities)
{
    AtomFunctions::removeSymbols();

    for (QHash<IDTriple,Expression>::iterator it = potentials_by_atoms.begin();
         it != potentials_by_atoms.end();
         ++it)
    {
        it.value() = it.value().substitute(identities);
        AtomFunctions::addSymbols(it.value().symbols());
    }
}

/** Return whether or not this is empty (has no potentials for any internals) */
bool ThreeAtomFunctions::isEmpty() const
{
    return potentials_by_atoms.isEmpty();
}

/** Return the function acting between the atoms 'atom0' to 'atom2'.
    This returns an empty expression if there is no expression between
    these atoms 
    
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
Expression ThreeAtomFunctions::potential(AtomIdx atom0, AtomIdx atom1,
                                         AtomIdx atom2) const
{
    quint32 atm0 = atom0.map( info().nAtoms() );
    quint32 atm1 = atom1.map( info().nAtoms() );
    quint32 atm2 = atom2.map( info().nAtoms() );
    
    if (atm0 == atm1 or atm0 == atm2 or atm1 == atm2)
        throw SireMol::duplicate_atom( QObject::tr(
            "There is no potential that acts between the same atoms! "
            "(%1-%2-%3)")
                .arg(atm0).arg(atm1).arg(atm2), CODELOC );
                
    return potentials_by_atoms.value( IDTriple(atm0,atm1,atm2) );
}

/** Return the function acting between the atoms 'atom0' to 'atom2'.
    This returns an empty expression if there is no expression between
    these atoms 
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression ThreeAtomFunctions::potential(const AtomID &atom0, 
                                         const AtomID &atom1,
                                         const AtomID &atom2) const
{
    return this->potential( info().atomIdx(atom0), info().atomIdx(atom1),
                            info().atomIdx(atom2) );
}

/** Return the function acting on the angle identified by 'angleid'.
    This returns an empty expression if there is no expression on
    this angle
    
    This search first for the function for 1-2-3, but if that
    is not found, then it returns the function for 3-2-1
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression ThreeAtomFunctions::potential(const AngleID &angleid) const
{
    AtomIdx atom0 = info().atomIdx( angleid.atom0() );
    AtomIdx atom1 = info().atomIdx( angleid.atom1() );
    AtomIdx atom2 = info().atomIdx( angleid.atom2() );
    
    quint32 atm0 = atom0.map( info().nAtoms() );
    quint32 atm1 = atom1.map( info().nAtoms() );
    quint32 atm2 = atom2.map( info().nAtoms() );

    if ( potentials_by_atoms.contains( IDTriple(atm0,atm1,atm2) ) )
    {
        return potentials_by_atoms.value( IDTriple(atm0,atm1,atm2) );
    }
    else
        return potentials_by_atoms.value( IDTriple(atm2,atm1,atm0) );
}

/** Return the force (derivative of the potential with respect to 'symbol')
    between the atoms 'atom0' to 'atom2'
    
    \throw SireError::invalid_index
*/
Expression ThreeAtomFunctions::force(AtomIdx atom0, AtomIdx atom1, 
                                     AtomIdx atom2,
                                     const Symbol &symbol) const
{
    return -(this->potential(atom0,atom1,atom2).differentiate(symbol));
}
                                 
/** Return the force (derivative of the potential with respect to 'symbol')
    between the atoms 'atom0' to 'atom2'
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression ThreeAtomFunctions::force(const AtomID &atom0, const AtomID &atom1,
                                     const AtomID &atom2,
                                     const Symbol &symbol) const
{
    return -(this->potential(atom0,atom1,atom2).differentiate(symbol));
}

/** Return the force (derivative of the potential with respect to 'symbol')
    on the angle identified by 'angleid'
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression ThreeAtomFunctions::force(const AngleID &angleid, 
                                     const Symbol &symbol) const
{
    return -(this->potential(angleid).differentiate(symbol));
}

/** Return the potential energy functions acting between the identified
    triples of atoms */
QVector<ThreeAtomFunction> ThreeAtomFunctions::potentials() const
{
    QVector<ThreeAtomFunction> funcs( potentials_by_atoms.count() );
    
    ThreeAtomFunction *funcs_array = funcs.data();
    
    int i = 0;
    
    for (QHash<IDTriple,Expression>::const_iterator
                                    it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd();
         ++it)
    {
        funcs_array[i] = ThreeAtomFunction( info().cgAtomIdx( AtomIdx(it.key().atom0) ),
                                            info().cgAtomIdx( AtomIdx(it.key().atom1) ),
                                            info().cgAtomIdx( AtomIdx(it.key().atom2) ),
                                            it.value() );
    
        ++i;
    }
    
    return funcs;
}

/** Return the force functions acting between the identified
    triples of atoms, for the given symbol */
QVector<ThreeAtomFunction> ThreeAtomFunctions::forces(const Symbol &symbol) const
{
    QVector<ThreeAtomFunction> forces;
    forces.reserve(potentials_by_atoms.count());
    
    for (QHash<IDTriple,Expression>::const_iterator
                                    it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd();
         ++it)
    {
        Expression force = it.value().differentiate(symbol);
        
        if (not force.isZero())
        {
            forces.append( ThreeAtomFunction( 
                              info().cgAtomIdx( AtomIdx(it.key().atom0) ),
                              info().cgAtomIdx( AtomIdx(it.key().atom1) ),
                              info().cgAtomIdx( AtomIdx(it.key().atom2) ),
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
ThreeAtomFunctions ThreeAtomFunctions::includeOnly(const AtomSelection &selected_atoms,
                                                   bool isstrict) const
{
    ThreeAtomFunctions ret(*this);

    QMutableHashIterator<IDTriple,Expression> it(ret.potentials_by_atoms);

    if (isstrict)
    {
        while (it.hasNext())
        {
            it.next();
            
            if (not (selected_atoms.selected(AtomIdx(it.key().atom0)) and
                     selected_atoms.selected(AtomIdx(it.key().atom1)) and
                     selected_atoms.selected(AtomIdx(it.key().atom2)) ) )
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
                     selected_atoms.selected(AtomIdx(it.key().atom2)) ) )
            {
                it.remove();
            }
        }
    }
    
    return ret;
}

/** Return the number of functions in this set */
int ThreeAtomFunctions::nFunctions() const
{
    return potentials_by_atoms.count();
}

/** Return a copy of this property that has been made to be compatible
    with the molecule layout in 'molinfo' - this uses the atom matching
    functions in 'atommatcher' to match atoms from the current molecule
    to the atoms in the molecule whose layout is in 'molinfo'

    This only copies the ThreeAtomFunction for pairs of atoms that
    are successfully matched - it does not copy functions for atoms
    that are not matched. Use ThreeAtomFunctions::nFunctions() to check
    if the number of functions in the returned set is the same as
    the number in this set, if you want to ensure that all of the 
    functions have been copied.
    
    \throw SireError::incompatible_error
*/
PropertyPtr 
ThreeAtomFunctions::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                            const AtomMatcher &atommatcher) const
{
    if (not atommatcher.changesOrder(this->info(), molinfo))
    {
        //the order of the atoms remains the same - this means that the 
        //AtomIdx indicies are still valid
        ThreeAtomFunctions ret(molinfo);
        ret.potentials_by_atoms = this->potentials_by_atoms;
        return ret;
    }

    QHash<AtomIdx,AtomIdx> matched_atoms = atommatcher.match(this->info(), molinfo);
    
    ThreeAtomFunctions ret(molinfo);
    
    for (QHash<IDTriple,Expression>::const_iterator it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd();
         ++it)
    {
        AtomIdx new_atom0 = matched_atoms.value( AtomIdx(it.key().atom0), AtomIdx(-1) );
        AtomIdx new_atom1 = matched_atoms.value( AtomIdx(it.key().atom1), AtomIdx(-1) );
        AtomIdx new_atom2 = matched_atoms.value( AtomIdx(it.key().atom2), AtomIdx(-1) );

        if (new_atom0 == -1 or new_atom1 == -1 or new_atom2 == -1)
            continue;
        
        ret.set( new_atom0, new_atom1, new_atom2, it.value() );
    }
    
    return ret;
}

const char* ThreeAtomFunctions::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ThreeAtomFunctions>() );
}
