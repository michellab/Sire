/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include <QElapsedTimer>

#include "cljatoms.h"

#include "atomljs.h"

#include "SireMol/atomcoords.h"
#include "SireMol/atomcharges.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/molecules.h"
#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/selector.hpp"
#include "SireMol/atom.h"
#include "SireMol/molidx.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireStream;

///////
/////// Implementation of CLJAtom
///////

static const RegisterMetaType<CLJAtom> r_cljatom(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJAtom &cljatom)
{
    writeHeader(ds, r_cljatom, 1);
    
    ds << cljatom.x << cljatom.y << cljatom.z
       << cljatom.chg << cljatom.sig << cljatom.eps
       << cljatom.idnum;
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJAtom &cljatom)
{
    VersionID v = readHeader(ds, r_cljatom);
    
    if (v == 1)
    {
        ds >> cljatom.x >> cljatom.y >> cljatom.z
           >> cljatom.chg >> cljatom.sig >> cljatom.eps
           >> cljatom.idnum;
    }
    else
        throw version_error(v, "1", r_cljatom, CODELOC);
    
    return ds;
}

/** Null constructor */
CLJAtom::CLJAtom() : x(0), y(0), z(0), chg(0), sig(0), eps(0), idnum(0)
{}

/** Return whether or not this is a dummy atom (has no ID number) */
bool CLJAtom::isDummy() const
{
    return idnum == 0;
}

/** Return whether or not this atom is null (no information) */
bool CLJAtom::isNull() const
{
    return this->operator==( CLJAtom() );
}

/** Construct from the passed coordinates, charge and LJ parameters */
CLJAtom::CLJAtom(Vector coords, Charge charge, LJParameter ljparam, qint32 atomid)
        : x(coords.x()), y(coords.y()), z(coords.z()),
          chg( charge.value() * std::sqrt(SireUnits::one_over_four_pi_eps0) ),
          sig( sqrt(ljparam.sigma()) ), eps( sqrt(4.0 * ljparam.epsilon()) ),
          idnum(atomid)
{}

/** Copy constructor */
CLJAtom::CLJAtom(const CLJAtom &other)
        : x(other.x), y(other.y), z(other.z), chg(other.chg),
          sig(other.sig), eps(other.eps), idnum(other.idnum)
{}

/** Destructor */
CLJAtom::~CLJAtom()
{}

/** Copy assignment operator */
CLJAtom& CLJAtom::operator=(const CLJAtom &other)
{
    if (this != &other)
    {
        x = other.x;
        y = other.y;
        z = other.z;
        chg = other.chg;
        sig = other.sig;
        eps = other.eps;
        idnum = other.idnum;
    }
    
    return *this;
}

/** Comparison operator */
bool CLJAtom::operator==(const CLJAtom &other) const
{
    return x == other.x and
           y == other.y and
           z == other.z and
           chg == other.chg and
           sig == other.sig and
           eps == other.eps and
           idnum == other.idnum;
}

/** Comparison operator */
bool CLJAtom::operator!=(const CLJAtom &other) const
{
    return not operator==(other);
}

QString CLJAtom::toString() const
{
    return QObject::tr("CLJAtom( %1, %2, %3, %4 )")
           .arg(coordinates().toString())
           .arg(charge().toString())
           .arg(ljParameter().toString())
           .arg(ID());
}

const char* CLJAtom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJAtom>() );
}

const char* CLJAtom::what() const
{
    return CLJAtom::typeName();
}

/** Construct an array of CLJAtom atoms from the passed molecule view */
QVector<CLJAtom> CLJAtom::buildFrom(const MoleculeView &molecule, const PropertyMap &map)
{
    //QElapsedTimer t;
    //t.start();
    
    const PropertyName coords_property = map["coordinates"];
    const PropertyName chg_property = map["charge"];
    const PropertyName lj_property = map["LJ"];
    
    QVector<CLJAtom> cljatoms;
    
    if (molecule.selectedAll())
    {
        Molecule mol = molecule.molecule();
        
        const int nats = mol.nAtoms();
        
        if (nats == 0)
            return QVector<CLJAtom>();
        
        //reserve space for the data
        cljatoms = QVector<CLJAtom>(nats);
        CLJAtom *cljatms = cljatoms.data();

        const AtomCoords &coords = mol.property(coords_property).asA<AtomCoords>();
        const AtomCharges &chgs = mol.property(chg_property).asA<AtomCharges>();
        const AtomLJs &ljs = mol.property(lj_property).asA<AtomLJs>();

        const quint32 molid = mol.number().value();
        const qint32 s_molid = *(reinterpret_cast<const qint32*>(&molid));
        
        int idx = 0;
        
        for (int i=0; i<coords.nCutGroups(); ++i)
        {
            const CGIdx cgidx(i);
        
            const Vector *icoords = coords.constData(cgidx);
            const Charge *ichg = chgs.constData(cgidx);
            const LJParameter *ilj = ljs.constData(cgidx);
            
            for (int j=0; j<coords.nAtoms(cgidx); ++j)
            {
                CLJAtom &cljatm = cljatms[idx];
            
                cljatm.x = icoords[j].x();
                cljatm.y = icoords[j].y();
                cljatm.z = icoords[j].z();
                
                cljatm.chg = ichg[j].value() * std::sqrt(SireUnits::one_over_four_pi_eps0);

                cljatm.sig = std::sqrt(ilj[j].sigma());
                cljatm.eps = std::sqrt(4.0 * ilj[j].epsilon());
                
                cljatm.idnum = s_molid;
                
                idx += 1;
            }
        }
    }
    else
    {
        Selector<Atom> atoms = molecule.atoms();
        
        QList<Vector> coords = atoms.property<Vector>(coords_property);
        QList<Charge> chgs = atoms.property<Charge>(chg_property);
        QList<LJParameter> ljs = atoms.property<LJParameter>(lj_property);
        
        if (coords.count() != chgs.count() or
            coords.count() != ljs.count())
        {
            throw SireError::program_bug( QObject::tr(
                    "It should not be possible for the number of coordinates (%1) "
                    "to be different to the number of charges (%2) or LJs (%3) "
                    "for molecule %4.")
                        .arg(coords.count())
                        .arg(chgs.count())
                        .arg(ljs.count())
                        .arg(atoms.molecule().toString()), CODELOC );
        }

        const quint32 molid = atoms.data().number().value();
        const qint32 s_molid = *(reinterpret_cast<const qint32*>(&molid));
        
        const int nats = coords.count();
        
        if (nats == 0)
            return QVector<CLJAtom>();
        
        //reserve space for the data
        cljatoms = QVector<CLJAtom>(nats);
        CLJAtom *cljatms = cljatoms.data();
        
        int idx = 0;
        
        for (int i=0; i<nats; ++i)
        {
            CLJAtom &cljatm = cljatms[idx];
        
            cljatm.x = coords[i].x();
            cljatm.y = coords[i].y();
            cljatm.z = coords[i].z();
            
            cljatm.chg = chgs[i].value() * std::sqrt(SireUnits::one_over_four_pi_eps0);

            cljatm.sig = std::sqrt(ljs[i].sigma());
            cljatm.eps = std::sqrt(4.0 * ljs[i].epsilon());
            
            cljatm.idnum = s_molid;
            
            idx += 1;
        }
    }
    
    //quint64 ns = t.nsecsElapsed();
    
    //qDebug() << "Converting" << cljatoms.count() << "atoms took"
    //         << (0.000001*ns) << "ms";
    
    return cljatoms;
}

/** Return the coordinates of the atom */
Vector CLJAtom::coordinates() const
{
    return Vector(x,y,z);
}

/** Return the partial charge of the atom */
Charge CLJAtom::charge() const
{
    return Charge(chg / std::sqrt(SireUnits::one_over_four_pi_eps0) );
}

/** Return the LJ parameters of the atom */
LJParameter CLJAtom::ljParameter() const
{
    double s = sig * sig;
    double e = eps * eps;
    e /= 4.0;
    
    return LJParameter( SireUnits::Dimension::Length(s), SireUnits::Dimension::MolarEnergy(e) );
}

/** Return the ID number for the atom */
qint32 CLJAtom::ID() const
{
    return idnum;
}

/** Return the negative of this atom - this returns a copy where
    the reduced charge and reduced epsilon values have been negated */
CLJAtom CLJAtom::negate() const
{
    CLJAtom ret(*this);
    ret.eps *= -1;
    ret.chg *= -1;
    return ret;
}

///////
/////// Implementation of CLJAtoms
///////

static const RegisterMetaType<CLJAtoms> r_cljatoms(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJAtoms &cljatoms)
{
    writeHeader(ds, r_cljatoms, 1);
    
    SharedDataStream sds(ds);
    sds << cljatoms.atoms();
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJAtoms &cljatoms)
{
    VersionID v = readHeader(ds, r_cljatoms);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        QVector<CLJAtom> atoms;
        
        sds >> atoms;
        
        cljatoms = CLJAtoms(atoms);
    }
    else
        throw version_error(v, "1", r_cljatoms, CODELOC);
    
    return ds;
}

qint32 CLJAtoms::id_of_dummy = 0;

/** Return a MultiFloat of the ID of a dummy atom */
MultiInt CLJAtoms::idOfDummy()
{
    return MultiInt(id_of_dummy);
}

/** Null constructor */
CLJAtoms::CLJAtoms()
{}

/** Constructor allowing implicit conversion from a single CLJAtom */
CLJAtoms::CLJAtoms(const CLJAtom &cljatom)
{
    if (cljatom.isNull())
        return;
/*
            xa[idx] = atm.x;
            ya[idx] = atm.y;
            za[idx] = atm.z;
            ca[idx] = atm.chg;
            sa[idx] = atm.sig;
            ea[idx] = atm.eps;
            ida[idx] = atm.idnum;
*/
    _x = MultiFloat::fromArray( &(cljatom.x), 1 );
    _y = MultiFloat::fromArray( &(cljatom.y), 1 );
    _z = MultiFloat::fromArray( &(cljatom.z), 1 );
    _q = MultiFloat::fromArray( &(cljatom.chg), 1 );
    _sig = MultiFloat::fromArray( &(cljatom.sig), 1 );
    _eps = MultiFloat::fromArray( &(cljatom.eps), 1 );
    _id = MultiInt::fromArray( &(cljatom.idnum), 1 );
}

/** Construct from the passed array of CLJAtom atoms */
CLJAtoms::CLJAtoms(const QVector<CLJAtom> &atoms)
{
    if (atoms.isEmpty())
        return;
    
    else if (atoms.count() == 1)
    {
        this->operator=( CLJAtoms(atoms.at(0)) );
        return;
    }
    
    /*QElapsedTimer t;
    t.start();*/
    
    //vectorise all of the parameters
    QVector<float> xf(atoms.count());
    QVector<float> yf(atoms.count());
    QVector<float> zf(atoms.count());
    
    QVector<float> cf(atoms.count());
    QVector<float> sf(atoms.count());
    QVector<float> ef(atoms.count());
    
    QVector<qint32> idf(atoms.count());
    
    float *xa = xf.data();
    float *ya = yf.data();
    float *za = zf.data();
    
    float *ca = cf.data();
    float *sa = sf.data();
    float *ea = ef.data();
    
    qint32 *ida = idf.data();
    
    const CLJAtom *atms = atoms.constData();
    
    int idx = 0;
    
    for (int i=0; i<atoms.count(); ++i)
    {
        const CLJAtom &atm = atms[i];
    
        if (atm.chg != 0 or atm.eps != 0)
        {
            xa[idx] = atm.x;
            ya[idx] = atm.y;
            za[idx] = atm.z;
            ca[idx] = atm.chg;
            sa[idx] = atm.sig;
            ea[idx] = atm.eps;
            ida[idx] = atm.idnum;
            
            idx += 1;
        }
    }
    
    if (idx > 0)
    {
        _x = MultiFloat::fromArray(xf.constData(), idx);
        _y = MultiFloat::fromArray(yf.constData(), idx);
        _z = MultiFloat::fromArray(zf.constData(), idx);
        
        _q = MultiFloat::fromArray(cf.constData(), idx);
        _sig = MultiFloat::fromArray(sf.constData(), idx);
        _eps = MultiFloat::fromArray(ef.constData(), idx);
        
        _id = MultiInt::fromArray(idf.constData(), idx);
    }
    
    /*quint64 ns = t.nsecsElapsed();

    qDebug() << "Converting" << (_q.count() * MultiFloat::count()) << "atoms took"
             << (0.000001*ns) << "ms";*/
}

/** Construct from the passed array of CLJAtom atoms */
CLJAtoms::CLJAtoms(const CLJAtom *atoms, int natoms)
{
    if (atoms == 0 or natoms <= 0)
        return;
    
    else if (natoms == 1)
    {
        this->operator=( CLJAtoms(atoms[0]) );
        return;
    }
    
    /*QElapsedTimer t;
    t.start();*/
    
    //vectorise all of the parameters
    QVector<float> xf(natoms);
    QVector<float> yf(natoms);
    QVector<float> zf(natoms);
    
    QVector<float> cf(natoms);
    QVector<float> sf(natoms);
    QVector<float> ef(natoms);
    
    QVector<qint32> idf(natoms);
    
    float *xa = xf.data();
    float *ya = yf.data();
    float *za = zf.data();
    
    float *ca = cf.data();
    float *sa = sf.data();
    float *ea = ef.data();
    
    qint32 *ida = idf.data();
    
    int idx = 0;
    
    for (int i=0; i<natoms; ++i)
    {
        const CLJAtom &atm = atoms[i];
    
        if (atm.chg != 0 or atm.eps != 0)
        {
            xa[idx] = atm.x;
            ya[idx] = atm.y;
            za[idx] = atm.z;
            ca[idx] = atm.chg;
            sa[idx] = atm.sig;
            ea[idx] = atm.eps;
            ida[idx] = atm.idnum;
            
            idx += 1;
        }
    }
    
    if (idx > 0)
    {
        _x = MultiFloat::fromArray(xf.constData(), idx);
        _y = MultiFloat::fromArray(yf.constData(), idx);
        _z = MultiFloat::fromArray(zf.constData(), idx);
        
        _q = MultiFloat::fromArray(cf.constData(), idx);
        _sig = MultiFloat::fromArray(sf.constData(), idx);
        _eps = MultiFloat::fromArray(ef.constData(), idx);
        
        _id = MultiInt::fromArray(idf.constData(), idx);
    }
    
    /*quint64 ns = t.nsecsElapsed();

    qDebug() << "Converting" << (_q.count() * MultiFloat::count()) << "atoms took"
             << (0.000001*ns) << "ms";*/
}

/** Construct from the passed list of CLJAtom atoms */
CLJAtoms::CLJAtoms(const QList<CLJAtom> &atoms)
{
    if (atoms.isEmpty())
        return;
    
    /*QElapsedTimer t;
    t.start();*/
    
    //vectorise all of the parameters
    QVector<float> xf(atoms.count());
    QVector<float> yf(atoms.count());
    QVector<float> zf(atoms.count());
    
    QVector<float> cf(atoms.count());
    QVector<float> sf(atoms.count());
    QVector<float> ef(atoms.count());
    
    QVector<qint32> idf(atoms.count());
    
    float *xa = xf.data();
    float *ya = yf.data();
    float *za = zf.data();
    
    float *ca = cf.data();
    float *sa = sf.data();
    float *ea = ef.data();
    
    qint32 *ida = idf.data();
    
    int idx = 0;
    
    for (int i=0; i<atoms.count(); ++i)
    {
        const CLJAtom &atm = atoms.at(i);
    
        if (atm.chg != 0 or atm.eps != 0)
        {
            xa[idx] = atm.x;
            ya[idx] = atm.y;
            za[idx] = atm.z;
            ca[idx] = atm.chg;
            sa[idx] = atm.sig;
            ea[idx] = atm.eps;
            ida[idx] = atm.idnum;
            
            idx += 1;
        }
    }
    
    if (idx > 0)
    {
        _x = MultiFloat::fromArray(xf.constData(), idx);
        _y = MultiFloat::fromArray(yf.constData(), idx);
        _z = MultiFloat::fromArray(zf.constData(), idx);
        
        _q = MultiFloat::fromArray(cf.constData(), idx);
        _sig = MultiFloat::fromArray(sf.constData(), idx);
        _eps = MultiFloat::fromArray(ef.constData(), idx);
        
        _id = MultiInt::fromArray(idf.constData(), idx);
    }
    
    /*quint64 ns = t.nsecsElapsed();

    qDebug() << "Converting" << (_q.count() * MultiFloat::count()) << "atoms took"
             << (0.000001*ns) << "ms";*/
}

/** Construct from the passed set of coordinates, partial charges and LJ parameters.
    Each atom is assumed to be part of the same molecule, with atom ID 'atomid' */
CLJAtoms::CLJAtoms(const QVector<Vector> &coordinates,
                   const QVector<Charge> &charges,
                   const QVector<LJParameter> &ljparams,
                   qint32 atomid)
{
    if (coordinates.count() != charges.count() or
        coordinates.count() != ljparams.count())
    {
        throw SireError::incompatible_error( QObject::tr(
                "You cannot construct a set of CLJAtoms where the number of coordinates "
                "(%1), charges (%2) and Lennard Jones parameters (%3) are different!")
                    .arg(coordinates.count())
                    .arg(charges.count())
                    .arg(ljparams.count()), CODELOC );
    }

    //vectorise the parameters, and convert to float as we can control
    //precision in the energy calculation
    {
        QVector<float> xf(coordinates.count());
        QVector<float> yf(coordinates.count());
        QVector<float> zf(coordinates.count());
        
        QVector<float> cf(charges.count());
        QVector<float> sf(ljparams.count());
        QVector<float> ef(ljparams.count());
        
        QVector<qint32> idf(coordinates.count());
        
        float *xa = xf.data();
        float *ya = yf.data();
        float *za = zf.data();
        
        float *ca = cf.data();
        float *sa = sf.data();
        float *ea = ef.data();
        
        qint32 *ida = idf.data();
        
        const Vector *coords = coordinates.constData();
        const Charge *chgs = charges.constData();
        const LJParameter *ljs = ljparams.constData();
        
        int idx = 0;
        
        for (int i=0; i<coordinates.count(); ++i)
        {
            if (chgs[i].value() != 0 or (not ljs[i].isDummy()))
            {
                xa[idx] = coords[i].x();
                ya[idx] = coords[i].y();
                za[idx] = coords[i].z();
                
                ca[idx] = chgs[i].value();
                sa[idx] = ljs[i].sigma();
                ea[idx] = ljs[i].epsilon();
                
                ida[idx] = atomid;
            
                idx += 1;
            }
        }
        
        if (idx > 0)
        {
            _x = MultiFloat::fromArray(xf.constData(), idx);
            _y = MultiFloat::fromArray(yf.constData(), idx);
            _z = MultiFloat::fromArray(zf.constData(), idx);
            
            _q = MultiFloat::fromArray(cf.constData(), idx);
            _sig = MultiFloat::fromArray(sf.constData(), idx);
            _eps = MultiFloat::fromArray(ef.constData(), idx);
            
            _id = MultiInt::fromArray(idf.constData(), idx);
        }
    }
    
    MultiFloat *q = _q.data();
    MultiFloat *s = _sig.data();
    MultiFloat *e = _eps.data();
    
    const MultiFloat four(4.0);
    const MultiFloat one_over_4_pi_eps_0( std::sqrt(SireUnits::one_over_four_pi_eps0) );
    
    //now reduce the charge, sigma and epsilon parameters
    for (int i=0; i<_q.count(); ++i)
    {
        q[i] = q[i] * one_over_4_pi_eps_0;
        s[i] = s[i].sqrt();
        e[i] = (e[i] * four).sqrt();
    }
}

/** Construct from the passed set of coordinates, partial charges, LJ parameters
    and atom IDs */
CLJAtoms::CLJAtoms(const QVector<Vector> &coordinates,
                   const QVector<Charge> &charges,
                   const QVector<LJParameter> &ljparams,
                   const QVector<qint32> &atomids)
{
    //QElapsedTimer t;
    //t.start();

    if (coordinates.count() != charges.count() or
        coordinates.count() != ljparams.count() or
        coordinates.count() != atomids.count())
    {
        throw SireError::incompatible_error( QObject::tr(
                "You cannot construct a set of CLJAtoms where the number of coordinates "
                "(%1), charges (%2), Lennard Jones parameters (%3) or "
                "atom IDs (%4) are different!")
                    .arg(coordinates.count())
                    .arg(charges.count())
                    .arg(ljparams.count())
                    .arg(atomids.count()), CODELOC );
    }

    //vectorise the parameters, and convert to float as we can control
    //precision in the energy calculation
    {
        QVector<float> xf(coordinates.count());
        QVector<float> yf(coordinates.count());
        QVector<float> zf(coordinates.count());
        
        QVector<float> cf(charges.count());
        QVector<float> sf(ljparams.count());
        QVector<float> ef(ljparams.count());
        
        QVector<qint32> idf(atomids.count());
        
        float *xa = xf.data();
        float *ya = yf.data();
        float *za = zf.data();
        
        float *ca = cf.data();
        float *sa = sf.data();
        float *ea = ef.data();
        
        qint32 *ida = idf.data();
        
        const Vector *coords = coordinates.constData();
        const Charge *chgs = charges.constData();
        const LJParameter *ljs = ljparams.constData();
        const qint32 *atmids = atomids.constData();
        
        int idx = 0;
        
        for (int i=0; i<coordinates.count(); ++i)
        {
            if (chgs[i].value() != 0 or (not ljs[i].isDummy()))
            {
                xa[idx] = coords[i].x();
                ya[idx] = coords[i].y();
                za[idx] = coords[i].z();
                
                ca[idx] = chgs[i].value();
                sa[idx] = ljs[i].sigma();
                ea[idx] = ljs[i].epsilon();
                
                ida[idx] = atmids[i];
            
                idx += 1;
            }
        }
        
        if (idx > 0)
        {
            _x = MultiFloat::fromArray(xf.constData(), idx);
            _y = MultiFloat::fromArray(yf.constData(), idx);
            _z = MultiFloat::fromArray(zf.constData(), idx);
            
            _q = MultiFloat::fromArray(cf.constData(), idx);
            _sig = MultiFloat::fromArray(sf.constData(), idx);
            _eps = MultiFloat::fromArray(ef.constData(), idx);
            
            _id = MultiInt::fromArray(idf.constData(), idx);
        }
    }
    
    MultiFloat *q = _q.data();
    MultiFloat *s = _sig.data();
    MultiFloat *e = _eps.data();
    
    const MultiFloat four(4.0);
    const MultiFloat one_over_4_pi_eps_0( std::sqrt(SireUnits::one_over_four_pi_eps0) );
    
    //now reduce the charge, sigma and epsilon parameters
    for (int i=0; i<_q.count(); ++i)
    {
        q[i] = q[i] * one_over_4_pi_eps_0;
        s[i] = s[i].sqrt();
        e[i] = (e[i] * four).sqrt();
    }
    
    //quint64 ns = t.nsecsElapsed();
    
    //qDebug() << "Converting" << (_q.count() * MultiFloat::count()) << "atoms took"
    //         << (0.000001*ns) << "ms";
}

/** Construct from the passed MoleculeView */
void CLJAtoms::constructFrom(const MoleculeView &molecule,
                             const ID_SOURCE id_source, const PropertyMap &map)
{
    if (molecule.isEmpty())
        return;

    //QElapsedTimer t;
    //t.start();
    
    //extract all of the data from the passed molecules
    {
        const PropertyName coords_property = map["coordinates"];
        const PropertyName chg_property = map["charge"];
        const PropertyName lj_property = map["LJ"];
        
        if (molecule.selectedAll())
        {
            Molecule mol = molecule.molecule();
            
            const int nats = mol.nAtoms();
            
            if (nats == 0)
                return;
            
            //reserve space for the data
            QVector<float> xf(nats);
            QVector<float> yf(nats);
            QVector<float> zf(nats);
            
            QVector<float> qf(nats);
            QVector<float> sigf(nats);
            QVector<float> epsf(nats);
            
            QVector<qint32> idf(nats);
            
            float *xa = xf.data();
            float *ya = yf.data();
            float *za = zf.data();
            
            float *qa = qf.data();
            float *siga = sigf.data();
            float *epsa = epsf.data();
            
            qint32 *ida = idf.data();
            
            int idx = 0;

            const AtomCoords &coords = mol.property(coords_property).asA<AtomCoords>();
            const AtomCharges &chgs = mol.property(chg_property).asA<AtomCharges>();
            const AtomLJs &ljs = mol.property(lj_property).asA<AtomLJs>();

            const qint32 s_molid = mol.number().value();
            
            const MoleculeInfoData &molinfo = mol.data().info();
            
            for (int i=0; i<coords.nCutGroups(); ++i)
            {
                const CGIdx cgidx(i);
            
                const Vector *icoords = coords.constData(cgidx);
                const Charge *ichg = chgs.constData(cgidx);
                const LJParameter *ilj = ljs.constData(cgidx);
                
                for (int j=0; j<coords.nAtoms(cgidx); ++j)
                {
                    if (ichg[j].value() != 0 or (not ilj[j].isDummy()))
                    {
                        xa[idx] = icoords[j].x();
                        ya[idx] = icoords[j].y();
                        za[idx] = icoords[j].z();

                        qa[idx] = ichg[j].value();

                        siga[idx] = ilj[j].sigma();
                        epsa[idx] = ilj[j].epsilon();
                        
                        if (id_source == USE_MOLNUM)
                        {
                            ida[idx] = s_molid;
                        }
                        else if (id_source == USE_ATOMIDX)
                        {
                            const AtomIdx atomidx = molinfo.atomIdx( CGAtomIdx(CGIdx(i),Index(j)) );
                            ida[idx] = atomidx.value() + 1;
                        }
                        else
                        {
                            throw SireError::program_bug( QObject::tr(
                                    "Unknown source used for ID (%1)")
                                        .arg(id_source), CODELOC );
                        }
                    
                        idx += 1;
                    }
                }
            }

            if (idx > 0)
            {
                _x = MultiFloat::fromArray(xf.constData(), idx);
                _y = MultiFloat::fromArray(yf.constData(), idx);
                _z = MultiFloat::fromArray(zf.constData(), idx);
                _q = MultiFloat::fromArray(qf.constData(), idx);
                _sig = MultiFloat::fromArray(sigf.constData(), idx);
                _eps = MultiFloat::fromArray(epsf.constData(), idx);
                _id = MultiInt::fromArray(idf.constData(), idx);
            }
        }
        else
        {
            Selector<Atom> atoms = molecule.atoms();
            
            QList<Vector> coords = atoms.property<Vector>(coords_property);
            QList<Charge> chgs = atoms.property<Charge>(chg_property);
            QList<LJParameter> ljs = atoms.property<LJParameter>(lj_property);
            
            if (coords.count() != chgs.count() or
                coords.count() != ljs.count())
            {
                throw SireError::program_bug( QObject::tr(
                        "It should not be possible for the number of coordinates (%1) "
                        "to be different to the number of charges (%2) or LJs (%3) "
                        "for molecule %4.")
                            .arg(coords.count())
                            .arg(chgs.count())
                            .arg(ljs.count())
                            .arg(atoms.molecule().toString()), CODELOC );
            }

            const qint32 s_molid = atoms.data().number().value();
            
            const int nats = coords.count();
            
            if (nats == 0)
                return;
            
            //reserve space for the data
            QVector<float> xf(nats);
            QVector<float> yf(nats);
            QVector<float> zf(nats);
            
            QVector<float> qf(nats);
            QVector<float> sigf(nats);
            QVector<float> epsf(nats);
            
            QVector<qint32> idf(nats);
            
            float *xa = xf.data();
            float *ya = yf.data();
            float *za = zf.data();
            
            float *qa = qf.data();
            float *siga = sigf.data();
            float *epsa = epsf.data();
            
            qint32 *ida = idf.data();
            
            int idx = 0;
            
            for (int i=0; i<nats; ++i)
            {
                if (chgs[i].value() != 0 or (not ljs[i].isDummy()))
                {
                    xa[idx] = coords[i].x();
                    ya[idx] = coords[i].y();
                    za[idx] = coords[i].z();
                    
                    qa[idx] = chgs[i].value();
                    siga[idx] = ljs[i].sigma();
                    epsa[idx] = ljs[i].epsilon();

                    if (id_source == USE_MOLNUM)
                    {
                        ida[idx] = s_molid;
                    }
                    else if (id_source == USE_ATOMIDX)
                    {
                        const AtomIdx atomidx = atoms[i].index();
                        ida[idx] = atomidx.value() + 1;
                    }
                    else
                    {
                        throw SireError::program_bug( QObject::tr(
                                "Unknown source used for ID (%1)")
                                    .arg(id_source), CODELOC );
                    }
                    
                    idx += 1;
                }
            }
        
            if (idx > 0)
            {
                _x = MultiFloat::fromArray(xf.constData(), idx);
                _y = MultiFloat::fromArray(yf.constData(), idx);
                _z = MultiFloat::fromArray(zf.constData(), idx);
                _q = MultiFloat::fromArray(qf.constData(), idx);
                _sig = MultiFloat::fromArray(sigf.constData(), idx);
                _eps = MultiFloat::fromArray(epsf.constData(), idx);
                _id = MultiInt::fromArray(idf.constData(), idx);
            }
        }
    }
    
    MultiFloat *q = _q.data();
    MultiFloat *s = _sig.data();
    MultiFloat *e = _eps.data();
    
    const MultiFloat four(4.0);
    const MultiFloat one_over_4_pi_eps_0( std::sqrt(SireUnits::one_over_four_pi_eps0) );
    
    //now reduce the charge, sigma and epsilon parameters
    for (int i=0; i<_q.count(); ++i)
    {
        q[i] = q[i] * one_over_4_pi_eps_0;
        s[i] = s[i].sqrt();
        e[i] = (e[i] * four).sqrt();
    }
    
    //quint64 ns = t.nsecsElapsed();
    
    //qDebug() << "Converting" << (_q.count() * MultiFloat::count()) << "atoms took"
    //         << (0.000001*ns) << "ms";
}

void CLJAtoms::reconstruct(const MoleculeView &molecule, const PropertyMap &map)
{
    this->constructFrom(molecule, USE_MOLNUM, map);
}

void CLJAtoms::reconstruct(const MoleculeView &molecule, ID_SOURCE source,
                           const PropertyMap &map)
{
    this->constructFrom(molecule, source, map);
}

/** Construct from the parameters in the passed set of Molecules */
void CLJAtoms::constructFrom(const Molecules &molecules,
                             ID_SOURCE id_source, const PropertyMap &map)
{
    if (molecules.isEmpty())
        return;
   
    //QElapsedTimer t;
    //t.start();
    
    //extract all of the data from the passed molecules
    {
        const PropertyName coords_property = map["coordinates"];
        const PropertyName chg_property = map["charge"];
        const PropertyName lj_property = map["LJ"];
        
        //calculate the number of atoms...
        int nats = 0;
        
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            nats += it.value().selection().nSelected();
        }
        
        //reserve space for the data
        QVector<float> xf(nats);
        QVector<float> yf(nats);
        QVector<float> zf(nats);
        
        QVector<float> qf(nats);
        QVector<float> sigf(nats);
        QVector<float> epsf(nats);
        
        QVector<qint32> idf(nats);
        
        float *xa = xf.data();
        float *ya = yf.data();
        float *za = zf.data();
        
        float *qa = qf.data();
        float *siga = sigf.data();
        float *epsa = epsf.data();
        
        qint32 *ida = idf.data();
        
        int idx = 0;
        
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            if (it.value().selectedAll())
            {
                Molecule mol = it.value().molecule();
                
                const AtomCoords &coords = mol.property(coords_property).asA<AtomCoords>();
                const AtomCharges &chgs = mol.property(chg_property).asA<AtomCharges>();
                const AtomLJs &ljs = mol.property(lj_property).asA<AtomLJs>();

                const qint32 s_molid = it.key().value();
                const MoleculeInfoData &molinfo = mol.data().info();
                
                for (int i=0; i<coords.nCutGroups(); ++i)
                {
                    const CGIdx cgidx(i);
                
                    const Vector *icoords = coords.constData(cgidx);
                    const Charge *ichg = chgs.constData(cgidx);
                    const LJParameter *ilj = ljs.constData(cgidx);
                    
                    for (int j=0; j<coords.nAtoms(cgidx); ++j)
                    {
                        if (ichg[j].value() != 0 or (not ilj[j].isDummy()))
                        {
                            xa[idx] = icoords[j].x();
                            ya[idx] = icoords[j].y();
                            za[idx] = icoords[j].z();

                            qa[idx] = ichg[j].value();

                            siga[idx] = ilj[j].sigma();
                            epsa[idx] = ilj[j].epsilon();
                            
                            if (id_source == USE_MOLNUM)
                            {
                                ida[idx] = s_molid;
                            }
                            else if (id_source == USE_ATOMIDX)
                            {
                                const AtomIdx atomidx
                                                = molinfo.atomIdx( CGAtomIdx(CGIdx(i),Index(j)) );
                                
                                ida[idx] = atomidx.value() + 1;
                            }
                            else
                            {
                                throw SireError::program_bug( QObject::tr(
                                        "Unknown source used for ID (%1)")
                                            .arg(id_source), CODELOC );
                            }
                        
                            idx += 1;
                        }
                    }
                }
            }
            else
            {
                Selector<Atom> atoms = it.value().atoms();
                
                QList<Vector> coords = atoms.property<Vector>(coords_property);
                QList<Charge> chgs = atoms.property<Charge>(chg_property);
                QList<LJParameter> ljs = atoms.property<LJParameter>(lj_property);
                
                if (coords.count() != chgs.count() or
                    coords.count() != ljs.count())
                {
                    throw SireError::program_bug( QObject::tr(
                            "It should not be possible for the number of coordinates (%1) "
                            "to be different to the number of charges (%2) or LJs (%3) "
                            "for molecule %4.")
                                .arg(coords.count())
                                .arg(chgs.count())
                                .arg(ljs.count())
                                .arg(atoms.molecule().toString()), CODELOC );
                }

                const qint32 s_molid = it.key().value();
                
                for (int i=0; i<coords.count(); ++i)
                {
                    if (chgs[i].value() != 0 or (not ljs[i].isDummy()))
                    {
                        xa[idx] = coords[i].x();
                        ya[idx] = coords[i].y();
                        za[idx] = coords[i].z();
                        
                        qa[idx] = chgs[i].value();
                        siga[idx] = ljs[i].sigma();
                        epsa[idx] = ljs[i].epsilon();

                        if (id_source == USE_MOLNUM)
                        {
                            ida[idx] = s_molid;
                        }
                        else if (id_source == USE_ATOMIDX)
                        {
                            ida[idx] = atoms[i].index().value() + 1;
                        }
                        else
                        {
                            throw SireError::program_bug( QObject::tr(
                                    "Unknown source used for ID (%1)")
                                        .arg(id_source), CODELOC );
                        }
                        
                        idx += 1;
                    }
                }
            }
        }
        
        if (idx > 0)
        {
            _x = MultiFloat::fromArray(xf.constData(), idx);
            _y = MultiFloat::fromArray(yf.constData(), idx);
            _z = MultiFloat::fromArray(zf.constData(), idx);
            _q = MultiFloat::fromArray(qf.constData(), idx);
            _sig = MultiFloat::fromArray(sigf.constData(), idx);
            _eps = MultiFloat::fromArray(epsf.constData(), idx);
            _id = MultiInt::fromArray(idf.constData(), idx);
        }
    }
    
    MultiFloat *q = _q.data();
    MultiFloat *s = _sig.data();
    MultiFloat *e = _eps.data();
    
    const MultiFloat four(4.0);
    const MultiFloat one_over_4_pi_eps_0( std::sqrt(SireUnits::one_over_four_pi_eps0) );
    
    //now reduce the charge, sigma and epsilon parameters
    for (int i=0; i<_q.count(); ++i)
    {
        q[i] = q[i] * one_over_4_pi_eps_0;
        s[i] = s[i].sqrt();
        e[i] = (e[i] * four).sqrt();
    }
    
    //quint64 ns = t.nsecsElapsed();
    
    //qDebug() << "Converting" << (_q.count() * MultiFloat::count()) << "atoms took"
    //         << (0.000001*ns) << "ms";
}

/** Construct from the parameters in the passed set of Molecules */
void CLJAtoms::constructFrom(const MoleculeGroup &molecules,
                             ID_SOURCE id_source, const PropertyMap &map)
{
    if (molecules.isEmpty())
        return;
   
    //QElapsedTimer t;
    //t.start();
    
    //extract all of the data from the passed molecules
    {
        const PropertyName coords_property = map["coordinates"];
        const PropertyName chg_property = map["charge"];
        const PropertyName lj_property = map["LJ"];
        
        //calculate the number of atoms...
        int nats = 0;
        
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            nats += it.value().selection().nSelected();
        }
        
        //reserve space for the data
        QVector<float> xf(nats);
        QVector<float> yf(nats);
        QVector<float> zf(nats);
        
        QVector<float> qf(nats);
        QVector<float> sigf(nats);
        QVector<float> epsf(nats);
        
        QVector<qint32> idf(nats);
        
        float *xa = xf.data();
        float *ya = yf.data();
        float *za = zf.data();
        
        float *qa = qf.data();
        float *siga = sigf.data();
        float *epsa = epsf.data();
        
        qint32 *ida = idf.data();
        
        int idx = 0;
        
        for (int i=0; i<molecules.nMolecules(); ++i)
        {
            const MoleculeView &view = molecules[MolIdx(i)];
        
            if (view.selectedAll())
            {
                Molecule mol = view.molecule();
                
                const AtomCoords &coords = mol.property(coords_property).asA<AtomCoords>();
                const AtomCharges &chgs = mol.property(chg_property).asA<AtomCharges>();
                const AtomLJs &ljs = mol.property(lj_property).asA<AtomLJs>();

                const qint32 s_molid = mol.number();
                const MoleculeInfoData &molinfo = mol.data().info();
                
                for (int i=0; i<coords.nCutGroups(); ++i)
                {
                    const CGIdx cgidx(i);
                
                    const Vector *icoords = coords.constData(cgidx);
                    const Charge *ichg = chgs.constData(cgidx);
                    const LJParameter *ilj = ljs.constData(cgidx);
                    
                    for (int j=0; j<coords.nAtoms(cgidx); ++j)
                    {
                        if (ichg[j].value() != 0 or (not ilj[j].isDummy()))
                        {
                            xa[idx] = icoords[j].x();
                            ya[idx] = icoords[j].y();
                            za[idx] = icoords[j].z();

                            qa[idx] = ichg[j].value();

                            siga[idx] = ilj[j].sigma();
                            epsa[idx] = ilj[j].epsilon();
                            
                            if (id_source == USE_MOLNUM)
                            {
                                ida[idx] = s_molid;
                            }
                            else if (id_source == USE_ATOMIDX)
                            {
                                const AtomIdx atomidx
                                                = molinfo.atomIdx( CGAtomIdx(CGIdx(i),Index(j)) );
                                
                                ida[idx] = atomidx.value() + 1;
                            }
                            else
                            {
                                throw SireError::program_bug( QObject::tr(
                                        "Unknown source used for ID (%1)")
                                            .arg(id_source), CODELOC );
                            }
                        
                            idx += 1;
                        }
                    }
                }
            }
            else
            {
                Selector<Atom> atoms = view.atoms();
                
                QList<Vector> coords = atoms.property<Vector>(coords_property);
                QList<Charge> chgs = atoms.property<Charge>(chg_property);
                QList<LJParameter> ljs = atoms.property<LJParameter>(lj_property);
                
                if (coords.count() != chgs.count() or
                    coords.count() != ljs.count())
                {
                    throw SireError::program_bug( QObject::tr(
                            "It should not be possible for the number of coordinates (%1) "
                            "to be different to the number of charges (%2) or LJs (%3) "
                            "for molecule %4.")
                                .arg(coords.count())
                                .arg(chgs.count())
                                .arg(ljs.count())
                                .arg(atoms.molecule().toString()), CODELOC );
                }

                const qint32 s_molid = view.data().number();
                
                for (int i=0; i<coords.count(); ++i)
                {
                    if (chgs[i].value() != 0 or (not ljs[i].isDummy()))
                    {
                        xa[idx] = coords[i].x();
                        ya[idx] = coords[i].y();
                        za[idx] = coords[i].z();
                        
                        qa[idx] = chgs[i].value();
                        siga[idx] = ljs[i].sigma();
                        epsa[idx] = ljs[i].epsilon();

                        if (id_source == USE_MOLNUM)
                        {
                            ida[idx] = s_molid;
                        }
                        else if (id_source == USE_ATOMIDX)
                        {
                            ida[idx] = atoms[i].index().value() + 1;
                        }
                        else
                        {
                            throw SireError::program_bug( QObject::tr(
                                    "Unknown source used for ID (%1)")
                                        .arg(id_source), CODELOC );
                        }
                        
                        idx += 1;
                    }
                }
            }
        }
        
        if (idx > 0)
        {
            _x = MultiFloat::fromArray(xf.constData(), idx);
            _y = MultiFloat::fromArray(yf.constData(), idx);
            _z = MultiFloat::fromArray(zf.constData(), idx);
            _q = MultiFloat::fromArray(qf.constData(), idx);
            _sig = MultiFloat::fromArray(sigf.constData(), idx);
            _eps = MultiFloat::fromArray(epsf.constData(), idx);
            _id = MultiInt::fromArray(idf.constData(), idx);
        }
    }
    
    MultiFloat *q = _q.data();
    MultiFloat *s = _sig.data();
    MultiFloat *e = _eps.data();
    
    const MultiFloat four(4.0);
    const MultiFloat one_over_4_pi_eps_0( std::sqrt(SireUnits::one_over_four_pi_eps0) );
    
    //now reduce the charge, sigma and epsilon parameters
    for (int i=0; i<_q.count(); ++i)
    {
        q[i] = q[i] * one_over_4_pi_eps_0;
        s[i] = s[i].sqrt();
        e[i] = (e[i] * four).sqrt();
    }
    
    //quint64 ns = t.nsecsElapsed();
    
    //qDebug() << "Converting" << (_q.count() * MultiFloat::count()) << "atoms took"
    //         << (0.000001*ns) << "ms";
}

/** Construct from the parameters in the passed molecule view */
CLJAtoms::CLJAtoms(const MoleculeView &view, const PropertyMap &map)
{
    constructFrom(view, USE_MOLNUM, map);
}

/** Construct from the parameters in the passed set of Molecules */
CLJAtoms::CLJAtoms(const Molecules &molecules, const PropertyMap &map)
{
    constructFrom(molecules, USE_MOLNUM, map);
}

/** Construct from the parameters in the passed set of Molecules */
CLJAtoms::CLJAtoms(const MoleculeGroup &molecules, const PropertyMap &map)
{
    constructFrom(molecules, USE_MOLNUM, map);
}

/** Construct from the parameters in the passed molecule view, specifying
    how the ID number should be obtained for each atom */
CLJAtoms::CLJAtoms(const MoleculeView &view, ID_SOURCE id_source, const PropertyMap &map)
{
    constructFrom(view, id_source, map);
}

/** Construct from the parameters in the passed molecule view, specifying
    how the ID number should be obtained for each atom */
CLJAtoms::CLJAtoms(const Molecules &molecules, ID_SOURCE id_source, const PropertyMap &map)
{
    constructFrom(molecules, id_source, map);
}

/** Construct from the parameters in the passed molecule view, specifying
    how the ID number should be obtained for each atom */
CLJAtoms::CLJAtoms(const MoleculeGroup &molecules, ID_SOURCE id_source, const PropertyMap &map)
{
    constructFrom(molecules, id_source, map);
}

/** Copy constructor */
CLJAtoms::CLJAtoms(const CLJAtoms &other)
         : _x(other._x), _y(other._y), _z(other._z),
           _q(other._q), _sig(other._sig), _eps(other._eps),
           _id(other._id)
{}

/** Destructor */
CLJAtoms::~CLJAtoms()
{}

/** Copy assignment operator */
CLJAtoms& CLJAtoms::operator=(const CLJAtoms &other)
{
    if (this != &other)
    {
        _x = other._x;
        _y = other._y;
        _z = other._z;
        _q = other._q;
        _sig = other._sig;
        _eps = other._eps;
        _id = other._id;
    }
    
    return *this;
}

/** Comparison operator */
bool CLJAtoms::operator==(const CLJAtoms &other) const
{
    return this == &other or
           (_x == other._x and _y == other._y and _z == other._z and
            _q == other._q and _sig == other._sig and _eps == other._eps and
            _id == other._id);
}

/** Comparison operator */
bool CLJAtoms::operator!=(const CLJAtoms &other) const
{
    return not operator==(other);
}

const char* CLJAtoms::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJAtoms>() );
}

const char* CLJAtoms::what() const
{
    return CLJAtoms::typeName();
}

/** Append the passed atom onto this array. Note that this is not space
    efficient as the passed atom will be padded up to the size of the vector */
CLJAtoms& CLJAtoms::operator+=(const CLJAtom &atom)
{
    const int idx = _x.count() * MultiFloat::count();

    _x.append( MultiFloat(0) );
    _y.append( MultiFloat(0) );
    _z.append( MultiFloat(0) );
    _q.append( MultiFloat(0) );
    _sig.append( MultiFloat(0) );
    _eps.append( MultiFloat(0) );
    _id.append( MultiFloat(0) );
    
    if (not atom.isDummy())
    {
        this->set(idx, atom);
    }
    
    return *this;
}

/** Append the passed array of CLJAtoms onto this array */
CLJAtoms& CLJAtoms::operator+=(const CLJAtoms &other)
{
    if (_x.isEmpty())
    {
        this->operator=(other);
    }
    else if (not other._x.isEmpty())
    {
        _x += other._x;
        _y += other._y;
        _z += other._z;
        _q += other._q;
        _sig += other._sig;
        _eps += other._eps;
        _id += other._id;
    }
    
    return *this;
}

/** Append the passed array of CLJAtom objects onto this array */
CLJAtoms& CLJAtoms::operator+=(const QVector<CLJAtom> &atoms)
{
    return this->operator+=( CLJAtoms(atoms) );
}

/** Return the combination of this CLJAtoms with 'other' pasted
    onto the end */
CLJAtoms CLJAtoms::operator+(const CLJAtoms &other) const
{
    CLJAtoms ret(*this);
    ret += other;
    return ret;
}

/** Return the combination of this CLJAtoms with 'atom' pasted onto the end */
CLJAtoms CLJAtoms::operator+(const CLJAtom &atom) const
{
    CLJAtoms ret(*this);
    ret += atom;
    return ret;
}

/** Return the combination of this CLJAToms with the array of CLJAtoms in 'atoms'
    pasted onto the end */
CLJAtoms CLJAtoms::operator+(const QVector<CLJAtom> &atoms) const
{
    CLJAtoms ret(*this);
    ret += atoms;
    return ret;
}

/** Append the passed atom onto the end of this set. Note that this is not space efficient
    as the atom will be added with a lot of padding */
void CLJAtoms::append(const CLJAtom &atom)
{
    this->operator+=(atom);
}

/** Append the passed set of atoms onto the end of this set. If n is greater than or
    equal to zero, then only 'n' of the passed atoms will be added onto this set */
void CLJAtoms::append(const CLJAtoms &other, int n)
{
    if (n == 0)
        return;
    
    else if (n < 0 or n >= other.count())
    {
        this->operator+=(other);
        return;
    }
    else
    {
        //how many partial vectors need to be added?
        const int n_partial = n % MultiFloat::count();
    
        //how many whole vectors need to be added?
        const int n_whole = n / MultiFloat::count();
        
        //start index of added vectors
        const int start_idx = _x.count();
        
        //what will be the new size of the vector?
        const int new_size = start_idx + n_whole + (n_partial == 0 ? 0 : 1);
        
        _x.resize(new_size);
        _y.resize(new_size);
        _z.resize(new_size);
        _q.resize(new_size);
        _sig.resize(new_size);
        _eps.resize(new_size);
        _id.resize(new_size);
        
        if (n_whole > 0)
        {
            std::memcpy( &(_x[start_idx]), &(other._x[0]), n_whole * sizeof(MultiFloat) );
            std::memcpy( &(_y[start_idx]), &(other._y[0]), n_whole * sizeof(MultiFloat) );
            std::memcpy( &(_z[start_idx]), &(other._z[0]), n_whole * sizeof(MultiFloat) );
            std::memcpy( &(_q[start_idx]), &(other._q[0]), n_whole * sizeof(MultiFloat) );
            std::memcpy( &(_sig[start_idx]), &(other._sig[0]), n_whole * sizeof(MultiFloat) );
            std::memcpy( &(_eps[start_idx]), &(other._eps[0]), n_whole * sizeof(MultiFloat) );
            std::memcpy( &(_id[start_idx]), &(other._id[0]), n_whole * sizeof(MultiInt) );
        }
        
        if (n_partial > 0)
        {
            for (int i=0; i<n_partial; ++i)
            {
                _x.last().set(i, other._x[n_whole][i]);
                _y.last().set(i, other._y[n_whole][i]);
                _z.last().set(i, other._z[n_whole][i]);
                _q.last().set(i, other._q[n_whole][i]);
                _sig.last().set(i, other._sig[n_whole][i]);
                _eps.last().set(i, other._eps[n_whole][i]);
                _id.last().set(i, other._id[n_whole][i]);
            }
        }
    }
}

/** Return whether or not this array is empty */
bool CLJAtoms::isEmpty() const
{
    return _x.isEmpty();
}

/** Return the number of atoms in this set. Note that vectorisation
    may mean that the array of atoms has been padded with dummy atoms */
int CLJAtoms::count() const
{
    return MultiFloat::size() * _x.count();
}

/** Return the number of atoms in this set. Note that vectorisation
    may mean that the array of atoms has been padded with dummy atoms */
int CLJAtoms::size() const
{
    return count();
}

QString CLJAtoms::toString() const
{
    QStringList lines;

    foreach( CLJAtom atom, this->atoms() )
    {
        lines.append( atom.toString() );
    }

    if (lines.isEmpty())
        return "CLJAtoms()";

    return QObject::tr("CLJAtoms( %1 )").arg(lines.join("\n"));
}

/** Return the ith atom in the vector */
CLJAtom CLJAtoms::operator[](int i) const
{
    i = SireID::Index(i).map(count());
    
    int idx = i / MultiFloat::count();
    int sub_idx = i % MultiFloat::count();
    
    CLJAtom atom;
    atom.x = _x[idx][sub_idx];
    atom.y = _y[idx][sub_idx];
    atom.z = _z[idx][sub_idx];
    atom.chg = _q[idx][sub_idx];
    atom.sig = _sig[idx][sub_idx];
    atom.eps = _eps[idx][sub_idx];
    atom.idnum = _id[idx][sub_idx];
    
    return atom;
}

/** Return the ith atom in the vector */
CLJAtom CLJAtoms::at(int i) const
{
    return operator[](i);
}

/** Return the ith atom in the vector */
CLJAtom CLJAtoms::getitem(int i) const
{
    return operator[](i);
}

/** Overwrite the atom at index i with the data in 'atom' */
void CLJAtoms::set(int i, const CLJAtom &atom)
{
    i = SireID::Index(i).map(count());
    
    int idx = i / MultiFloat::count();
    int sub_idx = i % MultiFloat::count();
    
    _x[idx].set(sub_idx, atom.x);
    _y[idx].set(sub_idx, atom.y);
    _z[idx].set(sub_idx, atom.z);
    _q[idx].set(sub_idx, atom.chg);
    _sig[idx].set(sub_idx, atom.sig);
    _eps[idx].set(sub_idx, atom.eps);
    _id[idx].set(sub_idx, atom.idnum);
}

/** Set the coordinates of the ith atom to 'coords' */
void CLJAtoms::setCoordinates(int i, Vector coords)
{
    i = SireID::Index(i).map(count());
    
    int idx = i / MultiFloat::count();
    int sub_idx = i % MultiFloat::count();
    
    _x[idx].set(sub_idx, coords.x());
    _y[idx].set(sub_idx, coords.y());
    _z[idx].set(sub_idx, coords.z());
}

/** Set the charge of the ith atom to 'charge' */
void CLJAtoms::setCharge(int i, Charge charge)
{
    float c = charge.value() * std::sqrt(SireUnits::one_over_four_pi_eps0);
    
    int idx = i / MultiFloat::count();
    int sub_idx = i % MultiFloat::count();

    _q[idx].set(sub_idx, c);
}

/** Set the LJ parameter of the ith atom to 'ljparam' */
void CLJAtoms::setLJParameter(int i, LJParameter ljparam)
{
    float s = ljparam.sigma();
    float e = ljparam.epsilon() * 4.0;
    
    int idx = i / MultiFloat::count();
    int sub_idx = i % MultiFloat::count();

    _sig[idx].set(sub_idx, sqrt(s));
    _eps[idx].set(sub_idx, sqrt(e));
}

/** Set the ID number for the ith atom to 'idnum' */
void CLJAtoms::setID(int i, qint32 idnum)
{
    int idx = i / MultiFloat::count();
    int sub_idx = i % MultiFloat::count();

    _id[idx].set(sub_idx, idnum);
}

/** Make the ith atom into a dummy atom (set the atom ID to 0) */
void CLJAtoms::makeDummy(int i)
{
    int idx = i / MultiFloat::count();
    int sub_idx = i % MultiFloat::count();
    _id[idx].set(sub_idx, 0);
}

/** Set the ID number of all (non-dummy) atoms to 'idnum' */
void CLJAtoms::setAllID(qint32 idnum)
{
    qint32 dummy_id = idOfDummy()[0];

    for (int i=0; i<_id.count(); ++i)
    {
        for (int j=0; j<MultiInt::count(); ++j)
        {
            if (_id[i][j] != dummy_id)
            {
                _id[i].set(j, idnum);
            }
        }
    }
}

/** Resize this collection to hold 'n' atoms. This will add dummy atoms (and padding)
    if necessary, or will delete elements (but keeping dummy padding) if needed */
void CLJAtoms::resize(int new_size)
{
    if (new_size == 0)
    {
        _x.resize(0);
        _y.resize(0);
        _z.resize(0);
        _q.resize(0);
        _sig.resize(0);
        _eps.resize(0);
        _id.resize(0);
        return;
    }

    const int old_natoms = this->nAtoms();

    if (new_size == old_natoms)
        //nothing needs to be done
        return;

    const int new_remainder = new_size % MultiFloat::count();
    const int new_nvectors = new_size / MultiFloat::count() + (new_remainder == 0 ? 0 : 1);

    //resize the actual vectors
    if (new_nvectors != _x.count())
    {
        _x.resize(new_nvectors);
        _y.resize(new_nvectors);
        _z.resize(new_nvectors);
        _q.resize(new_nvectors);
        _sig.resize(new_nvectors);
        _eps.resize(new_nvectors);
        _id.resize(new_nvectors);
    }
    
    //padd with dummies (if needed)
    if (new_remainder != 0)
    {
        for (int i=new_remainder; i<MultiFloat::count(); ++i)
        {
            _id.last().set(i, id_of_dummy);
        }
    }
    
    if (this->count() < new_size)
    {
        throw SireError::program_bug( QObject::tr(
                "Something went wrong resizing CLJAtoms from size %1 to size %2. "
                "New size is %3.")
                    .arg(old_natoms).arg(new_size).arg(this->count()), CODELOC );
    }
}

/** Copy the contents of 'other' into this vectors, starting at index 0,
    and copying all elements of 'other' */
void CLJAtoms::copyIn(const CLJAtoms &other)
{
    //straight copy
    if (other._x.count() >= _x.count())
    {
        this->operator=(other);
        return;
    }
    
    const int nelements = other._x.count();
    
    std::memcpy( _x.data(), other._x.constData(), nelements * sizeof(MultiFloat) );
    std::memcpy( _y.data(), other._y.constData(), nelements * sizeof(MultiFloat) );
    std::memcpy( _z.data(), other._z.constData(), nelements * sizeof(MultiFloat) );
    std::memcpy( _q.data(), other._q.constData(), nelements * sizeof(MultiFloat) );
    std::memcpy( _sig.data(), other._sig.constData(), nelements * sizeof(MultiFloat) );
    std::memcpy( _eps.data(), other._eps.constData(), nelements * sizeof(MultiFloat) );
    std::memcpy( _id.data(), other._id.constData(), nelements * sizeof(MultiInt) );
}

/** Return a copy of these CLJAtoms where the charge and LJ epsilon parameters
    are negated. This will mean that the negative of the energy of these CLJAtoms
    will be calculated by the CLJFunctions (useful for calculating energy differences) */
CLJAtoms CLJAtoms::negate() const
{
    CLJAtoms ret(*this);
    
    for (int i=0; i<_q.count(); ++i)
    {
        ret._q[i] = -(_q[i]);
        ret._eps[i] = -(_eps[i]);
    }
    
    return ret;
}

/** Return a squeezed copy of these CLJAtoms whereby all of the 
    dummy atoms are removed and atoms squeezed into a single, contiguous space */
CLJAtoms CLJAtoms::squeeze() const
{
    QVector<CLJAtom> atms = this->atoms();
    
    QMutableVectorIterator<CLJAtom> it(atms);
    
    qint32 dummy_id = CLJAtoms::idOfDummy()[0];
    
    while (it.hasNext())
    {
        if (it.next().ID() == dummy_id)
        {
            it.remove();
        }
    }
    
    return CLJAtoms(atms);
}

/** Return an array of all of the atoms */
QVector<CLJAtom> CLJAtoms::atoms() const
{
    if (this->isEmpty())
        return QVector<CLJAtom>();

    QVector<CLJAtom> atms( _x.count() * MultiFloat::count() );
    CLJAtom *a = atms.data();
    
    int idx = 0;
    
    for (int i=0; i<_x.count(); ++i)
    {
        const MultiFloat &xf = _x[i];
        const MultiFloat &yf = _y[i];
        const MultiFloat &zf = _z[i];
        const MultiFloat &qf = _q[i];
        const MultiFloat &sigf = _sig[i];
        const MultiFloat &epsf = _eps[i];
        const MultiInt &idf = _id[i];
        
        for (int j=0; j<MultiFloat::count(); ++j)
        {
            CLJAtom &atom = a[idx];
            ++idx;
            
            atom.x = xf[j];
            atom.y = yf[j];
            atom.z = zf[j];
            atom.chg = qf[j];
            atom.sig = sigf[j];
            atom.eps = epsf[j];
            atom.idnum = idf[j];
        }
    }
    
    return atms;
}

/** Return the coordinates of all of the atoms */
QVector<Vector> CLJAtoms::coordinates() const
{
    if (this->isEmpty())
        return QVector<Vector>();

    QVector<Vector> coords;
    coords.reserve( _x.count() * MultiFloat::count() );
    
    for (int i=0; i<_x.count(); ++i)
    {
        const MultiFloat &xf = _x[i];
        const MultiFloat &yf = _y[i];
        const MultiFloat &zf = _z[i];

        if (i < (_x.count() -1))
        {
            for (int j=0; j<MultiFloat::count(); ++j)
            {
                coords.append( Vector(xf[j], yf[j], zf[j]) );
            }
        }
        else
        {
            //make sure we don't add any padded atoms
            int npadding = 0;
            
            for (int j=MultiFloat::count()-1; j>=0; --j)
            {
                if (_id[i][j] == 0)
                {
                    npadding += 1;
                }
            }
            
            for (int j=0; j<MultiFloat::count()-npadding; ++j)
            {
                coords.append( Vector(xf[j], yf[j], zf[j]) );
            }
        }
    }
    
    return coords;
}

/** Return the charges of all of the atoms */
QVector<Charge> CLJAtoms::charges() const
{
    if (this->isEmpty())
        return QVector<Charge>();

    QVector<Charge> chgs;
    chgs.reserve( _x.count() * MultiFloat::count() );
    
    for (int i=0; i<_x.count(); ++i)
    {
        const MultiFloat &qf = _q[i];

        if (i < (_x.count() - 1))
        {
            for (int j=0; j<MultiFloat::count(); ++j)
            {
                chgs.append( Charge( qf[j] / std::sqrt(SireUnits::one_over_four_pi_eps0) ) );
            }
        }
        else
        {
            //make sure we don't add any padded atoms
            int npadding = 0;
            
            for (int j=MultiFloat::count()-1; j>=0; --j)
            {
                if (_id[i][j] == 0)
                {
                    npadding += 1;
                }
            }
            
            for (int j=0; j<MultiFloat::count()-npadding; ++j)
            {
                chgs.append( Charge( qf[j] / std::sqrt(SireUnits::one_over_four_pi_eps0) ) );
            }
        }
    }
    
    return chgs;
}

/** Return the charges of all of the atoms */
QVector<LJParameter> CLJAtoms::ljParameters() const
{
    if (this->isEmpty())
        return QVector<LJParameter>();

    QVector<LJParameter> ljs;
    ljs.reserve( _x.count() * MultiFloat::count() );
    
    for (int i=0; i<_x.count(); ++i)
    {
        const MultiFloat &sigf = _sig[i];
        const MultiFloat &epsf = _eps[i];

        if (i < (_x.count() - 1))
        {
            for (int j=0; j<MultiFloat::count(); ++j)
            {
                double sig = sigf[j] * sigf[j];
                double eps = epsf[j] * epsf[j];

                ljs.append( LJParameter(SireUnits::Dimension::Length(sig),
                                        SireUnits::Dimension::MolarEnergy(eps / 4.0) ) );
            }
        }
        else
        {
            //make sure we don't add any padded atoms
            int npadding = 0;
            
            for (int j=MultiFloat::count()-1; j>=0; --j)
            {
                if (_id[i][j] == 0)
                {
                    npadding += 1;
                }
            }
            
            for (int j=0; j<MultiFloat::count()-npadding; ++j)
            {
                double sig = sigf[j] * sigf[j];
                double eps = epsf[j] * epsf[j];

                ljs.append( LJParameter(SireUnits::Dimension::Length(sig),
                                        SireUnits::Dimension::MolarEnergy(eps / 4.0) ) );
            }
        }
    }
    
    return ljs;
}

/** Return the IDs of all of the atoms */
QVector<qint32> CLJAtoms::IDs() const
{
    if (this->isEmpty())
        return QVector<qint32>();

    QVector<qint32> ids( _id.count() * MultiFloat::count() );
    qint32 *idval = ids.data();

    int idx = 0;

    for (int i=0; i<_id.count(); ++i)
    {
        const MultiInt &idf = _id[i];

        for (int j=0; j<MultiInt::count(); ++j)
        {
            idval[idx] = idf[j];
            ++idx;
        }
    }

    return ids;
}

/** Return whether or not there are any dummy (or padded) atoms in this set */
bool CLJAtoms::hasDummies() const
{
    if (this->isEmpty())
        return false;
    
    for (int i=0; i<_id.count(); ++i)
    {
        const MultiInt &idf = _id[i];
        
        for (int j=0; j<MultiInt::count(); ++j)
        {
            if (idf[j] == id_of_dummy)
                return true;
        }
    }
    
    return false;
}

/** Return the number of dummy (or padded) atoms in this set */
int CLJAtoms::nDummies() const
{
    if (this->isEmpty())
        return 0;
    
    int ndummies = 0;
    
    for (int i=0; i<_id.count(); ++i)
    {
        const MultiInt &idf = _id[i];
        
        for (int j=0; j<MultiInt::count(); ++j)
        {
            if (idf[j] == id_of_dummy)
                ndummies += 1;
        }
    }
    
    return ndummies;
}

/** Return the number of non-dummy atoms in this set. This is equal to
    count() - nDummies() */
int CLJAtoms::nAtoms() const
{
    return count() - nDummies();
}

/** Return the minimum coordinates of these atoms (ignoring dummies) */
Vector CLJAtoms::minCoords() const
{
    Vector mincoords( std::numeric_limits<double>::max() );
    
    for (int i=0; i<_id.count(); ++i)
    {
        const MultiFloat &xf = _x[i];
        const MultiFloat &yf = _y[i];
        const MultiFloat &zf = _z[i];
        const MultiInt &idf = _id[i];
        
        for (int j=0; j<MultiInt::count(); ++j)
        {
            if (idf[j] != id_of_dummy)
            {
                mincoords.setMin( Vector(xf[j],yf[j],zf[j]) );
            }
        }
    }
    
    return mincoords;
}

/** Return the maximum coordinates of these atoms (ignoring dummies) */
Vector CLJAtoms::maxCoords() const
{
    Vector maxcoords( -std::numeric_limits<double>::max() );
    
    for (int i=0; i<_id.count(); ++i)
    {
        const MultiFloat &xf = _x[i];
        const MultiFloat &yf = _y[i];
        const MultiFloat &zf = _z[i];
        const MultiInt &idf = _id[i];
        
        for (int j=0; j<MultiInt::count(); ++j)
        {
            if (idf[j] != id_of_dummy)
            {
                maxcoords.setMax( Vector(xf[j],yf[j],zf[j]) );
            }
        }
    }
    
    return maxcoords;
}
