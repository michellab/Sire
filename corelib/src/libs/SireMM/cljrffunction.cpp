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

#include "cljrffunction.h"

#include "SireMaths/multifloat.h"
#include "SireMaths/multidouble.h"
#include "SireMaths/multiint.h"

#include "SireBase/numberproperty.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireVol/gridinfo.h"

#include <QElapsedTimer>
#include <QDebug>

using namespace SireMM;
using namespace SireMaths;
using namespace SireVol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireStream;

const float default_dielectric = 1.0;

/////////
///////// Implementation of CLJRFFunction
/////////

static const RegisterMetaType<CLJRFFunction> r_shift;

QDataStream &operator<<(QDataStream &ds, const CLJRFFunction &func)
{
    writeHeader(ds, r_shift, 1);
    
    ds << func.diel << static_cast<const CLJCutoffFunction&>(func);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJRFFunction &func)
{
    VersionID v = readHeader(ds, r_shift);
    
    if (v == 1)
    {
        ds >> func.diel >> static_cast<CLJCutoffFunction&>(func);
    }
    else
        throw version_error(v, "1", r_shift, CODELOC);
    
    return ds;
}

CLJRFFunction::CLJRFFunction()
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(),
                   diel(default_dielectric)
{}

CLJFunctionPtr CLJRFFunction::defaultRFFunction()
{
    static CLJFunctionPtr ptr( new CLJRFFunction() );
    return ptr;
}

CLJRFFunction::CLJRFFunction(Length cutoff)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(cutoff),
                   diel(default_dielectric)
{}

CLJRFFunction::CLJRFFunction(Length coul_cutoff, Length lj_cutoff)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(coul_cutoff, lj_cutoff),
                   diel(default_dielectric)
{}

CLJRFFunction::CLJRFFunction(const Space &space, Length cutoff)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(space, cutoff),
                   diel(default_dielectric)
{}

CLJRFFunction::CLJRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(space, coul_cutoff,
                                                                        lj_cutoff),
                   diel(default_dielectric)
{}

CLJRFFunction::CLJRFFunction(Length cutoff, COMBINING_RULES combining_rules)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(cutoff, combining_rules),
                   diel(default_dielectric)
{}

CLJRFFunction::CLJRFFunction(Length coul_cutoff, Length lj_cutoff,
                                   COMBINING_RULES combining_rules)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(
                                   coul_cutoff, lj_cutoff, combining_rules),
                   diel(default_dielectric)
{}

CLJRFFunction::CLJRFFunction(const Space &space, COMBINING_RULES combining_rules)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(space, combining_rules),
                   diel(default_dielectric)
{}

CLJRFFunction::CLJRFFunction(const Space &space, Length cutoff,
                                   COMBINING_RULES combining_rules)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(
                                   space, cutoff, combining_rules),
                   diel(default_dielectric)
{}

CLJRFFunction::CLJRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff,
                                   COMBINING_RULES combining_rules)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(
                                   space, coul_cutoff, lj_cutoff, combining_rules),
                   diel(default_dielectric)
{}

/** Copy constructor */
CLJRFFunction::CLJRFFunction(const CLJRFFunction &other)
                 : ConcreteProperty<CLJRFFunction,CLJCutoffFunction>(other),
                   diel(other.diel)
{}

/** Destructor */
CLJRFFunction::~CLJRFFunction()
{}

/** Copy assignment operator */
CLJRFFunction& CLJRFFunction::operator=(const CLJRFFunction &other)
{
    diel = other.diel;
    CLJCutoffFunction::operator=(other);
    return *this;
}

/** Comparison operator */
bool CLJRFFunction::operator==(const CLJRFFunction &other) const
{
    return diel == other.diel and CLJCutoffFunction::operator==(other);
}

/** Comparison operator */
bool CLJRFFunction::operator!=(const CLJRFFunction &other) const
{
    return not operator==(other);
}

const char* CLJRFFunction::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJRFFunction>() );
}

const char* CLJRFFunction::what() const
{
    return CLJRFFunction::typeName();
}

CLJRFFunction* CLJRFFunction::clone() const
{
    return new CLJRFFunction(*this);
}

QString CLJRFFunction::toString() const
{
    if (this->hasCutoff())
        return QObject::tr("CLJRFFunction( dielectric() == %4, coulombCutoff() == %1 A, "
                           "ljCutoff() == %2 A, space() == %3 )")
            .arg(coulombCutoff().to(angstrom))
            .arg(ljCutoff().to(angstrom))
            .arg(space().toString())
            .arg(dielectric());
    else
        return QObject::tr("CLJRFFunction( dielectric() == %2, no cutoff, space() == %1 )")
                    .arg(space().toString()).arg(dielectric());
}

/** Return the properties of this function */
Properties CLJRFFunction::properties() const
{
    Properties props = CLJCutoffFunction::properties();
    props.setProperty("dielectric", NumberProperty(dielectric()));
    return props;
}

/** Return a copy of this function where the property 'name' has been set to the
    value 'value' */
CLJFunctionPtr CLJRFFunction::setProperty(const QString &name, const Property &value) const
{
    if (name == "dielectric")
    {
        CLJFunctionPtr ret(*this);
        ret.edit().asA<CLJRFFunction>().setDielectric( value.asA<NumberProperty>().value() );
        return ret;
    }
    else
        return CLJCutoffFunction::setProperty(name, value);
}

/** Return the value of the property with name 'name' */
PropertyPtr CLJRFFunction::property(const QString &name) const
{
    if (name == "dielectric")
    {
        return NumberProperty(dielectric());
    }
    else
    {
        return CLJCutoffFunction::property(name);
    }
}

/** Return whether or not this function contains a property called 'name' */
bool CLJRFFunction::containsProperty(const QString &name) const
{
    return (name == "dielectric") or CLJCutoffFunction::containsProperty(name);
}

/** Set the dielectric constant to 'dielectric' */
void CLJRFFunction::setDielectric(float d)
{
    diel = d;
}

/** Return the value of the dielectric constant */
float CLJRFFunction::dielectric() const
{
    return diel;
}

/** Calculate the coulomb and LJ intermolecular energy of all of the atoms in 'atoms',
    returning the results in the arguments 'cnrg' and 'ljnrg' */
void CLJRFFunction::calcVacEnergyGeo(const CLJAtoms &atoms,
                                     double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                
                if (qa[i][ii] != 0)
                {
                    const MultiFloat q( qa[i][ii] );
                    
                    if (epsa[i][ii] == 0)
                    {
                        //coulomb calculation only
                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance between the fixed and mobile atoms
                            tmp = xa[j] - x;
                            r2 = tmp * tmp;
                            tmp = ya[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = za[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            r = r2.sqrt();

                            one_over_r = r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_r + (k_rf * r2) - c_rf;
                            tmp *= q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= r.compareLess(Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        //calculate both coulomb and LJ
                        const MultiFloat sig( siga[i][ii] );
                        const MultiFloat eps( epsa[i][ii] );

                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance between the fixed and mobile atoms
                            tmp = xa[j] - x;
                            r2 = tmp * tmp;
                            tmp = ya[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = za[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            r = r2.sqrt();

                            one_over_r = r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_r + k_rf * r2 - c_rf;
                            tmp *= q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= r.compareLess(Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);

                            //now the LJ energy
                            sig2_over_r2 = sig * siga[j] * one_over_r;
                            sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                            sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                            sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                            tmp = sig6_over_r6 * sig6_over_r6;
                            tmp -= sig6_over_r6;
                            tmp *= eps;
                            tmp *= epsa[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r.compareLess(Rlj);
                            iljnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else
                {
                    //LJ calculation only
                    const MultiFloat sig( siga[i][ii] );
                    const MultiFloat eps( epsa[i][ii] );

                    for (int j=i; j<n; ++j)
                    {
                        // if i == j then we double-calculate the energies, so must
                        // scale them by 0.5
                        const MultiFloat scale( i == j ? 0.5 : 1.0 );
                    
                        //calculate the distance between the fixed and mobile atoms
                        tmp = xa[j] - x;
                        r = tmp * tmp;
                        tmp = ya[j] - y;
                        r.multiplyAdd(tmp, tmp);
                        tmp = za[j] - z;
                        r.multiplyAdd(tmp, tmp);
                        r = r.sqrt();

                        one_over_r = r.reciprocal();

                        sig2_over_r2 = sig * siga[j] * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= epsa[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += scale * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intermolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', returning the result in the arguments 'cnrg' and 'ljnrg' */
void CLJRFFunction::calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                     double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int)
            {
                const MultiInt id(id0[i][ii]);
            
                if (q0[i][ii] != 0)
                {
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    if (eps0[i][ii] == 0)
                    {
                        //coulomb energy only
                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance between the fixed and mobile atoms
                            tmp = x1[j] - x;
                            r2 = tmp * tmp;
                            tmp = y1[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = z1[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            r = r2.sqrt();

                            one_over_r = r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_r + k_rf * r2 - c_rf;
                            tmp *= q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= r.compareLess(Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        const MultiFloat sig(sig0[i][ii]);
                        const MultiFloat eps(eps0[i][ii]);

                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance between the fixed and mobile atoms
                            tmp = x1[j] - x;
                            r2 = tmp * tmp;
                            tmp = y1[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = z1[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            r = r2.sqrt();

                            one_over_r = r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_r + k_rf * r2 - c_rf;
                            tmp *= q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= r.compareLess(Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            //logical and will remove all energies where id1 == 0 or id0 == id1
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);

                            icnrg += tmp.logicalAndNot(itmp);
                            
                            //Now do the LJ energy

                            //arithmetic combining rules
                            sig2_over_r2 = sig * sig1[j] * one_over_r;
                            sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                            sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                            sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                            tmp = sig6_over_r6 * sig6_over_r6;
                            tmp -= sig6_over_r6;
                            tmp *= eps;
                            tmp *= eps1[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r.compareLess(Rlj);
                            iljnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else if (eps0[i][ii] != 0)
                {
                    //LJ energy only
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r = tmp * tmp;
                        tmp = y1[j] - y;
                        r.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r.multiplyAdd(tmp, tmp);
                        r = r.sqrt();

                        one_over_r = r.reciprocal();
                
                        //geometric combining rules
                        sig2_over_r2 = sig * sig1[j] * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb and LJ intermolecular energy of all of the atoms in 'atoms',
    assuming periodic boundary conditions in a cubic box of size 'box_dimensions',
    returning the results in 'cnrg' and 'ljnrg' */
void CLJRFFunction::calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box_dimensions,
                                     double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    const MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    tmp = xa[j] - x;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                    r2 = tmp * tmp;

                    tmp = ya[j] - y;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    tmp = za[j] - z;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * qa[j];
                
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0.
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);
                    
                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    sig2_over_r2 = sig * siga[j] * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= epsa[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intermolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', assuming periodic boundary conditions in a cubic box
    of size 'box_dimensions, returning the result in the arguments 'cnrg' and 'ljnrg' */
void CLJRFFunction::calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                     const Vector &box_dimensions,
                                     double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int)
            {
                const MultiInt id(id0[i][ii]);
                const MultiFloat x(x0[i][ii]);
                const MultiFloat y(y0[i][ii]);
                const MultiFloat z(z0[i][ii]);
                const MultiFloat q(q0[i][ii]);

                const MultiFloat sig(sig0[i][ii]);
                const MultiFloat eps(eps0[i][ii]);

                for (int j=0; j<n1; ++j)
                {
                    tmp = x1[j] - x;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                    r2 = tmp * tmp;

                    tmp = y1[j] - y;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    tmp = z1[j] - z;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);
                    
                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * q1[j];
                    
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0.
                    //logical and will remove all energies where id1 == 0 or id0 == id1
                    itmp = id1[j].compareEqual(dummy_id);
                    itmp |= id1[j].compareEqual(id);

                    icnrg += tmp.logicalAndNot(itmp);
                    
                    //Now do the LJ energy

                    //arithmetic combining rules
                    sig2_over_r2 = sig * sig1[j] * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= eps1[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    iljnrg += tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb and LJ intermolecular energy of all of the atoms in 'atoms',
    returning the results in the arguments 'cnrg' and 'ljnrg' */
void CLJRFFunction::calcVacEnergyAri(const CLJAtoms &atoms,
                                     double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                
                if (qa[i][ii] != 0)
                {
                    const MultiFloat q( qa[i][ii] );
                    
                    if (epsa[i][ii] == 0)
                    {
                        //coulomb calculation only
                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance between the fixed and mobile atoms
                            tmp = xa[j] - x;
                            r2 = tmp * tmp;
                            tmp = ya[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = za[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            r = r2.sqrt();

                            one_over_r = r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_r + k_rf * r2 - c_rf;
                            tmp *= q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= r.compareLess(Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        //calculate both coulomb and LJ
                        const MultiFloat sig( siga[i][ii] * siga[i][ii] );
                        const MultiFloat eps( epsa[i][ii] );

                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance between the fixed and mobile atoms
                            tmp = xa[j] - x;
                            r2 = tmp * tmp;
                            tmp = ya[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = za[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            r = r2.sqrt();

                            one_over_r = r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_r + k_rf * r2 - c_rf;
                            tmp *= q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= r.compareLess(Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);

                            //now the LJ energy
                            tmp = sig + (siga[j]*siga[j]);
                            tmp *= half;
                        
                            sig2_over_r2 = tmp * one_over_r;
                            sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                            sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                            sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                            tmp = sig6_over_r6 * sig6_over_r6;
                            tmp -= sig6_over_r6;
                            tmp *= eps;
                            tmp *= epsa[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r.compareLess(Rlj);
                            iljnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else
                {
                    //LJ calculation only
                    const MultiFloat sig( siga[i][ii] * siga[i][ii] );
                    const MultiFloat eps( epsa[i][ii] );

                    for (int j=i; j<n; ++j)
                    {
                        // if i == j then we double-calculate the energies, so must
                        // scale them by 0.5
                        const MultiFloat scale( i == j ? 0.5 : 1.0 );
                    
                        //calculate the distance between the fixed and mobile atoms
                        tmp = xa[j] - x;
                        r = tmp * tmp;
                        tmp = ya[j] - y;
                        r.multiplyAdd(tmp, tmp);
                        tmp = za[j] - z;
                        r.multiplyAdd(tmp, tmp);
                        r = r.sqrt();

                        one_over_r = r.reciprocal();

                        tmp = sig + (siga[j]*siga[j]);
                        tmp *= half;
                    
                        sig2_over_r2 = tmp * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= epsa[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += scale * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intermolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', returning the result in the arguments 'cnrg' and 'ljnrg' */
void CLJRFFunction::calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                     double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int)
            {
                const MultiInt id(id0[i][ii]);
            
                if (q0[i][ii] != 0)
                {
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    if (eps0[i][ii] == 0)
                    {
                        //coulomb energy only
                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance between the fixed and mobile atoms
                            tmp = x1[j] - x;
                            r2 = tmp * tmp;
                            tmp = y1[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = z1[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            r = r2.sqrt();

                            one_over_r = r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_r + k_rf * r2 - c_rf;
                            tmp *= q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= r.compareLess(Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                        const MultiFloat eps(eps0[i][ii]);

                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance between the fixed and mobile atoms
                            tmp = x1[j] - x;
                            r2 = tmp * tmp;
                            tmp = y1[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = z1[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            r = r2.sqrt();

                            one_over_r = r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_r + k_rf * r2 - c_rf;
                            tmp *= q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= r.compareLess(Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            //logical and will remove all energies where id1 == 0 or id0 == id1
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);

                            icnrg += tmp.logicalAndNot(itmp);
                            
                            //Now do the LJ energy

                            //arithmetic combining rules
                            tmp = sig + (sig1[j]*sig1[j]);
                            tmp *= half;
                        
                            sig2_over_r2 = tmp * one_over_r;
                            sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                            sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                            sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                            tmp = sig6_over_r6 * sig6_over_r6;
                            tmp -= sig6_over_r6;
                            tmp *= eps;
                            tmp *= eps1[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r.compareLess(Rlj);
                            iljnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else if (eps0[i][ii] != 0)
                {
                    //LJ energy only
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r = tmp * tmp;
                        tmp = y1[j] - y;
                        r.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r.multiplyAdd(tmp, tmp);
                        r = r.sqrt();

                        one_over_r = r.reciprocal();
                
                        //arithmetic combining rules
                        tmp = sig + (sig1[j]*sig1[j]);
                        tmp *= half;
                    
                        sig2_over_r2 = tmp * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb intermolecular energy of all atoms in 'atoms' */
double CLJRFFunction::calcVacCoulombEnergyAri(const CLJAtoms &atoms) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r;
    MultiDouble icnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int and qa[i][ii] != 0)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                
                //coulomb calculation only
                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    const MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //calculate the distance between the fixed and mobile atoms
                    tmp = xa[j] - x;
                    r2 = tmp * tmp;
                    tmp = ya[j] - y;
                    r2.multiplyAdd(tmp, tmp);
                    tmp = za[j] - z;
                    r2.multiplyAdd(tmp, tmp);
                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * qa[j];
                
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0.
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);
                    
                    icnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    return icnrg.sum();
}

/** Calculate the coulomb intermolecular energy between all atoms in 'atoms0'
    and 'atoms1' */
double CLJRFFunction::calcVacCoulombEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                              float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rc(coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r;
    MultiDouble icnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int and q0[i][ii] != 0)
            {
                const MultiInt id(id0[i][ii]);
                const MultiFloat x(x0[i][ii]);
                const MultiFloat y(y0[i][ii]);
                const MultiFloat z(z0[i][ii]);
                const MultiFloat q(q0[i][ii]);

                //coulomb energy only
                for (int j=0; j<n1; ++j)
                {
                    //calculate the distance between the fixed and mobile atoms
                    tmp = x1[j] - x;
                    r2 = tmp * tmp;
                    tmp = y1[j] - y;
                    r2.multiplyAdd(tmp, tmp);
                    tmp = z1[j] - z;
                    r2.multiplyAdd(tmp, tmp);
                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * q1[j];
                
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0.
                    itmp = id1[j].compareEqual(dummy_id);
                    itmp |= id1[j].compareEqual(id);
                    
                    icnrg += tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    return icnrg.sum();
}
    
/** Calculate the LJ intermolecular energy of all atoms in 'atoms' */
double CLJRFFunction::calcVacLJEnergyAri(const CLJAtoms &atoms) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int and epsa[i][ii] != 0)
            {
                //LJ calculation only
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat sig( siga[i][ii] * siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    const MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //calculate the distance between the fixed and mobile atoms
                    tmp = xa[j] - x;
                    r = tmp * tmp;
                    tmp = ya[j] - y;
                    r.multiplyAdd(tmp, tmp);
                    tmp = za[j] - z;
                    r.multiplyAdd(tmp, tmp);
                    r = r.sqrt();

                    one_over_r = r.reciprocal();

                    tmp = sig + (siga[j]*siga[j]);
                    tmp *= half;
                
                    sig2_over_r2 = tmp * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= epsa[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    return iljnrg.sum();
}

/** Calculate the LJ intermolecular energy between all atoms in 'atoms0'
    and 'atoms1' */
double CLJRFFunction::calcVacLJEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                         float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rlj(lj_cutoff);
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble iljnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int and eps0[i][ii] != 0)
            {
                //LJ energy only
                const MultiFloat x(x0[i][ii]);
                const MultiFloat y(y0[i][ii]);
                const MultiFloat z(z0[i][ii]);
                const MultiInt id(id0[i][ii]);
                const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                const MultiFloat eps(eps0[i][ii]);

                for (int j=0; j<n1; ++j)
                {
                    //calculate the distance between the fixed and mobile atoms
                    tmp = x1[j] - x;
                    r = tmp * tmp;
                    tmp = y1[j] - y;
                    r.multiplyAdd(tmp, tmp);
                    tmp = z1[j] - z;
                    r.multiplyAdd(tmp, tmp);
                    r = r.sqrt();

                    one_over_r = r.reciprocal();
            
                    //arithmetic combining rules
                    tmp = sig + (sig1[j]*sig1[j]);
                    tmp *= half;
                
                    sig2_over_r2 = tmp * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= eps1[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    itmp = id1[j].compareEqual(dummy_id);
                    itmp |= id1[j].compareEqual(id);

                    iljnrg += tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    return iljnrg.sum();
}

/** Calculate the coulomb and LJ intermolecular energy of all of the atoms in 'atoms',
    assuming periodic boundary conditions in a cubic box of size 'box_dimensions',
    returning the results in 'cnrg' and 'ljnrg' */
void CLJRFFunction::calcBoxEnergyAri(const CLJAtoms &atoms, const Vector &box_dimensions,
                                     double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii] * siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    const MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    tmp = xa[j] - x;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                    r2 = tmp * tmp;

                    tmp = ya[j] - y;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    tmp = za[j] - z;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * qa[j];
                
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0.
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);
                    
                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    tmp = sig + (siga[j]*siga[j]);
                    tmp *= half;
                
                    sig2_over_r2 = tmp * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= epsa[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intermolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', assuming periodic boundary conditions in a cubic box
    of size 'box_dimensions, returning the result in the arguments 'cnrg' and 'ljnrg' */
void CLJRFFunction::calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                     const Vector &box_dimensions,
                                     double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int)
            {
                const MultiInt id(id0[i][ii]);
                const MultiFloat x(x0[i][ii]);
                const MultiFloat y(y0[i][ii]);
                const MultiFloat z(z0[i][ii]);
                const MultiFloat q(q0[i][ii]);

                const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                const MultiFloat eps(eps0[i][ii]);

                for (int j=0; j<n1; ++j)
                {
                    tmp = x1[j] - x;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                    r2 = tmp * tmp;

                    tmp = y1[j] - y;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    tmp = z1[j] - z;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);
                    
                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * q1[j];
                    
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0.
                    //logical and will remove all energies where id1 == 0 or id0 == id1
                    itmp = id1[j].compareEqual(dummy_id);
                    itmp |= id1[j].compareEqual(id);

                    icnrg += tmp.logicalAndNot(itmp);
                    
                    //Now do the LJ energy

                    //arithmetic combining rules
                    tmp = sig + (sig1[j]*sig1[j]);
                    tmp *= half;
                
                    sig2_over_r2 = tmp * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= eps1[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    iljnrg += tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** This function does support calculations using a grid */
bool CLJRFFunction::supportsGridCalculation() const
{
    return true;
}

/** Calculate the energy on the grid from the passed atoms using vacuum boundary conditions */
void CLJRFFunction::calcVacGrid(const CLJAtoms &atoms, const GridInfo &grid_info,
                                const int start, const int end, float *gridpot_array) const
{
    const MultiFloat* const x = atoms.x().constData();
    const MultiFloat* const y = atoms.y().constData();
    const MultiFloat* const z = atoms.z().constData();
    const MultiFloat* const q = atoms.q().constData();
    const MultiInt* const id = atoms.ID().constData();
    
    const MultiFloat Rc( coul_cutoff );

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiInt dummy_id = CLJAtoms::idOfDummy();

    MultiFloat tmp, r, r2, one_over_r, itmp;

    const int nats = atoms.x().count();

    for (int i = start; i < end; ++i)
    {
        const Vector grid_point = grid_info.point(i);
        
        const MultiFloat px(grid_point.x());
        const MultiFloat py(grid_point.y());
        const MultiFloat pz(grid_point.z());
        
        MultiDouble pot(0);
        
        for (int j=0; j<nats; ++j)
        {
            //calculate the distance between the atom and grid point
            tmp = px - x[j];
            r2 = tmp * tmp;
            tmp = py - y[j];
            r2.multiplyAdd(tmp, tmp);
            tmp = pz - z[j];
            r2.multiplyAdd(tmp, tmp);

            r = r2.sqrt();

            one_over_r = r.reciprocal();
    
            // calculate the coulomb energy using
            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
            // c = (1/r_c) * (3 eps)/(2 eps + 1)
            tmp = one_over_r + k_rf * r2 - c_rf;
            
            //exclude dummy atoms when building the grid
            tmp *= q[j].logicalAndNot( id[j].compareEqual(dummy_id) );
        
            //apply the cutoff - compare r against Rc. This will
            //return 1 if r is less than Rc, or 0 otherwise. Logical
            //and will then remove all energies where r >= Rc
            pot += tmp.logicalAnd( r.compareLess(Rc) );
        }
        
        *gridpot_array = pot.sum();
        ++gridpot_array;
    }
}

/** Calculate the energy on the grid from the passed atoms using vacuum boundary conditions */
void CLJRFFunction::calcBoxGrid(const CLJAtoms &atoms, const GridInfo &grid_info,
                                const Vector &box_dimensions,
                                const int start, const int end, float *gridpot_array) const
{
    const MultiFloat* const x = atoms.x().constData();
    const MultiFloat* const y = atoms.y().constData();
    const MultiFloat* const z = atoms.z().constData();
    const MultiFloat* const q = atoms.q().constData();
    const MultiInt* const id = atoms.ID().constData();
    
    const MultiFloat Rc( coul_cutoff );

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );

    const MultiInt dummy_id = CLJAtoms::idOfDummy();

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    MultiFloat tmp, r, r2, one_over_r, itmp;

    const int nats = atoms.x().count();

    for (int i = start; i < end; ++i)
    {
        const Vector grid_point = grid_info.point(i);
        
        const MultiFloat px(grid_point.x());
        const MultiFloat py(grid_point.y());
        const MultiFloat pz(grid_point.z());
        
        MultiDouble pot(0);
        
        for (int j=0; j<nats; ++j)
        {
            //calculate the distance between the atom and grid point
            tmp = px - x[j];
            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
            tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
            r2 = tmp * tmp;

            tmp = py - y[j];
            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
            tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
            r2.multiplyAdd(tmp, tmp);

            tmp = pz - z[j];
            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
            tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
            r2.multiplyAdd(tmp, tmp);

            r = r2.sqrt();

            one_over_r = r.reciprocal();
    
            // calculate the coulomb energy using
            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
            // c = (1/r_c) * (3 eps)/(2 eps + 1)
            tmp = one_over_r + k_rf * r2 - c_rf;
            
            //exclude dummy atoms when building the grid
            tmp *= q[j].logicalAndNot( id[j].compareEqual(dummy_id) );
        
            //apply the cutoff - compare r against Rc. This will
            //return 1 if r is less than Rc, or 0 otherwise. Logical
            //and will then remove all energies where r >= Rc
            pot += tmp.logicalAnd( r.compareLess(Rc) );
        }
        
        *gridpot_array = pot.sum();
        ++gridpot_array;
    }
}


/////////
///////// Implementation of CLJSoftRFFunction
/////////

static const RegisterMetaType<CLJSoftRFFunction> r_softshift;

QDataStream &operator<<(QDataStream &ds, const CLJSoftRFFunction &func)
{
    writeHeader(ds, r_softshift, 1);
    
    ds << func.diel << static_cast<const CLJSoftFunction&>(func);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJSoftRFFunction &func)
{
    VersionID v = readHeader(ds, r_softshift);
    
    if (v == 1)
    {
        ds >> func.diel >> static_cast<CLJSoftFunction&>(func);
    }
    else
        throw version_error(v, "1", r_softshift, CODELOC);
    
    return ds;
}

CLJSoftRFFunction::CLJSoftRFFunction()
                  : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(),
                    diel(default_dielectric)
{}

CLJFunctionPtr CLJSoftRFFunction::defaultRFFunction()
{
    static CLJFunctionPtr ptr( new CLJSoftRFFunction() );
    return ptr;
}

CLJSoftRFFunction::CLJSoftRFFunction(Length cutoff)
                     : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(cutoff),
                       diel(default_dielectric)
{}

CLJSoftRFFunction::CLJSoftRFFunction(Length coul_cutoff, Length lj_cutoff)
    : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(coul_cutoff, lj_cutoff),
                       diel(default_dielectric)
{}

CLJSoftRFFunction::CLJSoftRFFunction(const Space &space, Length cutoff)
    : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(space, cutoff),
                       diel(default_dielectric)
{}

CLJSoftRFFunction::CLJSoftRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff)
    : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(space, coul_cutoff, lj_cutoff),
                       diel(default_dielectric)
{}

CLJSoftRFFunction::CLJSoftRFFunction(Length cutoff, COMBINING_RULES combining_rules)
    : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(cutoff, combining_rules),
                       diel(default_dielectric)
{}

CLJSoftRFFunction::CLJSoftRFFunction(Length coul_cutoff, Length lj_cutoff,
                                           COMBINING_RULES combining_rules)
    : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(
                                   coul_cutoff, lj_cutoff, combining_rules),
                       diel(default_dielectric)
{}

CLJSoftRFFunction::CLJSoftRFFunction(const Space &space, COMBINING_RULES combining_rules)
    : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(space, combining_rules),
                       diel(default_dielectric)
{}

CLJSoftRFFunction::CLJSoftRFFunction(const Space &space, Length cutoff,
                                           COMBINING_RULES combining_rules)
    : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(
                                   space, cutoff, combining_rules),
                       diel(default_dielectric)
{}

CLJSoftRFFunction::CLJSoftRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff,
                                           COMBINING_RULES combining_rules)
    : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(
                                   space, coul_cutoff, lj_cutoff, combining_rules),
                       diel(default_dielectric)
{}

/** Copy constructor */
CLJSoftRFFunction::CLJSoftRFFunction(const CLJSoftRFFunction &other)
    : ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>(other),
                       diel(other.diel)
{}

/** Destructor */
CLJSoftRFFunction::~CLJSoftRFFunction()
{}

/** Copy assignment operator */
CLJSoftRFFunction& CLJSoftRFFunction::operator=(const CLJSoftRFFunction &other)
{
    diel = other.diel;
    CLJSoftFunction::operator=(other);
    return *this;
}

/** Comparison operator */
bool CLJSoftRFFunction::operator==(const CLJSoftRFFunction &other) const
{
    return diel == other.diel and CLJSoftFunction::operator==(other);
}

/** Comparison operator */
bool CLJSoftRFFunction::operator!=(const CLJSoftRFFunction &other) const
{
    return not operator==(other);
}

const char* CLJSoftRFFunction::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJSoftRFFunction>() );
}

const char* CLJSoftRFFunction::what() const
{
    return CLJSoftRFFunction::typeName();
}

CLJSoftRFFunction* CLJSoftRFFunction::clone() const
{
    return new CLJSoftRFFunction(*this);
}

/** Return the properties of this function */
Properties CLJSoftRFFunction::properties() const
{
    Properties props = CLJSoftFunction::properties();
    props.setProperty("dielectric", NumberProperty(dielectric()));
    return props;
}

/** Return a copy of this function where the property 'name' has been set to the
    value 'value' */
CLJFunctionPtr CLJSoftRFFunction::setProperty(const QString &name, const Property &value) const
{
    if (name == "dielectric")
    {
        CLJFunctionPtr ret(*this);
        ret.edit().asA<CLJSoftRFFunction>().setDielectric( value.asA<NumberProperty>().value() );
        return ret;
    }
    else
        return CLJSoftFunction::setProperty(name, value);
}

/** Return the value of the property with name 'name' */
PropertyPtr CLJSoftRFFunction::property(const QString &name) const
{
    if (name == "dielectric")
    {
        return NumberProperty(dielectric());
    }
    else
    {
        return CLJSoftFunction::property(name);
    }
}

/** Return whether or not this function contains a property called 'name' */
bool CLJSoftRFFunction::containsProperty(const QString &name) const
{
    return (name == "dielectric") or CLJSoftFunction::containsProperty(name);
}

/** Set the dielectric constant to 'dielectric' */
void CLJSoftRFFunction::setDielectric(float d)
{
    diel = d;
}

/** Return the value of the dielectric constant */
float CLJSoftRFFunction::dielectric() const
{
    return diel;
}

QString CLJSoftRFFunction::toString() const
{
    if (this->hasCutoff())
        return QObject::tr("CLJSoftRFFunction( dielectric() == %4, coulombCutoff() == %1 A, "
                           "ljCutoff() == %2 A, space() == %3 )")
            .arg(coulombCutoff().to(angstrom))
            .arg(ljCutoff().to(angstrom))
            .arg(space().toString())
            .arg(dielectric());
    else
        return QObject::tr("CLJSoftRFFunction( dielectric() == %2, no cutoff, space() == %1 )")
                    .arg(space().toString()).arg(dielectric());
}

/** Calculate the coulomb and LJ intermolecular energy of all of the atoms in 'atoms',
    returning the results in the arguments 'cnrg' and 'ljnrg' */
void CLJSoftRFFunction::calcVacEnergyGeo(const CLJAtoms &atoms,
                                         double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );

    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                
                if (qa[i][ii] != 0)
                {
                    const MultiFloat q( qa[i][ii] );
                    
                    if (epsa[i][ii] == 0)
                    {
                        //coulomb calculation only
                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance^2 between the fixed and mobile atoms
                            tmp = xa[j] - x;
                            r2 = tmp * tmp;
                            tmp = ya[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = za[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            
                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();

                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        //calculate both coulomb and LJ
                        const MultiFloat sig( siga[i][ii] );
                        const MultiFloat eps( epsa[i][ii] );

                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance between the fixed and mobile atoms
                            tmp = xa[j] - x;
                            r2 = tmp * tmp;
                            tmp = ya[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = za[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            
                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);

                            //now the LJ energy
                            sigma = sig * siga[j];
                            delta_sigma_r2 = delta * sigma + r2;
                            
                            sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                            sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                            tmp = sig6_over_delta3 * sig6_over_delta3;
                            tmp -= sig6_over_delta3;
                            tmp *= eps;
                            tmp *= epsa[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r2.compareLess(Rlj2);
                            iljnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else
                {
                    //LJ calculation only
                    const MultiFloat sig( siga[i][ii] );
                    const MultiFloat eps( epsa[i][ii] );

                    for (int j=i; j<n; ++j)
                    {
                        // if i == j then we double-calculate the energies, so must
                        // scale them by 0.5
                        const MultiFloat scale( i == j ? 0.5 : 1.0 );
                    
                        //calculate the distance between the fixed and mobile atoms
                        tmp = xa[j] - x;
                        r2 = tmp * tmp;
                        tmp = ya[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = za[j] - z;
                        r2.multiplyAdd(tmp, tmp);

                        sigma = sig * siga[j];
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= epsa[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += scale * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intermolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', returning the result in the arguments 'cnrg' and 'ljnrg' */
void CLJSoftRFFunction::calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                         double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);
    
    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );

    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );

    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int)
            {
                const MultiInt id(id0[i][ii]);
            
                if (q0[i][ii] != 0)
                {
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    if (eps0[i][ii] == 0)
                    {
                        //coulomb energy only
                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance between atoms
                            tmp = x1[j] - x;
                            r2 = tmp * tmp;
                            tmp = y1[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = z1[j] - z;
                            r2.multiplyAdd(tmp, tmp);

                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        const MultiFloat sig(sig0[i][ii]);
                        const MultiFloat eps(eps0[i][ii]);

                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance between atoms
                            tmp = x1[j] - x;
                            r2 = tmp * tmp;
                            tmp = y1[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = z1[j] - z;
                            r2.multiplyAdd(tmp, tmp);

                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                            
                            //now the LJ energy
                            sigma = sig * sig1[j];
                            delta_sigma_r2 = delta * sigma + r2;
                            
                            sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                            sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                            tmp = sig6_over_delta3 * sig6_over_delta3;
                            tmp -= sig6_over_delta3;
                            tmp *= eps;
                            tmp *= eps1[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r2.compareLess(Rlj2);
                            iljnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else if (eps0[i][ii] != 0)
                {
                    //LJ energy only
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance between atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);

                        sigma = sig * sig1[j];
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb and LJ intermolecular energy of all of the atoms in 'atoms',
    assuming periodic boundary conditions in a cubic box of size 'box_dimensions',
    returning the results in 'cnrg' and 'ljnrg' */
void CLJSoftRFFunction::calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box_dimensions,
                                         double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                
                if (qa[i][ii] != 0)
                {
                    const MultiFloat q( qa[i][ii] );
                    
                    if (epsa[i][ii] == 0)
                    {
                        //coulomb calculation only
                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance^2
                            tmp = xa[j] - x;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                            r2 = tmp * tmp;

                            tmp = ya[j] - y;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            tmp = za[j] - z;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);
                            
                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        //calculate both coulomb and LJ
                        const MultiFloat sig( siga[i][ii] );
                        const MultiFloat eps( epsa[i][ii] );

                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance^2
                            tmp = xa[j] - x;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                            r2 = tmp * tmp;

                            tmp = ya[j] - y;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            tmp = za[j] - z;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);
                            
                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);

                            //now the LJ energy
                            sigma = sig * siga[j];
                            delta_sigma_r2 = delta * sigma + r2;
                            
                            sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                            sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                            tmp = sig6_over_delta3 * sig6_over_delta3;
                            tmp -= sig6_over_delta3;
                            tmp *= eps;
                            tmp *= epsa[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r2.compareLess(Rlj2);
                            iljnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else
                {
                    //LJ calculation only
                    const MultiFloat sig( siga[i][ii] );
                    const MultiFloat eps( epsa[i][ii] );

                    for (int j=i; j<n; ++j)
                    {
                        // if i == j then we double-calculate the energies, so must
                        // scale them by 0.5
                        const MultiFloat scale( i == j ? 0.5 : 1.0 );
                    
                        //calculate the distance^2
                        tmp = xa[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = ya[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = za[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        sigma = sig * siga[j];
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= epsa[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += scale * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intermolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', assuming periodic boundary conditions in a cubic box
    of size 'box_dimensions, returning the result in the arguments 'cnrg' and 'ljnrg' */
void CLJSoftRFFunction::calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                         const Vector &box_dimensions,
                                         double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);
    
    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );

    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int)
            {
                const MultiInt id(id0[i][ii]);
            
                if (q0[i][ii] != 0)
                {
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    if (eps0[i][ii] == 0)
                    {
                        //coulomb energy only
                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance^2
                            tmp = x1[j] - x;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                            r2 = tmp * tmp;

                            tmp = y1[j] - y;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            tmp = z1[j] - z;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        const MultiFloat sig(sig0[i][ii]);
                        const MultiFloat eps(eps0[i][ii]);

                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance^2
                            tmp = x1[j] - x;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                            r2 = tmp * tmp;

                            tmp = y1[j] - y;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            tmp = z1[j] - z;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                            
                            //now the LJ energy
                            sigma = sig * sig1[j];
                            delta_sigma_r2 = delta * sigma + r2;
                            
                            sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                            sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                            tmp = sig6_over_delta3 * sig6_over_delta3;
                            tmp -= sig6_over_delta3;
                            tmp *= eps;
                            tmp *= eps1[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r2.compareLess(Rlj2);
                            iljnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else if (eps0[i][ii] != 0)
                {
                    //LJ energy only
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance^2
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        sigma = sig * sig1[j];
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb and LJ intermolecular energy of all of the atoms in 'atoms',
    returning the results in the arguments 'cnrg' and 'ljnrg' */
void CLJSoftRFFunction::calcVacEnergyAri(const CLJAtoms &atoms,
                                         double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );

    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                
                if (qa[i][ii] != 0)
                {
                    const MultiFloat q( qa[i][ii] );
                    
                    if (epsa[i][ii] == 0)
                    {
                        //coulomb calculation only
                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance^2 between the fixed and mobile atoms
                            tmp = xa[j] - x;
                            r2 = tmp * tmp;
                            tmp = ya[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = za[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            
                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        //calculate both coulomb and LJ
                        const MultiFloat sig( siga[i][ii] * siga[i][ii] );
                        const MultiFloat eps( epsa[i][ii] );

                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance between the fixed and mobile atoms
                            tmp = xa[j] - x;
                            r2 = tmp * tmp;
                            tmp = ya[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = za[j] - z;
                            r2.multiplyAdd(tmp, tmp);
                            
                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);

                            //now the LJ energy

                            //arithmetic combining rules
                            sigma = sig + (siga[j]*siga[j]);
                            sigma *= half;

                            delta_sigma_r2 = delta * sigma + r2;
                            
                            sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                            sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                            tmp = sig6_over_delta3 * sig6_over_delta3;
                            tmp -= sig6_over_delta3;
                            tmp *= eps;
                            tmp *= epsa[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r2.compareLess(Rlj2);
                            iljnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else
                {
                    //LJ calculation only
                    const MultiFloat sig( siga[i][ii] );
                    const MultiFloat eps( epsa[i][ii] );

                    for (int j=i; j<n; ++j)
                    {
                        // if i == j then we double-calculate the energies, so must
                        // scale them by 0.5
                        const MultiFloat scale( i == j ? 0.5 : 1.0 );
                    
                        //calculate the distance between the fixed and mobile atoms
                        tmp = xa[j] - x;
                        r2 = tmp * tmp;
                        tmp = ya[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = za[j] - z;
                        r2.multiplyAdd(tmp, tmp);

                        //arithmetic combining rules
                        sigma = sig + (siga[j]*siga[j]);
                        sigma *= half;

                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= epsa[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += scale * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intermolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', returning the result in the arguments 'cnrg' and 'ljnrg' */
void CLJSoftRFFunction::calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                         double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);
    
    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );

    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int)
            {
                const MultiInt id(id0[i][ii]);
            
                if (q0[i][ii] != 0)
                {
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    if (eps0[i][ii] == 0)
                    {
                        //coulomb energy only
                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance between atoms
                            tmp = x1[j] - x;
                            r2 = tmp * tmp;
                            tmp = y1[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = z1[j] - z;
                            r2.multiplyAdd(tmp, tmp);

                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                        const MultiFloat eps(eps0[i][ii]);

                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance between atoms
                            tmp = x1[j] - x;
                            r2 = tmp * tmp;
                            tmp = y1[j] - y;
                            r2.multiplyAdd(tmp, tmp);
                            tmp = z1[j] - z;
                            r2.multiplyAdd(tmp, tmp);

                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                            
                            //now the LJ energy

                            //arithmetic combining rules
                            sigma = sig + (sig1[j]*sig1[j]);
                            sigma *= half;

                            delta_sigma_r2 = delta * sigma + r2;
                            
                            sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                            sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                            tmp = sig6_over_delta3 * sig6_over_delta3;
                            tmp -= sig6_over_delta3;
                            tmp *= eps;
                            tmp *= eps1[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r2.compareLess(Rlj2);
                            iljnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else if (eps0[i][ii] != 0)
                {
                    //LJ energy only
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance between atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);

                        //arithmetic combining rules
                        sigma = sig + (sig1[j]*sig1[j]);
                        sigma *= half;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intermolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', assuming periodic boundary conditions in a cubic box
    of size 'box_dimensions, returning the result in the arguments 'cnrg' and 'ljnrg' */
void CLJSoftRFFunction::calcBoxEnergyAri(const CLJAtoms &atoms,
                                         const Vector &box_dimensions,
                                         double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );

    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                
                if (qa[i][ii] != 0)
                {
                    const MultiFloat q( qa[i][ii] );
                    
                    if (epsa[i][ii] == 0)
                    {
                        //coulomb calculation only
                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance^2
                            tmp = xa[j] - x;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                            r2 = tmp * tmp;

                            tmp = ya[j] - y;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            tmp = za[j] - z;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);
                            
                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        //calculate both coulomb and LJ
                        const MultiFloat sig( siga[i][ii] * siga[i][ii] );
                        const MultiFloat eps( epsa[i][ii] );

                        for (int j=i; j<n; ++j)
                        {
                            // if i == j then we double-calculate the energies, so must
                            // scale them by 0.5
                            const MultiFloat scale( i == j ? 0.5 : 1.0 );
                        
                            //calculate the distance^2
                            tmp = xa[j] - x;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                            r2 = tmp * tmp;

                            tmp = ya[j] - y;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            tmp = za[j] - z;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);
                            
                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * qa[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = ida[j].compareEqual(dummy_id);
                            itmp |= ida[j].compareEqual(id);
                            
                            icnrg += scale * tmp.logicalAndNot(itmp);

                            //now the LJ energy

                            //arithmetic combining rules
                            sigma = sig + (siga[j]*siga[j]);
                            sigma *= half;

                            delta_sigma_r2 = delta * sigma + r2;
                            
                            sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                            sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                            tmp = sig6_over_delta3 * sig6_over_delta3;
                            tmp -= sig6_over_delta3;
                            tmp *= eps;
                            tmp *= epsa[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r2.compareLess(Rlj2);
                            iljnrg += scale * tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else
                {
                    //LJ calculation only
                    const MultiFloat sig( siga[i][ii] );
                    const MultiFloat eps( epsa[i][ii] );

                    for (int j=i; j<n; ++j)
                    {
                        // if i == j then we double-calculate the energies, so must
                        // scale them by 0.5
                        const MultiFloat scale( i == j ? 0.5 : 1.0 );
                    
                        //calculate the distance^2
                        tmp = xa[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = ya[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = za[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        //arithmetic combining rules
                        sigma = sig + (siga[j]*siga[j]);
                        sigma *= half;

                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= epsa[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += scale * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intermolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', assuming periodic boundary conditions in a cubic box
    of size 'box_dimensions, returning the result in the arguments 'cnrg' and 'ljnrg' */
void CLJSoftRFFunction::calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                         const Vector &box_dimensions,
                                         double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);
    
    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );

    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();

    for (int i=0; i<n0; ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (id0[i][ii] != dummy_int)
            {
                const MultiInt id(id0[i][ii]);
            
                if (q0[i][ii] != 0)
                {
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    if (eps0[i][ii] == 0)
                    {
                        //coulomb energy only
                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance^2
                            tmp = x1[j] - x;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                            r2 = tmp * tmp;

                            tmp = y1[j] - y;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            tmp = z1[j] - z;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                    else
                    {
                        const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                        const MultiFloat eps(eps0[i][ii]);

                        for (int j=0; j<n1; ++j)
                        {
                            //calculate the distance^2
                            tmp = x1[j] - x;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                            r2 = tmp * tmp;

                            tmp = y1[j] - y;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            tmp = z1[j] - z;
                            tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                            tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                            r2.multiplyAdd(tmp, tmp);

                            soft_r2 = r2 + alfa;
                            soft_r = soft_r2.sqrt();

                            one_over_soft_r = soft_r.reciprocal();
                    
                            // calculate the coulomb energy using
                            // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                            // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                            // c = (1/r_c) * (3 eps)/(2 eps + 1)
                            tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                            tmp *= one_minus_alpha_to_n * q * q1[j];
                        
                            //apply the cutoff - compare r against Rc. This will
                            //return 1 if r is less than Rc, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rc
                            tmp &= soft_r.compareLess(soft_Rc);

                            //make sure that the ID of atoms1 is not zero, and is
                            //also not the same as the atoms0.
                            itmp = id1[j].compareEqual(dummy_id);
                            itmp |= id1[j].compareEqual(id);
                            
                            icnrg += tmp.logicalAndNot(itmp);
                            
                            //now the LJ energy

                            //arithmetic combining rules
                            sigma = sig + (sig1[j]*sig1[j]);
                            sigma *= half;

                            delta_sigma_r2 = delta * sigma + r2;
                            
                            sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                            sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                            tmp = sig6_over_delta3 * sig6_over_delta3;
                            tmp -= sig6_over_delta3;
                            tmp *= eps;
                            tmp *= eps1[j];
                        
                            //apply the cutoff - compare r against Rlj. This will
                            //return 1 if r is less than Rlj, or 0 otherwise. Logical
                            //and will then remove all energies where r >= Rlj
                            tmp &= r2.compareLess(Rlj2);
                            iljnrg += tmp.logicalAndNot(itmp);
                        }
                    }
                }
                else if (eps0[i][ii] != 0)
                {
                    //LJ energy only
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance^2
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        //arithmetic combining rules
                        sigma = sig + (sig1[j]*sig1[j]);
                        sigma *= half;

                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/////////
///////// Implementation of CLJIntraRFFunction
/////////

static const RegisterMetaType<CLJIntraRFFunction> r_intra;

QDataStream &operator<<(QDataStream &ds, const CLJIntraRFFunction &intra)
{
    writeHeader(ds, r_intra, 1);
    
    ds << intra.diel << static_cast<const CLJIntraFunction&>(intra);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJIntraRFFunction &intra)
{
    VersionID v = readHeader(ds, r_intra);
    
    if (v == 1)
    {
        ds >> intra.diel >> static_cast<CLJIntraFunction&>(intra);
    }
    else
        throw version_error(v, "1", r_intra, CODELOC);
    
    return ds;
}

CLJIntraRFFunction::CLJIntraRFFunction()
                   : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(),
                     diel(default_dielectric)
{}

CLJFunctionPtr CLJIntraRFFunction::defaultRFFunction()
{
    static CLJFunctionPtr ptr( new CLJIntraRFFunction() );
    return ptr;
}

CLJIntraRFFunction::CLJIntraRFFunction(Length cutoff)
                      : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(cutoff),
                        diel(default_dielectric)
{}

CLJIntraRFFunction::CLJIntraRFFunction(Length coul_cutoff, Length lj_cutoff)
                      : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(
                            coul_cutoff, lj_cutoff),
                       diel(default_dielectric)
{}

CLJIntraRFFunction::CLJIntraRFFunction(const Space &space, Length cutoff)
                      : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(
                            space, cutoff),
                       diel(default_dielectric)
{}

CLJIntraRFFunction::CLJIntraRFFunction(const Space &space,
                                             Length coul_cutoff, Length lj_cutoff)
                      : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(
                            space, coul_cutoff, lj_cutoff),
                       diel(default_dielectric)
{}

CLJIntraRFFunction::CLJIntraRFFunction(Length cutoff, COMBINING_RULES combining_rules)
                      : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(
                            cutoff, combining_rules),
                       diel(default_dielectric)
{}

CLJIntraRFFunction::CLJIntraRFFunction(Length coul_cutoff, Length lj_cutoff,
                                             COMBINING_RULES combining_rules)
                 : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(
                            coul_cutoff, lj_cutoff, combining_rules),
                       diel(default_dielectric)
{}

CLJIntraRFFunction::CLJIntraRFFunction(const Space &space, COMBINING_RULES combining_rules)
                 : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(
                            space, combining_rules),
                       diel(default_dielectric)
{}

CLJIntraRFFunction::CLJIntraRFFunction(const Space &space, Length cutoff,
                                             COMBINING_RULES combining_rules)
                 : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(
                            space, cutoff, combining_rules),
                       diel(default_dielectric)
{}

CLJIntraRFFunction::CLJIntraRFFunction(const Space &space, Length coul_cutoff,
                                             Length lj_cutoff,
                                             COMBINING_RULES combining_rules)
                 : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(
                            space, coul_cutoff, lj_cutoff, combining_rules),
                       diel(default_dielectric)
{}

/** Copy constructor */
CLJIntraRFFunction::CLJIntraRFFunction(const CLJIntraRFFunction &other)
                 : ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>(other),
                       diel(other.diel)
{}

/** Destructor */
CLJIntraRFFunction::~CLJIntraRFFunction()
{}

/** Copy assignment operator */
CLJIntraRFFunction& CLJIntraRFFunction::operator=(const CLJIntraRFFunction &other)
{
    diel = other.diel;
    CLJIntraFunction::operator=(other);
    return *this;
}

/** Comparison operator */
bool CLJIntraRFFunction::operator==(const CLJIntraRFFunction &other) const
{
    return diel == other.diel and CLJIntraFunction::operator==(other);
}

/** Comparison operator */
bool CLJIntraRFFunction::operator!=(const CLJIntraRFFunction &other) const
{
    return not operator==(other);
}

const char* CLJIntraRFFunction::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJIntraRFFunction>() );
}

const char* CLJIntraRFFunction::what() const
{
    return CLJIntraRFFunction::typeName();
}

CLJIntraRFFunction* CLJIntraRFFunction::clone() const
{
    return new CLJIntraRFFunction(*this);
}

/** Return the properties of this function */
Properties CLJIntraRFFunction::properties() const
{
    Properties props = CLJIntraFunction::properties();
    props.setProperty("dielectric", NumberProperty(dielectric()));
    return props;
}

/** Return a copy of this function where the property 'name' has been set to the
    value 'value' */
CLJFunctionPtr CLJIntraRFFunction::setProperty(const QString &name, const Property &value) const
{
    if (name == "dielectric")
    {
        CLJFunctionPtr ret(*this);
        ret.edit().asA<CLJIntraRFFunction>().setDielectric( value.asA<NumberProperty>().value() );
        return ret;
    }
    else
        return CLJIntraFunction::setProperty(name, value);
}

/** Return the value of the property with name 'name' */
PropertyPtr CLJIntraRFFunction::property(const QString &name) const
{
    if (name == "dielectric")
    {
        return NumberProperty(dielectric());
    }
    else
    {
        return CLJIntraFunction::property(name);
    }
}

/** Return whether or not this function contains a property called 'name' */
bool CLJIntraRFFunction::containsProperty(const QString &name) const
{
    return (name == "dielectric") or CLJIntraFunction::containsProperty(name);
}

/** Set the dielectric constant to 'dielectric' */
void CLJIntraRFFunction::setDielectric(float d)
{
    diel = d;
}

/** Return the value of the dielectric constant */
float CLJIntraRFFunction::dielectric() const
{
    return diel;
}

/** Calculate the coulomb and LJ intramolecular energy of all of the atoms in 'atoms',
    returning the results in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJIntraRFFunction::calcVacEnergyGeo(const CLJAtoms &atoms,
                                          double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                const bool *row = bondMatrix().constData()[id[0]].constData();

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //get the bond mask to screen out bonded interactions
                    for (int k=0; k<MultiFloat::count(); ++k)
                    {
                        bond_mask.quickSet(k, row[ida[j][k]] ? 0.0 : 1.0);
                    }
                
                    scale *= bond_mask;
                
                    //calculate the distance between the fixed and mobile atoms
                    tmp = xa[j] - x;
                    r2 = tmp * tmp;
                    tmp = ya[j] - y;
                    r2.multiplyAdd(tmp, tmp);
                    tmp = za[j] - z;
                    r2.multiplyAdd(tmp, tmp);
                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * qa[j];
                
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0 and that we are not
                    //including the energy of the atom with itself
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);

                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    sig2_over_r2 = sig * siga[j] * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= epsa[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intramolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', returning the result in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJIntraRFFunction::calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                          double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();
    
    bool not_bonded = true;
    
    if (min_distance < 5.5)
    {
        not_bonded = isNotBonded(atoms0.ID(), atoms1.ID());
    }
    
    if ( Q_LIKELY(not_bonded) )
    {
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);
                        r = r2.sqrt();

                        one_over_r = r.reciprocal();
                
                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_r + k_rf * r2 - c_rf;
                        tmp *= q * q1[j];
                    
                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= r.compareLess(Rc);

                        //make sure that the ID of atoms1 is not zero
                        itmp = id1[j].compareEqual(dummy_id);

                        icnrg += tmp.logicalAndNot(itmp);
                        
                        //Now do the LJ energy

                        sig2_over_r2 = sig * sig1[j] * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    else
    {
        MultiFloat bonded_mask(0);
    
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    const bool *row = bondMatrix().constData()[id[0]].constData();

                    for (int j=0; j<n1; ++j)
                    {
                        //create a mask to cancel out calculations of bonded pairs
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            bonded_mask.quickSet(k, row[id1[j][k]] ? 0.0 : 1.0);
                        }

                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);
                        r = r2.sqrt();

                        one_over_r = r.reciprocal();
                
                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_r + k_rf * r2 - c_rf;
                        tmp *= q * q1[j];
                    
                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= r.compareLess(Rc);

                        //make sure that the ID of atoms1 is not zero
                        itmp = id1[j].compareEqual(dummy_id);

                        icnrg += bonded_mask * tmp.logicalAndNot(itmp);
                        
                        //Now do the LJ energy
                        sig2_over_r2 = sig * sig1[j] * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += bonded_mask * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb and LJ intramolecular energy of all of the atoms in 'atoms',
    assuming periodic boundary conditions in a cubic box of size 'box_dimensions',
    returning the results in 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJIntraRFFunction::calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box_dimensions,
                                          double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                const bool *row = bondMatrix().constData()[id[0]].constData();

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //get the bond mask to screen out bonded interactions
                    for (int k=0; k<MultiFloat::count(); ++k)
                    {
                        bond_mask.quickSet(k, row[ida[j][k]] ? 0.0 : 1.0);
                    }
                
                    scale *= bond_mask;
                
                    //calculate the distance between the fixed and mobile atoms
                    tmp = xa[j] - x;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                    r2 = tmp * tmp;

                    tmp = ya[j] - y;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    tmp = za[j] - z;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * qa[j];
                
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0 and that we are not
                    //including the energy of the atom with itself
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);

                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    sig2_over_r2 = sig * siga[j] * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= epsa[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intramolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', assuming periodic boundary conditions in a cubic box
    of size 'box_dimensions, returning the result in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJIntraRFFunction::calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                          const Vector &box_dimensions,
                                          double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();
    
    bool not_bonded = true;
    
    if (min_distance < 5.5)
    {
        not_bonded = isNotBonded(atoms0.ID(), atoms1.ID());
    }
    
    if ( Q_LIKELY(not_bonded) )
    {
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);
                        
                        r = r2.sqrt();

                        one_over_r = r.reciprocal();
                
                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_r + k_rf * r2 - c_rf;
                        tmp *= q * q1[j];
                    
                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= r.compareLess(Rc);

                        //make sure that the ID of atoms1 is not zero
                        itmp = id1[j].compareEqual(dummy_id);

                        icnrg += tmp.logicalAndNot(itmp);
                        
                        //Now do the LJ energy
                        sig2_over_r2 = sig * sig1[j] * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    else
    {
        MultiFloat bonded_mask(0);
    
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    const bool *row = bondMatrix().constData()[id[0]].constData();

                    for (int j=0; j<n1; ++j)
                    {
                        //create a mask to cancel out calculations of bonded pairs
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            bonded_mask.quickSet(k, row[id1[j][k]] ? 0.0 : 1.0);
                        }

                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);
                        
                        r = r2.sqrt();

                        one_over_r = r.reciprocal();
                
                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_r + k_rf * r2 - c_rf;
                        tmp *= q * q1[j];
                    
                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= r.compareLess(Rc);

                        //make sure that the ID of atoms1 is not zero
                        itmp = id1[j].compareEqual(dummy_id);

                        icnrg += bonded_mask * tmp.logicalAndNot(itmp);
                        
                        //Now do the LJ energy
                        sig2_over_r2 = sig * sig1[j] * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += bonded_mask * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb and LJ intramolecular energy of all of the atoms in 'atoms',
    returning the results in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJIntraRFFunction::calcVacEnergyAri(const CLJAtoms &atoms,
                                          double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii] * siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                const bool *row = bondMatrix().constData()[id[0]].constData();

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //get the bond mask to screen out bonded interactions
                    for (int k=0; k<MultiFloat::count(); ++k)
                    {
                        bond_mask.quickSet(k, row[ida[j][k]] ? 0.0 : 1.0);
                    }
                
                    scale *= bond_mask;
                
                    //calculate the distance between the fixed and mobile atoms
                    tmp = xa[j] - x;
                    r2 = tmp * tmp;
                    tmp = ya[j] - y;
                    r2.multiplyAdd(tmp, tmp);
                    tmp = za[j] - z;
                    r2.multiplyAdd(tmp, tmp);
                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * qa[j];
                
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0 and that we are not
                    //including the energy of the atom with itself
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);

                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    tmp = sig + (siga[j]*siga[j]);
                    tmp *= half;
                
                    sig2_over_r2 = tmp * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= epsa[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intramolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', returning the result in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJIntraRFFunction::calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                          double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();
    
    bool not_bonded = true;
    
    if (min_distance < 5.5)
    {
        not_bonded = isNotBonded(atoms0.ID(), atoms1.ID());
    }
    
    if ( Q_LIKELY(not_bonded) )
    {
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);
                        r = r2.sqrt();

                        one_over_r = r.reciprocal();
                
                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_r + k_rf * r2 - c_rf;
                        tmp *= q * q1[j];
                    
                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= r.compareLess(Rc);

                        //make sure that the ID of atoms1 is not zero
                        itmp = id1[j].compareEqual(dummy_id);

                        icnrg += tmp.logicalAndNot(itmp);
                        
                        //Now do the LJ energy

                        //arithmetic combining rules
                        tmp = sig + (sig1[j]*sig1[j]);
                        tmp *= half;
                    
                        sig2_over_r2 = tmp * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    else
    {
        MultiFloat bonded_mask(0);
    
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    const bool *row = bondMatrix().constData()[id[0]].constData();

                    for (int j=0; j<n1; ++j)
                    {
                        //create a mask to cancel out calculations of bonded pairs
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            bonded_mask.quickSet(k, row[id1[j][k]] ? 0.0 : 1.0);
                        }

                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);
                        r = r2.sqrt();

                        one_over_r = r.reciprocal();
                
                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_r + k_rf * r2 - c_rf;
                        tmp *= q * q1[j];
                    
                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= r.compareLess(Rc);

                        //make sure that the ID of atoms1 is not zero
                        itmp = id1[j].compareEqual(dummy_id);

                        icnrg += bonded_mask * tmp.logicalAndNot(itmp);
                        
                        //Now do the LJ energy

                        //arithmetic combining rules
                        tmp = sig + (sig1[j]*sig1[j]);
                        tmp *= half;
                    
                        sig2_over_r2 = tmp * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += bonded_mask * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb and LJ intramolecular energy of all of the atoms in 'atoms',
    assuming periodic boundary conditions in a cubic box of size 'box_dimensions',
    returning the results in 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJIntraRFFunction::calcBoxEnergyAri(const CLJAtoms &atoms, const Vector &box_dimensions,
                                          double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii] * siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                const bool *row = bondMatrix().constData()[id[0]].constData();

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //get the bond mask to screen out bonded interactions
                    for (int k=0; k<MultiFloat::count(); ++k)
                    {
                        bond_mask.quickSet(k, row[ida[j][k]] ? 0.0 : 1.0);
                    }
                
                    scale *= bond_mask;
                
                    //calculate the distance between the fixed and mobile atoms
                    tmp = xa[j] - x;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                    r2 = tmp * tmp;

                    tmp = ya[j] - y;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    tmp = za[j] - z;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    r = r2.sqrt();

                    one_over_r = r.reciprocal();
            
                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_r + k_rf * r2 - c_rf;
                    tmp *= q * qa[j];
                
                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= r.compareLess(Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0 and that we are not
                    //including the energy of the atom with itself
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);

                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    tmp = sig + (siga[j]*siga[j]);
                    tmp *= half;
                
                    sig2_over_r2 = tmp * one_over_r;
                    sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                    sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                    tmp = sig6_over_r6 * sig6_over_r6;
                    tmp -= sig6_over_r6;
                    tmp *= eps;
                    tmp *= epsa[j];
                
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r.compareLess(Rlj);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intramolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', assuming periodic boundary conditions in a cubic box
    of size 'box_dimensions, returning the result in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJIntraRFFunction::calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                          const Vector &box_dimensions,
                                          double &cnrg, double &ljnrg, float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();
    
    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj(lj_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(coul_cutoff)) * ( (dielectric()-1) /
                                                          (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / coul_cutoff ) * ( (3*dielectric()) /
                                                    (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];

    MultiFloat tmp, r, r2, one_over_r, sig2_over_r2, sig6_over_r6;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();
    
    bool not_bonded = true;
    
    if (min_distance < 5.5)
    {
        not_bonded = isNotBonded(atoms0.ID(), atoms1.ID());
    }
    
    if ( Q_LIKELY(not_bonded) )
    {
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);
                        
                        r = r2.sqrt();

                        one_over_r = r.reciprocal();
                
                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_r + k_rf * r2 - c_rf;
                        tmp *= q * q1[j];
                    
                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= r.compareLess(Rc);

                        //make sure that the ID of atoms1 is not zero
                        itmp = id1[j].compareEqual(dummy_id);

                        icnrg += tmp.logicalAndNot(itmp);
                        
                        //Now do the LJ energy

                        //arithmetic combining rules
                        tmp = sig + (sig1[j]*sig1[j]);
                        tmp *= half;
                    
                        sig2_over_r2 = tmp * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    else
    {
        MultiFloat bonded_mask(0);
    
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii] * sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    const bool *row = bondMatrix().constData()[id[0]].constData();

                    for (int j=0; j<n1; ++j)
                    {
                        //create a mask to cancel out calculations of bonded pairs
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            bonded_mask.quickSet(k, row[id1[j][k]] ? 0.0 : 1.0);
                        }

                        //calculate the distance between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);
                        
                        r = r2.sqrt();

                        one_over_r = r.reciprocal();
                
                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_r + k_rf * r2 - c_rf;
                        tmp *= q * q1[j];
                    
                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= r.compareLess(Rc);

                        //make sure that the ID of atoms1 is not zero
                        itmp = id1[j].compareEqual(dummy_id);

                        icnrg += bonded_mask * tmp.logicalAndNot(itmp);
                        
                        //Now do the LJ energy

                        //arithmetic combining rules
                        tmp = sig + (sig1[j]*sig1[j]);
                        tmp *= half;
                    
                        sig2_over_r2 = tmp * one_over_r;
                        sig2_over_r2 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig2_over_r2*sig2_over_r2;
                        sig6_over_r6 = sig6_over_r6*sig2_over_r2;

                        tmp = sig6_over_r6 * sig6_over_r6;
                        tmp -= sig6_over_r6;
                        tmp *= eps;
                        tmp *= eps1[j];
                    
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r.compareLess(Rlj);
                        iljnrg += bonded_mask * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/////////
///////// Implementation of CLJSoftIntraRFFunction
/////////

static const RegisterMetaType<CLJSoftIntraRFFunction> r_softintra;

QDataStream &operator<<(QDataStream &ds, const CLJSoftIntraRFFunction &intra)
{
    writeHeader(ds, r_softintra, 1);
    
    ds << intra.diel << static_cast<const CLJSoftIntraFunction&>(intra);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJSoftIntraRFFunction &intra)
{
    VersionID v = readHeader(ds, r_softintra);
    
    if (v == 1)
    {
        ds >> intra.diel >> static_cast<CLJSoftIntraFunction&>(intra);
    }
    else
        throw version_error(v, "1", r_softintra, CODELOC);
    
    return ds;
}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction()
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(),
                       diel(default_dielectric)
{}

CLJFunctionPtr CLJSoftIntraRFFunction::defaultRFFunction()
{
    static CLJFunctionPtr ptr( new CLJSoftIntraRFFunction() );
    return ptr;
}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(Length cutoff)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(cutoff),
                       diel(default_dielectric)
{}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(Length coul_cutoff, Length lj_cutoff)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(
                            coul_cutoff, lj_cutoff),
                       diel(default_dielectric)
{}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(const Space &space, Length cutoff)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(
                            space, cutoff),
                       diel(default_dielectric)
{}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(const Space &space,
                                             Length coul_cutoff, Length lj_cutoff)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(
                            space, coul_cutoff, lj_cutoff),
                       diel(default_dielectric)
{}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(Length cutoff, COMBINING_RULES combining_rules)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(
                            cutoff, combining_rules),
                       diel(default_dielectric)
{}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(Length coul_cutoff, Length lj_cutoff,
                                             COMBINING_RULES combining_rules)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(
                            coul_cutoff, lj_cutoff, combining_rules),
                       diel(default_dielectric)
{}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(const Space &space,
                                                     COMBINING_RULES combining_rules)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(
                            space, combining_rules),
                       diel(default_dielectric)
{}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(const Space &space, Length cutoff,
                                             COMBINING_RULES combining_rules)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(
                            space, cutoff, combining_rules),
                       diel(default_dielectric)
{}

CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(const Space &space, Length coul_cutoff,
                                             Length lj_cutoff,
                                             COMBINING_RULES combining_rules)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(
                            space, coul_cutoff, lj_cutoff, combining_rules),
                       diel(default_dielectric)
{}

/** Copy constructor */
CLJSoftIntraRFFunction::CLJSoftIntraRFFunction(const CLJSoftIntraRFFunction &other)
        : ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>(other),
                       diel(other.diel)
{}

/** Destructor */
CLJSoftIntraRFFunction::~CLJSoftIntraRFFunction()
{}

/** Copy assignment operator */
CLJSoftIntraRFFunction&
CLJSoftIntraRFFunction::operator=(const CLJSoftIntraRFFunction &other)
{
    diel = other.diel;
    CLJSoftIntraFunction::operator=(other);
    return *this;
}

/** Comparison operator */
bool CLJSoftIntraRFFunction::operator==(const CLJSoftIntraRFFunction &other) const
{
    return diel == other.diel and CLJSoftIntraFunction::operator==(other);
}

/** Comparison operator */
bool CLJSoftIntraRFFunction::operator!=(const CLJSoftIntraRFFunction &other) const
{
    return not operator==(other);
}

const char* CLJSoftIntraRFFunction::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJSoftIntraRFFunction>() );
}

const char* CLJSoftIntraRFFunction::what() const
{
    return CLJSoftIntraRFFunction::typeName();
}

CLJSoftIntraRFFunction* CLJSoftIntraRFFunction::clone() const
{
    return new CLJSoftIntraRFFunction(*this);
}

/** Return the properties of this function */
Properties CLJSoftIntraRFFunction::properties() const
{
    Properties props = CLJSoftIntraFunction::properties();
    props.setProperty("dielectric", NumberProperty(dielectric()));
    return props;
}

/** Return a copy of this function where the property 'name' has been set to the
    value 'value' */
CLJFunctionPtr CLJSoftIntraRFFunction::setProperty(const QString &name, const Property &value) const
{
    if (name == "dielectric")
    {
        CLJFunctionPtr ret(*this);
        ret.edit().asA<CLJSoftIntraRFFunction>()
                  .setDielectric( value.asA<NumberProperty>().value() );
        return ret;
    }
    else
        return CLJSoftIntraFunction::setProperty(name, value);
}

/** Return the value of the property with name 'name' */
PropertyPtr CLJSoftIntraRFFunction::property(const QString &name) const
{
    if (name == "dielectric")
    {
        return NumberProperty(dielectric());
    }
    else
    {
        return CLJSoftIntraFunction::property(name);
    }
}

/** Return whether or not this function contains a property called 'name' */
bool CLJSoftIntraRFFunction::containsProperty(const QString &name) const
{
    return (name == "dielectric") or CLJSoftIntraFunction::containsProperty(name);
}

/** Set the dielectric constant to 'dielectric' */
void CLJSoftIntraRFFunction::setDielectric(float d)
{
    diel = d;
}

/** Return the value of the dielectric constant */
float CLJSoftIntraRFFunction::dielectric() const
{
    return diel;
}

/** Calculate the coulomb and LJ intramolecular energy of all of the atoms in 'atoms',
    returning the results in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJSoftIntraRFFunction::calcVacEnergyGeo(const CLJAtoms &atoms,
                                              double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );
    
    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                const bool *row = bondMatrix().constData()[id[0]].constData();

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //get the bond mask to screen out bonded interactions
                    for (int k=0; k<MultiFloat::count(); ++k)
                    {
                        bond_mask.quickSet(k, row[ida[j][k]] ? 0.0 : 1.0);
                    }
                
                    scale *= bond_mask;
                
                    //calculate the distance^2 between the fixed and mobile atoms
                    tmp = xa[j] - x;
                    r2 = tmp * tmp;
                    tmp = ya[j] - y;
                    r2.multiplyAdd(tmp, tmp);
                    tmp = za[j] - z;
                    r2.multiplyAdd(tmp, tmp);
                    
                    soft_r2 = r2 + alfa;
                    soft_r = soft_r2.sqrt();

                    one_over_soft_r = soft_r.reciprocal();

                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                    tmp *= one_minus_alpha_to_n * q * qa[j];

                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= soft_r.compareLess(soft_Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0 and that we are not
                    //including the energy of the atom with itself
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);

                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    sigma = sig * siga[j];
                    delta_sigma_r2 = delta * sigma + r2;
                    
                    sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                    sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                    tmp = sig6_over_delta3 * sig6_over_delta3;
                    tmp -= sig6_over_delta3;
                    tmp *= eps;
                    tmp *= epsa[j];
                    
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r2.compareLess(Rlj2);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intramolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', returning the result in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJSoftIntraRFFunction::calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                              double &cnrg, double &ljnrg,
                                              float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );
    
    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();
    
    bool not_bonded = true;
    
    if (min_distance < 5.5)
    {
        not_bonded = isNotBonded(atoms0.ID(), atoms1.ID());
    }
    
    if ( Q_LIKELY(not_bonded) )
    {
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance^2 between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);
                        
                        soft_r2 = r2 + alfa;
                        soft_r = soft_r2.sqrt();

                        one_over_soft_r = soft_r.reciprocal();

                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                        tmp *= one_minus_alpha_to_n * q * q1[j];

                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= soft_r.compareLess(soft_Rc);

                        //make sure that the ID of atoms1 is not zero, and is
                        //also not the same as the atoms0 and that we are not
                        //including the energy of the atom with itself
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        icnrg += tmp.logicalAndNot(itmp);

                        //now the LJ energy
                        sigma = sig * sig1[j];
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                        
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    else
    {
        MultiFloat bonded_mask(0);
    
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    const bool *row = bondMatrix().constData()[id[0]].constData();

                    for (int j=0; j<n1; ++j)
                    {
                        //create a mask to cancel out calculations of bonded pairs
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            bonded_mask.quickSet(k, row[id1[j][k]] ? 0.0 : 1.0);
                        }

                        //calculate the distance^2 between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);
                        
                        soft_r2 = r2 + alfa;
                        soft_r = soft_r2.sqrt();

                        one_over_soft_r = soft_r.reciprocal();

                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                        tmp *= one_minus_alpha_to_n * q * q1[j];

                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= soft_r.compareLess(soft_Rc);

                        //make sure that the ID of atoms1 is not zero, and is
                        //also not the same as the atoms0 and that we are not
                        //including the energy of the atom with itself
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        icnrg += bonded_mask * tmp.logicalAndNot(itmp);

                        //now the LJ energy
                        sigma = sig * sig1[j];
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                        
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += bonded_mask * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb and LJ intramolecular energy of all of the atoms in 'atoms',
    assuming periodic boundary conditions in a cubic box of size 'box_dimensions',
    returning the results in 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJSoftIntraRFFunction::calcBoxEnergyGeo(const CLJAtoms &atoms,
                                              const Vector &box_dimensions,
                                              double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );
    
    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                const bool *row = bondMatrix().constData()[id[0]].constData();

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //get the bond mask to screen out bonded interactions
                    for (int k=0; k<MultiFloat::count(); ++k)
                    {
                        bond_mask.quickSet(k, row[ida[j][k]] ? 0.0 : 1.0);
                    }
                
                    scale *= bond_mask;
                
                    //calculate the distance^2
                    tmp = xa[j] - x;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                    r2 = tmp * tmp;

                    tmp = ya[j] - y;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    tmp = za[j] - z;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);
                    
                    soft_r2 = r2 + alfa;
                    soft_r = soft_r2.sqrt();

                    one_over_soft_r = soft_r.reciprocal();

                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                    tmp *= one_minus_alpha_to_n * q * qa[j];

                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= soft_r.compareLess(soft_Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0 and that we are not
                    //including the energy of the atom with itself
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);

                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    sigma = sig * siga[j];
                    delta_sigma_r2 = delta * sigma + r2;
                    
                    sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                    sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                    tmp = sig6_over_delta3 * sig6_over_delta3;
                    tmp -= sig6_over_delta3;
                    tmp *= eps;
                    tmp *= epsa[j];
                    
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r2.compareLess(Rlj2);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intramolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', assuming periodic boundary conditions in a cubic box
    of size 'box_dimensions, returning the result in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJSoftIntraRFFunction::calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                              const Vector &box_dimensions,
                                              double &cnrg, double &ljnrg,
                                              float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );
    
    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();
    
    bool not_bonded = true;
    
    if (min_distance < 5.5)
    {
        not_bonded = isNotBonded(atoms0.ID(), atoms1.ID());
    }
    
    if ( Q_LIKELY(not_bonded) )
    {
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance^2
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);
                        
                        soft_r2 = r2 + alfa;
                        soft_r = soft_r2.sqrt();

                        one_over_soft_r = soft_r.reciprocal();

                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                        tmp *= one_minus_alpha_to_n * q * q1[j];

                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= soft_r.compareLess(soft_Rc);

                        //make sure that the ID of atoms1 is not zero, and is
                        //also not the same as the atoms0 and that we are not
                        //including the energy of the atom with itself
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        icnrg += tmp.logicalAndNot(itmp);

                        //now the LJ energy
                        sigma = sig * sig1[j];
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                        
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    else
    {
        MultiFloat bonded_mask(0);
    
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    const bool *row = bondMatrix().constData()[id[0]].constData();

                    for (int j=0; j<n1; ++j)
                    {
                        //create a mask to cancel out calculations of bonded pairs
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            bonded_mask.quickSet(k, row[id1[j][k]] ? 0.0 : 1.0);
                        }

                        //calculate the distance^2
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);
                        
                        soft_r2 = r2 + alfa;
                        soft_r = soft_r2.sqrt();

                        one_over_soft_r = soft_r.reciprocal();

                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                        tmp *= one_minus_alpha_to_n * q * q1[j];

                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= soft_r.compareLess(soft_Rc);

                        //make sure that the ID of atoms1 is not zero, and is
                        //also not the same as the atoms0 and that we are not
                        //including the energy of the atom with itself
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        icnrg += bonded_mask * tmp.logicalAndNot(itmp);

                        //now the LJ energy
                        sigma = sig * sig1[j];
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                        
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += bonded_mask * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}


/** Calculate the coulomb and LJ intramolecular energy of all of the atoms in 'atoms',
    returning the results in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJSoftIntraRFFunction::calcVacEnergyAri(const CLJAtoms &atoms,
                                              double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );
    
    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii] * siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                const bool *row = bondMatrix().constData()[id[0]].constData();

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //get the bond mask to screen out bonded interactions
                    for (int k=0; k<MultiFloat::count(); ++k)
                    {
                        bond_mask.quickSet(k, row[ida[j][k]] ? 0.0 : 1.0);
                    }
                
                    scale *= bond_mask;
                
                    //calculate the distance^2 between the fixed and mobile atoms
                    tmp = xa[j] - x;
                    r2 = tmp * tmp;
                    tmp = ya[j] - y;
                    r2.multiplyAdd(tmp, tmp);
                    tmp = za[j] - z;
                    r2.multiplyAdd(tmp, tmp);
                    
                    soft_r2 = r2 + alfa;
                    soft_r = soft_r2.sqrt();

                    one_over_soft_r = soft_r.reciprocal();

                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                    tmp *= one_minus_alpha_to_n * q * qa[j];

                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= soft_r.compareLess(soft_Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0 and that we are not
                    //including the energy of the atom with itself
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);

                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    sigma = half * (sig + (siga[j]*siga[j]));
                    delta_sigma_r2 = delta * sigma + r2;
                    
                    sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                    sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                    tmp = sig6_over_delta3 * sig6_over_delta3;
                    tmp -= sig6_over_delta3;
                    tmp *= eps;
                    tmp *= epsa[j];
                    
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r2.compareLess(Rlj2);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intramolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', returning the result in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJSoftIntraRFFunction::calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                              double &cnrg, double &ljnrg,
                                              float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );
    
    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();
    
    bool not_bonded = true;
    
    if (min_distance < 5.5)
    {
        not_bonded = isNotBonded(atoms0.ID(), atoms1.ID());
    }
    
    if ( Q_LIKELY(not_bonded) )
    {
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]*sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance^2 between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);
                        
                        soft_r2 = r2 + alfa;
                        soft_r = soft_r2.sqrt();

                        one_over_soft_r = soft_r.reciprocal();

                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                        tmp *= one_minus_alpha_to_n * q * q1[j];

                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= soft_r.compareLess(soft_Rc);

                        //make sure that the ID of atoms1 is not zero, and is
                        //also not the same as the atoms0 and that we are not
                        //including the energy of the atom with itself
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        icnrg += tmp.logicalAndNot(itmp);

                        //now the LJ energy
                        sigma = half * (sig + (sig1[j]*sig1[j]));
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                        
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    else
    {
        MultiFloat bonded_mask(0);
    
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]*sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    const bool *row = bondMatrix().constData()[id[0]].constData();

                    for (int j=0; j<n1; ++j)
                    {
                        //create a mask to cancel out calculations of bonded pairs
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            bonded_mask.quickSet(k, row[id1[j][k]] ? 0.0 : 1.0);
                        }

                        //calculate the distance^2 between the fixed and mobile atoms
                        tmp = x1[j] - x;
                        r2 = tmp * tmp;
                        tmp = y1[j] - y;
                        r2.multiplyAdd(tmp, tmp);
                        tmp = z1[j] - z;
                        r2.multiplyAdd(tmp, tmp);
                        
                        soft_r2 = r2 + alfa;
                        soft_r = soft_r2.sqrt();

                        one_over_soft_r = soft_r.reciprocal();

                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                        tmp *= one_minus_alpha_to_n * q * q1[j];

                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= soft_r.compareLess(soft_Rc);

                        //make sure that the ID of atoms1 is not zero, and is
                        //also not the same as the atoms0 and that we are not
                        //including the energy of the atom with itself
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        icnrg += bonded_mask * tmp.logicalAndNot(itmp);

                        //now the LJ energy
                        sigma = half * (sig + (sig1[j]*sig1[j]));
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                        
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += bonded_mask * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the coulomb and LJ intramolecular energy of all of the atoms in 'atoms',
    assuming periodic boundary conditions in a cubic box of size 'box_dimensions',
    returning the results in 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJSoftIntraRFFunction::calcBoxEnergyAri(const CLJAtoms &atoms,
                                              const Vector &box_dimensions,
                                              double &cnrg, double &ljnrg) const
{
    const MultiFloat *xa = atoms.x().constData();
    const MultiFloat *ya = atoms.y().constData();
    const MultiFloat *za = atoms.z().constData();
    const MultiFloat *qa = atoms.q().constData();
    const MultiFloat *siga = atoms.sigma().constData();
    const MultiFloat *epsa = atoms.epsilon().constData();
    const MultiInt *ida = atoms.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );
    
    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    int n = atoms.x().count();
    
    for (int i=0; i<n; ++i)
    {
        for (int ii=0; ii<MultiFloat::size(); ++ii)
        {
            if (ida[i][ii] != dummy_int)
            {
                const MultiInt id( ida[i][ii] );
                const MultiFloat x( xa[i][ii] );
                const MultiFloat y( ya[i][ii] );
                const MultiFloat z( za[i][ii] );
                const MultiFloat q( qa[i][ii] );
                const MultiFloat sig( siga[i][ii]*siga[i][ii] );
                const MultiFloat eps( epsa[i][ii] );

                const bool *row = bondMatrix().constData()[id[0]].constData();

                for (int j=i; j<n; ++j)
                {
                    // if i == j then we double-calculate the energies, so must
                    // scale them by 0.5
                    MultiFloat scale( i == j ? 0.5 : 1.0 );
                
                    //get the bond mask to screen out bonded interactions
                    for (int k=0; k<MultiFloat::count(); ++k)
                    {
                        bond_mask.quickSet(k, row[ida[j][k]] ? 0.0 : 1.0);
                    }
                
                    scale *= bond_mask;
                
                    //calculate the distance^2
                    tmp = xa[j] - x;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                    r2 = tmp * tmp;

                    tmp = ya[j] - y;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);

                    tmp = za[j] - z;
                    tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                    tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                    r2.multiplyAdd(tmp, tmp);
                    
                    soft_r2 = r2 + alfa;
                    soft_r = soft_r2.sqrt();

                    one_over_soft_r = soft_r.reciprocal();

                    // calculate the coulomb energy using
                    // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                    // c = (1/r_c) * (3 eps)/(2 eps + 1)
                    tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                    tmp *= one_minus_alpha_to_n * q * qa[j];

                    //apply the cutoff - compare r against Rc. This will
                    //return 1 if r is less than Rc, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rc
                    tmp &= soft_r.compareLess(soft_Rc);

                    //make sure that the ID of atoms1 is not zero, and is
                    //also not the same as the atoms0 and that we are not
                    //including the energy of the atom with itself
                    itmp = ida[j].compareEqual(dummy_id);
                    itmp |= ida[j].compareEqual(id);

                    icnrg += scale * tmp.logicalAndNot(itmp);

                    //now the LJ energy
                    sigma = half * (sig + (siga[j]*siga[j]));
                    delta_sigma_r2 = delta * sigma + r2;
                    
                    sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                    sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                    tmp = sig6_over_delta3 * sig6_over_delta3;
                    tmp -= sig6_over_delta3;
                    tmp *= eps;
                    tmp *= epsa[j];
                    
                    //apply the cutoff - compare r against Rlj. This will
                    //return 1 if r is less than Rlj, or 0 otherwise. Logical
                    //and will then remove all energies where r >= Rlj
                    tmp &= r2.compareLess(Rlj2);
                    iljnrg += scale * tmp.logicalAndNot(itmp);
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}

/** Calculate the intramolecular energy between all atoms in 'atoms0' and all
    atoms in 'atoms1', assuming periodic boundary conditions in a cubic box
    of size 'box_dimensions, returning the result in the arguments 'cnrg' and 'ljnrg'.
    Note that all of the atoms must be part of the same molecule, and must
    have their intramolecular non-bonded scale factors loaded into this function */
void CLJSoftIntraRFFunction::calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                              const Vector &box_dimensions,
                                              double &cnrg, double &ljnrg,
                                              float min_distance) const
{
    const MultiFloat *x0 = atoms0.x().constData();
    const MultiFloat *y0 = atoms0.y().constData();
    const MultiFloat *z0 = atoms0.z().constData();
    const MultiFloat *q0 = atoms0.q().constData();
    const MultiFloat *sig0 = atoms0.sigma().constData();
    const MultiFloat *eps0 = atoms0.epsilon().constData();
    const MultiInt *id0 = atoms0.ID().constData();

    const MultiFloat *x1 = atoms1.x().constData();
    const MultiFloat *y1 = atoms1.y().constData();
    const MultiFloat *z1 = atoms1.z().constData();
    const MultiFloat *q1 = atoms1.q().constData();
    const MultiFloat *sig1 = atoms1.sigma().constData();
    const MultiFloat *eps1 = atoms1.epsilon().constData();
    const MultiInt *id1 = atoms1.ID().constData();

    const MultiFloat Rc(coul_cutoff);
    const MultiFloat Rlj2(lj_cutoff*lj_cutoff);

    const float soft_coul_cutoff = std::sqrt(alpha() + coul_cutoff*coul_cutoff);

    const MultiFloat soft_Rc(soft_coul_cutoff);

    const MultiFloat k_rf( (1.0 / pow_3(soft_coul_cutoff)) * ( (dielectric()-1) /
                                                               (2*dielectric() + 1) ) );
    const MultiFloat c_rf( (1.0 / soft_coul_cutoff ) * ( (3*dielectric()) /
                                                         (2*dielectric() + 1) ) );
    const MultiFloat half(0.5);
    const MultiInt dummy_id = CLJAtoms::idOfDummy();
    const qint32 dummy_int = dummy_id[0];
    const MultiFloat one_minus_alpha_to_n( this->oneMinusAlphaToN() );
    const MultiFloat delta( this->alphaTimesShiftDelta() );
    const MultiFloat alfa( this->alpha() );
    
    MultiFloat tmp, r2, soft_r, soft_r2, one_over_soft_r, sigma, delta_sigma_r2;
    MultiFloat sig2_over_delta, sig6_over_delta3;
    MultiFloat bond_mask;
    MultiDouble icnrg(0), iljnrg(0);
    MultiInt itmp;

    const MultiFloat box_x( box_dimensions.x() );
    const MultiFloat box_y( box_dimensions.y() );
    const MultiFloat box_z( box_dimensions.z() );
    
    const MultiFloat half_box_x( 0.5 * box_dimensions.x() );
    const MultiFloat half_box_y( 0.5 * box_dimensions.y() );
    const MultiFloat half_box_z( 0.5 * box_dimensions.z() );

    const int n0 = atoms0.x().count();
    const int n1 = atoms1.x().count();
    
    bool not_bonded = true;
    
    if (min_distance < 5.5)
    {
        not_bonded = isNotBonded(atoms0.ID(), atoms1.ID());
    }
    
    if ( Q_LIKELY(not_bonded) )
    {
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]*sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    for (int j=0; j<n1; ++j)
                    {
                        //calculate the distance^2
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);
                        
                        soft_r2 = r2 + alfa;
                        soft_r = soft_r2.sqrt();

                        one_over_soft_r = soft_r.reciprocal();

                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                        tmp *= one_minus_alpha_to_n * q * q1[j];

                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= soft_r.compareLess(soft_Rc);

                        //make sure that the ID of atoms1 is not zero, and is
                        //also not the same as the atoms0 and that we are not
                        //including the energy of the atom with itself
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        icnrg += tmp.logicalAndNot(itmp);

                        //now the LJ energy
                        sigma = half * (sig + (sig1[j]*sig1[j]));
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                        
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    else
    {
        MultiFloat bonded_mask(0);
    
        for (int i=0; i<n0; ++i)
        {
            for (int ii=0; ii<MultiFloat::count(); ++ii)
            {
                if (id0[i][ii] != dummy_int)
                {
                    const MultiInt id(id0[i][ii]);
                    const MultiFloat x(x0[i][ii]);
                    const MultiFloat y(y0[i][ii]);
                    const MultiFloat z(z0[i][ii]);
                    const MultiFloat q(q0[i][ii]);

                    const MultiFloat sig(sig0[i][ii]*sig0[i][ii]);
                    const MultiFloat eps(eps0[i][ii]);

                    const bool *row = bondMatrix().constData()[id[0]].constData();

                    for (int j=0; j<n1; ++j)
                    {
                        //create a mask to cancel out calculations of bonded pairs
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            bonded_mask.quickSet(k, row[id1[j][k]] ? 0.0 : 1.0);
                        }

                        //calculate the distance^2
                        tmp = x1[j] - x;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_x.logicalAnd( half_box_x.compareLess(tmp) );
                        r2 = tmp * tmp;

                        tmp = y1[j] - y;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_y.logicalAnd( half_box_y.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);

                        tmp = z1[j] - z;
                        tmp &= MULTIFLOAT_POS_MASK;  // this creates the absolute value :-)
                        tmp -= box_z.logicalAnd( half_box_z.compareLess(tmp) );
                        r2.multiplyAdd(tmp, tmp);
                        
                        soft_r2 = r2 + alfa;
                        soft_r = soft_r2.sqrt();

                        one_over_soft_r = soft_r.reciprocal();

                        // calculate the coulomb energy using
                        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
                        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
                        // c = (1/r_c) * (3 eps)/(2 eps + 1)
                        tmp = one_over_soft_r + k_rf * soft_r2 - c_rf;
                        tmp *= one_minus_alpha_to_n * q * q1[j];

                        //apply the cutoff - compare r against Rc. This will
                        //return 1 if r is less than Rc, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rc
                        tmp &= soft_r.compareLess(soft_Rc);

                        //make sure that the ID of atoms1 is not zero, and is
                        //also not the same as the atoms0 and that we are not
                        //including the energy of the atom with itself
                        itmp = id1[j].compareEqual(dummy_id);
                        itmp |= id1[j].compareEqual(id);

                        icnrg += bonded_mask * tmp.logicalAndNot(itmp);

                        //now the LJ energy
                        sigma = half * (sig + (sig1[j]*sig1[j]));
                        delta_sigma_r2 = delta * sigma + r2;
                        
                        sig2_over_delta = (sigma*sigma) / delta_sigma_r2;
                        sig6_over_delta3 = sig2_over_delta * sig2_over_delta * sig2_over_delta;

                        tmp = sig6_over_delta3 * sig6_over_delta3;
                        tmp -= sig6_over_delta3;
                        tmp *= eps;
                        tmp *= eps1[j];
                        
                        //apply the cutoff - compare r against Rlj. This will
                        //return 1 if r is less than Rlj, or 0 otherwise. Logical
                        //and will then remove all energies where r >= Rlj
                        tmp &= r2.compareLess(Rlj2);
                        iljnrg += bonded_mask * tmp.logicalAndNot(itmp);
                    }
                }
            }
        }
    }
    
    cnrg = icnrg.sum();
    ljnrg = iljnrg.sum();
}
