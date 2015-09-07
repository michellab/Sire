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

#ifndef SIREMM_CLJRSFUNCTION_H
#define SIREMM_CLJRSFUNCTION_H

#include "cljfunction.h"

#include <QHash>

namespace SireMM
{
class CLJRFFunction;
class CLJIntraRFFunction;
class CLJSoftRFFunction;
class CLJSoftIntraRFFunction;
}

QDataStream& operator<<(QDataStream&, const SireMM::CLJRFFunction&);
QDataStream& operator>>(QDataStream&, SireMM::CLJRFFunction&);

QDataStream& operator<<(QDataStream&, const SireMM::CLJIntraRFFunction&);
QDataStream& operator>>(QDataStream&, SireMM::CLJIntraRFFunction&);

QDataStream& operator<<(QDataStream&, const SireMM::CLJSoftRFFunction&);
QDataStream& operator>>(QDataStream&, SireMM::CLJSoftRFFunction&);

QDataStream& operator<<(QDataStream&, const SireMM::CLJSoftIntraRFFunction&);
QDataStream& operator>>(QDataStream&, SireMM::CLJSoftIntraRFFunction&);

namespace SireMM
{

/** This CLJFunction calculates the intermolecular coulomb and LJ energy of the passed
    CLJAtoms using a reaction field cutoff
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJRFFunction
        : public SireBase::ConcreteProperty<CLJRFFunction,CLJCutoffFunction>
{

friend QDataStream& ::operator<<(QDataStream&, const CLJRFFunction&);
friend QDataStream& ::operator>>(QDataStream&, CLJRFFunction&);

public:
    CLJRFFunction();
    CLJRFFunction(Length cutoff);
    CLJRFFunction(Length coul_cutoff, Length lj_cutoff);
    
    CLJRFFunction(const Space &space, Length cutoff);
    CLJRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff);
    
    CLJRFFunction(Length cutoff, COMBINING_RULES combining_rules);
    CLJRFFunction(Length coul_cutoff, Length lj_cutoff, COMBINING_RULES combining_rules);
    
    CLJRFFunction(const Space &space, COMBINING_RULES combining_rules);
    CLJRFFunction(const Space &space, Length cutoff, COMBINING_RULES combining_rules);
    CLJRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff,
                     COMBINING_RULES combining_rules);
    
    CLJRFFunction(const CLJRFFunction &other);
    
    ~CLJRFFunction();
    
    CLJRFFunction& operator=(const CLJRFFunction &other);
    
    bool operator==(const CLJRFFunction &other) const;
    bool operator!=(const CLJRFFunction &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    QString toString() const;
    
    CLJRFFunction* clone() const;

    Properties properties() const;
    CLJFunctionPtr setProperty(const QString &name, const Property &value) const;
    PropertyPtr property(const QString &name) const;
    bool containsProperty(const QString &name) const;

    void setDielectric(float dielectric);
    float dielectric() const;

    bool supportsGridCalculation() const;

    static CLJFunctionPtr defaultRFFunction();

protected:
    void calcVacEnergyAri(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcVacEnergyGeo(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcBoxEnergyAri(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

    void calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

    double calcVacCoulombEnergyAri(const CLJAtoms &atoms) const;
    double calcVacCoulombEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                   float min_distance) const;
    
    double calcVacLJEnergyAri(const CLJAtoms &atoms) const;
    double calcVacLJEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                              float min_distance) const;

    void calcVacGrid(const CLJAtoms &atoms, const GridInfo &gridinfo,
                     const int start, const int end, float *gridpot) const;
    
    void calcBoxGrid(const CLJAtoms &atoms, const GridInfo &gridinfo,
                     const Vector &box_dimensions,
                     const int start, const int end, float *gridpot) const;

private:
    /** The dielectric constant */
    float diel;
};

/** This CLJFunction calculates the intramolecular coulomb and LJ energy of the passed
    CLJAtoms using a reaction field cutoff
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJIntraRFFunction
        : public SireBase::ConcreteProperty<CLJIntraRFFunction,CLJIntraFunction>
{

friend QDataStream& ::operator<<(QDataStream&, const CLJIntraRFFunction&);
friend QDataStream& ::operator>>(QDataStream&, CLJIntraRFFunction&);

public:
    CLJIntraRFFunction();
    CLJIntraRFFunction(Length cutoff);
    CLJIntraRFFunction(Length coul_cutoff, Length lj_cutoff);
    
    CLJIntraRFFunction(const Space &space, Length cutoff);
    CLJIntraRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff);
    
    CLJIntraRFFunction(Length cutoff, COMBINING_RULES combining_rules);
    CLJIntraRFFunction(Length coul_cutoff, Length lj_cutoff, COMBINING_RULES combining_rules);
    
    CLJIntraRFFunction(const Space &space, COMBINING_RULES combining_rules);
    CLJIntraRFFunction(const Space &space, Length cutoff, COMBINING_RULES combining_rules);
    CLJIntraRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff,
                          COMBINING_RULES combining_rules);
    
    CLJIntraRFFunction(const CLJIntraRFFunction &other);
    
    ~CLJIntraRFFunction();
    
    CLJIntraRFFunction& operator=(const CLJIntraRFFunction &other);
    
    bool operator==(const CLJIntraRFFunction &other) const;
    bool operator!=(const CLJIntraRFFunction &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    CLJIntraRFFunction* clone() const;

    Properties properties() const;
    CLJFunctionPtr setProperty(const QString &name, const Property &value) const;
    PropertyPtr property(const QString &name) const;
    bool containsProperty(const QString &name) const;

    void setDielectric(float dielectric);
    float dielectric() const;

    static CLJFunctionPtr defaultRFFunction();

protected:
    void calcVacEnergyAri(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcVacEnergyGeo(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcBoxEnergyAri(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg
                          , float min_distance) const;

    void calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

private:
    /** The dielectric constant */
    float diel;
};


/** This CLJFunction calculates the intermolecular coulomb and LJ energy of the passed
    CLJAtoms using a reaction field cutoff, and provides a soft-core
    potential that can soften molecules that are involved in free energy calculations
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJSoftRFFunction
        : public SireBase::ConcreteProperty<CLJSoftRFFunction,CLJSoftFunction>
{

friend QDataStream& ::operator<<(QDataStream&, const CLJSoftRFFunction&);
friend QDataStream& ::operator>>(QDataStream&, CLJSoftRFFunction&);

public:
    CLJSoftRFFunction();
    CLJSoftRFFunction(Length cutoff);
    CLJSoftRFFunction(Length coul_cutoff, Length lj_cutoff);
    
    CLJSoftRFFunction(const Space &space, Length cutoff);
    CLJSoftRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff);
    
    CLJSoftRFFunction(Length cutoff, COMBINING_RULES combining_rules);
    CLJSoftRFFunction(Length coul_cutoff, Length lj_cutoff, COMBINING_RULES combining_rules);
    
    CLJSoftRFFunction(const Space &space, COMBINING_RULES combining_rules);
    CLJSoftRFFunction(const Space &space, Length cutoff, COMBINING_RULES combining_rules);
    CLJSoftRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff,
                     COMBINING_RULES combining_rules);
    
    CLJSoftRFFunction(const CLJSoftRFFunction &other);
    
    ~CLJSoftRFFunction();
    
    CLJSoftRFFunction& operator=(const CLJSoftRFFunction &other);
    
    bool operator==(const CLJSoftRFFunction &other) const;
    bool operator!=(const CLJSoftRFFunction &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    QString toString() const;
    
    CLJSoftRFFunction* clone() const;

    Properties properties() const;
    CLJFunctionPtr setProperty(const QString &name, const Property &value) const;
    PropertyPtr property(const QString &name) const;
    bool containsProperty(const QString &name) const;

    void setDielectric(float dielectric);
    float dielectric() const;

    static CLJFunctionPtr defaultRFFunction();

protected:
    void calcVacEnergyAri(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcVacEnergyGeo(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcBoxEnergyAri(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

    void calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

private:
    /** The dielectric constant */
    float diel;
};

/** This CLJFunction calculates the intramolecular coulomb and LJ energy of the passed
    CLJAtoms using a reaction field cutoff, and provides a soft-core
    potential that can soften molecules that are involved in free energy calculations
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJSoftIntraRFFunction
        : public SireBase::ConcreteProperty<CLJSoftIntraRFFunction,CLJSoftIntraFunction>
{

friend QDataStream& ::operator<<(QDataStream&, const CLJSoftIntraRFFunction&);
friend QDataStream& ::operator>>(QDataStream&, CLJSoftIntraRFFunction&);

public:
    CLJSoftIntraRFFunction();
    CLJSoftIntraRFFunction(Length cutoff);
    CLJSoftIntraRFFunction(Length coul_cutoff, Length lj_cutoff);
    
    CLJSoftIntraRFFunction(const Space &space, Length cutoff);
    CLJSoftIntraRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff);
    
    CLJSoftIntraRFFunction(Length cutoff, COMBINING_RULES combining_rules);
    CLJSoftIntraRFFunction(Length coul_cutoff, Length lj_cutoff,
                              COMBINING_RULES combining_rules);
    
    CLJSoftIntraRFFunction(const Space &space, COMBINING_RULES combining_rules);
    CLJSoftIntraRFFunction(const Space &space, Length cutoff, COMBINING_RULES combining_rules);
    CLJSoftIntraRFFunction(const Space &space, Length coul_cutoff, Length lj_cutoff,
                              COMBINING_RULES combining_rules);
    
    CLJSoftIntraRFFunction(const CLJSoftIntraRFFunction &other);
    
    ~CLJSoftIntraRFFunction();
    
    CLJSoftIntraRFFunction& operator=(const CLJSoftIntraRFFunction &other);
    
    bool operator==(const CLJSoftIntraRFFunction &other) const;
    bool operator!=(const CLJSoftIntraRFFunction &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    CLJSoftIntraRFFunction* clone() const;

    Properties properties() const;
    CLJFunctionPtr setProperty(const QString &name, const Property &value) const;
    PropertyPtr property(const QString &name) const;
    bool containsProperty(const QString &name) const;

    void setDielectric(float dielectric);
    float dielectric() const;

    static CLJFunctionPtr defaultRFFunction();

protected:
    void calcVacEnergyAri(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcVacEnergyGeo(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcBoxEnergyAri(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

    void calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

private:
    /** The dielectric constant */
    float diel;
};

}

Q_DECLARE_METATYPE( SireMM::CLJRFFunction )
Q_DECLARE_METATYPE( SireMM::CLJIntraRFFunction )

Q_DECLARE_METATYPE( SireMM::CLJSoftRFFunction )
Q_DECLARE_METATYPE( SireMM::CLJSoftIntraRFFunction );

SIRE_EXPOSE_CLASS( SireMM::CLJRFFunction )
SIRE_EXPOSE_CLASS( SireMM::CLJIntraRFFunction )

SIRE_EXPOSE_CLASS( SireMM::CLJSoftRFFunction )
SIRE_EXPOSE_CLASS( SireMM::CLJSoftIntraRFFunction )

#endif
