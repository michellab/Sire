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

#include "cljfunction.h"
#include "cljboxes.h"
#include "switchingfunction.h"

#include "SireMaths/multidouble.h"

#include "SireVol/cartesian.h"
#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"
#include "SireVol/gridinfo.h"

#include "SireError/errors.h"

#include "SireBase/properties.h"
#include "SireBase/stringproperty.h"
#include "SireBase/numberproperty.h"
#include "SireBase/generalunitproperty.h"
#include "SireBase/lengthproperty.h"
#include "SireBase/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

using namespace SireMM;
using namespace SireVol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireID;
using namespace SireUnits::Dimension;
using namespace SireStream;

/////////
///////// Implementation of CLJFunction
/////////

static const RegisterMetaType<CLJFunction> r_cljfunc(MAGIC_ONLY, CLJFunction::typeName());

QDataStream &operator<<(QDataStream &ds, const CLJFunction &cljfunc)
{
    writeHeader(ds, r_cljfunc, 1);

    SharedDataStream sds(ds);

    sds << cljfunc.spce << cljfunc.use_arithmetic;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJFunction &cljfunc)
{
    VersionID v = readHeader(ds, r_cljfunc);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> cljfunc.spce >> cljfunc.use_arithmetic;
        cljfunc.extractDetailsFromSpace();
    }
    else
        throw version_error(v, "1", r_cljfunc, CODELOC);

    return ds;
}

/** Constructor. By default we will use vacuum boundary conditions with
    arithmetic combining rules */
CLJFunction::CLJFunction() : Property(), use_arithmetic(true), use_box(false)
{}

void CLJFunction::extractDetailsFromRules(SireMM::CLJFunction::COMBINING_RULES rules)
{
    switch (rules)
    {
        case ARITHMETIC:
            use_arithmetic = true;
            break;
        case GEOMETRIC:
            use_arithmetic = false;
            break;
        default:
            throw SireError::program_bug( QObject::tr(
                        "Cannot understand the passed combining rules (%1)")
                            .arg(rules), CODELOC );
    }
}

/** Construct, using vacuum boundary conditions, but specifying the combining rules */
CLJFunction::CLJFunction(COMBINING_RULES combining_rules)
            : Property(), use_arithmetic(true), use_box(false)
{
    extractDetailsFromRules(combining_rules);
}

void CLJFunction::extractDetailsFromSpace()
{
    if (spce.isNull())
    {
        use_box = false;
        box_dimensions = Vector(0);
    }
    else if (spce.read().isA<PeriodicBox>())
    {
        use_box = true;
        box_dimensions = spce.read().asA<PeriodicBox>().dimensions();
    }
    else if (spce.read().isA<TriclinicBox>())
    {
        // Support Cartesian triclinic boxes. Here we work out the appropriate
        // box dimensions by extracting the lattice vectors of the triclinic cell.
        if (spce.read().asA<TriclinicBox>().isCartesian())
        {
            use_box = true;
            auto v0 = spce.read().asA<TriclinicBox>().vector0();
            auto v1 = spce.read().asA<TriclinicBox>().vector1();
            auto v2 = spce.read().asA<TriclinicBox>().vector2();
            box_dimensions = Vector(v0.x(), v1.y(), v2.z());
        }
    }
    else if (spce.read().isA<Cartesian>())
    {
        use_box = false;
        box_dimensions = Vector(0);
    }
    else
        throw SireError::unsupported( QObject::tr(
                "CLJFunction-based forcefields currently only support using either "
                "periodic (cubic) boundary conditions, or vacuum boundary conditions. "
                "They are not compatible with the passed space \"%1\".")
                    .arg(spce.read().toString()), CODELOC );
}

/** Construct, using arithmetic combining rules, but specifying the space */
CLJFunction::CLJFunction(const Space &space)
            : Property(), spce(space), use_arithmetic(true), use_box(false)
{
    extractDetailsFromSpace();
}

/** Construct, specifying both the combining rules and simulation space */
CLJFunction::CLJFunction(const Space &space, COMBINING_RULES combining_rules)
            : Property(), spce(space), use_arithmetic(true), use_box(false)
{
    extractDetailsFromSpace();
    extractDetailsFromRules(combining_rules);
}

/** Copy constructor */
CLJFunction::CLJFunction(const CLJFunction &other)
            : Property(other), spce(other.spce), box_dimensions(other.box_dimensions),
              use_arithmetic(other.use_arithmetic), use_box(other.use_box)
{}

/** Destructor */
CLJFunction::~CLJFunction()
{}

/** Copy assignment operator */
CLJFunction& CLJFunction::operator=(const CLJFunction &other)
{
    if (this != &other)
    {
        spce = other.spce;
        box_dimensions = other.box_dimensions;
        use_arithmetic = other.use_arithmetic;
        use_box = other.use_box;
        Property::operator=(other);
    }

    return *this;
}

const char* CLJFunction::typeName()
{
    return "SireMM::CLJFunction";
}

/** Comparison operator */
bool CLJFunction::operator==(const CLJFunction &other) const
{
    return use_arithmetic == other.use_arithmetic and
           spce == other.spce and Property::operator==(other);
}

/** Return the null (do nothing) function */
const NullCLJFunction& CLJFunction::null()
{
    static NullCLJFunction nullfunc;
    return nullfunc;
}

/** Return all of the configurable properties of this function */
Properties CLJFunction::properties() const
{
    Properties props;
    props.setProperty( "space", spce );

    if (use_arithmetic)
        props.setProperty( "combiningRules", StringProperty("arithmetic") );
    else
        props.setProperty( "combiningRules", StringProperty("geometric") );

    return props;
}

/** Return a copy of this function where the property 'name' has been set to
    the value 'value'

    \throw SireBase::missing_property
*/
CLJFunctionPtr CLJFunction::setProperty(const QString &name, const Property &value) const
{
    CLJFunctionPtr ret(*this);

    if (name == "space")
    {
        ret.edit().setSpace( value.asA<Space>() );
    }
    else if (name == "combiningRules")
    {
        QString typ = value.asAString();

        if (typ.toLower() == "arithmetic")
            ret.edit().setCombiningRules( CLJFunction::ARITHMETIC );
        else if (typ.toLower() == "geometric")
            ret.edit().setCombiningRules( CLJFunction::GEOMETRIC );
        else
            throw SireError::invalid_arg( QObject::tr(
                    "Cannot interpret combining rules from value \"%1\". Valid "
                    "values are \"arithmetic\" or \"geometric\".")
                        .arg(typ), CODELOC );
    }
    else
    {
        throw SireBase::missing_property( QObject::tr(
                "There is no property with the name \"%1\" in function \"%2\". "
                "Available property names are %3.")
                    .arg(name).arg(this->toString())
                    .arg( Sire::toString(this->properties().propertyKeys()) ), CODELOC );
    }

    return ret;
}

/** Return the value of the property with name 'name' */
PropertyPtr CLJFunction::property(const QString &name) const
{
    if (name == "space")
        return spce;

    else if (name == "combiningRules")
    {
        static StringProperty arithmetic_property("arithmetic");
        static StringProperty geometric_property("geometric");

        if (use_arithmetic)
            return arithmetic_property;
        else
            return geometric_property;
    }
    else
    {
        throw SireBase::missing_property( QObject::tr(
                "There is no property with the name \"%1\" in function \"%2\". "
                "Available property names are %3.")
                    .arg(name).arg(this->toString())
                    .arg( Sire::toString(this->properties().propertyKeys()) ), CODELOC );

        // the code below is never executed
        static NullCLJFunction func;
        return func;
    }
}

/** Return whether or not this function contains a property with name 'name' */
bool CLJFunction::containsProperty(const QString &name) const
{
    return (name == "space") or (name == "combiningRules");
}

/** Tell the function to use arithmetic combining rules for LJ parameters */
void CLJFunction::setArithmeticCombiningRules(bool on)
{
    use_arithmetic = on;
}

/** Tell the function to use geometric combining rules for LJ parameters */
void CLJFunction::setGeometricCombiningRules(bool on)
{
    use_arithmetic = not on;
}

/** Return the combining rules used by the function */
SireMM::CLJFunction::COMBINING_RULES CLJFunction::combiningRules() const
{
    if (use_arithmetic)
        return ARITHMETIC;
    else
        return GEOMETRIC;
}

/** Return whether or not this function uses a softened (soft-core) potential */
bool CLJFunction::isSoftened() const
{
    return false;
}

/** Set the combining rules used by the function */
void CLJFunction::setCombiningRules(COMBINING_RULES rules)
{
    extractDetailsFromRules(rules);
}

/** Return whether or not this function uses a cutoff */
bool CLJFunction::hasCutoff() const
{
    return false;
}

/** Return the coulomb cutoff if this function uses one */
Length CLJFunction::coulombCutoff() const
{
    return Length( std::numeric_limits<double>::max() );
}

/** Return the LJ cutoff if this function uses one */
Length CLJFunction::ljCutoff() const
{
    return Length( std::numeric_limits<double>::max() );
}

/** Set the coulomb and LJ cutoff to 'distance', if this function
    has a cutoff */
void CLJFunction::setCutoff(Length distance)
{
    setCoulombCutoff(distance);
    setLJCutoff(distance);
}

/** Set the coulomb and LJ cutoff to the specified values, if this function
    has a cutoff */
void CLJFunction::setCutoff(Length coulomb_cutoff, Length lj_cutoff)
{
    setCoulombCutoff(coulomb_cutoff);
    setLJCutoff(lj_cutoff);
}

/** Set the coulomb cutoff if this function has a cutoff */
void CLJFunction::setCoulombCutoff(Length)
{
    throw SireError::unsupported( QObject::tr( "The CLJFunction \"%1\" does not "
                    "support use of a cutoff, so one cannot be set.")
                            .arg(this->toString()), CODELOC );
}

/** Set the LJ cutoff if this function has a cutoff */
void CLJFunction::setLJCutoff(Length)
{
    throw SireError::unsupported( QObject::tr( "The CLJFunction \"%1\" does not "
                    "support use of a cutoff, so one cannot be set.")
                            .arg(this->toString()), CODELOC );
}

/** Return the space represented by the function */
const Space& CLJFunction::space() const
{
    return spce.read();
}

/** Return whether or not the space of the function is periodic */
bool CLJFunction::isPeriodic() const
{
    return use_box;
}

/** Set the space used by the function */
void CLJFunction::setSpace(const Space &space)
{
    SpacePtr old_space = spce;

    try
    {
        spce = space;
        extractDetailsFromSpace();
    }
    catch(...)
    {
        spce = old_space;
        throw;
    }
}

/** Return whether or not arithmetic combining rules are used */
bool CLJFunction::usingArithmeticCombiningRules() const
{
    return use_arithmetic;
}

/** Return whether or not geometric combining rules are used */
bool CLJFunction::usingGeometricCombiningRules() const
{
    return not use_arithmetic;
}

/** Return whether or not this function supports calculating potentials on grids */
bool CLJFunction::supportsGridCalculation() const
{
    return false;
}

/** Dummy function that needs to be overridden to support grid calculations */
void CLJFunction::calcBoxGrid(const CLJAtoms &atoms, const GridInfo &gridinfo,
                              const Vector &box_dimensions,
                              const int start, const int end, float *gridpot) const
{
    for (int i=start; i<end; ++i)
    {
        gridpot[i] = 0;
    }
}

/** Dummy function that needs to be overridden to support grid calculations */
void CLJFunction::calcVacGrid(const CLJAtoms &atoms, const GridInfo &gridinfo,
                              const int start, const int end, float *gridpot) const
{
    for (int i=start; i<end; ++i)
    {
        gridpot[i] = 0;
    }
}

namespace SireMM
{
    namespace detail
    {
        class CLJGridCalculator
        {
        public:
            CLJGridCalculator()
                : atoms(0), grid_info(0), cljfunc(0), gridpot(0), use_box(false)
            {}

            CLJGridCalculator(const CLJAtoms &_atoms, const GridInfo &_grid_info,
                              const CLJFunction &_cljfunc, float *_gridpot)
                : atoms(&_atoms), grid_info(&_grid_info), cljfunc(&_cljfunc),
                  gridpot(_gridpot), use_box(_cljfunc.use_box)
            {}

            ~CLJGridCalculator()
            {}

            void operator()(const tbb::blocked_range<int> &range) const
            {
                if (use_box)
                {
                    cljfunc->calcBoxGrid(*atoms, *grid_info, cljfunc->box_dimensions,
                                         range.begin(), range.end(),
                                         &(gridpot[range.begin()]));
                }
                else
                {
                    cljfunc->calcVacGrid(*atoms, *grid_info,
                                         range.begin(), range.end(),
                                         &(gridpot[range.begin()]));
                }
            }

        private:
            const CLJAtoms* const atoms;
            const GridInfo* const grid_info;
            const CLJFunction* const cljfunc;
            float *gridpot;
            const bool use_box;
        };
    }
}

/** Return the potential on the described grid of the passed atoms using
    this function. This returns an empty grid if this function doesn't support
    grid calculations */
QVector<float> CLJFunction::calculate(const CLJAtoms &atoms, const GridInfo &gridinfo) const
{
    QVector<float> gridpot( gridinfo.nPoints(), 0.0 );

    if (this->supportsGridCalculation() and not gridinfo.isEmpty())
    {
        SireMM::detail::CLJGridCalculator calc(atoms, gridinfo, *this, gridpot.data());
        tbb::parallel_for(tbb::blocked_range<int>(0,gridinfo.nPoints(),4096), calc);
    }

    return gridpot;
}

/** Calculate the coulomb energy between all atoms in 'atoms' */
double CLJFunction::calcVacCoulombEnergyAri(const CLJAtoms &atoms) const
{
    double cnrg, ljnrg;
    calcVacEnergyAri(atoms, cnrg, ljnrg);
    return cnrg;
}

/** Calculate the coulomb energy between the atoms in 'atoms0' and the atoms in 'atoms1' */
double CLJFunction::calcVacCoulombEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                            float min_distance) const
{
    double cnrg, ljnrg;
    calcVacEnergyAri(atoms0, atoms1, cnrg, ljnrg, min_distance);
    return cnrg;
}

/** Calculate the coulomb energy between all atoms in 'atoms' */
double CLJFunction::calcVacCoulombEnergyGeo(const CLJAtoms &atoms) const
{
    double cnrg, ljnrg;
    calcVacEnergyGeo(atoms, cnrg, ljnrg);
    return cnrg;
}

/** Calculate the coulomb energy between the atoms in 'atoms0' and the atoms in 'atoms1' */
double CLJFunction::calcVacCoulombEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                            float min_distance) const
{
    double cnrg, ljnrg;
    calcVacEnergyGeo(atoms0, atoms1, cnrg, ljnrg, min_distance);
    return cnrg;
}

/** Calculate the LJ energy between all atoms in 'atoms' */
double CLJFunction::calcVacLJEnergyAri(const CLJAtoms &atoms) const
{
    double cnrg, ljnrg;
    calcVacEnergyAri(atoms, cnrg, ljnrg);
    return ljnrg;
}

/** Calculate the LJ energy between the atoms in 'atoms0' and the atoms in 'atoms1' */
double CLJFunction::calcVacLJEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                       float min_distance) const
{
    double cnrg, ljnrg;
    calcVacEnergyAri(atoms0, atoms1, cnrg, ljnrg, min_distance);
    return ljnrg;
}

/** Calculate the LJ energy between all atoms in 'atoms' */
double CLJFunction::calcVacLJEnergyGeo(const CLJAtoms &atoms) const
{
    double cnrg, ljnrg;
    calcVacEnergyGeo(atoms, cnrg, ljnrg);
    return ljnrg;
}

/** Calculate the LJ energy between the atoms in 'atoms0' and the atoms in 'atoms1' */
double CLJFunction::calcVacLJEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                       float min_distance) const
{
    double cnrg, ljnrg;
    calcVacEnergyGeo(atoms0, atoms1, cnrg, ljnrg, min_distance);
    return ljnrg;
}

/** Calculate the coulomb energy between all atoms in 'atoms' */
double CLJFunction::calcBoxCoulombEnergyAri(const CLJAtoms &atoms,
                                            const Vector &box_dimensions) const
{
    double cnrg, ljnrg;
    calcBoxEnergyAri(atoms, box_dimensions, cnrg, ljnrg);
    return cnrg;
}

/** Calculate the coulomb energy between the atoms in 'atoms0' and the atoms in 'atoms1' */
double CLJFunction::calcBoxCoulombEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                            const Vector &box_dimensions,
                                            float min_distance) const
{
    double cnrg, ljnrg;
    calcBoxEnergyAri(atoms0, atoms1, box_dimensions, cnrg, ljnrg, min_distance);
    return cnrg;
}

/** Calculate the coulomb energy between all atoms in 'atoms' */
double CLJFunction::calcBoxCoulombEnergyGeo(const CLJAtoms &atoms,
                                            const Vector &box_dimensions) const
{
    double cnrg, ljnrg;
    calcBoxEnergyGeo(atoms, box_dimensions, cnrg, ljnrg);
    return cnrg;
}

/** Calculate the coulomb energy between the atoms in 'atoms0' and the atoms in 'atoms1' */
double CLJFunction::calcBoxCoulombEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                            const Vector &box_dimensions,
                                            float min_distance) const
{
    double cnrg, ljnrg;
    calcBoxEnergyGeo(atoms0, atoms1, box_dimensions, cnrg, ljnrg, min_distance);
    return cnrg;
}

/** Calculate the LJ energy between all atoms in 'atoms' */
double CLJFunction::calcBoxLJEnergyAri(const CLJAtoms &atoms,
                                       const Vector &box_dimensions) const
{
    double cnrg, ljnrg;
    calcBoxEnergyAri(atoms, box_dimensions, cnrg, ljnrg);
    return ljnrg;
}

/** Calculate the LJ energy between the atoms in 'atoms0' and the atoms in 'atoms1' */
double CLJFunction::calcBoxLJEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                       const Vector &box_dimensions, float min_distance) const
{
    double cnrg, ljnrg;
    calcBoxEnergyAri(atoms0, atoms1, box_dimensions, cnrg, ljnrg, min_distance);
    return ljnrg;
}

/** Calculate the LJ energy between all atoms in 'atoms' */
double CLJFunction::calcBoxLJEnergyGeo(const CLJAtoms &atoms,
                                       const Vector &box_dimensions) const
{
    double cnrg, ljnrg;
    calcBoxEnergyGeo(atoms, box_dimensions, cnrg, ljnrg);
    return ljnrg;
}

/** Calculate the LJ energy between the atoms in 'atoms0' and the atoms in 'atoms1' */
double CLJFunction::calcBoxLJEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                       const Vector &box_dimensions,
                                       float min_distance) const
{
    double cnrg, ljnrg;
    calcBoxEnergyGeo(atoms0, atoms1, box_dimensions, cnrg, ljnrg, min_distance);
    return ljnrg;
}

/** Return the total energy between 'atoms', returning the coulomb part in 'cnrg'
    and the LJ part in 'ljnrg' */
void CLJFunction::operator()(const CLJAtoms &atoms,
                             double &cnrg, double &ljnrg) const
{
    if (atoms.isEmpty())
    {
        cnrg = 0;
        ljnrg = 0;
    }
    else
    {
        if (use_arithmetic)
        {
            if (use_box)
                this->calcBoxEnergyAri(atoms, box_dimensions, cnrg, ljnrg);
            else
                this->calcVacEnergyAri(atoms, cnrg, ljnrg);
        }
        else
        {
            if (use_box)
                this->calcBoxEnergyGeo(atoms, box_dimensions, cnrg, ljnrg);
            else
                this->calcVacEnergyGeo(atoms, cnrg, ljnrg);
        }
    }
}

/** Return the total energy between 'atoms0' and 'atoms1', returning the coulomb part in 'cnrg'
    and the LJ part in 'ljnrg' */
void CLJFunction::operator()(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                             double &cnrg, double &ljnrg, float min_distance) const
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
    {
        cnrg = 0;
        ljnrg = 0;
        return;
    }
    else if (atoms0.count() > atoms1.count())
    {
        if (use_arithmetic)
        {
            if (use_box)
                this->calcBoxEnergyAri(atoms1, atoms0, box_dimensions, cnrg, ljnrg, min_distance);
            else
                this->calcVacEnergyAri(atoms1, atoms0, cnrg, ljnrg, min_distance);
        }
        else
        {
            if (use_box)
                this->calcBoxEnergyGeo(atoms1, atoms0, box_dimensions, cnrg, ljnrg, min_distance);
            else
                this->calcVacEnergyGeo(atoms1, atoms0, cnrg, ljnrg, min_distance);
        }
    }
    else
    {
        if (use_arithmetic)
        {
            if (use_box)
                this->calcBoxEnergyAri(atoms0, atoms1, box_dimensions, cnrg, ljnrg, min_distance);
            else
                this->calcVacEnergyAri(atoms0, atoms1, cnrg, ljnrg, min_distance);
        }
        else
        {
            if (use_box)
                this->calcBoxEnergyGeo(atoms0, atoms1, box_dimensions, cnrg, ljnrg, min_distance);
            else
                this->calcVacEnergyGeo(atoms0, atoms1, cnrg, ljnrg, min_distance);
        }
    }
}

/** Return the total energy between 'atoms', returning the coulomb part in 'cnrg'
    and the LJ part in 'ljnrg' */
void CLJFunction::operator()(const CLJBoxes &atoms, double &cnrg, double &ljnrg) const
{
    cnrg = 0;
    ljnrg = 0;

    if (this->hasCutoff())
    {
        const float min_cutoff = qMax( this->coulombCutoff(), this->ljCutoff() );

        const CLJBoxes::Container &boxes = atoms.occupiedBoxes();

        for (CLJBoxes::const_iterator it0 = boxes.constBegin();
             it0 != boxes.constEnd();
             ++it0)
        {
            double icnrg(0), iljnrg(0);

            //calculate the self-energy of the box
            this->total(it0->read().atoms(), icnrg, iljnrg);

            cnrg += icnrg;
            ljnrg += iljnrg;

            //now calculate its interaction with all other boxes
            CLJBoxes::const_iterator it1 = it0;

            const CLJBoxIndex &idx0 = it0->read().index();

            for (++it1; it1 != boxes.constEnd(); ++it1)
            {
                const CLJBoxIndex &idx1 = it1->read().index();

                float mindist = atoms.getDistance(spce.read(), idx0, idx1);

                if (mindist < min_cutoff)
                {
                    this->total(it0->read().atoms(), it1->read().atoms(),
                                icnrg, iljnrg, mindist);

                    cnrg += icnrg;
                    ljnrg += iljnrg;
                }
            }
        }
    }
    else
    {
        const CLJBoxes::Container &boxes = atoms.occupiedBoxes();

        for (CLJBoxes::const_iterator it0 = boxes.constBegin();
             it0 != boxes.constEnd();
             ++it0)
        {
            double icnrg(0), iljnrg(0);

            //calculate the self-energy of the box
            this->total(it0->read().atoms(), icnrg, iljnrg);

            cnrg += icnrg;
            ljnrg += iljnrg;

            //now calculate its interaction with all other boxes
            CLJBoxes::const_iterator it1 = it0;

            for (++it1; it1 != boxes.constEnd(); ++it1)
            {
                this->total(it0->read().atoms(), it1->read().atoms(),
                            icnrg, iljnrg);

                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
}

/** Return the total energy between 'atoms0' and 'atoms1', returning the coulomb part in 'cnrg'
    and the LJ part in 'ljnrg' */
void CLJFunction::operator()(const CLJBoxes &atoms0, const CLJBoxes &atoms1,
                             double &cnrg, double &ljnrg) const
{
    cnrg = 0;
    ljnrg = 0;

    if (this->hasCutoff() and atoms0.length() == atoms1.length())
    {
        const CLJBoxes::Container &boxes0 = atoms0.occupiedBoxes();
        const CLJBoxes::Container &boxes1 = atoms1.occupiedBoxes();

        const float min_cutoff = qMax(this->coulombCutoff(), this->ljCutoff());

        for (CLJBoxes::const_iterator it0 = boxes0.constBegin();
             it0 != boxes0.constEnd();
             ++it0)
        {
            const CLJBoxIndex &idx0 = it0->read().index();

            for (CLJBoxes::const_iterator it1 = boxes1.constBegin();
                 it1 != boxes1.constEnd();
                 ++it1)
            {
                double icnrg(0), iljnrg(0);
                const CLJBoxIndex &idx1 = it1->read().index();

                const float mindist = atoms0.getDistance(spce.read(), idx0, idx1);

                if (mindist < min_cutoff)
                {
                    this->operator()(it0->read().atoms(), it1->read().atoms(),
                                     icnrg, iljnrg, mindist);
                }

                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
    else
    {
        const CLJBoxes::Container &boxes0 = atoms0.occupiedBoxes();
        const CLJBoxes::Container &boxes1 = atoms1.occupiedBoxes();

        for (CLJBoxes::const_iterator it0 = boxes0.constBegin();
             it0 != boxes0.constEnd();
             ++it0)
        {
            for (CLJBoxes::const_iterator it1 = boxes1.constBegin();
                 it1 != boxes1.constEnd();
                 ++it1)
            {
                double icnrg(0), iljnrg(0);

                this->operator()(it0->read().atoms(), it1->read().atoms(),
                                 icnrg, iljnrg);

                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
}

/** Return the total energy between 'atoms0' and 'atoms1', returning the coulomb part in 'cnrg'
    and the LJ part in 'ljnrg' */
void CLJFunction::operator()(const CLJAtoms &atoms0, const CLJBoxes &atoms1,
                             double &cnrg, double &ljnrg) const
{
    cnrg = 0;
    ljnrg = 0;

    if (this->hasCutoff())
    {
        this->operator()( CLJBoxes(atoms0, atoms1.length()), atoms1, cnrg, ljnrg );
    }
    else
    {
        const CLJBoxes::Container &boxes1 = atoms1.occupiedBoxes();

        for (CLJBoxes::const_iterator it1 = boxes1.constBegin();
             it1 != boxes1.constEnd();
             ++it1)
        {
            double icnrg(0), iljnrg(0);

            this->operator()(atoms0, it1->read().atoms(),
                             icnrg, iljnrg);

            cnrg += icnrg;
            ljnrg += iljnrg;
        }
    }
}

/** Return the total energy between 'atoms', returning the coulomb part in 'cnrg'
    and the LJ part in 'ljnrg' */
void CLJFunction::total(const CLJAtoms &atoms, double &cnrg, double &ljnrg) const
{
    this->operator()(atoms, cnrg, ljnrg);
}

/** Return the total energy between 'atoms0' and 'atoms1', returning the coulomb part in 'cnrg'
    and the LJ part in 'ljnrg' */
void CLJFunction::total(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                        double &cnrg, double &ljnrg, float min_distance) const
{
    this->operator()(atoms0, atoms1, cnrg, ljnrg, min_distance);
}

/** Return the total energy between 'atoms', returning the coulomb part in 'cnrg'
    and the LJ part in 'ljnrg' */
void CLJFunction::total(const CLJBoxes &atoms, double &cnrg, double &ljnrg) const
{
    this->operator()(atoms, cnrg, ljnrg);
}

/** Return the total energy between 'atoms0' and 'atoms1', returning the coulomb part in 'cnrg'
    and the LJ part in 'ljnrg' */
void CLJFunction::total(const CLJBoxes &atoms0, const CLJBoxes &atoms1,
                        double &cnrg, double &ljnrg) const
{
    this->operator()(atoms0, atoms1, cnrg, ljnrg);
}

/** Return the total energy between 'atoms', returning the coulomb part as the first
    element of the tuple and the LJ part as the second */
boost::tuple<double,double> CLJFunction::calculate(const CLJAtoms &atoms) const
{
    double cnrg, ljnrg;
    this->operator()(atoms, cnrg, ljnrg);
    return boost::tuple<double,double>(cnrg, ljnrg);
}

/** Return the total energy between 'atoms0' and 'atoms1', returning the coulomb part as the first
    element of the tuple and the LJ part as the second */
boost::tuple<double,double> CLJFunction::calculate(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                                   float min_distance) const
{
    double cnrg, ljnrg;
    this->operator()(atoms0, atoms1, cnrg, ljnrg, min_distance);
    return boost::tuple<double,double>(cnrg, ljnrg);
}

/** Return the total energy between 'atoms', returning the coulomb part as the first
    element of the tuple and the LJ part as the second */
boost::tuple<double,double> CLJFunction::calculate(const CLJBoxes &atoms) const
{
    double cnrg, ljnrg;
    this->operator()(atoms, cnrg, ljnrg);
    return boost::tuple<double,double>(cnrg, ljnrg);
}

/** Return the total energy between 'atoms0' and 'atoms1', returning the coulomb part as the first
    element of the tuple and the LJ part as the second */
boost::tuple<double,double> CLJFunction::calculate(const CLJBoxes &atoms0,
                                                   const CLJBoxes &atoms1) const
{
    double cnrg, ljnrg;
    this->operator()(atoms0, atoms1, cnrg, ljnrg);
    return boost::tuple<double,double>(cnrg, ljnrg);
}

/** Return the total energy between 'atoms0' and 'atoms1', returning the coulomb part as the first
    element of the tuple and the LJ part as the second */
boost::tuple<double,double> CLJFunction::calculate(const CLJAtoms &atoms0,
                                                   const CLJBoxes &atoms1) const
{
    double cnrg, ljnrg;
    this->operator()(atoms0, atoms1, cnrg, ljnrg);
    return boost::tuple<double,double>(cnrg, ljnrg);
}

/** Return the coulomb energy between the atoms in 'atoms' */
double CLJFunction::coulomb(const CLJAtoms &atoms) const
{
    if (atoms.isEmpty())
    {
        return 0;
    }
    else
    {
        if (use_arithmetic)
        {
            if (use_box)
                return this->calcBoxCoulombEnergyAri(atoms, box_dimensions);
            else
                return this->calcVacCoulombEnergyAri(atoms);
        }
        else
        {
            if (use_box)
                return this->calcBoxCoulombEnergyGeo(atoms, box_dimensions);
            else
                return this->calcVacCoulombEnergyGeo(atoms);
        }
    }
}

/** Return the coulomb energy between the atoms in 'atoms' */
double CLJFunction::coulomb(const CLJBoxes &atoms) const
{
    return this->calculate(atoms).get<0>();
}

/** Return the coulomb energy between the atoms in 'atoms0' and in 'atoms1' */
double CLJFunction::coulomb(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                            float min_distance) const
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
    {
        return 0;
    }
    else if (atoms0.count() > atoms1.count())
    {
        if (use_arithmetic)
        {
            if (use_box)
                return this->calcBoxCoulombEnergyAri(atoms1, atoms0, box_dimensions, min_distance);
            else
                return this->calcVacCoulombEnergyAri(atoms1, atoms0, min_distance);
        }
        else
        {
            if (use_box)
                return this->calcBoxCoulombEnergyGeo(atoms1, atoms0, box_dimensions, min_distance);
            else
                return this->calcVacCoulombEnergyGeo(atoms1, atoms0, min_distance);
        }
    }
    else
    {
        if (use_arithmetic)
        {
            if (use_box)
                return this->calcBoxCoulombEnergyAri(atoms0, atoms1, box_dimensions, min_distance);
            else
                return this->calcVacCoulombEnergyAri(atoms0, atoms1, min_distance);
        }
        else
        {
            if (use_box)
                return this->calcBoxCoulombEnergyGeo(atoms0, atoms1, box_dimensions, min_distance);
            else
                return this->calcVacCoulombEnergyGeo(atoms0, atoms1, min_distance);
        }
    }
}

/** Return the coulomb energy between the atoms in 'atoms0' and in 'atoms1' */
double CLJFunction::coulomb(const CLJBoxes &atoms0, const CLJBoxes &atoms1) const
{
    return this->calculate(atoms0, atoms1).get<0>();
}

/** Return the LJ energy between the atoms in 'atoms' */
double CLJFunction::lj(const CLJAtoms &atoms) const
{
    if (atoms.isEmpty())
    {
        return 0;
    }
    else
    {
        if (use_arithmetic)
        {
            if (use_box)
                return this->calcBoxLJEnergyAri(atoms, box_dimensions);
            else
                return this->calcVacLJEnergyAri(atoms);
        }
        else
        {
            if (use_box)
                return this->calcBoxLJEnergyGeo(atoms, box_dimensions);
            else
                return this->calcVacLJEnergyGeo(atoms);
        }
    }
}

/** Return the LJ energy between the atoms in 'atoms' */
double CLJFunction::lj(const CLJBoxes &atoms) const
{
    return this->calculate(atoms).get<1>();
}

/** Return the LJ energy between the atoms in 'atoms0' and in 'atoms1' */
double CLJFunction::lj(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                       float min_distance) const
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
    {
        return 0;
    }
    else if (atoms0.count() > atoms1.count())
    {
        if (use_arithmetic)
        {
            if (use_box)
                return this->calcBoxLJEnergyAri(atoms1, atoms0, box_dimensions, min_distance);
            else
                return this->calcVacLJEnergyAri(atoms1, atoms0, min_distance);
        }
        else
        {
            if (use_box)
                return this->calcBoxLJEnergyGeo(atoms1, atoms0, box_dimensions, min_distance);
            else
                return this->calcVacLJEnergyGeo(atoms1, atoms0, min_distance);
        }
    }
    else
    {
        if (use_arithmetic)
        {
            if (use_box)
                return this->calcBoxLJEnergyAri(atoms0, atoms1, box_dimensions, min_distance);
            else
                return this->calcVacLJEnergyAri(atoms0, atoms1, min_distance);
        }
        else
        {
            if (use_box)
                return this->calcBoxLJEnergyGeo(atoms0, atoms1, box_dimensions, min_distance);
            else
                return this->calcVacLJEnergyGeo(atoms0, atoms1, min_distance);
        }
    }
}

/** Return the LJ energy between the atoms in 'atoms0' and in 'atoms1' */
double CLJFunction::lj(const CLJBoxes &atoms0, const CLJBoxes &atoms1) const
{
    return this->calculate(atoms0, atoms1).get<1>();
}

tuple< QVector<double>,QVector<double> >
CLJFunction::multiCalculate(const QVector<CLJFunctionPtr> &funcs, const CLJAtoms &atoms)
{
    if (funcs.isEmpty())
        return tuple< QVector<double>,QVector<double> >();

    QVector<double> cnrgs(funcs.count()), ljnrgs(funcs.count());

    for (int i=0; i<funcs.count(); ++i)
    {
        tuple<double,double> nrgs = funcs.constData()[i].read().calculate(atoms);
        cnrgs[i] = nrgs.get<0>();
        ljnrgs[i] = nrgs.get<1>();
    }

    return tuple< QVector<double>,QVector<double> >(cnrgs, ljnrgs);
}

tuple< QVector<double>,QVector<double> >
CLJFunction::multiCalculate(const QVector<CLJFunctionPtr> &funcs,
                            const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                            float min_distance)
{
    if (funcs.isEmpty())
        return tuple< QVector<double>,QVector<double> >();

    QVector<double> cnrgs(funcs.count()), ljnrgs(funcs.count());

    for (int i=0; i<funcs.count(); ++i)
    {
        tuple<double,double> nrgs = funcs.constData()[i].read()
                                                        .calculate(atoms0, atoms1, min_distance);
        cnrgs[i] = nrgs.get<0>();
        ljnrgs[i] = nrgs.get<1>();
    }

    return tuple< QVector<double>,QVector<double> >(cnrgs, ljnrgs);
}

tuple< QVector<double>,QVector<double> >
CLJFunction::multiCalculate(const QVector<CLJFunctionPtr> &funcs, const CLJBoxes &atoms)
{
    if (funcs.isEmpty())
        return tuple< QVector<double>,QVector<double> >();

    QVector<double> cnrgs(funcs.count()), ljnrgs(funcs.count());

    for (int i=0; i<funcs.count(); ++i)
    {
        tuple<double,double> nrgs = funcs.constData()[i].read().calculate(atoms);
        cnrgs[i] = nrgs.get<0>();
        ljnrgs[i] = nrgs.get<1>();
    }

    return tuple< QVector<double>,QVector<double> >(cnrgs, ljnrgs);
}

tuple< QVector<double>,QVector<double> >
CLJFunction::multiCalculate(const QVector<CLJFunctionPtr> &funcs,
                            const CLJBoxes &atoms0, const CLJBoxes &atoms1)
{
    if (funcs.isEmpty())
        return tuple< QVector<double>,QVector<double> >();

    QVector<double> cnrgs(funcs.count()), ljnrgs(funcs.count());

    for (int i=0; i<funcs.count(); ++i)
    {
        tuple<double,double> nrgs = funcs.constData()[i].read().calculate(atoms0, atoms1);
        cnrgs[i] = nrgs.get<0>();
        ljnrgs[i] = nrgs.get<1>();
    }

    return tuple< QVector<double>,QVector<double> >(cnrgs, ljnrgs);
}

tuple< QVector<double>,QVector<double> >
CLJFunction::multiCalculate(const QVector<CLJFunctionPtr> &funcs,
                            const CLJAtoms &atoms0, const CLJBoxes &atoms1)
{
    if (funcs.isEmpty())
        return tuple< QVector<double>,QVector<double> >();

    QVector<double> cnrgs(funcs.count()), ljnrgs(funcs.count());

    for (int i=0; i<funcs.count(); ++i)
    {
        tuple<double,double> nrgs = funcs.constData()[i].read().calculate(atoms0, atoms1);
        cnrgs[i] = nrgs.get<0>();
        ljnrgs[i] = nrgs.get<1>();
    }

    return tuple< QVector<double>,QVector<double> >(cnrgs, ljnrgs);
}

/////////
///////// Implementation of NulCLJFunction
/////////

static const RegisterMetaType<NullCLJFunction> r_nullfunc;

QDataStream &operator<<(QDataStream &ds, const NullCLJFunction &nullfunc)
{
    writeHeader(ds, r_nullfunc, 1);

    ds << static_cast<const CLJFunction&>(nullfunc);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, NullCLJFunction &nullfunc)
{
    VersionID v = readHeader(ds, r_nullfunc);

    if (v == 1)
    {
        ds >> static_cast<CLJFunction&>(nullfunc);
    }
    else
        throw version_error(v, "1", r_nullfunc, CODELOC);

    return ds;
}

/** Constructor */
NullCLJFunction::NullCLJFunction() : ConcreteProperty<NullCLJFunction,CLJFunction>()
{}

/** Copy constructor */
NullCLJFunction::NullCLJFunction(const NullCLJFunction &other)
                : ConcreteProperty<NullCLJFunction,CLJFunction>(other)
{}

/** Destructor */
NullCLJFunction::~NullCLJFunction()
{}

/** Copy assignment operator */
NullCLJFunction& NullCLJFunction::operator=(const NullCLJFunction &other)
{
    CLJFunction::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullCLJFunction::operator==(const NullCLJFunction &other) const
{
    return CLJFunction::operator==(other);
}

/** Comparison operator */
bool NullCLJFunction::operator!=(const NullCLJFunction &other) const
{
    return not operator==(other);
}

const char* NullCLJFunction::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullCLJFunction>() );
}

const char* NullCLJFunction::what() const
{
    return NullCLJFunction::typeName();
}

void NullCLJFunction::calcVacEnergyAri(const CLJAtoms &atoms,
                                      double &cnrg, double &ljnrg) const
{
    cnrg = 0;
    ljnrg = 0;
}

void NullCLJFunction::calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                       double &cnrg, double &ljnrg, float min_distance) const
{
    cnrg = 0;
    ljnrg = 0;
}

void NullCLJFunction::calcVacEnergyGeo(const CLJAtoms &atoms,
                                       double &cnrg, double &ljnrg) const
{
    cnrg = 0;
    ljnrg = 0;
}

void NullCLJFunction::calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                       double &cnrg, double &ljnrg, float min_distance) const
{
    cnrg = 0;
    ljnrg = 0;
}

void NullCLJFunction::calcBoxEnergyAri(const CLJAtoms &atoms, const Vector &box,
                                       double &cnrg, double &ljnrg) const
{
    cnrg = 0;
    ljnrg = 0;
}

void NullCLJFunction::calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                       const Vector &box, double &cnrg, double &ljnrg,
                                       float min_distance) const
{
    cnrg = 0;
    ljnrg = 0;
}

void NullCLJFunction::calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box,
                                       double &cnrg, double &ljnrg) const
{
    cnrg = 0;
    ljnrg = 0;
}

void NullCLJFunction::calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                       const Vector &box, double &cnrg, double &ljnrg,
                                       float min_distance) const
{
    cnrg = 0;
    ljnrg = 0;
}

/////////
///////// Implementation of CLJCutoffFunction
/////////

static const RegisterMetaType<CLJCutoffFunction> r_cutoff(MAGIC_ONLY,
                                                          CLJCutoffFunction::typeName());

QDataStream &operator<<(QDataStream &ds, const CLJCutoffFunction &func)
{
    writeHeader(ds, r_cutoff, 1);

    ds << func.coul_cutoff << func.lj_cutoff
       << static_cast<const CLJFunction&>(func);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJCutoffFunction &func)
{
    VersionID v = readHeader(ds, r_cutoff);

    if (v == 1)
    {
        ds >> func.coul_cutoff >> func.lj_cutoff
           >> static_cast<CLJFunction&>(func);
    }
    else
        throw version_error(v, "1", r_cutoff, CODELOC);

    return ds;
}

/** Default constructor - uses a cutoff of 10 angstrom */
CLJCutoffFunction::CLJCutoffFunction()
                  : CLJFunction(), coul_cutoff(10), lj_cutoff(10)
{}

void CLJCutoffFunction::pvt_setCutoff(Length coulomb, Length lj)
{
    if (coulomb.value() < 0)
    {
        coul_cutoff = 0;
    }
    else if (coulomb.value() > std::numeric_limits<float>::max())
    {
        coul_cutoff = std::numeric_limits<float>::max();
    }
    else
    {
        coul_cutoff = coulomb.value();
    }

    if (lj.value() < 0)
    {
        lj_cutoff = 0;
    }
    else if (lj.value() > std::numeric_limits<float>::max())
    {
        lj_cutoff = std::numeric_limits<float>::max();
    }
    else
    {
        lj_cutoff = lj.value();
    }
}

/** Construct, specifying the cutoff */
CLJCutoffFunction::CLJCutoffFunction(Length cutoff)
                  : CLJFunction(), coul_cutoff(10), lj_cutoff(10)
{
    pvt_setCutoff(cutoff, cutoff);
}

/** Construct, specifying the coulomb and LJ cutoffs */
CLJCutoffFunction::CLJCutoffFunction(Length coulomb, Length lj)
                  : CLJFunction(), coul_cutoff(10), lj_cutoff(10)
{
    pvt_setCutoff(coulomb, lj);
}

CLJCutoffFunction::CLJCutoffFunction(const Space &space, Length cutoff)
                  : CLJFunction(space), coul_cutoff(10), lj_cutoff(10)
{
    pvt_setCutoff(cutoff, cutoff);
}

CLJCutoffFunction::CLJCutoffFunction(const Space &space, Length coulomb, Length lj)
                  : CLJFunction(space), coul_cutoff(10), lj_cutoff(10)
{
    pvt_setCutoff(coulomb, lj);
}

CLJCutoffFunction::CLJCutoffFunction(Length cutoff, COMBINING_RULES combining_rules)
                  : CLJFunction(combining_rules), coul_cutoff(10), lj_cutoff(10)
{
    pvt_setCutoff(cutoff, cutoff);
}

CLJCutoffFunction::CLJCutoffFunction(Length coulomb, Length lj, COMBINING_RULES combining_rules)
                  : CLJFunction(combining_rules), coul_cutoff(10), lj_cutoff(10)
{
    pvt_setCutoff(coulomb, lj);
}

CLJCutoffFunction::CLJCutoffFunction(const Space &space, COMBINING_RULES combining_rules)
                  : CLJFunction(space, combining_rules), coul_cutoff(10), lj_cutoff(10)
{}

CLJCutoffFunction::CLJCutoffFunction(const Space &space, Length cutoff,
                                     COMBINING_RULES combining_rules)
                  : CLJFunction(space, combining_rules), coul_cutoff(10), lj_cutoff(10)
{
    pvt_setCutoff(cutoff, cutoff);
}

CLJCutoffFunction::CLJCutoffFunction(const Space &space, Length coulomb, Length lj,
                                     COMBINING_RULES combining_rules)
                  : CLJFunction(space, combining_rules), coul_cutoff(10), lj_cutoff(10)
{
    pvt_setCutoff(coulomb, lj);
}

/** Copy constructor */
CLJCutoffFunction::CLJCutoffFunction(const CLJCutoffFunction &other)
                  : CLJFunction(other), coul_cutoff(other.coul_cutoff), lj_cutoff(other.lj_cutoff)
{}

/** Destructor */
CLJCutoffFunction::~CLJCutoffFunction()
{}

/** Copy assignment operator */
CLJCutoffFunction& CLJCutoffFunction::operator=(const CLJCutoffFunction &other)
{
    if (this != &other)
    {
        coul_cutoff = other.coul_cutoff;
        lj_cutoff = other.lj_cutoff;
        CLJFunction::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool CLJCutoffFunction::operator==(const CLJCutoffFunction &other) const
{
    return coul_cutoff == other.coul_cutoff and lj_cutoff == other.lj_cutoff and
           CLJFunction::operator==(other);
}

const char* CLJCutoffFunction::typeName()
{
    return "SireMM::CLJCutoffFunction";
}

QString CLJCutoffFunction::toString() const
{
    return QObject::tr("%1( coulombCutoff() = %2 A, ljCutoff() = %3 A, space() = %4 )")
                .arg(this->what())
                .arg(coul_cutoff).arg(lj_cutoff).arg(this->space().toString());
}

/** Return the properties that can be set in this function */
Properties CLJCutoffFunction::properties() const
{
    Properties props = CLJFunction::properties();

    props.setProperty( "coulombCutoff", GeneralUnitProperty(Length(coul_cutoff)) );
    props.setProperty( "ljCutoff", GeneralUnitProperty(Length(lj_cutoff)) );
    props.setProperty( "switchingFunction", HarmonicSwitchingFunction( Length(coul_cutoff),
                                                                       Length(lj_cutoff) ) );

    return props;
}

/** Set the property with name 'name' to value 'value' */
CLJFunctionPtr CLJCutoffFunction::setProperty(const QString &name, const Property &value) const
{
    CLJFunctionPtr ret(*this);

    if (name == "switchingFunction")
    {
        const SwitchingFunction &switchfunc = value.asA<SwitchingFunction>();
        ret.edit().setCoulombCutoff( Length(switchfunc.electrostaticCutoffDistance()) );
        ret.edit().setLJCutoff( Length(switchfunc.vdwCutoffDistance()) );
    }
    else if (name == "coulombCutoff")
    {
        if (value.isA<LengthProperty>())
            ret.edit().setCoulombCutoff( GeneralUnitProperty(value.asA<LengthProperty>().value()) );
        else
            ret.edit().setCoulombCutoff( value.asA<GeneralUnitProperty>() );
    }
    else if (name == "ljCutoff")
    {
        if (value.isA<LengthProperty>())
            ret.edit().setLJCutoff( GeneralUnitProperty(value.asA<LengthProperty>().value()) );
        else
            ret.edit().setLJCutoff( value.asA<GeneralUnitProperty>() );
    }
    else
    {
        ret = CLJFunction::setProperty(name, value);
    }

    return ret;
}

/** Return the value of the property with name 'name' */
PropertyPtr CLJCutoffFunction::property(const QString &name) const
{
    if (name == "switchingFunction")
    {
        return HarmonicSwitchingFunction( Length(coul_cutoff), Length(lj_cutoff) );
    }
    else if (name == "coulombCutoff")
    {
        return GeneralUnitProperty(Length(coul_cutoff));
    }
    else if (name == "ljCutoff")
    {
        return GeneralUnitProperty(Length(lj_cutoff));
    }
    else
    {
        return CLJFunction::property(name);
    }
}

/** Return whether or not this function contains a property called 'name' */
bool CLJCutoffFunction::containsProperty(const QString &name) const
{
    return (name == "switchingFunction") or (name == "coulombCutoff") or
           (name == "ljCutoff") or CLJFunction::containsProperty(name);
}

/** Return whether or not this function has a cutoff */
bool CLJCutoffFunction::hasCutoff() const
{
    return true;
}

/** Return the coulomb cutoff distance */
Length CLJCutoffFunction::coulombCutoff() const
{
    return Length(coul_cutoff);
}

/** Return the LJ cutoff distance */
Length CLJCutoffFunction::ljCutoff() const
{
    return Length(lj_cutoff);
}

/** Set the coulomb and LJ cutoff distances to 'distance' */
void CLJCutoffFunction::setCutoff(Length distance)
{
    pvt_setCutoff(distance, distance);
}

/** Set the coulomb and LJ cutoff distances to the specified values */
void CLJCutoffFunction::setCutoff(Length coulomb, Length lj)
{
    pvt_setCutoff(coulomb, lj);
}

/** Set the coulomb cutoff to the specified distance */
void CLJCutoffFunction::setCoulombCutoff(Length distance)
{
    pvt_setCutoff(distance, ljCutoff());
}

/** Set the LJ cutoff to the specified distance */
void CLJCutoffFunction::setLJCutoff(Length distance)
{
    pvt_setCutoff(coulombCutoff(), distance);
}

/////////
///////// Implementation of CLJIntraFunction
/////////

static const RegisterMetaType<CLJIntraFunction> r_intra(MAGIC_ONLY, CLJIntraFunction::typeName());

QDataStream &operator<<(QDataStream &ds, const CLJIntraFunction &func)
{
    writeHeader(ds, r_intra, 1);

    SharedDataStream sds(ds);

    sds << func.cty
        << static_cast<const CLJCutoffFunction&>(func);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJIntraFunction &func)
{
    VersionID v = readHeader(ds, r_intra);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        func.bond_matrix.clear();
        func.cty = Connectivity();

        Connectivity connectivity;

        sds >> connectivity
            >> static_cast<CLJCutoffFunction&>(func);

        func.setConnectivity(connectivity);
    }
    else
        throw version_error(v, "1", r_intra, CODELOC);

    return ds;
}


CLJIntraFunction::CLJIntraFunction() : CLJCutoffFunction()
{}

CLJIntraFunction::CLJIntraFunction(Length cutoff) : CLJCutoffFunction(cutoff)
{}

CLJIntraFunction::CLJIntraFunction(Length coul_cutoff, Length lj_cutoff)
                 : CLJCutoffFunction(coul_cutoff, lj_cutoff)
{}

CLJIntraFunction::CLJIntraFunction(const Space &space, Length cutoff)
                 : CLJCutoffFunction(space, cutoff)
{}

CLJIntraFunction::CLJIntraFunction(const Space &space,
                                   Length coul_cutoff, Length lj_cutoff)
                 : CLJCutoffFunction(space, coul_cutoff, lj_cutoff)
{}

CLJIntraFunction::CLJIntraFunction(Length cutoff, COMBINING_RULES combining_rules)
                 : CLJCutoffFunction(cutoff, combining_rules)
{}

CLJIntraFunction::CLJIntraFunction(Length coul_cutoff, Length lj_cutoff,
                                   COMBINING_RULES combining_rules)
                 : CLJCutoffFunction(coul_cutoff, lj_cutoff, combining_rules)
{}

CLJIntraFunction::CLJIntraFunction(const Space &space, COMBINING_RULES combining_rules)
                 : CLJCutoffFunction(space, combining_rules)
{}

CLJIntraFunction::CLJIntraFunction(const Space &space, Length cutoff,
                                   COMBINING_RULES combining_rules)
                 : CLJCutoffFunction(space, cutoff, combining_rules)
{}

CLJIntraFunction::CLJIntraFunction(const Space &space, Length coul_cutoff,
                                   Length lj_cutoff,
                                   COMBINING_RULES combining_rules)
                 : CLJCutoffFunction(space, coul_cutoff, lj_cutoff, combining_rules)
{}

/** Copy constructor */
CLJIntraFunction::CLJIntraFunction(const CLJIntraFunction &other)
                 : CLJCutoffFunction(other), cty(other.cty), bond_matrix(other.bond_matrix)
{}

/** Destructor */
CLJIntraFunction::~CLJIntraFunction()
{}

const char* CLJIntraFunction::typeName()
{
    return "SireMM::CLJIntraFunction";
}

/** Copy assignment operator */
CLJIntraFunction& CLJIntraFunction::operator=(const CLJIntraFunction &other)
{
    if (this != &other)
    {
        cty = other.cty;
        bond_matrix = other.bond_matrix;
        CLJCutoffFunction::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool CLJIntraFunction::operator==(const CLJIntraFunction &other) const
{
    return cty == other.cty and CLJCutoffFunction::operator==(other);
}

/** Return the properties that can be set in this function */
Properties CLJIntraFunction::properties() const
{
    Properties props = CLJCutoffFunction::properties();

    props.setProperty( "connectivity", cty );

    return props;
}

/** Set the property with name 'name' to value 'value' */
CLJFunctionPtr CLJIntraFunction::setProperty(const QString &name, const Property &value) const
{
    CLJFunctionPtr ret(*this);

    if (name == "connectivity")
    {
        ret.edit().asA<CLJIntraFunction>().setConnectivity( value.asA<Connectivity>() );
    }
    else
    {
        ret = CLJCutoffFunction::setProperty(name, value);
    }

    return ret;
}

/** Return the value of the property with name 'name' */
PropertyPtr CLJIntraFunction::property(const QString &name) const
{
    if (name == "connectivity")
    {
        return cty;
    }
    else
    {
        return CLJCutoffFunction::property(name);
    }
}

/** Return whether or not this function contains a property called 'name' */
bool CLJIntraFunction::containsProperty(const QString &name) const
{
    return (name == "connectivity") or CLJCutoffFunction::containsProperty(name);
}

/** Return the connectivity used to find the non-bonded pairs */
const Connectivity& CLJIntraFunction::connectivity() const
{
    return cty;
}

/** Set the connectivity used to find the non-bonded pairs */
void CLJIntraFunction::setConnectivity(const Connectivity &c)
{
    if (cty != c)
    {
        cty = c;
        bond_matrix = cty.getBondMatrix(1,4);

        if (bond_matrix.count() > 0)
        {
            //we now need to pad this with zeroes, as we use AtomIdx(0) to mean a dummy atom
            bond_matrix.prepend( QVector<bool>(bond_matrix.count() + 1, false) );

            bond_matrix[0].squeeze();

            for (int i=1; i<bond_matrix.count(); ++i)
            {
                bond_matrix[i].prepend(false);
                bond_matrix[i].squeeze();
            }

            bond_matrix.squeeze();
        }
    }
}

/** Set the connectivity by copying the specified property from the passed molecule */
void CLJIntraFunction::setConnectivity(const MoleculeView &molecule, const PropertyMap &map)
{
    setConnectivity( molecule.data().property( map["connectivity"] ).asA<Connectivity>() );
}

/** Return whether or not there are no bonded pairs between the atoms in 'ids0' and 'ids1' */
bool CLJIntraFunction::isNotBonded(const QVector<MultiInt> &ids0,
                                   const QVector<MultiInt> &ids1) const
{
    const int nats0 = ids0.count();
    const int nats1 = ids1.count();

    const MultiInt *aid0 = ids0.constData();
    const MultiInt *aid1 = ids1.constData();

    for (int i=0; i<nats0; ++i)
    {
        const MultiInt &id0 = aid0[i];

        for (int ii=0; ii<MultiInt::count(); ++ii)
        {
            const bool *row = bond_matrix.constData()[ id0[ii] ].constData();

            for (int j=0; j<nats1; ++j)
            {
                const MultiInt &id1 = aid1[j];

                for (int jj=0; jj<MultiInt::count(); ++jj)
                {
                    if (row[id1[jj]])
                        return false;
                }
            }
        }

        /*for (int j=0; j<nats1; ++j)
        {
            const MultiInt &id1 = aid1[j];

            for (int ii=0; ii<MultiInt::count(); ++ii)
            {
                const QVector<bool> &row = bond_matrix.constData()[ id0[ii] ];

                for (int jj=0; jj<MultiInt::count(); ++jj)
                {
                    if (row.constData()[ id1[jj] ])
                    {
                        return false;
                    }
                }
            }
        }*/
    }

    return true;
}

/////////
///////// Implementation of CLJSoftFunction
/////////

static const RegisterMetaType<CLJSoftFunction> r_soft(MAGIC_ONLY, CLJSoftFunction::typeName());

QDataStream &operator<<(QDataStream &ds, const CLJSoftFunction &func)
{
    writeHeader(ds, r_soft, 1);

    ds << func.alpha_value << func.shift_delta << func.coulomb_power
       << static_cast<const CLJCutoffFunction&>(func);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJSoftFunction &func)
{
    VersionID v = readHeader(ds, r_soft);

    if (v == 1)
    {
        ds >> func.alpha_value >> func.shift_delta >> func.coulomb_power
           >> static_cast<CLJCutoffFunction&>(func);
    }
    else
        throw version_error(v, "1", r_soft, CODELOC);

    return ds;
}

/** Constructor */
CLJSoftFunction::CLJSoftFunction()
                : CLJCutoffFunction(), alpha_value(0), shift_delta(0), coulomb_power(1)
{}


/** Construct, specifying the cutoff */
CLJSoftFunction::CLJSoftFunction(Length cutoff)
                : CLJCutoffFunction(cutoff)
{}

/** Construct, specifying the coulomb and LJ cutoffs */
CLJSoftFunction::CLJSoftFunction(Length coulomb, Length lj)
                : CLJCutoffFunction(coulomb, lj)
{}

CLJSoftFunction::CLJSoftFunction(const Space &space, Length cutoff)
                : CLJCutoffFunction(space, cutoff)
{}

CLJSoftFunction::CLJSoftFunction(const Space &space, Length coulomb, Length lj)
                : CLJCutoffFunction(space, coulomb, lj)
{}

CLJSoftFunction::CLJSoftFunction(Length cutoff, COMBINING_RULES combining_rules)
                : CLJCutoffFunction(cutoff, combining_rules)
{}

CLJSoftFunction::CLJSoftFunction(Length coulomb, Length lj, COMBINING_RULES combining_rules)
                : CLJCutoffFunction(coulomb, lj, combining_rules)
{}

CLJSoftFunction::CLJSoftFunction(const Space &space, COMBINING_RULES combining_rules)
                : CLJCutoffFunction(space, combining_rules)
{}

CLJSoftFunction::CLJSoftFunction(const Space &space, Length cutoff,
                                 COMBINING_RULES combining_rules)
                : CLJCutoffFunction(space, cutoff, combining_rules)
{}

CLJSoftFunction::CLJSoftFunction(const Space &space, Length coulomb, Length lj,
                                 COMBINING_RULES combining_rules)
                : CLJCutoffFunction(space, coulomb, lj, combining_rules)
{}

/** Copy constructor */
CLJSoftFunction::CLJSoftFunction(const CLJSoftFunction &other)
                : CLJCutoffFunction(other), alpha_value(other.alpha_value),
                  shift_delta(other.shift_delta), coulomb_power(other.coulomb_power)
{}

/** Destructor */
CLJSoftFunction::~CLJSoftFunction()
{}

/** Copy assignment operator */
CLJSoftFunction& CLJSoftFunction::operator=(const CLJSoftFunction &other)
{
    if (this != &other)
    {
        alpha_value = other.alpha_value;
        shift_delta = other.shift_delta;
        coulomb_power = other.coulomb_power;
        CLJCutoffFunction::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool CLJSoftFunction::operator==(const CLJSoftFunction &other) const
{
    return alpha_value == other.alpha_value and
           shift_delta == other.shift_delta and
           coulomb_power == other.coulomb_power and
           CLJCutoffFunction::operator==(other);
}

const char* CLJSoftFunction::typeName()
{
    return "SireMM::CLJSoftFunction";
}

/** Return whether or not this is a softened function */
bool CLJSoftFunction::isSoftened() const
{
    return true;
}

/** Return the properties that can be set in this function */
Properties CLJSoftFunction::properties() const
{
    Properties props = CLJCutoffFunction::properties();

    props.setProperty( "alpha", NumberProperty(alpha_value) );
    props.setProperty( "shiftDelta", NumberProperty(shift_delta) );
    props.setProperty( "coulombPower", NumberProperty(coulomb_power) );

    return props;
}

/** Set the property with name 'name' to value 'value' */
CLJFunctionPtr CLJSoftFunction::setProperty(const QString &name, const Property &value) const
{
    CLJFunctionPtr ret(*this);

    if (name == "alpha")
    {
        ret.edit().asA<CLJSoftFunction>()
                  .setAlpha( value.asADouble() );
    }
    else if (name == "shiftDelta")
    {
        ret.edit().asA<CLJSoftFunction>()
                  .setShiftDelta( value.asADouble() );
    }
    else if (name == "coulombPower")
    {
        ret.edit().asA<CLJSoftFunction>()
                  .setCoulombPower( value.asADouble() );
    }
    else
    {
        ret = CLJCutoffFunction::setProperty(name, value);
    }

    return ret;
}

/** Return the value of the property with name 'name' */
PropertyPtr CLJSoftFunction::property(const QString &name) const
{
    if (name == "alpha")
    {
        return NumberProperty(alpha_value);
    }
    else if (name == "shiftDelta")
    {
        return NumberProperty(shift_delta);
    }
    else if (name == "coulombPower")
    {
        return NumberProperty(coulomb_power);
    }
    else
    {
        return CLJCutoffFunction::property(name);
    }
}

/** Return whether or not this function contains a property called 'name' */
bool CLJSoftFunction::containsProperty(const QString &name) const
{
    return (name == "alpha") or (name == "shiftDelta") or
           (name == "coulombPower") or CLJCutoffFunction::containsProperty(name);
}

/** Return the soft-core alpha value. A value of 0 is a completely hard
    potential, while increasing values of alpha will increasingly soften
    the potential */
float CLJSoftFunction::alpha() const
{
    return alpha_value;
}

/** Return the soft-core shift_delta parameter. This is used to soften
    the LJ interactions */
float CLJSoftFunction::shiftDelta() const
{
    return shift_delta;
}

/** Return the soft-core coulomb_power parameter. This is used to soften
    the electrostatic interactions */
float CLJSoftFunction::coulombPower() const
{
    return coulomb_power;
}

void CLJSoftFunction::pvt_set(float alpha, float shift, float power)
{
    if (alpha < 0)
        alpha = 0;

    if (shift < 0)
        shift = 0;

    if (alpha > 1)
    {
        if ( power != int(power) )
            throw SireError::incompatible_error( QObject::tr(
                    "You cannot have a soft-core function where alpha > 1 and "
                    "you have a non-integer coulomb power (alpha = %1, coulomb power = %2)")
                        .arg(alpha).arg(power), CODELOC );
    }

    alpha_value = alpha;
    shift_delta = shift;
    coulomb_power = power;
}

/** Set the soft-core alpha parameter */
void CLJSoftFunction::setAlpha(float alp)
{
    pvt_set(alp, shiftDelta(), coulombPower());
}

/** Set the soft-core shift delta parameter */
void CLJSoftFunction::setShiftDelta(float shift)
{
    pvt_set(alpha(), shift, coulombPower());
}

/** Set the soft-core coulomb power parameter */
void CLJSoftFunction::setCoulombPower(float power)
{
    pvt_set(alpha(), shiftDelta(), power);
}

/** Return (1-alpha)^(coulomb_power) */
float CLJSoftFunction::oneMinusAlphaToN() const
{
    if (coulombPower() == 0)
        return 1.0;
    else if (coulombPower() == 1)
        return (1.0 - alpha());
    else
        return std::pow( float(1.0 - alpha()), coulombPower() );
}

/** Return alpha * shift_delta */
float CLJSoftFunction::alphaTimesShiftDelta() const
{
    return alpha() * shiftDelta();
}

/////////
///////// Implementation of CLJSoftIntraFunction
/////////

static const RegisterMetaType<CLJSoftIntraFunction> r_softintra(MAGIC_ONLY,
                                                                CLJSoftIntraFunction::typeName());

QDataStream &operator<<(QDataStream &ds, const CLJSoftIntraFunction &func)
{
    writeHeader(ds, r_softintra, 1);

    ds << func.alpha_value << func.shift_delta << func.coulomb_power
       << static_cast<const CLJIntraFunction&>(func);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJSoftIntraFunction &func)
{
    VersionID v = readHeader(ds, r_softintra);

    if (v == 1)
    {
        ds >> func.alpha_value >> func.shift_delta >> func.coulomb_power
           >> static_cast<CLJIntraFunction&>(func);
    }
    else
        throw version_error(v, "1", r_softintra, CODELOC);

    return ds;
}

/** Constructor */
CLJSoftIntraFunction::CLJSoftIntraFunction()
                     : CLJIntraFunction(), alpha_value(0), shift_delta(0), coulomb_power(1)
{}

CLJSoftIntraFunction::CLJSoftIntraFunction(Length cutoff) : CLJIntraFunction(cutoff)
{}

CLJSoftIntraFunction::CLJSoftIntraFunction(Length coul_cutoff, Length lj_cutoff)
                     : CLJIntraFunction(coul_cutoff, lj_cutoff)
{}

CLJSoftIntraFunction::CLJSoftIntraFunction(const Space &space, Length cutoff)
                     : CLJIntraFunction(space, cutoff)
{}

CLJSoftIntraFunction::CLJSoftIntraFunction(const Space &space,
                                           Length coul_cutoff, Length lj_cutoff)
                     : CLJIntraFunction(space, coul_cutoff, lj_cutoff)
{}

CLJSoftIntraFunction::CLJSoftIntraFunction(Length cutoff, COMBINING_RULES combining_rules)
                     : CLJIntraFunction(cutoff, combining_rules)
{}

CLJSoftIntraFunction::CLJSoftIntraFunction(Length coul_cutoff, Length lj_cutoff,
                                           COMBINING_RULES combining_rules)
                     : CLJIntraFunction(coul_cutoff, lj_cutoff, combining_rules)
{}

CLJSoftIntraFunction::CLJSoftIntraFunction(const Space &space, COMBINING_RULES combining_rules)
                     : CLJIntraFunction(space, combining_rules)
{}

CLJSoftIntraFunction::CLJSoftIntraFunction(const Space &space, Length cutoff,
                                           COMBINING_RULES combining_rules)
                     : CLJIntraFunction(space, cutoff, combining_rules)
{}

CLJSoftIntraFunction::CLJSoftIntraFunction(const Space &space, Length coul_cutoff,
                                           Length lj_cutoff,
                                           COMBINING_RULES combining_rules)
                     : CLJIntraFunction(space, coul_cutoff, lj_cutoff, combining_rules)
{}

/** Copy constructor */
CLJSoftIntraFunction::CLJSoftIntraFunction(const CLJSoftIntraFunction &other)
                     : CLJIntraFunction(other), alpha_value(other.alpha_value),
                       shift_delta(other.shift_delta), coulomb_power(other.coulomb_power)
{}

/** Destructor */
CLJSoftIntraFunction::~CLJSoftIntraFunction()
{}

/** Copy assignment operator */
CLJSoftIntraFunction& CLJSoftIntraFunction::operator=(const CLJSoftIntraFunction &other)
{
    if (this != &other)
    {
        alpha_value = other.alpha_value;
        shift_delta = other.shift_delta;
        coulomb_power = other.coulomb_power;
        CLJIntraFunction::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool CLJSoftIntraFunction::operator==(const CLJSoftIntraFunction &other) const
{
    return alpha_value == other.alpha_value and
           shift_delta == other.shift_delta and
           coulomb_power == other.coulomb_power and
           CLJIntraFunction::operator==(other);
}

const char* CLJSoftIntraFunction::typeName()
{
    return "SireMM::CLJSoftIntraFunction";
}

/** Return whether or not this is a softened function */
bool CLJSoftIntraFunction::isSoftened() const
{
    return true;
}

/** Return the properties that can be set in this function */
Properties CLJSoftIntraFunction::properties() const
{
    Properties props = CLJIntraFunction::properties();

    props.setProperty( "alpha", NumberProperty(alpha_value) );
    props.setProperty( "shiftDelta", NumberProperty(shift_delta) );
    props.setProperty( "coulombPower", NumberProperty(coulomb_power) );

    return props;
}

/** Set the property with name 'name' to value 'value' */
CLJFunctionPtr CLJSoftIntraFunction::setProperty(const QString &name, const Property &value) const
{
    CLJFunctionPtr ret(*this);

    if (name == "alpha")
    {
        ret.edit().asA<CLJSoftIntraFunction>()
                  .setAlpha( value.asADouble() );
    }
    else if (name == "shiftDelta")
    {
        ret.edit().asA<CLJSoftIntraFunction>()
                  .setShiftDelta( value.asADouble() );
    }
    else if (name == "coulombPower")
    {
        ret.edit().asA<CLJSoftIntraFunction>()
                  .setCoulombPower( value.asADouble() );
    }
    else
    {
        ret = CLJIntraFunction::setProperty(name, value);
    }

    return ret;
}

/** Return the value of the property with name 'name' */
PropertyPtr CLJSoftIntraFunction::property(const QString &name) const
{
    if (name == "alpha")
    {
        return NumberProperty(alpha_value);
    }
    else if (name == "shiftDelta")
    {
        return NumberProperty(shift_delta);
    }
    else if (name == "coulombPower")
    {
        return NumberProperty(coulomb_power);
    }
    else
    {
        return CLJIntraFunction::property(name);
    }
}

/** Return whether or not this function contains a property called 'name' */
bool CLJSoftIntraFunction::containsProperty(const QString &name) const
{
    return (name == "alpha") or (name == "shiftDelta") or
           (name == "coulombPower") or CLJIntraFunction::containsProperty(name);
}

/** Return the soft-core alpha value. A value of 0 is a completely hard
    potential, while increasing values of alpha will increasingly soften
    the potential */
float CLJSoftIntraFunction::alpha() const
{
    return alpha_value;
}

/** Return the soft-core shift_delta parameter. This is used to soften
    the LJ interactions */
float CLJSoftIntraFunction::shiftDelta() const
{
    return shift_delta;
}

/** Return the soft-core coulomb_power parameter. This is used to soften
    the electrostatic interactions */
float CLJSoftIntraFunction::coulombPower() const
{
    return coulomb_power;
}

/** Return (1-alpha)^(coulomb_power) */
float CLJSoftIntraFunction::oneMinusAlphaToN() const
{
    if (coulombPower() == 0)
        return 1.0;
    else if (coulombPower() == 1)
        return (1.0 - alpha());
    else
        return std::pow( float(1.0 - alpha()), coulombPower() );
}

/** Return alpha * shift_delta */
float CLJSoftIntraFunction::alphaTimesShiftDelta() const
{
    return alpha() * shiftDelta();
}

void CLJSoftIntraFunction::pvt_set(float alpha, float shift, float power)
{
    if (alpha < 0)
        alpha = 0;

    if (shift < 0)
        shift = 0;

    if (alpha > 1)
    {
        if ( power != int(power) )
            throw SireError::incompatible_error( QObject::tr(
                    "You cannot have a soft-core function where alpha > 1 and "
                    "you have a non-integer coulomb power (alpha = %1, coulomb power = %2)")
                        .arg(alpha).arg(power), CODELOC );
    }

    alpha_value = alpha;
    shift_delta = shift;
    coulomb_power = power;
}

/** Set the soft-core alpha parameter */
void CLJSoftIntraFunction::setAlpha(float alp)
{
    pvt_set(alp, shiftDelta(), coulombPower());
}

/** Set the soft-core shift delta parameter */
void CLJSoftIntraFunction::setShiftDelta(float shift)
{
    pvt_set(alpha(), shift, coulombPower());
}

/** Set the soft-core coulomb power parameter */
void CLJSoftIntraFunction::setCoulombPower(float power)
{
    pvt_set(alpha(), shiftDelta(), power);
}
