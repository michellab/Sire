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

#include "testff.h"
#include "cljcalculator.h"
#include "cljshiftfunction.h"

#include "SireVol/cartesian.h"

#include "SireUnits/units.h"

#include <QDebug>
#include <QElapsedTimer>

#include "tbb/task_scheduler_init.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

#include "tostring.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireMaths;
using namespace SireUnits;
using namespace SireVol;

using boost::tuple;

TestFF::TestFF() : cljfunc(new CLJShiftFunction(15*angstrom, 15*angstrom))
{}

TestFF::TestFF(const TestFF &other)
       : atoms0(other.atoms0), atoms1(other.atoms1), cljfunc(other.cljfunc)
{}

TestFF::~TestFF()
{}

TestFF& TestFF::operator=(const TestFF &other)
{
    atoms0 = other.atoms0;
    atoms1 = other.atoms1;
    cljfunc = other.cljfunc;
    return *this;
}

void TestFF::add(const Molecules &molecules)
{
    atoms0 = CLJAtoms(molecules);
    cljboxes0 = CLJBoxes(atoms0);
}

void TestFF::addFixedAtoms(const Molecules &molecules)
{
    atoms1 = CLJAtoms(molecules);
    cljboxes1 = CLJBoxes(atoms1);
}

void TestFF::setCutoff(Length coul_cutoff, Length lj_cutoff)
{
    cljfunc.reset( new CLJShiftFunction(coul_cutoff, lj_cutoff) );
}

class CLJCalculator
{
public:
    CLJCalculator();
    CLJCalculator(const CLJFunction* const function,
                  const CLJBoxDistance* const distances,
                  const CLJBoxes::Container &cljboxes,
                  const float coulomb_cutoff, const float lenj_cutoff,
                  double *coulomb_energy, double *lj_energy)
        : func(function), dists(distances), boxes(&cljboxes),
          coul_nrg(coulomb_energy), lj_nrg(lj_energy),
          coul_cutoff(coulomb_cutoff), lj_cutoff(lenj_cutoff)
    {}
    
    ~CLJCalculator()
    {}
    
    void operator()(const tbb::blocked_range<int> &range) const
    {
        const CLJBoxDistance* ptr = dists + range.begin();
        const CLJBoxPtr* const b = boxes->constData();
    
        for (int i = range.begin(); i != range.end(); ++i)
        {
            if (ptr->box0() == ptr->box1())
            {
                (*func)(b[ptr->box0()].read().atoms(),
                        coul_nrg[i], lj_nrg[i]);
            }
            else
            {
                (*func)(b[ptr->box0()].read().atoms(),
                        b[ptr->box1()].read().atoms(),
                        coul_nrg[i], lj_nrg[i]);
            }
            
            ptr += 1;
        }
    }
    
private:
    const CLJFunction* const func;
    const CLJBoxDistance* const dists;
    const CLJBoxes::Container* boxes;
    
    double *coul_nrg;
    double *lj_nrg;
    
    const float coul_cutoff;
    const float lj_cutoff;
};

void TestFF::calculateEnergy()
{
    QElapsedTimer t;

    double cnrg;
    double ljnrg;

    qDebug() << "inter group energy";

    t.start();

    (*cljfunc)(atoms0, atoms1, cnrg, ljnrg);

    quint64 ns = t.nsecsElapsed();

    qDebug() << "TestFF" << (cnrg+ljnrg) << cnrg << ljnrg << "took" << (0.000001*ns) << "ms";

    qDebug() << "\nintra group energy";
    
    t.start();
    
    (*cljfunc)(atoms1, cnrg, ljnrg);
    
    ns = t.nsecsElapsed();
    
    qDebug() << "TestFF" << (cnrg+ljnrg) << cnrg << ljnrg << "took" << (0.000001*ns) << "ms";

    qDebug() << "\nUsing CLJBoxes to accelerate the calculation";
    
    Cartesian space;
    
    Length coul_cutoff = cljfunc->coulombCutoff();
    Length lj_cutoff = cljfunc->ljCutoff();
    
    QVector<CLJBoxDistance> dists;

    const CLJBoxes::Container &boxes0 = cljboxes0.occupiedBoxes();
    const CLJBoxes::Container &boxes1 = cljboxes1.occupiedBoxes();
    
    qDebug() << "inter energy";

    double icnrg, iljnrg;

    dists = CLJBoxes::getDistances(space, cljboxes0, cljboxes1, coul_cutoff);
    
    cnrg = 0;
    ljnrg = 0;
    
    t.start();
    
    const CLJBoxDistance *ptr = dists.constData();
    
    for (int i=0; i<dists.count(); ++i)
    {
        (*cljfunc)(boxes0[ptr->box0()].read().atoms(), boxes1[ptr->box1()].read().atoms(),
                   icnrg, iljnrg);
        
        cnrg += icnrg;
        ljnrg += iljnrg;
        
        ptr += 1;
    }
    
    ns = t.nsecsElapsed();
    
    qDebug() << "Boxed" << (cnrg+ljnrg) << cnrg << ljnrg;
    qDebug() << "Took" << (0.000001*ns) << "ms";

    qDebug() << "\nintra energy";
    
    dists = CLJBoxes::getDistances(space, cljboxes1, coul_cutoff);
    
    cnrg = 0;
    ljnrg = 0;
    
    ptr = dists.constData();
    
    t.start();
    
    for (int i=0; i<dists.count(); ++i)
    {
        if (ptr->box0() == ptr->box1())
        {
            (*cljfunc)(boxes1[ptr->box0()].read().atoms(), icnrg, iljnrg);
        }
        else
        {
            (*cljfunc)(boxes1[ptr->box0()].read().atoms(), boxes1[ptr->box1()].read().atoms(),
                       icnrg, iljnrg);
        }

        cnrg += icnrg;
        ljnrg += iljnrg;
        
        ptr += 1;
    }
    
    ns = t.nsecsElapsed();
    
    qDebug() << "Boxed" << (cnrg+ljnrg) << cnrg << ljnrg;
    qDebug() << "Took" << (0.000001*ns) << "ms";

    QVector<double> coul_nrgs(dists.count());
    QVector<double> lj_nrgs(dists.count());
    
    ::CLJCalculator calc(cljfunc.get(), dists.constData(), boxes1, coul_cutoff.value(),
                         lj_cutoff.value(), coul_nrgs.data(), lj_nrgs.data());
    
    t.start();
    tbb::parallel_for(tbb::blocked_range<int>(0,dists.count(),25), calc);
    
    cnrg = 0;
    ljnrg = 0;
    
    for (int i=0; i<dists.count(); ++i)
    {
        cnrg += coul_nrgs.constData()[i];
        ljnrg += lj_nrgs.constData()[i];
    }
    
    ns = t.nsecsElapsed();
    
    qDebug() << "\nParallel" << (cnrg+ljnrg) << cnrg << ljnrg;
    qDebug() << "Parallel version took" << (0.000001*ns) << "ms";
    
    qDebug() << "\nCLJCalculator (reproducible) version";
    SireMM::CLJCalculator calculator(true);
    t.start();
    tuple<double,double> nrgs = calculator.calculate(*cljfunc, cljboxes1);
    ns = t.nsecsElapsed();
    
    cnrg = nrgs.get<0>();
    ljnrg = nrgs.get<1>();
    
    qDebug() << "\nCLJCalculator" << (cnrg+ljnrg) << cnrg << ljnrg;
    qDebug() << "Took" << (0.000001*ns) << "ms";

    qDebug() << "\nCLJCalculator (non-reproducible) version";
    calculator = SireMM::CLJCalculator(false);
    t.start();
    nrgs = calculator.calculate(*cljfunc, cljboxes1);
    ns = t.nsecsElapsed();
    
    cnrg = nrgs.get<0>();
    ljnrg = nrgs.get<1>();
    
    qDebug() << "\nCLJCalculator" << (cnrg+ljnrg) << cnrg << ljnrg;
    qDebug() << "Took" << (0.000001*ns) << "ms";

    qDebug() << "\nCLJCalculator inter (reproducible) version";
    calculator = SireMM::CLJCalculator(true);
    t.start();
    nrgs = calculator.calculate(*cljfunc, cljboxes0, cljboxes1);
    ns = t.nsecsElapsed();
    
    cnrg = nrgs.get<0>();
    ljnrg = nrgs.get<1>();
    
    qDebug() << "\nCLJCalculator" << (cnrg+ljnrg) << cnrg << ljnrg;
    qDebug() << "Took" << (0.000001*ns) << "ms";

    qDebug() << "\nCLJCalculator inter (non-reproducible) version";
    calculator = SireMM::CLJCalculator(false);
    t.start();
    nrgs = calculator.calculate(*cljfunc, cljboxes0, cljboxes1);
    ns = t.nsecsElapsed();
    
    cnrg = nrgs.get<0>();
    ljnrg = nrgs.get<1>();
    
    qDebug() << "\nCLJCalculator" << (cnrg+ljnrg) << cnrg << ljnrg;
    qDebug() << "Took" << (0.000001*ns) << "ms";
}
