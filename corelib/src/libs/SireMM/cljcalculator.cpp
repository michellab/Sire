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

#include "cljcalculator.h"
#include "cljfunction.h"
#include "cljatoms.h"
#include "cljdelta.h"
#include "cljboxes.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

#include <QElapsedTimer>

using namespace SireMM;
using namespace SireBase;
using namespace SireStream;

using boost::tuple;

static const RegisterMetaType<CLJCalculator> r_cljcalc(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJCalculator &calc)
{
    writeHeader(ds, r_cljcalc, 1);
    ds << calc.reproducible_sum;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJCalculator &calc)
{
    VersionID v = readHeader(ds, r_cljcalc);
    
    if (v == 1)
    {
        ds >> calc.reproducible_sum;
    }
    else
        throw version_error(v, "1", r_cljcalc, CODELOC);
    
    return ds;
}

/** Constructor */
CLJCalculator::CLJCalculator(bool repro_sum) : reproducible_sum(repro_sum)
{}

/** Copy constructor */
CLJCalculator::CLJCalculator(const CLJCalculator &other) : reproducible_sum(other.reproducible_sum)
{}

/** Destructor */
CLJCalculator::~CLJCalculator()
{}

/** Copy assignment operator */
CLJCalculator& CLJCalculator::operator=(const CLJCalculator &other)
{
    reproducible_sum = other.reproducible_sum;
    return *this;
}

/** Comarison operator */
bool CLJCalculator::operator==(const CLJCalculator &other) const
{
    return reproducible_sum == other.reproducible_sum;
}

/** Comparison operator */
bool CLJCalculator::operator!=(const CLJCalculator&) const
{
    return false;
}

const char* CLJCalculator::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJCalculator>() );
}

const char* CLJCalculator::what() const
{
    return CLJCalculator::typeName();
}

QString CLJCalculator::toString() const
{
    return QObject::tr("CLJCalculator()");
}

namespace SireMM
{
    namespace detail
    {
        /** This is a private helper class that is used to calculate the
            coulomb and LJ energy in parallel using Intel TBB */
        class TotalWithCutoff
        {
        public:
            TotalWithCutoff() : func(0), dists(0), boxes(0), coul_cutoff(0), lj_cutoff(0)
            {}
        
            TotalWithCutoff(const CLJFunction* const function,
                            const CLJBoxDistance* const distances,
                            const CLJBoxes::Container &cljboxes,
                            const float coulomb_cutoff, const float lenj_cutoff,
                            double *coulomb_energy, double *lj_energy)
                : func(function), dists(distances), boxes(&cljboxes),
                  coul_nrg(coulomb_energy), lj_nrg(lj_energy),
                  coul_cutoff(coulomb_cutoff), lj_cutoff(lenj_cutoff)
            {}
            
            ~TotalWithCutoff()
            {}
            
            void operator()(const tbb::blocked_range<int> &range) const
            {
                const CLJBoxDistance* ptr = dists + range.begin();
                const CLJBoxPtr* const b = boxes->constData();
            
                for (int i = range.begin(); i != range.end(); ++i)
                {
                    if (ptr->box0() == ptr->box1())
                    {
                        func->total(b[ptr->box0()].read().atoms(),
                                    coul_nrg[i], lj_nrg[i]);
                    }
                    else if (ptr->distance() < coul_cutoff and ptr->distance() < lj_cutoff)
                    {
                        func->total(b[ptr->box0()].read().atoms(),
                                    b[ptr->box1()].read().atoms(),
                                    coul_nrg[i], lj_nrg[i], ptr->distance());
                    }
                    else if (ptr->distance() < coul_cutoff)
                    {
                        coul_nrg[i]
                            = func->coulomb(b[ptr->box0()].read().atoms(),
                                            b[ptr->box1()].read().atoms(),
                                            ptr->distance());
                        
                        lj_nrg[i] = 0;
                    }
                    else if (ptr->distance() < lj_cutoff)
                    {
                        lj_nrg[i]
                            = func->lj(b[ptr->box0()].read().atoms(),
                                       b[ptr->box1()].read().atoms(),
                                       ptr->distance());
                        
                        coul_nrg[i] = 0;
                    }
                    else
                    {
                        coul_nrg[i] = 0;
                        lj_nrg[i] = 0;
                    }
                    
                    ptr += 1;
                }
            }
            
        private:
            const CLJFunction* const func;
            const CLJBoxDistance* const dists;
            const CLJBoxes::Container* const boxes;
            
            double *coul_nrg;
            double *lj_nrg;
            
            const float coul_cutoff;
            const float lj_cutoff;
        };

        /** This is a private helper class that is used to calculate the
            coulomb and LJ energy in parallel using Intel TBB */
        class TotalWithoutCutoff
        {
        public:
            TotalWithoutCutoff() : func(0), dists(0), boxes(0)
            {}
        
            TotalWithoutCutoff(const CLJFunction* const function,
                               const CLJBoxDistance* const distances,
                               const CLJBoxes::Container &cljboxes,
                               double *coulomb_energy, double *lj_energy)
                : func(function), dists(distances), boxes(&cljboxes),
                  coul_nrg(coulomb_energy), lj_nrg(lj_energy)
            {}
            
            ~TotalWithoutCutoff()
            {}
            
            void operator()(const tbb::blocked_range<int> &range) const
            {
                const CLJBoxDistance* ptr = dists + range.begin();
                const CLJBoxPtr* const b = boxes->constData();
            
                for (int i = range.begin(); i != range.end(); ++i)
                {
                    if (ptr->box0() == ptr->box1())
                    {
                        func->total(b[ptr->box0()].read().atoms(),
                                    coul_nrg[i], lj_nrg[i]);
                    }
                    else
                    {
                        func->total(b[ptr->box0()].read().atoms(),
                                    b[ptr->box1()].read().atoms(),
                                    coul_nrg[i], lj_nrg[i], ptr->distance());
                    }
                    
                    ptr += 1;
                }
            }
            
        private:
            const CLJFunction* const func;
            const CLJBoxDistance* const dists;
            const CLJBoxes::Container* const boxes;
            
            double *coul_nrg;
            double *lj_nrg;
        };

        /** This is a private helper class that is used to calculate the
            coulomb and LJ energy in parallel using Intel TBB */
        class TotalWithCutoff2
        {
        public:
            TotalWithCutoff2() : func(0), dists(0), boxes0(0), boxes1(0),
                                 coul_cutoff(0), lj_cutoff(0)
            {}
        
            TotalWithCutoff2(const CLJFunction* const function,
                             const CLJBoxDistance* const distances,
                             const CLJBoxes::Container &cljboxes0,
                             const CLJBoxes::Container &cljboxes1,
                             const float coulomb_cutoff, const float lenj_cutoff,
                             double *coulomb_energy, double *lj_energy)
                : func(function), dists(distances), boxes0(&cljboxes0),
                  boxes1(&cljboxes1), coul_nrg(coulomb_energy), lj_nrg(lj_energy),
                  coul_cutoff(coulomb_cutoff), lj_cutoff(lenj_cutoff)
            {}
            
            ~TotalWithCutoff2()
            {}
            
            void operator()(const tbb::blocked_range<int> &range) const
            {
                const CLJBoxDistance* ptr = dists + range.begin();
                const CLJBoxPtr* const b0 = boxes0->constData();
                const CLJBoxPtr* const b1 = boxes1->constData();
            
                for (int i = range.begin(); i != range.end(); ++i)
                {
                    if (ptr->distance() < coul_cutoff and ptr->distance() < lj_cutoff)
                    {
                        func->total(b0[ptr->box0()].read().atoms(),
                                    b1[ptr->box1()].read().atoms(),
                                    coul_nrg[i], lj_nrg[i], ptr->distance());
                    }
                    else if (ptr->distance() < coul_cutoff)
                    {
                        coul_nrg[i]
                            = func->coulomb(b0[ptr->box0()].read().atoms(),
                                            b1[ptr->box1()].read().atoms(),
                                            ptr->distance());
                        
                        lj_nrg[i] = 0;
                    }
                    else if (ptr->distance() < lj_cutoff)
                    {
                        lj_nrg[i]
                            = func->lj(b0[ptr->box0()].read().atoms(),
                                       b1[ptr->box1()].read().atoms(),
                                       ptr->distance());
                        
                        coul_nrg[i] = 0;
                    }
                    else
                    {
                        coul_nrg[i] = 0;
                        lj_nrg[i] = 0;
                    }
                    
                    ptr += 1;
                }
            }
            
        private:
            const CLJFunction* const func;
            const CLJBoxDistance* const dists;
            const CLJBoxes::Container* const boxes0;
            const CLJBoxes::Container* const boxes1;
            
            double *coul_nrg;
            double *lj_nrg;
            
            const float coul_cutoff;
            const float lj_cutoff;
        };

        /** This is a private helper class that is used to calculate the
            coulomb and LJ energy in parallel using Intel TBB */
        class TotalWithoutCutoff2
        {
        public:
            TotalWithoutCutoff2() : func(0), dists(0), boxes0(0), boxes1(0)
            {}
        
            TotalWithoutCutoff2(const CLJFunction* const function,
                                const CLJBoxDistance* const distances,
                                const CLJBoxes::Container &cljboxes0,
                                const CLJBoxes::Container &cljboxes1,
                                double *coulomb_energy, double *lj_energy)
                : func(function), dists(distances), boxes0(&cljboxes0),
                  boxes1(&cljboxes1), coul_nrg(coulomb_energy), lj_nrg(lj_energy)
            {}
            
            ~TotalWithoutCutoff2()
            {}
            
            void operator()(const tbb::blocked_range<int> &range) const
            {
                const CLJBoxDistance* ptr = dists + range.begin();
                const CLJBoxPtr* const b0 = boxes0->constData();
                const CLJBoxPtr* const b1 = boxes1->constData();
            
                for (int i = range.begin(); i != range.end(); ++i)
                {
                    func->total(b0[ptr->box0()].read().atoms(),
                                b1[ptr->box1()].read().atoms(),
                                coul_nrg[i], lj_nrg[i], ptr->distance());
                    
                    ptr += 1;
                }
            }
            
        private:
            const CLJFunction* const func;
            const CLJBoxDistance* const dists;
            const CLJBoxes::Container* const boxes0;
            const CLJBoxes::Container* const boxes1;
            
            double *coul_nrg;
            double *lj_nrg;
        };

        /** This is a private helper class that is used to calculate the
            coulomb and LJ energy in parallel using Intel TBB */
        class DeltaWithCutoff
        {
        public:
            DeltaWithCutoff() : func(0), dists(0), atoms0(0), boxes1(0),
                                coul_cutoff(0), lj_cutoff(0)
            {}
        
            DeltaWithCutoff(const CLJFunction* const function,
                            const CLJBoxDistance* const distances,
                            const CLJAtoms &cljatoms0,
                            const CLJBoxes::Container &cljboxes1,
                            const float coulomb_cutoff, const float lenj_cutoff,
                            double *coulomb_energy, double *lj_energy)
                : func(function), dists(distances), atoms0(&cljatoms0),
                  boxes1(&cljboxes1), coul_nrg(coulomb_energy), lj_nrg(lj_energy),
                  coul_cutoff(coulomb_cutoff), lj_cutoff(lenj_cutoff)
            {}
            
            ~DeltaWithCutoff()
            {}
            
            void operator()(const tbb::blocked_range<int> &range) const
            {
                const CLJBoxDistance* ptr = dists + range.begin();
                const CLJBoxPtr* const b1 = boxes1->constData();
            
                for (int i = range.begin(); i != range.end(); ++i)
                {
                    if (ptr->distance() < coul_cutoff and ptr->distance() < lj_cutoff)
                    {
                        func->total(*atoms0,
                                    b1[ptr->box1()].read().atoms(),
                                    coul_nrg[i], lj_nrg[i], ptr->distance());
                    }
                    else if (ptr->distance() < coul_cutoff)
                    {
                        coul_nrg[i]
                            = func->coulomb(*atoms0,
                                            b1[ptr->box1()].read().atoms(),
                                            ptr->distance());
                        
                        lj_nrg[i] = 0;
                    }
                    else if (ptr->distance() < lj_cutoff)
                    {
                        lj_nrg[i]
                            = func->lj(*atoms0,
                                       b1[ptr->box1()].read().atoms(),
                                       ptr->distance());
                        
                        coul_nrg[i] = 0;
                    }
                    
                    ptr += 1;
                }
            }
            
        private:
            const CLJFunction* const func;
            const CLJBoxDistance* const dists;
            const CLJAtoms* const atoms0;
            const CLJBoxes::Container* const boxes1;
            
            double *coul_nrg;
            double *lj_nrg;
            
            const float coul_cutoff;
            const float lj_cutoff;
        };

        /** This is a private helper class that is used to calculate the
            coulomb and LJ energy in parallel using Intel TBB */
        class DeltaWithoutCutoff
        {
        public:
            DeltaWithoutCutoff() : func(0), atoms0(0), boxes1(0),
                                   coul_nrg(0), lj_nrg(0)
            {}
        
            DeltaWithoutCutoff(const CLJFunction* const function,
                               const CLJAtoms &cljatoms0,
                               const CLJBoxes::Container &cljboxes1,
                               double *coulomb_energy, double *lj_energy)
                : func(function), atoms0(&cljatoms0),
                  boxes1(&cljboxes1), coul_nrg(coulomb_energy), lj_nrg(lj_energy)
            {}
            
            ~DeltaWithoutCutoff()
            {}
            
            void operator()(const tbb::blocked_range<int> &range) const
            {
                const CLJBoxPtr* const b1 = boxes1->constData();
            
                for (int i = range.begin(); i != range.end(); ++i)
                {
                    func->total(*atoms0,
                                b1[i].read().atoms(),
                                coul_nrg[i], lj_nrg[i]);
                }
            }
            
        private:
            const CLJFunction* const func;
            const CLJAtoms* const atoms0;
            const CLJBoxes::Container* const boxes1;
            
            double *coul_nrg;
            double *lj_nrg;
        };
    
    } // end of namespace detail
} // end of namespace SireMM

/** Calculate the energy between all of the atoms in the passed CLJBoxes
    using the passed CLJFunction, returning
    the coulomb and LJ energy as a tuple (coulomb,lj) */
tuple<double,double> CLJCalculator::calculate(const CLJFunction &func, const CLJBoxes &boxes) const
{
    //get the cutoffs for the function
    if (func.hasCutoff())
    {
        Length coul_cutoff = func.coulombCutoff();
        Length lj_cutoff = func.ljCutoff();
        
        Length max_cutoff( coul_cutoff.value() >= lj_cutoff.value() ?
                           coul_cutoff : lj_cutoff );
        
        //get the list of box pairs that are within the cutoff distance
        QVector<CLJBoxDistance> dists = CLJBoxes::getDistances(func.space(), boxes, max_cutoff);

        //now create the space to hold the calculated energies
        QVarLengthArray<double> coul_nrgs( dists.count() );
        QVarLengthArray<double> lj_nrgs( dists.count() );
        
        //now create the object that will be used by TBB to calculate the energies
        detail::TotalWithCutoff helper(&func, dists.constData(), boxes.occupiedBoxes(),
                                       coul_cutoff.value(), lj_cutoff.value(),
                                       coul_nrgs.data(), lj_nrgs.data());

        //now perform the calculation in parallel
        tbb::parallel_for(tbb::blocked_range<int>(0,dists.count()), helper);
        
        if (reproducible_sum)
        {
            //do a sorted sum of energies so that we get the same result no matter the order
            //of calculation
            qSort(coul_nrgs);
            qSort(lj_nrgs);
        }
        
        double cnrg = 0;
        double ljnrg = 0;
        
        const double *coul_nrgs_array = coul_nrgs.constData();
        const double *lj_nrgs_array = lj_nrgs.constData();
        
        for (int i=0; i<coul_nrgs.count(); ++i)
        {
            cnrg += *coul_nrgs_array;
            ljnrg += *lj_nrgs_array;
            
            ++coul_nrgs_array;
            ++lj_nrgs_array;
        }
        
        return tuple<double,double>(cnrg, ljnrg);
    }
    else
    {
        //get the list of box pairs that are within the cutoff distance
        QVector<CLJBoxDistance> dists = CLJBoxes::getDistances(func.space(), boxes);

        //now create the space to hold the calculated energies
        QVarLengthArray<double> coul_nrgs( dists.count() );
        QVarLengthArray<double> lj_nrgs( dists.count() );
        
        //now create the object that will be used by TBB to calculate the energies
        detail::TotalWithoutCutoff helper(&func, dists.constData(), boxes.occupiedBoxes(),
                                          coul_nrgs.data(), lj_nrgs.data());

        //now perform the calculation in parallel
        tbb::parallel_for(tbb::blocked_range<int>(0,dists.count()), helper);
        
        if (reproducible_sum)
        {
            //do a sorted sum of energies so that we get the same result no matter the order
            //of calculation
            qSort(coul_nrgs);
            qSort(lj_nrgs);
        }
        
        double cnrg = 0;
        double ljnrg = 0;
        
        const double *coul_nrgs_array = coul_nrgs.constData();
        const double *lj_nrgs_array = lj_nrgs.constData();
        
        for (int i=0; i<coul_nrgs.count(); ++i)
        {
            cnrg += *coul_nrgs_array;
            ljnrg += *lj_nrgs_array;
            
            ++coul_nrgs_array;
            ++lj_nrgs_array;
        }
        
        return tuple<double,double>(cnrg, ljnrg);
    }
}

/** Calculate the energy between all of the atoms in the passed two CLJBoxes
    using the passed CLJFunction, returning
    the coulomb and LJ energy as a tuple (coulomb,lj) */
tuple<double,double> CLJCalculator::calculate(const CLJFunction &func,
                                              const CLJBoxes &boxes0, const CLJBoxes &boxes1) const
{
    if (boxes0.nOccupiedBoxes() > boxes1.nOccupiedBoxes())
        return this->calculate(func, boxes1, boxes0);

    if (boxes0.nOccupiedBoxes() == 1 and boxes0.length() == boxes1.length())
    {
        //there is only a single box, so use a different parallelisation strategy
        const CLJBoxPtr* const b0 = boxes0.constData();
        const CLJBoxPtr* const b1 = boxes1.constData();
        
        const CLJBoxIndex &idx0 = b0[0].read().index();
        const CLJAtoms &atoms0 = b0[0].read().atoms();
        const int n1 = boxes1.nOccupiedBoxes();
        
        const float coul_cutoff = func.coulombCutoff().value();
        const float lj_cutoff = func.ljCutoff().value();
        
        const float max_cutoff = qMax(coul_cutoff, lj_cutoff);
        
        double cnrg(0), ljnrg(0);
        
        for (int j=0; j<n1; ++j)
        {
            const float mindist = boxes0.getDistance(func.space(), idx0, b1[j].read().index());
            
            if (mindist < max_cutoff)
            {
                double icnrg(0), iljnrg(0);
                func(atoms0, b1[j].read().atoms(), icnrg, iljnrg, mindist);
                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
        
        return tuple<double,double>(cnrg,ljnrg);
    }
    else
    {
        //get the cutoffs for the function
        if (func.hasCutoff())
        {
            Length coul_cutoff = func.coulombCutoff();
            Length lj_cutoff = func.ljCutoff();
            
            Length max_cutoff( coul_cutoff.value() >= lj_cutoff.value() ?
                               coul_cutoff : lj_cutoff );
            
            //get the list of box pairs that are within the cutoff distance
            QVector<CLJBoxDistance> dists = CLJBoxes::getDistances(func.space(),
                                                                   boxes0, boxes1, max_cutoff);

            //now create the space to hold the calculated energies
            QVarLengthArray<double> coul_nrgs( dists.count() );
            QVarLengthArray<double> lj_nrgs( dists.count() );
            
            //now create the object that will be used by TBB to calculate the energies
            detail::TotalWithCutoff2 helper(&func, dists.constData(),
                                            boxes0.occupiedBoxes(),
                                            boxes1.occupiedBoxes(),
                                            coul_cutoff.value(), lj_cutoff.value(),
                                            coul_nrgs.data(), lj_nrgs.data());

            //now perform the calculation in parallel
            tbb::parallel_for(tbb::blocked_range<int>(0,dists.count()), helper);
            
            if (reproducible_sum)
            {
                //do a sorted sum of energies so that we get the same result no matter the order
                //of calculation
                qSort(coul_nrgs);
                qSort(lj_nrgs);
            }
            
            double cnrg = 0;
            double ljnrg = 0;
            
            const double *coul_nrgs_array = coul_nrgs.constData();
            const double *lj_nrgs_array = lj_nrgs.constData();
            
            for (int i=0; i<coul_nrgs.count(); ++i)
            {
                cnrg += *coul_nrgs_array;
                ljnrg += *lj_nrgs_array;
                
                ++coul_nrgs_array;
                ++lj_nrgs_array;
            }
            
            return tuple<double,double>(cnrg, ljnrg);
        }
        else
        {
            //get the list of box pairs that are within the cutoff distance
            QVector<CLJBoxDistance> dists = CLJBoxes::getDistances(func.space(), boxes0, boxes1);

            //now create the space to hold the calculated energies
            QVarLengthArray<double> coul_nrgs( dists.count() );
            QVarLengthArray<double> lj_nrgs( dists.count() );
            
            //now create the object that will be used by TBB to calculate the energies
            detail::TotalWithoutCutoff2 helper(&func, dists.constData(),
                                               boxes0.occupiedBoxes(),
                                               boxes1.occupiedBoxes(),
                                               coul_nrgs.data(), lj_nrgs.data());

            //now perform the calculation in parallel
            tbb::parallel_for(tbb::blocked_range<int>(0,dists.count()), helper);
            
            if (reproducible_sum)
            {
                //do a sorted sum of energies so that we get the same result no matter the order
                //of calculation
                qSort(coul_nrgs);
                qSort(lj_nrgs);
            }
            
            double cnrg = 0;
            double ljnrg = 0;
            
            const double *coul_nrgs_array = coul_nrgs.constData();
            const double *lj_nrgs_array = lj_nrgs.constData();
            
            for (int i=0; i<coul_nrgs.count(); ++i)
            {
                cnrg += *coul_nrgs_array;
                ljnrg += *lj_nrgs_array;
                
                ++coul_nrgs_array;
                ++lj_nrgs_array;
            }
            
            return tuple<double,double>(cnrg, ljnrg);
        }
    }
}

/** Calculate the energy between all of the atoms in the 'atoms0' and 'atoms1'
    using the passed CLJFunction, returning
    the coulomb and LJ energy as a tuple (coulomb,lj) */
tuple<double,double> CLJCalculator::calculate(const CLJFunction &func,
                                              const CLJAtoms &atoms0, const CLJBoxes &boxes1) const
{
    const int n1 = boxes1.nOccupiedBoxes();

    if (n1 == 0 or atoms0.count() == 0)
        return tuple<double,double>(0,0);
    
    else if (n1 == 1)
    {
        return func.calculate(atoms0, boxes1);
    }
    else if (func.hasCutoff())
    {
        //we have a cutoff, so need to find all boxes that are
        //within range of these atoms
        Length coul_cutoff = func.coulombCutoff();
        Length lj_cutoff = func.ljCutoff();
        
        Length max_cutoff( coul_cutoff.value() >= lj_cutoff.value() ?
                           coul_cutoff : lj_cutoff );
        
        //get the list of box pairs that are within the cutoff distance
        QVector<CLJBoxDistance> dists = CLJBoxes::getDistances(func.space(),
                                                               atoms0, boxes1, max_cutoff);

        //first, create the space to hold the calculated energies
        QVarLengthArray<double> coul_nrgs(dists.count());
        QVarLengthArray<double> lj_nrgs(dists.count());
        
        //now create the object that will be used by TBB to calculate the energies
        detail::DeltaWithCutoff helper(&func, dists.constData(),
                                       atoms0, boxes1.occupiedBoxes(),
                                       coul_cutoff.value(), lj_cutoff.value(),
                                       coul_nrgs.data(), lj_nrgs.data());

        //now perform the calculation in parallel
        tbb::parallel_for(tbb::blocked_range<int>(0,dists.count()), helper);
        
        if (reproducible_sum)
        {
            //do a sorted sum of energies so that we get the same result no matter the order
            //of calculation
            qSort(coul_nrgs);
            qSort(lj_nrgs);
        }
        
        double cnrg = 0;
        double ljnrg = 0;
        
        const double *coul_nrgs_array = coul_nrgs.constData();
        const double *lj_nrgs_array = lj_nrgs.constData();
        
        for (int i=0; i<coul_nrgs.count(); ++i)
        {
            cnrg += *coul_nrgs_array;
            ljnrg += *lj_nrgs_array;
            
            ++coul_nrgs_array;
            ++lj_nrgs_array;
        }

        return tuple<double,double>(cnrg,ljnrg);
    }
    else
    {
        //no cutoff, so just loop over all boxes in parallel
        //first, create the space to hold the calculated energies
        QVarLengthArray<double> coul_nrgs(n1);
        QVarLengthArray<double> lj_nrgs(n1);
        
        //now create the object that will be used by TBB to calculate the energies
        detail::DeltaWithoutCutoff helper(&func, atoms0, boxes1.occupiedBoxes(),
                                          coul_nrgs.data(), lj_nrgs.data());

        //now perform the calculation in parallel
        tbb::parallel_for(tbb::blocked_range<int>(0,n1), helper);
        
        if (reproducible_sum)
        {
            //do a sorted sum of energies so that we get the same result no matter the order
            //of calculation
            qSort(coul_nrgs);
            qSort(lj_nrgs);
        }
        
        double cnrg = 0;
        double ljnrg = 0;
        
        const double *coul_nrgs_array = coul_nrgs.constData();
        const double *lj_nrgs_array = lj_nrgs.constData();
        
        for (int i=0; i<coul_nrgs.count(); ++i)
        {
            cnrg += *coul_nrgs_array;
            ljnrg += *lj_nrgs_array;
            
            ++coul_nrgs_array;
            ++lj_nrgs_array;
        }

        return tuple<double,double>(cnrg,ljnrg);
    }
}

/** Calculate the energy between all of the atoms in the passed CLJBoxes
    using the passed array of CLJFunctions, returning the energies as
    a tuple of arrays of the coulomb and LJ energy (coulomb,lj) */
tuple< QVector<double>, QVector<double> > CLJCalculator::calculate(
                                const QVector<CLJFunctionPtr> &funcs,
                                const CLJBoxes &boxes) const
{
    QVector<double> cnrgs;
    QVector<double> ljnrgs;
    
    if (funcs.count() == 1)
    {
        tuple<double,double> nrgs = calculate( *(funcs.at(0)), boxes );
        cnrgs = QVector<double>(1, nrgs.get<0>());
        ljnrgs = QVector<double>(1, nrgs.get<1>());
    }
    else if (funcs.count() > 1)
    {
        for (int i=0; i<funcs.count(); ++i)
        {
            tuple<double,double> nrgs = calculate( *(funcs.at(i)), boxes );
            cnrgs.append( nrgs.get<0>() );
            ljnrgs.append( nrgs.get<1>() );
        }
    }
    
    return tuple< QVector<double>,QVector<double> >(cnrgs, ljnrgs);
}

/** Calculate the energy between all of the atoms in the passed two CLJBoxes
    using the passed array of CLJFunctions, returning the energies as
    a tuple of arrays of the coulomb and LJ energy (coulomb,lj) */
tuple< QVector<double>, QVector<double> > CLJCalculator::calculate(
                                const QVector<CLJFunctionPtr> &funcs,
                                const CLJBoxes &boxes0, const CLJBoxes &boxes1) const
{
    QVector<double> cnrgs;
    QVector<double> ljnrgs;
    
    if (funcs.count() == 1)
    {
        tuple<double,double> nrgs = calculate( *(funcs.at(0)), boxes0, boxes1 );
        cnrgs = QVector<double>(1, nrgs.get<0>());
        ljnrgs = QVector<double>(1, nrgs.get<1>());
    }
    else if (funcs.count() > 1)
    {
        for (int i=0; i<funcs.count(); ++i)
        {
            tuple<double,double> nrgs = calculate( *(funcs.at(i)), boxes0, boxes1 );
            cnrgs.append( nrgs.get<0>() );
            ljnrgs.append( nrgs.get<1>() );
        }
    }
    
    return tuple< QVector<double>,QVector<double> >(cnrgs, ljnrgs);
}

/** Calculate the energy between all of the atoms in the passed atoms0 and atoms1
    using the passed array of CLJFunctions, returning the energies as
    a tuple of arrays of the coulomb and LJ energy (coulomb,lj) */
tuple< QVector<double>, QVector<double> > CLJCalculator::calculate(
                                const QVector<CLJFunctionPtr> &funcs,
                                const CLJAtoms &atoms0, const CLJBoxes &boxes1) const
{
    QVector<double> cnrgs;
    QVector<double> ljnrgs;
    
    if (funcs.count() == 1)
    {
        tuple<double,double> nrgs = calculate( *(funcs.at(0)), atoms0, boxes1 );
        cnrgs = QVector<double>(1, nrgs.get<0>());
        ljnrgs = QVector<double>(1, nrgs.get<1>());
    }
    else if (funcs.count() > 1)
    {
        for (int i=0; i<funcs.count(); ++i)
        {
            tuple<double,double> nrgs = calculate( *(funcs.at(i)), atoms0, boxes1 );
            cnrgs.append( nrgs.get<0>() );
            ljnrgs.append( nrgs.get<1>() );
        }
    }
    
    return tuple< QVector<double>,QVector<double> >(cnrgs, ljnrgs);
}
