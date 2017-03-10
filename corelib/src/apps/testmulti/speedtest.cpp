
#include "SireMaths/multifloat.h"
#include "SireMaths/multidouble.h"
#include "SireMaths/rangenerator.h"

#include "SireBase/parallel.h"

#include <QElapsedTimer>
#include <QVector>

#include <QDebug>

static const qint64 NTEST = 16 * 2000;

using namespace SireMaths;

int main(int argc, const char **argv)
{
    QElapsedTimer t;

    qDebug() << "Calculating the coulomb energy between two groups of"
             << NTEST << "particles. This involves calculating and summing "
             << (NTEST*NTEST) << "coulomb interaction energies.";

    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        float *x0a = (float*)_mm_malloc( NTEST*sizeof(float), 32 );
        float *y0a = (float*)_mm_malloc( NTEST*sizeof(float), 32 );
        float *z0a = (float*)_mm_malloc( NTEST*sizeof(float), 32 );
        float *x1a = (float*)_mm_malloc( NTEST*sizeof(float), 32 );
        float *y1a = (float*)_mm_malloc( NTEST*sizeof(float), 32 );
        float *z1a = (float*)_mm_malloc( NTEST*sizeof(float), 32 );

        float *q0a = (float*)_mm_malloc( NTEST*sizeof(float), 32 );
        float *q1a = (float*)_mm_malloc( NTEST*sizeof(float), 32 );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        float *x0a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *y0a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *z0a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *x1a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *y1a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *z1a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );

        float *q0a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *q1a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
    #else
        float *x0a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *y0a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *z0a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *x1a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *y1a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *z1a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );

        float *q0a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
        float *q1a = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
    #endif
    #endif

    RanGenerator rand;
    rand.seed(1234561);

    qDebug() << "Generating the random numbers...";
    {
        t.start();

        for (int i=0; i<NTEST; ++i)
        {
            x0a[i] = rand.rand(-5,-3);
            y0a[i] = rand.rand(-5,-3);
            z0a[i] = rand.rand(-5,-3);
            q0a[i] = rand.rand(-1, 1);

            x1a[i] = rand.rand(3,5);
            y1a[i] = rand.rand(3,5);
            z1a[i] = rand.rand(3,5);
            q1a[i] = rand.rand(1,1);
        }

        qint64 ns = t.nsecsElapsed();
        qDebug() << "Generating the random numbers took" << (ns*1e-6) << "ms";
    }

    qDebug() << "Converting to MultiFloat...";
    t.start();
    QVector<MultiFloat> x0f = MultiFloat::fromArray(x0a, NTEST);
    QVector<MultiFloat> y0f = MultiFloat::fromArray(y0a, NTEST);
    QVector<MultiFloat> z0f = MultiFloat::fromArray(z0a, NTEST);
    QVector<MultiFloat> q0f = MultiFloat::fromArray(q0a, NTEST);

    QVector<MultiFloat> x1f = MultiFloat::fromArray(x1a, NTEST);
    QVector<MultiFloat> y1f = MultiFloat::fromArray(y1a, NTEST);
    QVector<MultiFloat> z1f = MultiFloat::fromArray(z1a, NTEST);
    QVector<MultiFloat> q1f = MultiFloat::fromArray(q1a, NTEST);
    qDebug() << "Took" << (t.nsecsElapsed()*1e-6) << "ms";

    qDebug() << "Converting to MultiDouble...";
    t.start();
    QVector<MultiDouble> x0d = MultiDouble::fromArray(x0a, NTEST);
    QVector<MultiDouble> y0d = MultiDouble::fromArray(y0a, NTEST);
    QVector<MultiDouble> z0d = MultiDouble::fromArray(z0a, NTEST);
    QVector<MultiDouble> q0d = MultiDouble::fromArray(q0a, NTEST);

    QVector<MultiDouble> x1d = MultiDouble::fromArray(x1a, NTEST);
    QVector<MultiDouble> y1d = MultiDouble::fromArray(y1a, NTEST);
    QVector<MultiDouble> z1d = MultiDouble::fromArray(z1a, NTEST);
    QVector<MultiDouble> q1d = MultiDouble::fromArray(q1a, NTEST);
    qDebug() << "Took" << (t.nsecsElapsed()*1e-6) << "ms";

    qDebug() << "Calculating the coulomb energy using floats...";
    for (int i=0; i<5; ++i)
    {
        t.start();
        const MultiFloat *x0 = x0f.constData();
        const MultiFloat *y0 = y0f.constData();
        const MultiFloat *z0 = z0f.constData();
        const MultiFloat *q0 = q0f.constData();

        const MultiFloat *x1 = x1f.constData();
        const MultiFloat *y1 = y1f.constData();
        const MultiFloat *z1 = z1f.constData();
        const MultiFloat *q1 = q1f.constData();

        MultiDouble coul_nrg = tbb::parallel_reduce( tbb::blocked_range<int>(0,x0f.count()), 
                               MultiDouble(0),
                               [&](tbb::blocked_range<int> r, MultiDouble my_coul_nrg )
        {
            MultiFloat delta;
            MultiFloat dist2;
            MultiFloat ox, oy, oz, oq;

            for (int i=r.begin(); i<r.end(); ++i)
            {
                const MultiFloat &x = x0[i];
                const MultiFloat &y = y0[i];
                const MultiFloat &z = z0[i];
                const MultiFloat &q = q0[i];

                for (int j=0; j<x1f.count(); ++j)
                {
                    ox = x1[j];
                    oy = y1[j];
                    oz = z1[j];
                    oq = q1[j];

                    for (int k=0; k<MultiFloat::count(); ++k)
                    {
                        delta = x - ox;
                        dist2 = delta*delta;
                        delta = y - oy;
                        dist2.multiplyAdd(delta, delta);
                        delta = z - oz;
                        dist2.multiplyAdd(delta, delta);

                        my_coul_nrg += q * oq * dist2.rsqrt_approx_nr();

                        //rotate the multifloat to process the other distances
                        ox = ox.rotate();
                        oy = oy.rotate();
                        oz = oz.rotate();
                        oq = oq.rotate();
                    }
                }
            }

            return my_coul_nrg;

        }, std::plus<MultiFloat>() );

        qint64 ns = t.nsecsElapsed();
        qDebug() << "Energies" << coul_nrg.toString() << "TOTAL" << coul_nrg.doubleSum();
        qDebug() << "took" << (ns*1e-6) << "ms";
        float gflops = (NTEST * NTEST * 17) / (1.f * ns);   // approx_nr has 6 ops
        qDebug() << "Speed is" << gflops << "GFLOPS";
    }

    qDebug() << "Calculating the coulomb energy using doubles...";
    for (int i=0; i<5; ++i)
    {
        t.start();
        const MultiDouble *x0 = x0d.constData();
        const MultiDouble *y0 = y0d.constData();
        const MultiDouble *z0 = z0d.constData();
        const MultiDouble *q0 = q0d.constData();

        const MultiDouble *x1 = x1d.constData();
        const MultiDouble *y1 = y1d.constData();
        const MultiDouble *z1 = z1d.constData();
        const MultiDouble *q1 = q1d.constData();

        MultiDouble coul_nrg = tbb::parallel_reduce( tbb::blocked_range<int>(0,x0d.count()),
                                                     MultiDouble(0),
                                                     [&](tbb::blocked_range<int> r, MultiDouble my_coul_nrg)
        {
            MultiDouble delta;
            MultiDouble dist2;
            MultiDouble cnrg;
            MultiDouble ox, oy, oz, oq;

            for (int i=r.begin(); i<r.end(); ++i)
            {
                const MultiDouble &x = x0[i];
                const MultiDouble &y = y0[i];
                const MultiDouble &z = z0[i];
                const MultiDouble &q = q0[i];

                for (int j=0; j<x1d.count(); ++j)
                {
                    ox = x1[j];
                    oy = y1[j];
                    oz = z1[j];
                    oq = q1[j];

                    for (int k=0; k<MultiDouble::count(); ++k)
                    {
                        delta = x - ox;
                        dist2 = delta*delta;
                        delta = y - oy;
                        dist2.multiplyAdd(delta, delta);
                        delta = z - oz;
                        dist2.multiplyAdd(delta, delta);

                        my_coul_nrg += q * oq * dist2.rsqrt_approx_nr();
    
                        ox = ox.rotate();
                        oy = oy.rotate();
                        oz = oz.rotate();
                        oq = oq.rotate();
                    }
                }
            }

            return my_coul_nrg;
        }, std::plus<MultiDouble>() );

        qint64 ns = t.nsecsElapsed();
        qDebug() << "Energies" << coul_nrg.toString() << "TOTAL" << coul_nrg.sum();
        qDebug() << "took" << (1e-6*ns) << "ms";
        float gflops = (NTEST * NTEST * 17) / (1.f*ns);  // approx_nr has 7 ops
        qDebug() << "Speed is" << gflops << "GFLOPS";
    }

    _mm_free(x0a);
    _mm_free(x1a);
    _mm_free(y0a);
    _mm_free(y1a);
    _mm_free(z0a);
    _mm_free(z1a);
    _mm_free(q0a);
    _mm_free(q1a);
 
    return 0;
}

