
#include "SireMaths/multifloat.h"
#include "SireMaths/multidouble.h"
#include "SireMaths/rangenerator.h"

#include <QElapsedTimer>
#include <QVector>

#include <QDebug>

static const int NTEST = 1024;
static const int NLOOPS = 2000000;

using namespace SireMaths;

void compareErrors(const QVector<MultiFloat> &exact, const QVector<MultiFloat> &approx,
                   const QVector<MultiFloat> &approx_nr)
{
    MultiFloat err(0), err_nr(0);
    for (int i=0; i<exact.count(); ++i)
    {
        MultiFloat del = exact[i] - approx[i];
        err += del*del;
        del = exact[i] - approx_nr[i];
        err_nr += del*del;
    }

    double rmsd_err = err.sum();
    double rmsd_err_nr = err_nr.sum();
    long count = MultiFloat::count() * exact.count();
    qDebug() << "Approx error = " << std::sqrt(rmsd_err/count);
    qDebug() << "Approx_nr error = " << std::sqrt(rmsd_err_nr/count);
}

int main(int argc, const char **argv)
{
    QElapsedTimer t;

    qDebug() << "Calculating lots of square roots and approximations....";

    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        float *x = (float*)_mm_malloc( NTEST*sizeof(float), 32 );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        float *x = (float*)_mm_malloc( NTEST*sizeof(float), 16 );
    #else
        float *x = (float*)malloc( NTEST*sizeof(float) );
    #endif
    #endif

    RanGenerator rand;
    rand.seed(1234561);

    qDebug() << "Generating the random numbers...";
    {
        t.start();

        for (int i=0; i<NTEST; ++i)
        {
            x[i] = rand.rand(0.0001, 10.0);
        }

        int ms = t.elapsed();
        qDebug() << "Generating the random numbers took" << ms << "ms";
    }

    qDebug() << "Converting to MultiFloat...";
    t.start();
    QVector<MultiFloat> xf = MultiFloat::fromArray(x, NTEST);
    qDebug() << "Took" << t.elapsed() << "ms";

    QVector<MultiFloat> exact = xf;
    QVector<MultiFloat> approx = xf;
    QVector<MultiFloat> approx_nr = xf;
    exact[0] = MultiFloat(0);
    approx[0] = MultiFloat(0);
    approx_nr[0] = MultiFloat(0);

    qDebug() << "\nTESTING THE SQUARE ROOT AND APPROXIMATIONS";
    qDebug() << "Calculating exact answers...";
    t.start();
    {
        MultiFloat sum(0);

        for (int j=0; j<NLOOPS; ++j)
        {
            for (int i=0; i<xf.count(); ++i)
            {
                exact[i] = xf[i].sqrt();
            }

            sum += exact[0];
        }

        qDebug() << "Sum" << sum.sum();
        int ms = t.elapsed();
        qDebug() << "Took" << ms << "ms";
        qDebug() << "Nanoseconds per operation:" << ((1000000.0 * ms) / (NLOOPS*NTEST));
    }

    qDebug() << "Calculating approx answers...";
    t.start();
    {
        MultiFloat sum(0);

        for (int j=0; j<NLOOPS; ++j)
        {
            for (int i=0; i<xf.count(); ++i)
            {
                approx[i] = xf[i].sqrt_approx();
            }

            sum += approx[0];
        }

        qDebug() << "Sum" << sum.sum();
        int ms = t.elapsed();
        qDebug() << "Took" << ms << "ms";
        qDebug() << "Nanoseconds per operation:" << ((1000000.0 * ms) / (NLOOPS*NTEST));
    }

    qDebug() << "Calculating approx_nr answers...";
    t.start();
    {
        MultiFloat sum(0);

        for (int j=0; j<NLOOPS; ++j)
        {
            for (int i=0; i<xf.count(); ++i)
            {
                approx_nr[i] = xf[i].sqrt_approx_nr();
            }

            sum += approx_nr[0];
        }

        qDebug() << "Sum" << sum.sum();
        int ms = t.elapsed();
        qDebug() << "Took" << ms << "ms";
        qDebug() << "Nanoseconds per operation:" << ((1000000.0 * ms) / (NLOOPS*NTEST));
    }

    qDebug() << xf[0].toString();
    qDebug() << xf[0].sqrt().toString();
    qDebug() << xf[0].sqrt_approx().toString();
    qDebug() << xf[0].sqrt_approx_nr().toString();

    compareErrors(exact, approx, approx_nr);

    qDebug() << "\nTESTING THE INVERSE SQUARE ROOT AND APPROXIMATIONS";
    qDebug() << "Calculating exact answers...";
    t.start();
    {
        MultiFloat sum(0);

        for (int j=0; j<NLOOPS; ++j)
        {
            for (int i=0; i<xf.count(); ++i)
            {
                exact[i] = xf[i].rsqrt();
            }

            sum += exact[0];
        }

        qDebug() << "Sum" << sum.sum();
        int ms = t.elapsed();
        qDebug() << "Took" << ms << "ms";
        qDebug() << "Nanoseconds per operation:" << ((1000000.0 * ms) / (NLOOPS*NTEST));
    }

    qDebug() << "Calculating approx answers...";
    t.start();
    {
        MultiFloat sum(0);

        for (int j=0; j<NLOOPS; ++j)
        {
            for (int i=0; i<xf.count(); ++i)
            {
                approx[i] = xf[i].rsqrt_approx();
            }

            sum += approx[0];
        }

        qDebug() << "Sum" << sum.sum();
        int ms = t.elapsed();
        qDebug() << "Took" << ms << "ms";
        qDebug() << "Nanoseconds per operation:" << ((1000000.0 * ms) / (NLOOPS*NTEST));
    }

    qDebug() << "Calculating approx_nr answers...";
    t.start();
    {
        MultiFloat sum(0);

        for (int j=0; j<NLOOPS; ++j)
        {
            for (int i=0; i<xf.count(); ++i)
            {
                approx_nr[i] = xf[i].rsqrt_approx_nr();
            }

            sum += approx_nr[0];
        }

        qDebug() << "Sum" << sum.sum();
        int ms = t.elapsed();
        qDebug() << "Took" << ms << "ms";
        qDebug() << "Nanoseconds per operation:" << ((1000000.0 * ms) / (NLOOPS*NTEST));
    }

    qDebug() << xf[0].toString();
    qDebug() << xf[0].rsqrt().toString();
    qDebug() << xf[0].rsqrt_approx().toString();
    qDebug() << xf[0].rsqrt_approx_nr().toString();

    compareErrors(exact, approx, approx_nr);

    qDebug() << "\nTESTING THE RECIPROCAL AND APPROXIMATIONS";
    qDebug() << "Calculating exact answers...";
    t.start();
    {
        MultiFloat sum(0);

        for (int j=0; j<NLOOPS; ++j)
        {
            for (int i=0; i<xf.count(); ++i)
            {
                exact[i] = xf[i].reciprocal();
            }

            sum += exact[0];
        }

        qDebug() << "Sum" << sum.sum();
        int ms = t.elapsed();
        qDebug() << "Took" << ms << "ms";
        qDebug() << "Nanoseconds per operation:" << ((1000000.0 * ms) / (NLOOPS*NTEST));
    }

    qDebug() << "Calculating approx answers...";
    t.start();
    {
        MultiFloat sum(0);

        for (int j=0; j<NLOOPS; ++j)
        {
            for (int i=0; i<xf.count(); ++i)
            {
                approx[i] = xf[i].reciprocal_approx();
            }

            sum += approx[0];
        }

        qDebug() << "Sum" << sum.sum();
        int ms = t.elapsed();
        qDebug() << "Took" << ms << "ms";
        qDebug() << "Nanoseconds per operation:" << ((1000000.0 * ms) / (NLOOPS*NTEST));
    }

    qDebug() << "Calculating approx_nr answers...";
    t.start();
    {
        MultiFloat sum(0);

        for (int j=0; j<NLOOPS; ++j)
        {
            for (int i=0; i<xf.count(); ++i)
            {
                approx_nr[i] = xf[i].reciprocal_approx_nr();
            }

            sum += approx_nr[0];
        }

        qDebug() << "Sum" << sum.sum();
        int ms = t.elapsed();
        qDebug() << "Took" << ms << "ms";
        qDebug() << "Nanoseconds per operation:" << ((1000000.0 * ms) / (NLOOPS*NTEST));
    }

    qDebug() << xf[0].toString();
    qDebug() << xf[0].reciprocal().toString();
    qDebug() << xf[0].reciprocal_approx().toString();
    qDebug() << xf[0].reciprocal_approx_nr().toString();

    compareErrors(exact, approx, approx_nr);

    free(x);

    return 0;
}

