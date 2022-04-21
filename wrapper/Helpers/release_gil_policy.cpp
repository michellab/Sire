
#include "release_gil_policy.hpp"

#include "SireError/getbacktrace.h"

#include <QDebug>

//#define SIRE_DISABLE_GIL_POLICY 1

boost::python::detail::GilHolder::GilHolder() : thread_state(0)
{
    #ifndef SIRE_DISABLE_GIL_POLICY
        qDebug() << "RELEASE GIL";
        thread_state = PyEval_SaveThread();
    #endif
}

boost::python::detail::GilHolder::~GilHolder()
{
    #ifndef SIRE_DISABLE_GIL_POLICY
        if (thread_state)
        {
            if (_Py_IsFinalizing())
            {
                qDebug() << "FINALIZING!";
            }
            else
            {
                qDebug() << "ACQUIRE GIL";
                PyEval_RestoreThread(thread_state);
            }
        }
    #endif
}

boost::python::detail::GilRaiiData::GilRaiiData()
{}

boost::python::detail::GilRaiiData::~GilRaiiData()
{
    #ifndef SIRE_DISABLE_GIL_POLICY
        if (boost::python::release_gil_policy::gil.hasLocalData())
        {
            qDebug() << "WARNING - DOUBLE HOLD GIL - POTENTIAL FOR DEADLOCK!";
        }

        boost::python::release_gil_policy::gil.setLocalData(new boost::python::detail::GilHolder());
    #endif
}

boost::python::GilRaii::GilRaii(bool acquired)
{
    #ifndef SIRE_DISABLE_GIL_POLICY
        if (acquired)
        {
            d.reset(new boost::python::detail::GilRaiiData());
        }
    #endif
}

boost::python::GilRaii::~GilRaii()
{}

/** Acquire the GIL, returning a RAII object which will release
 *  the GIL when it is destroyed
 */
boost::python::GilRaii boost::python::release_gil_policy::acquire_gil()
{
    #ifdef SIRE_DISABLE_GIL_POLICY
        return boost::python::GilRaii(false);
    #else
        if (gil.hasLocalData())
        {
            gil.setLocalData(0);
            return GilRaii(true);
        }
        else
        {
            return GilRaii(false);
        }
    #endif
}

QThreadStorage<boost::python::detail::GilHolder*> boost::python::release_gil_policy::gil;
