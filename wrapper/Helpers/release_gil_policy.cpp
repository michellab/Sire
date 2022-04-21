
#include "release_gil_policy.hpp"

#include "SireError/getbacktrace.h"

#include <QDebug>

boost::python::detail::GilHolder::GilHolder() : thread_state(0)
{
    qDebug() << "RELEASE GIL";
    thread_state = PyEval_SaveThread();
}

boost::python::detail::GilHolder::~GilHolder()
{
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
}

boost::python::detail::GilRaiiData::GilRaiiData()
{}

boost::python::detail::GilRaiiData::~GilRaiiData()
{
    if (boost::python::release_gil_policy::gil.hasLocalData())
    {
        qDebug() << "WARNING - DOUBLE HOLD GIL - POTENTIAL FOR DEADLOCK!";
    }

    boost::python::release_gil_policy::gil.setLocalData(new boost::python::detail::GilHolder());
}

boost::python::GilRaii::GilRaii(bool acquired)
{
    if (acquired)
    {
        d.reset(new boost::python::detail::GilRaiiData());
    }
}

boost::python::GilRaii::~GilRaii()
{}

/** Acquire the GIL, returning a RAII object which will release
 *  the GIL when it is destroyed
 */
boost::python::GilRaii boost::python::release_gil_policy::acquire_gil()
{
    if (gil.hasLocalData())
    {
        gil.setLocalData(0);
        return GilRaii(true);
    }
    else
    {
        return GilRaii(false);
    }
}

QThreadStorage<boost::python::detail::GilHolder*> boost::python::release_gil_policy::gil;
