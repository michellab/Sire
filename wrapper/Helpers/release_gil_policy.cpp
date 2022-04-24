
#include "release_gil_policy.hpp"

#include "SireError/getbacktrace.h"

#include <QDebug>

//#define SIRE_DISABLE_GIL_POLICY 1
//#define SIRE_PRINT_GIL_STATUS 1

boost::python::detail::GilHolder::GilHolder() : thread_state(0)
{
    #ifndef SIRE_DISABLE_GIL_POLICY
        #ifdef SIRE_PRINT_GIL_STATUS
            qDebug() << "--RELEASE GIL";
            //for (const auto &bt : SireError::getBackTrace())
            //{
            //    qDebug() << bt;
            //}
        #endif
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
                #ifdef SIRE_PRINT_GIL_STATUS
                    qDebug() << "--ACQUIRE GIL";
                #endif
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
        else
        {
            boost::python::release_gil_policy::gil.setLocalData(new boost::python::detail::GilHolder());
        }
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

/** Acquire the GIL without handling via a RAII object */
void boost::python::release_gil_policy::acquire_gil_no_raii()
{
    #ifndef SIRE_DISABLE_GIL_POLICY
        if (gil.hasLocalData())
        {
            gil.setLocalData(0);
        }
    #endif
}

/** Release the GIL without handling via a RAII object. You will
    need to match this with a call to re-acquire the GIL before
    you go back to Python
*/
void boost::python::release_gil_policy::release_gil_no_raii()
{
    #ifndef SIRE_DISABLE_GIL_POLICY
        if (boost::python::release_gil_policy::gil.hasLocalData())
        {
            qDebug() << "WARNING - DOUBLE HOLD GIL - POTENTIAL FOR DEADLOCK!";
            return;
        }

        boost::python::release_gil_policy::gil.setLocalData(new boost::python::detail::GilHolder());
    #endif
}

QThreadStorage<boost::python::detail::GilHolder*> boost::python::release_gil_policy::gil;
