#ifndef _HELPERS_RELEASE_GIL_POLICY_HPP_
#define _HELPERS_RELEASE_GIL_POLICY_HPP_

#include "boost/python.hpp"

#include <QThreadStorage>
#include <QDebug>

namespace boost {
    namespace python {

        namespace detail {
            class GilHolder
            {
            public:
                GilHolder() : thread_state(0)
                {
                    qDebug() << "release gil";
                    thread_state = PyEval_SaveThread();
                }

                ~GilHolder()
                {
                    if (thread_state)
                    {
                        qDebug() << "acquire gil";
                        PyEval_RestoreThread(thread_state);
                    }
                }

            private:
                PyThreadState *thread_state;
            };
        }

        struct release_gil_policy
        {
            template <class ArgumentPackage>
            static bool precall(ArgumentPackage const&)
            {
                gil.setLocalData(new detail::GilHolder());
                return true;
            }

            template <class ArgumentPackage>
            static PyObject* postcall(ArgumentPackage const&, PyObject* result)
            {
                if (gil.hasLocalData())
                {
                    gil.setLocalData(0);
                }

                return result;
            }

            typedef default_result_converter result_converter;
            typedef PyObject* argument_package;

            template <class Sig>
            struct extract_return_type : mpl::front<Sig>
            {};

        private:
            static QThreadStorage<detail::GilHolder*> gil;
        };
    }
}

#endif
