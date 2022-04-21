#ifndef _HELPERS_RELEASE_GIL_POLICY_HPP_
#define _HELPERS_RELEASE_GIL_POLICY_HPP_

#include "boost/python.hpp"

#include <QThreadStorage>
#include <QDebug>

#include <memory>

#include "sireglobal.h"

namespace boost {
    namespace python {

        namespace detail {
            class SIRE_EXPORT GilHolder
            {
            public:
                GilHolder();

                ~GilHolder();

            private:
                PyThreadState *thread_state;
            };

            class SIRE_EXPORT GilRaiiData
            {
            public:
                GilRaiiData();
                ~GilRaiiData();
            };
        }

        class SIRE_EXPORT GilRaii
        {
        public:
            GilRaii(bool acquired=true);
            ~GilRaii();
        private:
            std::shared_ptr<detail::GilRaiiData> d;
        };

        class SIRE_EXPORT release_gil_policy
        {
        friend class detail::GilRaiiData;

        public:
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

            static GilRaii acquire_gil();

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
