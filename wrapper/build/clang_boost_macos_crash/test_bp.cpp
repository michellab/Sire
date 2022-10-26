
#include "boost/python.hpp"

namespace bp = boost::python;

class TestClass
{
public:
    TestClass(int y, int m)
    {
    }

    ~TestClass()
    {
    }
};

void register_class()
{
    bp::class_<TestClass> c("TestClass",
                            bp::init<int, int>(
                                (bp::arg("y"), bp::arg("m"))));
}

BOOST_PYTHON_MODULE(_test_bp)
{
    register_class();
}
