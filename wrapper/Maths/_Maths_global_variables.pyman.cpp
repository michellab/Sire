
#include "_Maths_global_variables.pyman.hpp"
#include <boost/python.hpp>
#include "SireMaths/constants.h"

using namespace boost::python;
using namespace SireMaths;

void register_man_global_variables()
{
    scope().attr("e") = e;

    scope().attr("log2_e") = log2_e;

    scope().attr("log10_e") = log10_e;

    scope().attr("sqrt_two") = sqrt_two;

    scope().attr("sqrt_half") = sqrt_half;

    scope().attr("sqrt_three") = sqrt_three;

    scope().attr("pi") = pi;

    scope().attr("two_pi") = two_pi;

    scope().attr("pi_over_two") = pi_over_two;

    scope().attr("pi_4") = pi_4;

    scope().attr("sqrtpi") = sqrtpi;

    scope().attr("two_sqrtpi") = two_sqrtpi;

    scope().attr("one_over_pi") = one_over_pi;

    scope().attr("two_over_pi") = two_over_pi;

    scope().attr("ln_ten") = ln_ten;

    scope().attr("ln_two") = ln_two;

    scope().attr("ln_pi") = ln_pi;

    scope().attr("euler") = euler;

    scope().attr("small") = small;

    scope().attr("tiny") = tiny;

    scope().attr("smallest") = smallest;

    scope().attr("large") = large;

    scope().attr("huge") = huge;

    scope().attr("largest") = largest;

    scope().attr("smallest_int") = smallest_int;

    scope().attr("largest_int") = largest_int;

}

