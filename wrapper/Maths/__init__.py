"""
.. currentmodule:: Sire.Maths

Classes
=======

.. autosummary::
    :toctree: generated/

    Accumulator
    Average
    AverageAndStddev
    AxisSet
    BennettsFreeEnergyAverage
    Complex
    DistVector
    ExpAverage
    FreeEnergyAverage
    Histogram
    HistogramBin
    HistogramValue
    Line
    Matrix
    Median
    MultiDouble
    MultiFixed
    MultiFloat
    MultiInt
    MultiUInt
    MultiVector
    N4Matrix
    NMatrix
    NVector
    Plane
    Quaternion
    RanGenerator
    Rational
    RecordValues
    Sphere
    Torsion
    Transform
    Triangle
    TrigMatrix
    Vector
    VectorArrayProperty
    VectorProperty

Functions
=========

.. autosummary::
    :toctree: generated/

    gamma
    align
    boys
    brute_force_linear_assignment
    calculate_total_cost
    getAlignment
    getCentroid
    getRMSD
    get_lowest_total_cost
    incomplete_gamma_higher
    incomplete_gamma_lower
    kabasch
    kabaschFit
    multi_boys
    rotate
    solve_linear_assignment
    wrap

"""

import Sire.Qt
import Sire.Error
import Sire.Base

# Import all of the classes and functions from the C++ library
from Sire.Maths._Maths import *

__all__ = [ "Accumulator", "Average", "AverageAndStddev", "AxisSet",
            "BennettsFreeEnergyAverage", "Complex", "DistVector", "ExpAverage",
            "FreeEnergyAverage", "Histogram", "HistogramBin", "HistogramValue",
            "Line", "Matrix", "Median", "MultiDouble",
            "MultiFixed", "MultiFloat", "MultiInt", "MultiUInt",
            "MultiVector", "N4Matrix", "NMatrix", "NVector",
            "Plane", "Quaternion", "RanGenerator", "Rational",
            "RecordValues", "Sphere", "Torsion", "Transform",
            "Triangle", "TrigMatrix", "Vector", "VectorArrayProperty",
            "VectorProperty", "gamma", "align", "boys",
            "brute_force_linear_assignment", "calculate_total_cost", "getAlignment", "getCentroid",
            "getRMSD", "get_lowest_total_cost", "incomplete_gamma_higher", "incomplete_gamma_lower",
            "kabasch", "kabaschFit", "multi_boys", "rotate",
            "solve_linear_assignment", "wrap" ]

# Now define some pure Python functions and classes that are part of
# this library...

wrap = Sire.Base._add_wrap_function(wrap)

# No QVector<float> exposed (would have horrible casting bugs)
MultiFloat.toArray = staticmethod( MultiFloat.toDoubleArray )
