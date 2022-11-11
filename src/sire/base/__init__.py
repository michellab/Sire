"""
.. currentmodule:: sire.base

"""

from ..legacy import Base as _Base

from .. import use_new_api as _use_new_api
_use_new_api()

from ..legacy.Base import Property, Properties, PropertyMap, PropertyList, \
                          BooleanProperty, NumberProperty, StringProperty

wrap = _Base.wrap


def _fix_wrapped_property_base_type(CLASS):
    """Private function that adds in all of the operators so that,
       as much as possible, the base type PropertyWrappers (e.g.
       NumberProperty, StringProperty) behave like numbers or
       strings, and can auto-convert to those types as needed.
    """
    def __value(v):
        try:
            v = v.value()
        except Exception:
            pass

        try:
            if v.is_integer():
                return int(v)
        except Exception:
            pass

        return v

    def __add__(obj, other):
        return __value(obj) + __value(other)

    CLASS.__add__ = __add__
    CLASS.__radd__ = __add__

    def __sub__(obj, other):
        return __value(obj) - __value(other)

    def __rsub__(obj, other):
        return __value(other) - __value(obj)

    CLASS.__sub__ = __sub__
    CLASS.__rsub__ = __rsub__

    def __eq__(obj, other):
        return __value(obj) == __value(other)

    CLASS.__eq__ = __eq__

    def __ne__(obj, other):
        return __value(obj) != __value(other)

    CLASS.__ne__ = __ne__

    def __gt__(obj, other):
        return __value(obj) > __value(other)

    CLASS.__gt__ = __gt__

    def __ge__(obj, other):
        return __value(obj) >= __value(other)

    CLASS.__ge__ = __ge__

    def __lt__(obj, other):
        return __value(obj) < __value(other)

    CLASS.__lt__ = __lt__

    def __le__(obj, other):
        return __value(obj) <= __value(other)

    CLASS.__le__ = __le__

    def __float__(obj):
        return float(obj.value())

    CLASS.__float__ = __float__

    def __int__(obj):
        return int(obj.value())

    CLASS.__int__ = __int__

    def __str__(obj):
        return str(__value(obj))

    CLASS.__str__ = __str__


_fix_wrapped_property_base_type(BooleanProperty)
_fix_wrapped_property_base_type(NumberProperty)
_fix_wrapped_property_base_type(StringProperty)
