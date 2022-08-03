"""
.. currentmodule:: sire.base

"""

from ..legacy import Base as _Base

from .. import use_new_api as _use_new_api
_use_new_api()

Property = _Base.Property
Properties = _Base.Properties

wrap = _Base.wrap


def _add_property_operators(C):
    import operator

    def __property_op__(lhs, rhs, op):
        if rhs.what().endswith("Property"):
            return op(lhs.value(), rhs.value())
        else:
            return op(lhs.value(), rhs)

    C.__add__ = lambda x, y: __property_op__(x, y, operator.add)
    C.__sub__ = lambda x, y: __property_op__(x, y, operator.sub)
    C.__mul__ = lambda x, y: __property_op__(x, y, operator.mul)
    C.__div__ = lambda x, y: __property_op__(x, y, operator.div)

    C.__lt__ = lambda x, y: __property_op__(x, y, operator.lt)
    C.__le__ = lambda x, y: __property_op__(x, y, operator.le)
    C.__ge__ = lambda x, y: __property_op__(x, y, operator.ge)
    C.__gt__ = lambda x, y: __property_op__(x, y, operator.gt)

    C.__eq__ = lambda x, y: __property_op__(x, y, operator.eq)
    C.__ne__ = lambda x, y: __property_op__(x, y, operator.ne)

    def __property_neg__(lhs):
        return -(lhs.value())

    C.__neg__ = __property_neg__


for C in [_Base.LengthProperty, _Base.TimeProperty]:
    _add_property_operators(C)
