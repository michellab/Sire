#############################
##
## The SireBase module
##
## This contains some base classes
## that don't fit anywhere else.
##
## This module may be removed
## as I find homes for some of
## these classes
##
## (C) Christopher Woods
##

from Sire.Base._Base import *

_wrap_functions = []

_base_wrap = wrap

def wrap(value):
    # First, try to wrap the python concrete classes
    if isinstance(value, bool):
        return BooleanProperty(value)

    elif isinstance(value, int) or isinstance(value, float):
        return NumberProperty(value)

    elif isinstance(value, str):
        return StringProperty(value)

    for func in _wrap_functions:
        try:
            return func(value)
        except:
            pass

    try:
        return _base_wrap(value)
    except:
        pass

    # if this is a dictionary, then wrap as a Properties object
    if type(value) is dict:
        p = Properties()

        for key, value in value.items():
            p[key] = value

        return p
    else:
        return PropertyList(value)

_original_wrap = wrap

def _add_wrap_function(func):
    _wrap_functions.append(func)
    return _original_wrap

# cludgy quick fix for an anaconda install
_getBundledLibDir = getBundledLibDir
def getBundledLibDir():
    try:
        return _getBundledLibDir()
    except:
        return "%s/lib" % getInstallDir()

def __set_property__(obj, key, property):
    try:
        return obj.__setProperty__(key, property)
    except Exception as e:
        if e.__class__.__name__ == "ArgumentError":
            return obj.__setProperty__(key, wrap(property))
        else:
            raise e

def __getitem__(props, i):
    try:
        return props.__orig_getitem__(i)
    except Exception as e:
        if e.__class__.__name__ == "ArgumentError":
            key = props.propertyKeys()[i]
            val = props[key]
            return (key, val)
        else:
            raise e

def __properties_values__(props):
    vals = []

    for key in props.propertyKeys():
        vals.append(props.property(key))

    return vals

def __properties_items__(props):
    items = []

    for key in props.propertyKeys():
        items.append((key, props.property(key)))

    return items

Properties.__setProperty__ = Properties.setProperty
Properties.setProperty = __set_property__
Properties.__orig_getitem__ = Properties.__getitem__
Properties.__getitem__ = __getitem__
Properties.__setitem__ = __set_property__
Properties.keys = Properties.propertyKeys
Properties.values = __properties_values__
Properties.items = __properties_items__
