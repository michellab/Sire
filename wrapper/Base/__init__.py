"""
.. currentmodule:: Sire.Base

Classes
=======

.. autosummary::
    :toctree: generated/

    Property
    Propeties
    PropertyList
    PropertyMap
    PropertyName
    Range
    SimpleRange
    TempDir
    TimeProperty
    TrimString
    UpperCaseString
    VariantProperty
    Version

Functions
=========

.. autosummary::
    :toctree: generated/

    findExe
    getBinDir
    getBundledLibDir
    getInstallDir
    getLibDir
    getReleaseVersion
    getRepositoryBranch
    getRepositoryURL
    getRepositoryVersion
    getRepositoryVersionIsClean
    getShareDir
    getSireDir
    increment
    wrap
"""

from Sire.Base._Base import *

__all__ = ["Property", "Properties", "PropertyList", "PropertyMap",
           "PropertyName", "Range", "SimpleRange", "TempDir",
           "TimeProperty", "TrimString", "UpperCaseString",
           "VariantProperty", "Version", "findExe",
           "getBinDir", "getBundledLibDir", "getInstallDir",
           "getLibDir", "getReleaseVersion", "getRepositoryBranch",
           "getRepositoryURL", "getRepositoryVersion",
           "getRepositoryVersionIsClean", "getShareDir",
           "getSireDir", "increment", "wrap"]

_wrap_functions = []

_base_wrap = wrap

def wrap(value):
    """Wrap the passed value into a :class:`~Sire.Base.Property`
       object. This works recursively, wrapping all items in
       a container, such that the returned value is derived
       from :class:`~Sire.Base.Property` and can be passed to
       the C++ code in Sire. Note that you normally don't
       need to call this yourself, as wrapping is handled
       automatically.
    """
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
