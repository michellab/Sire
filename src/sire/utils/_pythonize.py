
__all__ = ["pythonize", "pythonize_module", "disable_pythonize"]

_is_pythonize_enabled = True


def _upper_split(data):
    """Thanks to this stackoverflow post for the inspiration
       https://stackoverflow.com/questions/7322028/how-to-replace-uppercase-with-underscore
    """
    buff = ''
    for item in data:
        if item.isupper():
            if buff:
                if buff == "An" or buff == "A":
                    buff = ''
                else:
                    yield buff
                    buff = ''

        buff += item

    yield buff


def pythonize(C, delete_old: bool=True) -> None:
    """Pythonize the passed class. This will rename the functions
       so that they better match the Python style
       (changing mixedCase function names into underscore_case)

       This makes the change with full awareness of the naming
       convention used in sire. This includes changing
       cutGroup to cutgroup and typeName to typename,
       and nSomething to num_something. This ignores functions
       that start with an underscore.

       Note that this changes the API of the class globally.
       All objects created from this class (and to be created)
       will now have the new API. The old API is deleted.

       Args:
        C
            The Class type to be pythonized

        delete_old: bool (defaults to True)
            Whether or not to remove the old function name.

       Returns:
        None
    """
    global _is_pythonize_enabled

    if not _is_pythonize_enabled:
        return

    if type(C) is list:
        for CLS in C:
            pythonize(CLS)
        return

    import re

    for attr in dir(C):
        if attr.startswith("_"):
            continue

        new_attr = attr

        # change typeName into typename
        if attr == "typeName":
            new_attr = "typename"

        # change 'cutGroup' into 'cutgroup'
        new_attr = new_attr.replace("utGroup", "utgroup")

        if new_attr.startswith("asAn"):
            new_attr = new_attr.replace("asAn", "as")
        elif new_attr.startswith("asA"):
            new_attr = new_attr.replace("asA", "as")
        elif new_attr.startswith("isAn"):
            new_attr = new_attr.replace("isAn", "is")
        elif new_attr.startswith("isA"):
            new_attr = new_attr.replace("isA", "is")

        # change nSomething into num_somthing
        m = re.match("^n([A-Z])[a-z]", new_attr)

        if m:
            new_attr = f"num_{m.groups()[0].lower()}{new_attr[2:]}"

        # now change anyCapitalLetter into any_capital_letter
        new_attr = "_".join(_upper_split(new_attr)).lower()

        if new_attr != attr:
            try:
                setattr(C, new_attr, getattr(C, attr))

                if delete_old:
                    delattr(C, attr)
            except Exception:
                # this is a base-class function
                pass

            # need to do this in all of the bases too
            #for B in C.mro():
            #    try:
            #        setattr(B, new_attr, getattr(B, attr))
            #
            #        if delete_old:
            #            delattr(B, attr)
            #    except Exception:
            #        # this is not in the base class
            #        pass


_pythonized = {}


def pythonize_module(MOD):
    """Pythonize all classes in the passed module"""
    global _is_pythonize_enabled

    if not _is_pythonize_enabled:
        return

    global _pythonized

    if MOD in _pythonized:
        return

    import inspect

    for (key, cls) in inspect.getmembers(MOD, inspect.isclass):
        pythonize(cls)

    print(f"PYTHONIZED {MOD}")
    _pythonized[MOD] = 1


def disable_pythonize():
    """Call this function to disable pythonizing of classes. This
       should be called after importing sire.utils in legacy
       scripts to preserve the old API
    """
    global _is_pythonize_enabled

    _is_pythonize_enabled = False
