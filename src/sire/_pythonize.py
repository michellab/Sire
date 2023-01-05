__all__ = ["use_old_api", "use_new_api", "use_mixed_api"]


def _upper_split(data):
    """
    Thanks to this stackoverflow post for the inspiration
    https://stackoverflow.com/questions/7322028/how-to-replace-uppercase-with-underscore
    """
    buff = ""
    for item in data:
        if item.isupper():
            if buff:
                if buff == "An" or buff == "A":
                    buff = ""
                else:
                    yield buff
                    buff = ""

        buff += item

    yield buff


def _pythonize(C, delete_old: bool = True) -> None:
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
    if type(C) is list:
        for CLS in C:
            _pythonize(CLS)
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

        # change 'ID()' into 'id()'
        if new_attr == "ID":
            new_attr = "id"
        elif new_attr == "IDs":
            new_attr = "ids"

        # change all caps into lowercase
        if new_attr.isupper():
            new_attr = new_attr.lower()

        # change "MCSmatches" into "Mcs_matches" (it will then be
        # converted to _mcs_matches by the code below)
        new_attr = new_attr.replace("MCSmatches", "Mcs_matches")

        # change "MCS" into "Mcs" (it will then be converted to _mcs by
        # the code below)
        new_attr = new_attr.replace("MCS", "Mcs")

        # change "MC" into "Mc" (it will be converted to _mc by the code below)
        new_attr = new_attr.replace("MC", "Mc")

        # change "aaBox" into "aabox"
        new_attr = new_attr.replace("aaBox", "aabox")

        # change "CONECT" to "Conect"
        new_attr = new_attr.replace("CONECT", "Conect")

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
            # for B in C.mro():
            #    try:
            #        setattr(B, new_attr, getattr(B, attr))
            #
            #        if delete_old:
            #            delattr(B, attr)
            #    except Exception:
            #        # this is not in the base class
            #        pass


def _pythonize_modules(modules, delete_old: bool = True):
    """Pythonize all classes in the passed module"""

    for MOD in modules:
        import inspect

        try:
            for (key, cls) in inspect.getmembers(MOD, inspect.isclass):
                _pythonize(cls, delete_old=delete_old)
        except Exception as e:
            print(e)
            print(f"Failed to pythonize {MOD}")


_is_using_old_api = None
_is_using_new_api = None


def use_mixed_api(support_old_module_names: bool = False):
    """Load Sire using both the new (python-style) and old APIs. This
    is useful for migrating old scripts as a temporary porting option.
    You can start writing functions using the new API, safe in the
    knowledge that the old API functions will still work.

    Do aim to finish your port though, else you will forever have
    a duplicated API (e.g. have both X.nAtoms() and X.num_atoms() etc.)
    """
    global _is_using_new_api, _is_using_old_api

    if _is_using_old_api and _is_using_new_api:
        # don't need to do this twice
        return

    if _is_using_old_api or _is_using_new_api:
        msg = (
            "Cannot import sire using the mixed API as either the old "
            "or new APIs have already been activated."
        )
        print(msg)

        raise ImportError(msg)

    # First, bring in the old API
    if support_old_module_names:
        print("Loading Sire with support for old module names.")
        print(
            "Note that this can cause problems with classes importing twice."
        )
        use_old_api()
    else:
        _is_using_old_api = True

    # Now, bring in the new API
    _is_using_new_api = True

    # call Pythonize on all of the new modules
    from . import (
        move,
        io,
        system,
        squire,
        mm,
        ff,
        mol,
        analysis,
        base,
        cas,
        cluster,
        error,
        id,
        maths,
        qt,
        stream,
        units,
        vol,
    )

    _pythonize_modules(
        [
            analysis._Analysis,
            base._Base,
            cas._CAS,
            cluster._Cluster,
            error._Error,
            ff._FF,
            id._ID,
            io._IO,
            maths._Maths,
            mm._MM,
            mol._Mol,
            move._Move,
            qt._Qt,
            squire._Squire,
            stream._Stream,
            system._System,
            units._Units,
            vol._Vol,
        ],
        delete_old=False,
    )


def use_new_api():
    """Load Sire using the new (python-style) API. This will be called
    automatically when you load any of the new Python modules, so you
    shouldn't need to call this yourself.
    """
    global _is_using_new_api, _is_using_old_api

    if _is_using_new_api:
        # already done
        return

    if _is_using_old_api:
        msg = (
            "Cannot import sire using the new API as the old API has "
            "already been activated. Both APIs cannot be active at "
            "the same time."
        )
        print(msg)

        raise ImportError(msg)

    _is_using_new_api = True

    # call Pythonize on all of the new modules
    from .legacy import (
        Move,
        IO,
        System,
        Squire,
        MM,
        FF,
        Mol,
        Analysis,
        Base,
        CAS,
        Cluster,
        Error,
        ID,
        Maths,
        Qt,
        Stream,
        Units,
        Vol,
    )

    _pythonize_modules(
        [
            Analysis._Analysis,
            Base._Base,
            CAS._CAS,
            Cluster._Cluster,
            Error._Error,
            FF._FF,
            ID._ID,
            IO._IO,
            Maths._Maths,
            MM._MM,
            Mol._Mol,
            Move._Move,
            Qt._Qt,
            Squire._Squire,
            Stream._Stream,
            System._System,
            Units._Units,
            Vol._Vol,
        ]
    )


def use_old_api():
    """Load Sire using the old (C++-style) API. This is for
    compatibility reasons for old code only. This should
    not be used with new code
    """
    global _is_using_old_api, _is_using_new_api

    if _is_using_new_api:
        raise ImportError(
            "Cannot import Sire using the old API as the new API has "
            "already been activated. Both APIs cannot be active at "
            "the same time."
        )

    if _is_using_old_api:
        # already active
        return

    _is_using_old_api = True

    from . import legacy

    # set up the meta-importer with these modules - thanks to this post
    # for all of the info
    # https://dev.to/dangerontheranger/dependency-injection-with-import-hooks-in-python-3-5hap
    import importlib.abc
    import importlib.machinery
    import sys
    import types

    class DependencyInjectorLoader(importlib.abc.Loader):
        def __init__(self):
            self._services = {}
            self._dummy_module = types.ModuleType("Sire")
            self._dummy_module.__path__ = []

        def provide(self, service_name, module):
            """Register a service as provided via the given module
            A service is any Python object in this context - an imported
            module, a class, etc."""
            self._services[service_name] = module

        def provides(self, fullname):
            if fullname in self._services:
                return True
            else:
                return False

        def create_module(self, spec):
            """Create the given module from the supplied module spec
            Under the hood, this module returns a service or a dummy module,
            depending on whether Python is still importing one of the names
            listed in _COMMON_PREFIX.
            """
            service_name = spec.name

            if service_name not in self._services:
                # return our dummy module since at this point we're loading
                # *something* along the lines of "myapp.virtual" that's not
                # a service
                return self._dummy_module

            module = self._services[service_name]
            return module

        def exec_module(self, module):
            """Execute the given module in its own namespace
            This method is required to be present by importlib.abc.Loader,
            but since we know our module object is already fully-formed,
            this method merely no-ops.
            """
            pass

    class DependencyInjectorFinder(importlib.abc.MetaPathFinder):
        def __init__(self, loader):
            # we'll write the loader in a minute, hang tight
            self._loader = loader

        def find_spec(self, fullname, path, target=None):
            """Attempt to locate the requested module
            fullname is the fully-qualified name of the module,
            path is set to __path__ for sub-modules/packages, or None
            otherwise target can be a module object, but is unused in this
            example.
            """
            if self._loader.provides(fullname):
                return self._gen_spec(fullname)

        def _gen_spec(self, fullname):
            spec = importlib.machinery.ModuleSpec(fullname, self._loader)
            return spec

    class DependencyInjector:
        """
        Convenience wrapper for DependencyInjectorLoader and
        DependencyInjectorFinder.
        """

        def __init__(self):
            self._loader = DependencyInjectorLoader()
            self._finder = DependencyInjectorFinder(self._loader)

        def install(self):
            sys.meta_path.append(self._finder)

        def provide(self, service_name, module):
            self._loader.provide(service_name, module)

    injector = DependencyInjector()

    # Use 'legacy' as a stand-in for Sire
    injector.provide("Sire", legacy)
    injector.install()
