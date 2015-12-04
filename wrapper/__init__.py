#############################
##
## The Sire python module
##
## This contains the parts of the main Sire program
## that are exposed to Python.
##
## (C) Christopher Woods
##

_module_to_package = {}

def _install_package(name, package_registry):
    """Internal function used to install the module
       called 'name', using the passed 'package_registry'
       to find the package that contains the package
       that contains this module"""

    # get the directory containing the python executable,
    # we will assume this will also contain 'conda', 'pip'
    # or 'easy_install'
    from os.path import realpath, dirname
    from os import system
    from sys import executable

    binpath = dirname( realpath(executable) )

    # ensure that we have the root package name
    try:
        package = name.split(".")[0]
    except:
        package = name

    if package in package_registry:
        package = package_registry[name]

    try:
        print("\nTrying to install %s from package %s using %s/conda...\n" \
                  % (name,package,binpath))
        ok = system("%s/conda install %s" % (binpath,package))

        if ok == 0:
            # installed ok
            return
    except:
        pass

    try:
        print("\nTrying to install %s from package %s using %s/pip...\n" \
                  % (name,package,binpath))
        ok = system("%s/pip install %s" % (binpath,package))

        if ok == 0:
            # installed ok
            return
    except:
        pass

    try:
        print("\nTrying to install %s from package %s using %s/easy_install...\n" \
                  % (name,package,binpath))
        ok = system("%s/easy_install %s" % (binpath,package))

        if ok == 0:
            # installed ok
            return
    except:
        pass

    print("\nWARNING: Unable to install '%s' from package '%s'\n" \
               % (name,package))

    return

def try_import(name, package_registry=_module_to_package):
    """Try to import the module called 'name', returning
       the loaded module as an argument. If the module
       is not available, then it looks up the name of 
       the package to install using "package_registry"
       (or if this is not available, using just the name
       of the module). This will then be installed using
       "conda", then "pip" then "easy_install" (first one
       that works will return).

       For example, use this via

       sys = try_import("sys")
       mdtraj = try_import("mdtraj")

       Note that you can also rename modules, e.g. by using

       md = try_import("mdtraj")

       Note that you should use try_import_from if you
       want only specific symbols, e.g.

       (argv, stdout) = try_import_from("sys", ["argv","stdout"])
    """

    try:
        mod = __import__(name)
        return mod
    except:
        pass

    if not (package_registry is None):
        _install_package(name, package_registry)
        return try_import(name, package_registry=None)

    raise ImportError("Failed to install module %s" % name)

def try_import_from(name, fromlist, package_registry=_module_to_package):
    """Try to import from the module called 'name' the passed symbol
       (or list of symbols) contained in 'fromlist', returning
       the symbol (or list of symbols).

       If the module cannot be loaded, then the package containing
       the module is looked up in 'module_to_package' (or just guessed
       from the name if it does not exist in 'module_to_package'. 
       An attempt is made to load the package, using first conda,
       then pip, then easy_install.

       Example usage:

       Mol = try_import_from("Sire", "Mol")
       (argv,stdout = try_import_from("sys", ["argv", "stdout"])
       mapReduce = try_import_from("scoop.Futures", "mapReduce")
       ut = try_import_from("mdtraj", "utils")
    """

    if isinstance(fromlist, str):
        # we are importing only a single module - put
        # this string into a list for the user
        fromlist = [fromlist]

    try:
        nsyms = len(fromlist)
    except:
        return try_import(name, package_registry)

    if nsyms == 0:
        # just import the entire module
        return try_import(name, package_registry)

    is_loaded = False

    try:
        mod = __import__(name, globals(), locals(), fromlist)
        is_loaded = True
    except:
        is_loaded = False

    if not is_loaded:
        if not (package_registry is None):
            _install_package(name, package_registry)
            return try_import_from(name, fromlist, package_registry=None)
        else:
            raise ImportError("Failed to install module '%s'" % name)

    if nsyms == 1:
        try:
            return getattr(mod, fromlist[0])
        except:
            raise ImportError("Cannot find the symbol '%s' in module '%s'" \
                     % (fromlist[0],name) )
    else:
        ret = []
        missing_symbols = []

        for sym in fromlist:
            try:
                ret.append( getattr(mod, sym) )
            except:
                missing_symbols.append(sym)

        if len(missing_symbols) > 0:
            raise ImportError("Cannot find the following symbols in module '%s' : [ %s ]" \
                      % (name, ", ".join(missing_symbols)))

        return ret

#ensure that the SireQt and SireError libraries are loaded as
#these are vital for the rest of the module
import Sire.Qt
import Sire.Error
import Sire.Config

__version__ = Sire.Config.__version__

sent_usage_data = None

def _getOSInfo():
    import platform as _pf

    data = {}
    data["platform"] = _pf.system()

    if _pf.system().startswith("Darwin"):
        data["OS"] = _pf.mac_ver()[0]
    elif _pf.system().startswith("Linux"):
        ld = _pf.linux_distribution()
        data["OS"] = "%s (%s %s)" % (ld[0],ld[1],ld[2])
    else:
        data["OS"] = "unknown"

    u = _pf.uname()
    data["uname"] = "%s | %s | %s | %s" % (u.system,u.release,u.machine,u.processor)

    data["OS"] = "%s : %s"

# Now try to upload usage data to siremol.org
def _uploadUsageData():
    try:
        global sent_usage_data

        if not sent_usage_data is None:
            # don't send data twice
            return

        import os as _os
    
        if "SIRE_DONT_PHONEHOME" in _os.environ:
           # respect user wish to not phone home
           return

        from Sire.Base import CPUID as _CPUID

        id = _CPUID()

        data = {}

        # get information about the processor
        data["processor"] = id.brand()
        data["vendor"] = id.vendor()
        data["clockspeed"] = id.clockSpeed()
        data["numcores"] = id.numCores()

        # get information about the operating system
        import platform as _pf

        data["platform"] = _pf.system()

        if _pf.system().startswith("Darwin"):
            data["OS"] = _pf.mac_ver()[0]
        elif _pf.system().startswith("Linux"):
            ld = _pf.linux_distribution()
            data["OS"] = "%s (%s %s)" % (ld[0],ld[1],ld[2])
        else:
            data["OS"] = "unknown"

        u = _pf.uname()
        data["uname"] = "%s | %s | %s | %s" % (u.system,u.release,u.machine,u.processor)

        # get information about the version of Sire
        data["version"] = Sire.__version__
        data["repository"] = Sire.Config.sire_repository_url
        data["repository_version"] = Sire.Config.sire_repository_version

        # now get information about which Sire app is running
        import sys as _sys
        # get the executable name, but make sure we don't get the path
        # (as it may contain sensitive user information)
        data["executable"] = _os.path.basename( _sys.executable )

        import json as _json

        import http.client as _htc
        import urllib.parse as _parse

        params = _parse.urlencode({'data' : _json.dumps(data)})
        headers = {"Content-type": "application/x-www-form-urlencoded", 
                   "Accept": "text/plain"}

        #print(_parse.urlencode({'data' : _json.dumps(data)}))
        #print(headers)

        sent_usage_data = data

        conn = _htc.HTTPConnection("siremol.org")
        conn.request("POST", "/phonehome/postusagestats.php", params, headers)
        #r1 = conn.getresponse()
        #print(r1.status, r1.reason)
        #print(r1.read())
        
    except:
        # something went wrong - just ignore the error
        # and cancel the phone home
        return

sent_usage_data = None

if not sent_usage_data:
    import threading as _threading

    _thread = _threading.Thread(target=_uploadUsageData)
    _thread.daemon = True
    _thread.start()
