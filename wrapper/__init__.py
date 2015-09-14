#############################
##
## The Sire python module
##
## This contains the parts of the main Sire program
## that are exposed to Python.
##
## (C) Christopher Woods
##

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
        conn.request("POST", "/cgi-bin/postusagestats.php", params, headers)
        r1 = conn.getresponse()
        #print(r1.status, r1.reason)
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

