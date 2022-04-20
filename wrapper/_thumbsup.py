
__all__ = ["thumbs_up", "disable_thumbs_up",
           "get_thumbs_up_info"]

def _get_deny_filename():
    from pathlib import Path
    return Path.joinpath(Path.home(), ".sire_no_thumbs_up.txt")


def _is_thumbs_up_denied():
    """Return whether or not thumbs up are denied"""
    import os

    if os.path.exists(_get_deny_filename()):
        return True
    else:
        return False


def disable_thumbs_up():
    """Permanently disable thumbs up. This will write a file
       into your home directory called '.sire_no_thumbs_up.txt'.
       If this file exists, then thumbs up are denied and
       no data will be sent to siremol.org.

       The only way to re-enable thumbs up is to remove
       the $HOME/.sire_no_thumbs_up.txt file.
    """
    if not _is_thumbs_up_denied():
        with open(_get_deny_filename(), "w") as FILE:
            FILE.write("Sire's thumbs up is disabled.\n"
                       "Remove this file if you want to re-enable it.")

_thumbs_up_data = None


def get_thumbs_up_info():
    """Return the info that will be sent to siremol.org if you call
       the thumbs_up() function. This will return nothing if you
       have disabled thumbs up using the 'disable_thumbs_up()' function.
    """
    if disable_thumbs_up():
        return None

    global _thumbs_up_data

    if _thumbs_up_data is not None:
        return _thumbs_up_data

    import Sire
    from Sire.Base import CPUID
    import os
    import sys
    import platform

    id = CPUID()

    data = {}

    # get information about the processor
    data["processor"] = id.brand()
    data["vendor"] = id.vendor()
    data["clockspeed"] = id.clockSpeed()
    data["numcores"] = id.numCores()

    # get information about the operating system
    data["platform"] = platform.system()

    if platform.system().startswith("Darwin"):
        data["OS"] = platform.mac_ver()[0]
    elif platform.system().startswith("Linux"):
        ld = platform.linux_distribution()
        data["OS"] = "%s (%s %s)" % (ld[0],ld[1],ld[2])
    elif platform.system().startswith("Windows"):
        ld = platform.win32_ver()
        data["OS"] = "%s (%s %s)" % (ld[0],ld[1],ld[2])
    else:
        data["OS"] = "unknown"

    u = platform.uname()
    data["uname"] = "%s | %s | %s | %s" % (u.system,u.release,u.machine,u.processor)

    # get information about the version of Sire
    data["version"] = Sire.__version__
    data["repository"] = Sire.Config.sire_repository_url
    data["repository_version"] = Sire.Config.sire_repository_version

    # now get information about which Sire app is running
    import sys as _sys
    # get the executable name, but make sure we don't get the path
    # (as it may contain sensitive user information)
    data["executable"] = os.path.basename( _sys.executable )

    # Was Sire was imported as part of BioSimSpace?
    # If so, then rename the executable.
    if "BioSimSpace" in sys.modules:
        data["executable"] = "BioSimSpace"

    _thumbs_up_data = data

    return _thumbs_up_data


_sent_usage_data = None


def thumbs_up():
    """Give Sire a thumbs up! This will send a small amount of data
       to siremol.org to let us know that you like Sire, and what
       operating system and version of Sire you are using.

       (you can see all of the information that would be sent
        by calling "get_thumbs_up_info()")

       You can get more information about what is sent and why
       this is useful to us by visiting https://siremol.org/thumbs_up

       You can permanently disable thumbs_up() on your computer by
       calling "disable_thumbs_up()". This will write a small
       file to your home directory that will tell Sire to not
       allow any more thumbs_up data to be sent.
    """
    global _sent_usage_data

    if _sent_usage_data is not None:
        # silently return - we will only send one thumbs up per run
        return

    data = get_thumbs_up_info()

    if data is None:
        # we aren't allowed to send anything
        return

    try:
        import json
        import os

        import http.client as htc
        import urllib.parse as parse

        params = parse.urlencode({'data' : json.dumps(data)})
        headers = {"Content-type": "application/x-www-form-urlencoded",
                   "Accept": "text/plain"}

        _sent_usage_data = data

        conn = htc.HTTPSConnection("siremol.org")
        conn.request("POST", "/phonehome/postusagestats.php", params, headers)

        # Next time this breaks, remember to uncomment the below lines so that
        # we can inspect the response code and error from the server...

        #r1 = conn.getresponse()
        #print(r1.status, r1.reason)
        #print(r1.read())

    except Exception:
        # something went wrong - just ignore the error
        # and cancel the phone home
        return
