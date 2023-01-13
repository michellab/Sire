#############################
##
## The Sire Tools module
##
## This contains the python scripts that provide
## useful tools and shortcuts in Sire
##
## (C) Christopher Woods
##

#ensure that the SireQt and SireError libraries are loaded as
#these are vital for the rest of the module
from .. import Qt as _Qt
from .. import Error as _Error
from ..Units import *

import sys

# Create a "readParams" function that reads a text file and returns a dictionary
# of key - value pairs
def readParams( filename ):
    """Read a text file containing key-value pairs and return as a dictionary"""
    try:
        lines = open(filename, "r").readlines()
    except UnicodeDecodeError:
        import codecs
        lines = codecs.open(filename, "r", encoding="UTF-8").readlines()

    params = {}

    for line in lines:
        # format is "key = value", "#" is used for commenting
        s = line.find("=")
        if s != -1:
            key = line[0:s].lstrip().rstrip()
            value = line[s+1:].lstrip().rstrip()

            if key.find("#") != -1:
                #ignore this line, as it must be a weird comment!
                next

            if value.find("#") != -1:
                value = value[0:value.find("#")].lstrip().rstrip()

            params[key] = value

    return params

class Parameter:
    """This class is used to help manage all script-level user-definable
       parameters. You create a new parameter using;

       my_val = Parameter("my value", default_value, "A test variable")

       When you run the script call resolve to set all of the user variables;

       Parameter.push( user_params )

       Now you can get the user-supplied value using

       my_val.val

       When you have finished running the function, call

       Parameter.pop()

       to clear the set of parameters and restore the previous set.

       Note that the "resolveParameters" decorator does this for
       you automatically
    """

    _old_user_params = []
    _user_params = {}
    _all_params = {}

    @staticmethod
    def push(params):
        Parameter._old_user_params.append(Parameter._user_params)
        Parameter._user_params = params

    @staticmethod
    def pop():
        Parameter._user_params = Parameter._old_user_params.pop()

    @staticmethod
    def printAll(verbose=False):
        keys = list(Parameter._all_params.keys())
        keys.sort()

        for key in keys:
            print(Parameter._all_params[key])
            if verbose:
                print("%s\n" % Parameter._all_params[key].description)

    def __init__(self, key, default_value, description):
        """Create a new parameter with specified key, default value
           and variable description"""
        self._key = key
        self._default_value = default_value
        self._desc = " ".join(description.replace("\n", " ").split())
        Parameter._all_params[key] = self

    def __str__(self):
        return "%s = %s" % (self._key, self.val)

    @property
    def description(self):
        return self._desc

    @property
    def val(self):
        if self._key in Parameter._user_params:
            return Parameter._user_params[self._key]
        else:
            return self._default_value

def resolveParameters(func):
    """Decorator that automatically pushes the "params" user-supplied
       variables onto the Parameter stack before the wrapped function
       is called, and pops them off after"""
    def inner(params = {}):

        broken_parameters = []

        print("Using parameters:")
        keys = list(params.keys())
        keys.sort()
        print("===============")
        for key in keys:
            value = str(params[key]).lstrip().rstrip()

            # execute the parameter, so that it is parsed
            try:
                exec("_pvt_val = %s" % value, globals())
            except:
                _pvt_words = value.split()

                _pvt_ok = False

                if len(_pvt_words) == 2:
                    # try the form N * unit
                    try:
                        exec("_pvt_val = %s * %s" % (_pvt_words[0],_pvt_words[1]), globals())
                        _pvt_ok = True
                    except:
                        pass
                else:
                    # maybe this is a string?
                    try:
                        exec("_pvt_val = \"%s\"" % value.replace("\\","\\\\"), globals())
                        _pvt_ok = True
                    except:
                        pass

                if not _pvt_ok:
                    broken_parameters.append( (key, value) )

            params[key] = _pvt_val
            print("%s == %s" % (key,params[key]))

        print("===============")

        if len(broken_parameters) > 0:
            print("\n!!!FATAL!!!\nCould not understand the following parameters!\n")

            for p in broken_parameters:
                print("%s = %s" % p)

            print("\nCANNOT CONTINUE - PROGRAM EXITING")
            sys.exit(-1)

        Parameter.push(params)
        try:
            retval = func()
        except:
            retval = None
            sys.exc_info()[0]
            Parameter.pop()
            raise

        Parameter.pop()

        return retval

    return inner
