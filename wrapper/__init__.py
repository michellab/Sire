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
