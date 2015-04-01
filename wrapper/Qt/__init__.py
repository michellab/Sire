#############################
##
## The SireQt module
##
## This holds the wrapping of the small
## number of Qt classes that are necessary
## to use Sire from python
##
## I don't yet know how these wrappers interact
## with PyQt4 or any other Qt4 wrappers...
##
## (C) Christopher Woods
##

# Eventually put in a test here that will prevent this module
# from being loaded if another Qt wrapper library has already
# been loaded

from Sire.Qt._Qt import *
