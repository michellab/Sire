###############################################
#
# This file contains special code to help
# with the wrapping of SireVol classes
#
#

import sys
import pickle

from pyplusplus.module_builder import call_policies

def fixMB(mb):
    mb.add_declaration_code( "void register_SireStream_functions();" )

    mb.add_registration_code( "register_SireStream_functions();" )
