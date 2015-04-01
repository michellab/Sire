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
    mb.add_declaration_code("#include \"SireID/name.h\"")

    for enum in mb.enums():
        if str(enum).find("SireID::CaseSensitivity") != -1:
            enum.include()
            enum.include_files.append( "SireID/name.h" )


