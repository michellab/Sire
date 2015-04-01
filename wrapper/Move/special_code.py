###############################################
#
# This file contains special code to help
# with the wrapping of SireMove classes
#
#

import sys
import pickle

from pyplusplus.module_builder import call_policies

sys.path.append("../AutoGenerate")
from scanheaders import *

def fix_Moves(c):
    
    for decl in c.decls():
        if str(decl).find("SimController") != -1:
            decl.exclude()

special_code = { "SireMove::MovesBase" : fix_Moves,
                 "SireMove::SameMoves" : fix_Moves,
                 "SireMove::WeightedMoves" : fix_Moves
               }
