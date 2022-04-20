

import sys

sys.path.insert(0, "../AutoGenerate")

from create_wrappers import call_with_released_gil

def fix_MoleculeParser(c):
    call_with_released_gil(c, "load")

special_code = { "SireIO::MoleculeParser" : fix_MoleculeParser }
