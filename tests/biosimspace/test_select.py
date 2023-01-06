
# We keep this here as this is a test that
# the legacy API still works when mixed_api
# mode is enabled
import sire as sr
sr.use_mixed_api()

import sire.legacy.Mol

def test_selector_wrap(ala_mols):
    mols = ala_mols

    # These will fail if importing Sire under a mixed
    # api and then importing sire.legacy modules causes
    # an issue
    print(mols["atoms within 2 of element C"])

    res = mols["resnum 1"]
