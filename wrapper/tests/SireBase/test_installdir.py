
from Sire.Base import *

from nose.tools import assert_equal

import os

def test_dirs_exist():
    assert( os.path.exists( getInstallDir() ) )
    assert( os.path.exists( getBinDir() ) )
    assert( os.path.exists( getShareDir() ) )
    assert( os.path.exists( getLibDir() ) )
    assert( os.path.exists( getBundledLibDir() ) )

