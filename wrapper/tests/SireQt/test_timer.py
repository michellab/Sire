
from Sire.Qt import *

import time

def test_timer(verbose=True):
    t = QElapsedTimer()

    t.start()
    ns = t.nsecsElapsed()

    if verbose:
        print("Minimum time is %s ms" % (0.000001*ns))

    assert( ns < 1000000 )

    t.start()

    time.sleep(1)

    ns = t.nsecsElapsed()

    if verbose:
        print("Slept for 1 s == %s ms" % (0.000001*ns))

    assert( ns >= 998000000 )
    assert( ns <= 1002000000 )

if __name__ == "__main__":
    test_timer(True)

