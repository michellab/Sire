#!/bin/env python

from Sire.Vol import *
from Sire.Maths import *

import copy
import time
import sys

c = CoordGroup(30000)

for i in range(0,10000):
  
  print(i, file=sys.stderr)
  
  print("Copy...", file=sys.stderr)
  c2 = copy.deepcopy(c)
  
  print("Edit...", file=sys.stderr)
  editor = c.edit()
  
  print("Translate...", file=sys.stderr)
  editor.translate( Vector(1,0,0) )
  
  print("Commit...", file=sys.stderr)
  c = editor.commit()
  
  #time.sleep(0.5)

