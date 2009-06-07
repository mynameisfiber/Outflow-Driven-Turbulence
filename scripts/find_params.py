#!/usr/bin/env python

from __future__ import division

def findparams(l,v):
  return {"S":v/l**(4), "I":v*l**3}
while True:
  print findparams(float(input("l_merge> ")), \
                   float(input("v_merge> ")))
