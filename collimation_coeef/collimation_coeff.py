#!/bin/env python

import sys, numpy

def sum3d(f,n,theta):
  sum = 0.0
  for x in range(n):
    for y in range(n):
      for z in range(n):
        if all((x,y,z)):
          sum += f(theta,x,y,z)
  return sum

P = lambda theta,x,y,z : 1.0/(numpy.log(2/.1) * (1 + (theta)**2 - (x / numpy.sqrt(x**2+y**2+z**2) * 2 / numpy.pi - 1 )**2))

print "Getting Data"
#datar = []
#rr=range(5,15)
#for i in rr:
#  print "i=%d" % i
#  datar.append(1.0/sum3d(P,i,0.1))

datat = {8:[], 10:[], 6:[]}
rt=numpy.arange(.01,.3,.01)
for r in datat.keys():
  print "r=%d"%r
  for i in rt:
    datat[r].append(1.0/sum3d(P,r,i))

print "Making Graph"
from pylab import *
#clf()
#plot(rr,datar)
#xlabel("r")
#ylabel("K")
#title("Collimation correction with theta=.1")
#savefig("collimation_coeef_r.png")

clf()
for r in datat.keys():
  plot(rt,datat[r],label="R="+str(r))
xlabel("t")
ylabel("K")
legend()
title("Collimation correction with r=8")
savefig("collimation_coeef_t.png")
