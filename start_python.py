#!/bin/env python

import numpy, pylab as py
NNODES = 8
SIZE = 33
GHOST = 3
GSIZE = SIZE*NNODES**(1.0/3)

def load(i):
  global NNODES,SIZE,GSIZE,GHOST
  ggsize = GSIZE - 2*NNODES**(1/3.0)*GHOST
  rsize = SIZE-2*GHOST
  rho = numpy.zeros((ggsize,ggsize,ggsize))
  rhovx = numpy.zeros((ggsize,ggsize,ggsize))
  rhovy = numpy.zeros((ggsize,ggsize,ggsize))
  rhovz = numpy.zeros((ggsize,ggsize,ggsize))
  offset = numpy.zeros(3)
  isint = lambda x : x == int(x)
  for n in range(NNODES):
    offset[2] = n%2
    if (n%2 == 0 and n!=0): offset[1] += 1 - 2*isint(n/4.0)
    if (n%4 == 0 and n!=0): offset[0] += 1 - 2*isint(n/8.0)
    if (n%8 == 0 and n!=0): offset += 1
    #print offset, n
  
    fd = file("output-%.8d-%.3d"%(i,n))
    a, b,c,d = numpy.transpose([[float(num) for num in line.split()] for line in fd.readlines()])

    a = a.reshape(SIZE,SIZE,SIZE)
    b = b.reshape(SIZE,SIZE,SIZE)
    c = c.reshape(SIZE,SIZE,SIZE)
    d = d.reshape(SIZE,SIZE,SIZE)    
    
    rho[offset[0]*rsize:(offset[0]+1)*rsize,offset[1]*rsize:(offset[1]+1)*rsize,offset[2]*rsize:(offset[2]+1)*rsize] = a[GHOST:SIZE-GHOST,GHOST:SIZE-GHOST,GHOST:SIZE-GHOST]
    rhovx[offset[0]*rsize:(offset[0]+1)*rsize,offset[1]*rsize:(offset[1]+1)*rsize,offset[2]*rsize:(offset[2]+1)*rsize] = b[GHOST:SIZE-GHOST,GHOST:SIZE-GHOST,GHOST:SIZE-GHOST]
    rhovy[offset[0]*rsize:(offset[0]+1)*rsize,offset[1]*rsize:(offset[1]+1)*rsize,offset[2]*rsize:(offset[2]+1)*rsize] = c[GHOST:SIZE-GHOST,GHOST:SIZE-GHOST,GHOST:SIZE-GHOST]
    rhovz[offset[0]*rsize:(offset[0]+1)*rsize,offset[1]*rsize:(offset[1]+1)*rsize,offset[2]*rsize:(offset[2]+1)*rsize] = d[GHOST:SIZE-GHOST,GHOST:SIZE-GHOST,GHOST:SIZE-GHOST]
          
    fd.close()
    # py.clf()
    # py.pcolor(a[:][:][17])
    # py.colorbar()
    # py.savefig("test-%d.png"%n)
  return (rho,rhovx,rhovy,rhovz)

for t in range(5,101,5):
  print "t = %d" % t
  rho,rhovx,rhovy,rhovz = load(t)
  py.clf()
  py.pcolor(rho[:][:][33])
  #py.clim(.9,1.1)
  py.axhline(SIZE-2*GHOST)
  py.axvline(SIZE-2*GHOST)
  # py.axhline(SIZE-3)
  # py.axhline(SIZE+3)
  # py.axvline(SIZE-3)
  # py.axvline(SIZE+3)
  py.colorbar()
  py.savefig("analysis-x-%.8d.png"%t)
  
  py.clf()
  py.pcolor(rho[:][33][:])
  #py.clim(.9,1.1)
  py.colorbar()
  py.axhline(SIZE-2*GHOST)
  py.axvline(SIZE-2*GHOST)
  # py.axhline(SIZE-3)
  # py.axhline(SIZE+3)
  # py.axvline(SIZE-3)
  # py.axvline(SIZE-3)
  py.savefig("analysis-y-%.8d.png"%t)
