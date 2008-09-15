#!/bin/env python

import numpy, pylab as py,sys
NNODES = 8
SIZE = 116
GHOST = 0
GSIZE = SIZE*NNODES**(1.0/3)
times = numpy.transpose([[float(i) for i in x.strip().split()] for x in file("output-times").readlines()])
times = dict([(times[1][i],times[0][i]) for i in range(len(times[0]))])

def kinetic(rho, rhovx, rhovy, rhovz):
  return numpy.mean(1.0/2.0 * (rhovx**2 + rhovy**2 + rhovz**2) / rho)

def sigma(rho,rhovx,rhovy,rhovz,bin):
  bins = 100
  rhomin = rho.min()
  rhomax = rho.max()
  width = (rhomax - rhomin)/bins
  if width:
    X = numpy.arange(rhomin,rhomax+width/2,width)
    Y = numpy.zeros(bins+1)
    count = numpy.zeros(bins+1)
    for x in range(len(rho)):
      for y in range(len(rho[0])):
        for z in range(len(rho[0][0])):
          i = int((rho[x,y,z] - rhomin)/width)
          Y[i] += (rhovx[x,y,z]**2 + rhovy[x,y,z]**2 + rhovz[x,y,z]**2)**.5 / rho[x,y,z]
          count[i] += 1
    Y = Y / count
    return (Y,X,rhomin, rhomax, width)
  return None
  

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
    print "%f%%" % (n*100.0/NNODES),
    sys.stdout.flush()
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
    
    print "\r",
  return (rho,rhovx,rhovy,rhovz)


kin = []
for t in range(10,521,10):
  print "t = %d" % t
  rho,rhovx,rhovy,rhovz = load(t)
  py.clf()
  py.pcolor(rho[:][:][33])
  py.axhline(SIZE-2*GHOST)
  py.axvline(SIZE-2*GHOST)
  py.title(times[t])
  py.colorbar()
  py.savefig("analysis-rho-%.8d.png"%t)
  
  print "Calculating sigma",
  sys.stdout.flush()
  Y,X, rhomin, rhomax, width = sigma(rho,rhovx,rhovy,rhovz,100)
  py.clf()
  py.gca()
  py.clim(rhomin,rhomax)
  py.bar(X,numpy.nan_to_num(Y),width=width)
  py.axhline(1)
  py.title("Average velocity of density with %d bins @ t=%f" % (100,times[t]))
  py.xlabel("rho")
  py.ylabel("Average velocity of cells with rho")
  py.savefig("analysis-sigma-%.8d.png" % t)
  print "                      \r",
  
  print "Calculating kinetic energy",
  sys.stdout.flush()
  kin.append(kinetic(rho,rhovx,rhovy,rhovz))
  print "                           \r",
  
  sys.stdout.flush()

py.clf()
py.plot(sorted(times.values())[1:],kin)
py.title("Mean kinetic energy per timestep")
py.savefig("analysis-kinetic.png")