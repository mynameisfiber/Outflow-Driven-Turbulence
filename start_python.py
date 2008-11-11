#!/bin/env python

import numpy, pylab as py,sys, os
NNODES = 8
SIZE = 14
GSIZE = SIZE*NNODES**(1.0/3)
timesslice = numpy.transpose([[float(i) for i in x.strip().split()] for x in file("output-times-slice").readlines()])
timesslice = dict([(timesslice[1][i],timesslice[0][i]) for i in range(len(timesslice[0]))])
timesanal = numpy.transpose([[float(i) for i in x.strip().split()] for x in file("output-times-kinetic").readlines()])
timesanal = dict([(timesanal[1][i],timesanal[0][i]) for i in range(len(timesanal[0]))])

def load_slice(i):
  global NNODES,SIZE,GSIZE
  
  rho = numpy.zeros((GSIZE,GSIZE,2))
  rhovx = numpy.zeros((GSIZE,GSIZE,2))
  rhovy = numpy.zeros((GSIZE,GSIZE,2))
  rhovz = numpy.zeros((GSIZE,GSIZE,2))
  offset = numpy.zeros(3)
  isint = lambda x : x == int(x)
  for n in range(NNODES):
    print "%f%%" % ((n+1)*100.0/NNODES),
    sys.stdout.flush()
    offset[0] = n%2
    if (n%2 == 0 and n!=0): offset[1] += 1 - 2*isint(n/4.0)
    if (n%4 == 0 and n!=0): offset[2] += 1 - 2*isint(n/8.0)
    if (n%8 == 0 and n!=0): offset += 1
    #print offset, n

    fd = file("output-slice-%.8d-%.3d"%(i,n))
    a, b,c,d = numpy.transpose([[float(num) for num in line.split()] for line in fd.readlines()])

    a = a.reshape(SIZE,SIZE)
    b = b.reshape(SIZE,SIZE)
    c = c.reshape(SIZE,SIZE)
    d = d.reshape(SIZE,SIZE)    
    
    rho[offset[0]*SIZE:(offset[0]+1)*SIZE,offset[1]*SIZE:(offset[1]+1)*SIZE,offset[2]] = a
    rhovx[offset[0]*SIZE:(offset[0]+1)*SIZE,offset[1]*SIZE:(offset[1]+1)*SIZE,offset[2]] = b
    rhovy[offset[0]*SIZE:(offset[0]+1)*SIZE,offset[1]*SIZE:(offset[1]+1)*SIZE,offset[2]] = c
    rhovz[offset[0]*SIZE:(offset[0]+1)*SIZE,offset[1]*SIZE:(offset[1]+1)*SIZE,offset[2]] = d
          
    fd.close()
    
    print "\r",
  return (rho,rhovx,rhovy,rhovz)
  
def load_sigma(i):
  global NNODES
  
  for n in range(NNODES):
    print "%f%%" % ((n+1)*100.0/NNODES),
    fd = file("output-sigma-%.8d-%.3d"%(i,n))
    new = numpy.array([float(x) for x in fd.readlines()])
    if n == 0:
      sigma = new
    else:
      sigma += new
    fd.close()
    print "\r",
    
  return sigma/NNODES
  
def load_kinetic():
  global NNODES

  for n in range(NNODES):
    print "%f%%" % (n*100.0/NNODES),
    fd = file("output-kinetic-%.3d"%n)
    new = numpy.array([float(x) for x in fd.readlines()])
    if n == 0:
      kinetic = new
    else:
      kinetic += new
    fd.close()
    print "\r",

  return kinetic/NNODES
  
def mkdir(dir):
  try:
    os.mkdir(dir)
  except Exception, e:
    pass


print "Rendering %d Slices" % len(timesslice)
files = timesslice.keys()
files.sort()
count = 1
for t in files:
  print "t = %d [%f%%]" % (t,count*100.0/len(timesslice)) + " "*10
  rho,rhovx,rhovy,rhovz = load_slice(t)
  print "\rRendering and Saving" + " "*5 + "\r",
  sys.stdout.flush()
  
  mkdir("rho0")
  py.clf()
  py.pcolor(rho[:,:,0])
  py.axhline(SIZE)
  py.axvline(SIZE)
  py.title(r'$ \rho @ %f $' % timesslice[t])
  py.colorbar()
  py.savefig("rho0/analysis-rho0-%.8d.png"%t)
  
  mkdir("rho1")
  py.clf()
  py.pcolor(rho[:,:,1])
  py.axhline(SIZE)
  py.axvline(SIZE)
  py.title(r'$ \rho @ %f $' % timesslice[t])
  py.colorbar()
  py.savefig("rho1/analysis-rho1-%.8d.png"%t)
  
  mkdir("rhovx0")
  py.clf()
  py.pcolor(rhovx[:,:,0])
  py.axhline(SIZE)
  py.axvline(SIZE)
  py.title(r'$ \rho_x @ %f $' % timesslice[t])
  py.colorbar()
  py.savefig("rhovx0/analysis-rhovx0-%.8d.png"%t)
  
  mkdir("rhovx1")
  py.clf()
  py.pcolor(rhovx[:,:,1])
  py.axhline(SIZE)
  py.axvline(SIZE)
  py.title(r'$ \rho_x @ %f $' % timesslice[t])
  py.colorbar()
  py.savefig("rhovx1/analysis-rhovx1-%.8d.png"%t)
  
  count += 1
  
print "Rendering %d Sigma" % len(timesanal) + " "*10
files = timesanal.keys()
files.sort()
count = 1
for t in files:
  print "t = %d [%f%%]" % (t,count*100.0/len(timesanal)) + " "*10
  
  sigma = load_sigma(t)
  mkdir("sigma")
  py.clf()
  py.bar(range(200),numpy.nan_to_num(sigma)[:200])
  py.axhline(1)
  py.title(r'$ \rho(\|v\|) @ t=%f $' % (timesanal[t]))
  py.xlabel(r'$ \frac{<\|v\|>}{v_{char}} $')
  py.ylabel(r'$ <\rho>  $')
  py.savefig("sigma/analysis-sigma-%.8d.png" % t)
  
  count += 1
  
print "Rendering Kinetic"
atimes = timesanal.values()
atimes.sort()
kinetic = load_kinetic()
mkdir("kinetic")
py.clf()
py.plot(atimes,kinetic)
py.title("Mean kinetic energy")
py.xlabel("time (s)")
py.ylabel("<Kinetic>")
py.savefig("kinetic/analysis-kinetic.png")
  
