#!/bin/env python

from __future__ import division
import numpy, pylab as py,sys, os, re

def load_slice(i,dir,ovars):
  
  rho = numpy.zeros((ovars['gsize'],ovars['gsize'],2))
  rhovx = numpy.zeros((ovars['gsize'],ovars['gsize'],2))
  rhovy = numpy.zeros((ovars['gsize'],ovars['gsize'],2))
  rhovz = numpy.zeros((ovars['gsize'],ovars['gsize'],2))
  offset = numpy.zeros(3)
  isint = lambda x : x == int(x)
  for n in range(ovars['procs']):
    print "%f%%" % ((n+1)*100.0/ovars['procs']),
    sys.stdout.flush()
    offset[0] = n%2
    if (n%2 == 0 and n!=0): offset[1] += 1 - 2*isint(n/4.0)
    if (n%4 == 0 and n!=0): offset[2] += 1 - 2*isint(n/8.0)
    if (n%8 == 0 and n!=0): offset += 1
    #print offset, n

    fd = file(dir+"output-slice-%.8d-%.3d"%(i,n))
    a, b,c,d = numpy.transpose([[float(num) for num in line.split()] for line in fd.readlines()])

    a = a.reshape(ovars['size'],ovars['size'])
    b = b.reshape(ovars['size'],ovars['size'])
    c = c.reshape(ovars['size'],ovars['size'])
    d = d.reshape(ovars['size'],ovars['size'])    
    
    rho[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]] = a
    rhovx[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]] = b
    rhovy[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]] = c
    rhovz[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]] = d
          
    fd.close()
    
    print "\r",
  return (rho,rhovx,rhovy,rhovz)
  
def load_sigma(i,dir,ovars):
  for n in range(ovars['procs']):
    print "%f%%" % ((n+1)*100.0/ovars['procs']),
    fd = file(dir+"output-sigma-%.8d-%.3d"%(i,n))
    new = numpy.array([float(x) for x in fd.readlines()])
    if n == 0:
      sigma = new
    else:
      sigma += new
    fd.close()
    print "\r",
    
  return sigma/ovars['procs']
  
def load_kinetic(dir,ovars):
  for n in range(ovars['procs']):
    print "%f%%" % (n*100.0/ovars['procs']),
    fd = file(dir+"output-kinetic-%.3d"%n)
    new = numpy.array([float(x) for x in fd.readlines()])
    if n == 0:
      kinetic = new
    else:
      kinetic += new
    fd.close()
    print "\r",

  return kinetic/ovars['procs']
  
def loadparams(dir):
  ovars = {}
  try:
    for line in file(dir+"params.f90").readlines():
      vars = re.search("(REAL|INTEGER), [\w]* :: ([\w]*)[ ]?=[ ]?([.\-0-9e]*)", line)
      if vars is not None:
        if vars.group(1) == "REAL":
          ovars[vars.group(2)] = float(vars.group(3))
        elif vars.group(1) == "INTEGER":
          ovars[vars.group(2)] = int(vars.group(3))
        
    ovars['tmerge'] = ovars['op']**(3/7) / (ovars['oImp']**(3/7) * ovars['oSnorm0']**(4/7))
    ovars['lmerge'] = (ovars['oImp'] / (ovars['op'] * ovars['oSnorm0']))**(1/7)
    ovars['vchar'] = ovars['oImp']**(4/7) * ovars['oSnorm0']**(3/7) / ovars['op']**(4/7)
    ovars['size'] = ovars['n'] - 2*ovars['ghost']
    ovars['gsize'] = ovars['size'] * ovars['procs']**(1/3)
    
    return ovars
  except IOError:
    return None
  
def mkdir(dir):
  try:
    os.mkdir(dir)
  except Exception, e:
    pass
    
    
if __name__ == '__main__':
  
  print "Agregating Parameters"
  if len(sys.argv) == 3:
    dir = sys.argv[1] + "/"
  else:
    dir = "./"
  ovars = loadparams(dir)
    
  print "Getting Timing Information"
  timesslice = numpy.transpose([[float(i) for i in x.strip().split()] for x in file(dir+"output-times-slice").readlines()])
  timesslice = dict([(timesslice[1][i],timesslice[0][i]) for i in range(len(timesslice[0]))])
  timesanal = numpy.transpose([[float(i) for i in x.strip().split()] for x in file(dir+"output-times-kinetic").readlines()])
  timesanal = dict([(timesanal[1][i],timesanal[0][i]) for i in range(len(timesanal[0]))])

  print "Rendering %d Slices" % len(timesslice)
  files = timesslice.keys()
  files.sort()
  count = 1
  mkdir(dir+"rho0")
  mkdir(dir+"rho1")
  mkdir(dir+"rhovx0")
  mkdir(dir+"rhovx1")
  for t in files:
    print "t = %d [%f%%]" % (t,count*100.0/len(timesslice)) + " "*10
    rho,rhovx,rhovy,rhovz = load_slice(t,dir,ovars)
    print "\rRendering and Saving" + " "*5 + "\r",
    sys.stdout.flush()
  
    tnorm = timesslice[t]/ovars['tmerge']
    
    py.clf()
    py.pcolor(rho[:,:,0])
    py.axhline(ovars['gsize'])
    py.axvline(ovars['gsize'])
    py.title(r'$ \rho @ %f \cdot t_{merge} $' % tnorm)
    py.colorbar()
    py.savefig(dir+"rho0/analysis-rho0-%.8d.png"%t)
  
    py.clf()
    py.pcolor(rho[:,:,1])
    py.axhline(ovars['gsize'])
    py.axvline(ovars['gsize'])
    py.title(r'$ \rho @ %f \cdot t_{merge}$' % tnorm)
    py.colorbar()
    py.savefig(dir+"rho1/analysis-rho1-%.8d.png"%t)
  
    py.clf()
    py.pcolor(rhovx[:,:,0])
    py.axhline(ovars['gsize'])
    py.axvline(ovars['gsize'])
    py.title(r'$ \rho_x @ %f \cdot t_{merge}$' % tnorm)
    py.colorbar()
    py.savefig(dir+"rhovx0/analysis-rhovx0-%.8d.png"%t)
  
    py.clf()
    py.pcolor(rhovx[:,:,1])
    py.axhline(ovars['gsize'])
    py.axvline(ovars['gsize'])
    py.title(r'$ \rho_x @ %f \cdot t_{merge}$' % tnorm)
    py.colorbar()
    py.savefig(dir+"rhovx1/analysis-rhovx1-%.8d.png"%t)
  
    count += 1
  
  print "Rendering %d Sigma" % len(timesanal) + " "*10
  files = timesanal.keys()
  files.sort()
  count = 1
  mkdir(dir+"sigma")
  for t in files:
    print "t = %d [%f%%]" % (t,count*100.0/len(timesanal)) + " "*10
  
    tnorm = timesanal[int(t)] / ovars['tmerge']
    sigma = load_sigma(t,dir,ovars)
    py.clf()
    py.bar(range(200),numpy.nan_to_num(sigma)[:200])
    py.axhline(1)
    py.title(r'$ \rho(\|v\|) @ t=%f \cdot t_{merge} $' % tnorm)
    py.xlabel(r'$ \frac{<\|v\|>}{v_{char}} $')
    py.ylabel(r'$ <\rho>  $')
    py.savefig(dir+"sigma/analysis-sigma-%.8d.png" % t)
  
    count += 1
  
  
  print "Rendering Kinetic"
  atimes = timesanal.values()
  atimes.sort()
  atimes = numpy.array(atimes)
  kinetic = load_kinetic(dir,ovars)

  # find the timestep that corresponds to merging plus 5%
  mergestep = 0
  while atimes[mergestep] < ovars['tmerge']:
    mergestep += 1
 # mergestep += (len(atimes) - mergestep)
    
  #fit post-merge saturation to horizontal line:
  kinsat = kinetic[mergestep:].mean()

  mkdir(dir+"kinetic")
  py.clf()
  py.scatter(atimes/ovars['tmerge'],kinetic)
  py.axvline(1)
  py.hlines(kinsat,1,atimes[-1]/ovars['tmerge'])
  py.xlim(xmin=0)
  py.title("Mean kinetic energy")
  py.xlabel("time ($t_{merge}$)")
  py.ylabel("<Kinetic>")
  py.savefig(dir+"kinetic/analysis-kinetic.png")
  
