#!/bin/env python

from __future__ import division
import numpy, pylab as py,sys, os, re
from itertools import izip
reverse_enumerate = lambda l: izip(xrange(len(l)-1, -1, -1), reversed(l))


def load_slice(i,dir,ovars):
  rho = numpy.zeros((ovars['gsize'],ovars['gsize'],2))
  rhovx = numpy.zeros((ovars['gsize'],ovars['gsize'],2))
  rhovy = numpy.zeros((ovars['gsize'],ovars['gsize'],2))
  rhovz = numpy.zeros((ovars['gsize'],ovars['gsize'],2))
  offset = numpy.zeros(333)
  isint = lambda x : x == int(x)
  for n in range(ovars['procs']):
    print "%0.2f%%" % ((n+1)*100.0/ovars['procs']) + " "*10,
    sys.stdout.flush()
    offset[0] = n%2
    if (n%2 == 0 and n!=0): offset[1] += 1 - 2*isint(n/4.0)
    if (n%4 == 0 and n!=0): offset[2] += 1 - 2*isint(n/8.0)
    if (n%8 == 0 and n!=0): offset += 1
    #print offset, n

    #fd = file(dir+"output-slice-%.8d-%.3d"%(i,n))
    #a, b,c,d = numpy.transpose([[float(num) for num in line.split()] for line in fd.readlines()])
    fd = open(dir+"output-slice-%.8d-%.3d"%(i,n))
    a,b,c,d = numpy.fromfile(file=fd,dtype=numpy.float32).reshape((4,ovars['size'],ovars['size']))

    rho[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]] = a
    rhovx[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]] = b
    rhovy[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]] = c
    rhovz[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]] = d
          
    fd.close()
    
    print "\r",
  return (rho,rhovx,rhovy,rhovz)
  
def load_cube(i,dir,ovars):

  rho = numpy.zeros((ovars['gsize'],ovars['gsize'],ovars['gsize']))
  rhovx = numpy.zeros((ovars['gsize'],ovars['gsize'],ovars['gsize']))
  rhovy = numpy.zeros((ovars['gsize'],ovars['gsize'],ovars['gsize']))
  rhovz = numpy.zeros((ovars['gsize'],ovars['gsize'],ovars['gsize']))
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

#    try:
#      fd = file(dir+"output-cube-%.8d-%.3d"%(i,n))
#      a, b,c,d = numpy.transpose([[float(num) for num in line.split()] for line in fd.readlines()])
#    except Exception, e:
#      a,b,c,d =  numpy.zeros((4,ovars['size']**3))


    fd = open(dir+"output-cube-%.8d-%.3d"%(i,n))
    a,b,c,d = numpy.fromfile(file=fd, dtype=numpy.float32).reshape((4,ovars['size'],ovars['size'],ovars['size']))

#    a = a.reshape(ovars['size'],ovars['size'],ovars['size'])
#    b = b.reshape(ovars['size'],ovars['size'],ovars['size'])
#    c = c.reshape(ovars['size'],ovars['size'],ovars['size'])
#    d = d.reshape(ovars['size'],ovars['size'],ovars['size'])
    
    rho[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]*ovars['size']:(offset[2]+1)*ovars['size']] = a
    rhovx[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]*ovars['size']:(offset[2]+1)*ovars['size']] = b
    rhovy[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]*ovars['size']:(offset[2]+1)*ovars['size']] = c
    rhovz[offset[0]*ovars['size']:(offset[0]+1)*ovars['size'],offset[1]*ovars['size']:(offset[1]+1)*ovars['size'],offset[2]*ovars['size']:(offset[2]+1)*ovars['size']] = d

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
    print "%0.2f%%" % (n*100.0/ovars['procs'])+" "*10,
    fd = file(dir+"output-kinetic-%.3d"%n)
    new = numpy.array([float(x) for x in fd.readlines()])
    if n == 0:
      kinetic = new
    else:
      kinetic += new
    fd.close()
    print "\r",
  return kinetic/ovars['procs']
  
def load_compressional(dir,ovars):
  for n in range(ovars['procs']):
    print "%f%%" % (n*100.0/ovars['procs']),
    fd = file(dir+"output-compress-%.3d"%n)
    new = numpy.array([float(x) for x in fd.readlines()])
    if n == 0:
      compress = new
    else:
      compress += new
    fd.close()
    print "\r",
  return compress/ovars['procs']
  
def load_mass(dir,ovars):
  rho = []
  for n in range(ovars['procs']):
    print "%f%%" % (n*100.0/ovars['procs']),
    fd = file(dir+"output-totalrho-%.3d"%n)
    rho.append(numpy.array([float(x) for x in fd.readlines()]))
    fd.close()
    print "\r",
  return numpy.array(rho).sum(0)
  return rho
  
def loadparams(dir):
  ovars = {}
  try:
    for line in file(dir+"params.f90").readlines():
      vars = re.search("[^!]+[ \t]*(REAL|INTEGER), [\w]* :: ([\w]*)[ ]?=[ ]?([.\-0-9e]*)", line)
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
    print ovars

    return ovars
  except IOError:
    return None
  
def mkdir(outdir):
  try:
    os.mkdir(outdir)
  except Exception, e:
    pass
    
    
if __name__ == '__main__':
  
  print "Agregating Parameters"
  try:
    dir = sys.argv[1] + "/"
  except Exception,e:
    dir =  "./"
  try:
    outdir = sys.argv[2] + "/"
  except Exception,e:
    outdir = "./"
  try:
    methods = int(sys.argv[3])
  except Exception,e:
    methods = 0
  ovars = loadparams(dir)
  print "Dir: %s"%dir
  print "OutDir: %s"%outdir
  
  print "Getting Timing Information"
  try:
    timesslice = numpy.transpose([[float(i) for i in x.strip().split()] for x in file(dir+"output-times-slice").readlines()])
    timesslice = dict([(timesslice[1][i],timesslice[0][i]) for i in range(len(timesslice[0]))])
  except Exception, e:
    print "No slice timing information"
  try:
    timesanal = numpy.transpose([[float(i) for i in x.strip().split()] for x in file(dir+"output-times-kinetic").readlines()])
    timesanal = dict([(timesanal[1][i],timesanal[0][i]) for i in range(len(timesanal[0]))])
  except Exception, e:
      print "No Analysis timing information"
  try:
    timescube = numpy.transpose([[float(i) for i in x.strip().split()] for x in file(dir+"output-times-cube").readlines()])
    timescube = dict([(timescube[1][i],timescube[0][i]) for i in range(len(timescube[0]))])
  except Exception, e:
    print "No cube timing information"
  

  if methods % 2 == 0:
    print "Rendering %d Slices" % len(timesslice)
    files = timesslice.keys()
    files.sort()
    count = 1
    mkdir(outdir+"rho_slicez0")
    mkdir(outdir+"rho_slicez1")
    mkdir(outdir+"rhovx_slice0")
    mkdir(outdir+"rhovx_slice1")
    for t in files:
      print "t = %d [%f%%]" % (t,count*100.0/len(timesslice)) + " "*10
      rho,rhovx,rhovy,rhovz = load_slice(t,dir,ovars)
      rho = numpy.log10(rho)
      print "\rRendering and Saving" + " "*5 + "\r",
      sys.stdout.flush()
    
      tnorm = timesslice[t]/ovars['tmerge']
      
      py.clf()
      py.pcolor(rho[:,:,0])
      py.axhline(ovars['size'])
      py.axvline(ovars['size'])
      py.title(r'$ log_{10}(\rho)\ @\  %0.4f \cdot t_{merge} $' % tnorm)
      py.colorbar()
      py.savefig(outdir+"rho_slicez0/analysis-slice-rho0-%.8d.png"%t)
    
      py.clf()
      py.pcolor(rho[:,:,1])
      py.axhline(ovars['size'])
      py.axvline(ovars['size'])
      py.title(r'$ \log_{10}(\rho)\  @\  %0.4f \cdot t_{merge}$' % tnorm)
      py.colorbar()
      py.savefig(outdir+"rho_slicez1/analysis-slice-rho1-%.8d.png"%t)
    
      py.clf()
      py.pcolor(rhovx[:,:,0])
      py.axhline(ovars['size'])
      py.axvline(ovars['size'])
      py.title(r'$ \rho_x\ @\ %0.4f \cdot t_{merge}$' % tnorm)
      py.colorbar()
      py.savefig(outdir+"rhovx_slice0/analysis-slice-rhovx0-%.8d.png"%t)
    
      py.clf()
      py.pcolor(rhovx[:,:,1])
      py.axhline(ovars['size'])
      py.axvline(ovars['size'])
      py.title(r'$ \rho_x\ @\ %0.4f \cdot t_{merge}$' % tnorm)
      py.colorbar()
      py.savefig(outdir+"rhovx_slice1/analysis-slice-rhovx1-%.8d.png"%t)
    
      count += 1
      
  if False and methods % 3 == 0:
    print "Rendering %d Cubes" % len(timescube)
    files = timescube.keys()
    files.sort()
    count = 1
    mkdir(outdir+"rho_cubez")
    mkdir(outdir+"rho_cubey")
    mkdir(outdir+"rhovx_cube")
    for t in files:
      print "t = %d [%f%%]" % (t,count*100.0/len(timescube)) + " "*25
      rho,rhovx,rhovy,rhovz = load_cube(t,dir,ovars)
      tnorm = timescube[t]/ovars['tmerge']
      
      for i in xrange(0,int(ovars["gsize"]),5):
        print " "*25 + "\rRendering and Saving\t" + "z=%d [%f%%]"%(i,i*100.0/ovars["gsize"]) + " "*5 + "\r",
        sys.stdout.flush()
        
        py.clf()
        py.pcolor(rho[:,:,i])
        py.axhline(y=ovars['size'])
        py.axvline(x=ovars['size'])
        py.title(r'$ \rho_z @ %f \cdot t_{merge} $ & z=%d & min=%f' % (tnorm,i,rho[:,:,i].min()))
        py.colorbar()
        py.savefig(outdir+"rho_cubez/analysis-cube-rhoz-%.8d-%.3d.png"%(t,i))
        
        py.clf()
        py.pcolor(rho[:,i,:])
        py.axhline(y=ovars['size'])
        py.axvline(x=ovars['size'])
        py.axvline(x=ovars['ghost'])
        py.title(r'$ \rho_y @ %f \cdot t_{merge} $' % tnorm)
        py.colorbar()
        py.savefig(outdir+"rho_cubey/analysis-cube-rhoy-%.8d-%.3d.png"%(t,i))
    
        py.clf()
        py.pcolor(rhovx[:,:,i])
        py.axhline(y=ovars['size'])
        py.axvline(x=ovars['size'])
        py.axvline(x=ovars['ghost'])
        py.title(r'$ \rho_x @ %f \cdot t_{merge}$' % tnorm)
        py.colorbar()
        py.savefig(outdir+"rhovx_cube/analysis-cube-rhovx0-%.8d-%.3d.png"%(t,i))
    
      count += 1

  if methods % 5 == 0:  
    print "Rendering %d Sigma" % len(timesanal) + " "*10
    files = timesanal.keys()
    files.sort()
    count = 1
    mkdir(outdir+"sigma")
    for t in files:
      print "t = %d [%f%%]" % (t,count*100.0/len(timesanal)) + " "*10
    
      tnorm = timesanal[int(t)] / ovars['tmerge']
      sigma = load_sigma(t,dir,ovars)
      py.clf()
      py.bar(range(100),numpy.nan_to_num(sigma)[:100])
      py.axhline(1)
      py.title(r'$ \rho(\|v\|) @ t=%f \cdot t_{merge} $' % tnorm)
      py.xlabel(r'$ \frac{<\|v\|>}{v_{char}} $')
      py.ylabel(r'$ <\rho>  $')
      py.savefig(outdir+"sigma/analysis-sigma-%.8d.png" % t)
    
      count += 1
  
  if methods % 7 == 0:
    print "Rendering Kinetic"
    atimes = timesanal.values()
    atimes.sort()
    atimes = numpy.array(atimes)
    kinetic = load_kinetic(dir,ovars)
    mass = load_mass(dir,ovars)

    try:
      mergestep = cubetimes.values()
      mergestep.sort()
      mergestep = mergestep[0]
    except Exception:
      mergestep = 0

    print "\rSaving\r",
    sys.stdout.flush()
    mkdir(outdir+"kinetic")
    py.clf()
    py.scatter(atimes/ovars['tmerge'],kinetic)
    py.axvline(1)
    py.axhline(kinetic[mergestep:].mean())
    py.xlim(xmin=0)
    py.title("Mean kinetic energy ($<E_{kin}(t>t_{merge}> = %f$)"%kinetic[mergestep:].mean())
    py.xlabel("time ($t_{merge}$)")
    py.ylabel("$<E_{kin}>$")
    py.savefig(outdir+"kinetic/analysis-kinetic.png")

    massratio = ovars['op'] * ovars['gsize']**3 / mass
    corrkin = kinetic * massratio
    py.clf()
    py.scatter(atimes/ovars['tmerge'],corrkin)
    py.axvline(1)
    py.axhline(corrkin[mergestep:].mean())
    py.xlim(xmin=0)
    py.title("Mean kinetic energy ($<E_{kin}(t>t_{merge})> = %f$)"%corrkin[mergestep:].mean())
    py.xlabel("time ($t_{merge}$)")
    py.ylabel("$<E_{kin}>$")
    py.savefig(outdir+"kinetic/analysis-kinetic-masscorrected.png")


  if methods % 11 == 0:
    print "Rendering Mass"
    atimes = timesanal.values()
    atimes.sort()
    atimes = numpy.array(atimes)
    mass = load_mass(dir,ovars)

    mkdir(outdir+"mass")
    py.clf()
    py.scatter(atimes/ovars['tmerge'],mass)
    py.axvline(1)
    py.xlim(xmin=0)
    py.xlabel("time ($t_{merge}$)")
    py.ylabel("Mass")
    py.title("Total mass")
    py.savefig(outdir+"mass/analysis-mass.png")


  if methods % 13 == 0:
    print "Rendering Compressional"
    atimes = timesanal.values()
    atimes.sort()
    atimes = numpy.array(atimes)
    compress = load_compressional(dir,ovars)

    try:
      mergestep = cubetimes.values()
      mergestep.sort()
      mergestep = mergestep[0]
    except Exception:
      mergestep = 0

    print "\rSaving\r",
    sys.stdout.flush()
    mkdir(outdir+"compress")
    py.clf()
    py.scatter(atimes/ovars['tmerge'],compress)
    py.axvline(1)
    py.axhline(compress[mergestep:].mean())
    py.xlim(xmin=0)
    py.title("Mean Compressional energy ($<E_{compress}(t>t_{merge})> = %f$)"%compress[mergestep:].mean())
    py.xlabel("time ($t_{merge}$)")
    py.ylabel("$<E_{compress}>$")
    py.savefig(outdir+"compress/analysis-compress.png")

  def meanrhovelocity(i, slice, dir, ovars):
    if slice:
      slice = "slice"
      size = ovars['size']**2
    else:
      slice = "cube"
      size = ovars['size']**3
    print ovars['size']
    sum = numpy.zeros(4)
    for n in range(ovars['procs']):
      fd = open(dir+"output-%s-%.8d-%.3d"%(slice,i,n))
      data = numpy.fromfile(file=fd,dtype=numpy.float32).reshape((4,size))
      sum += data.sum(axis=1)
    sum[1:] /= sum[0]
    return sum

  def calculatesigma(i, slice, means, dir, ovars):
    if slice:
      slice = "slice"
      size = ovars['size']**2
    else:
      slice = "cube"
      size = ovars['size']**3
    sigma = numpy.zeros(3)
    for n in range(ovars['procs']):
      fd = open(dir+"output-%s-%.8d-%.3d"%(slice,i,n))
      data = numpy.fromfile(file=fd,dtype=numpy.float32).reshape((4,size))
      for i in range(size):
        line = data[:,i]
        sigma += (line[1:]/line[0] - means[1:])**2 * line[0]
    return sigma / means[0]
        
  if methods % 17 == 0:
    print "Calculating velocity dispersion"
    steps = timesslice.keys() #+ timescube.keys()
    steps.sort()
    steps = numpy.array(steps)
    times = timesslice.values() #+ timescube.values()
    times.sort()
    times = numpy.array(times)
   
    #steps = steps[::5]
    #times = times[::5]
   
    vrms2 = []
    means = []
    totsteps = len(steps)
    for progress,i in enumerate(steps):
      isslice = True #not timescube.has_key(i)
      print "\r%0.2f%% (%s)"%(progress/totsteps*100, "slice" if isslice else "cube") + " "*10,
      sys.stdout.flush()
      means.append(meanrhovelocity(i, isslice, dir, ovars))
      vrms2.append(numpy.sum(calculatesigma(i, isslice, means[-1],dir,ovars)**2))
    means = numpy.array(means)
    vrms2 = numpy.array(vrms2)
    
    mergestep = 0
    while times[mergestep] < ovars['tmerge'] and mergestep < len(times):
      mergestep += 1
    if mergestep == len(times):
      mergestep = 0

    mkdir(outdir+"meanvel")
    py.clf()
    py.scatter(times/ovars['tmerge'],means[:,1],label='$<v_x>$',c='r')
    py.scatter(times/ovars['tmerge'],means[:,2],label='$<v_y>$',c='b')
    py.scatter(times/ovars['tmerge'],means[:,3],label='$<v_z>$',c='g')
    py.legend()
    py.axvline(1)
    py.axhline(0)
    py.xlim(xmin=0)
    py.title("Mean Velocity ($<v_i(t)$)")
    py.xlabel("time ($t_{merge}$)")
    py.ylabel("$<v_i(t)>$")
    py.savefig(outdir+"meanvel/analysis-meanvel.png")
   
    mkdir(outdir+"veldist")
    py.clf()
    py.scatter(times/ovars['tmerge'],vrms2)
    py.axvline(1)
    py.axhline(y=vrms2[mergestep:].mean())
    py.xlim(xmin=0)
    py.title("3D Velocity Dispersion ($v_{rms}(t>t_{merge})^2 = %0.4f$)"%vrms2[mergestep:].mean())
    py.xlabel("time ($t_{merge}$)")
    py.ylabel("$v_{rms}(t)^2$")
    py.savefig(outdir+"veldist/analysis-3Dvelocitydispersion.png")

    disp3dmean = (numpy.sqrt(vrms2/3))[mergestep:].mean()
    py.clf()
    py.scatter(times/ovars['tmerge'],numpy.sqrt(vrms2/3))
    py.axvline(1)
    py.axhline(y=disp3dmean)
    py.xlim(xmin=0)
    py.title("1D Velocity Dispersion ($ < \sigma (t>t_{merge}) > = %0.4f $)"%disp3dmean)
    py.xlabel("time ($t_{merge}$)")
    py.ylabel("$ \sigma(t) $")
    py.savefig(outdir+"veldist/analysis-1Dvelocitydispersion.png")

def getegyspec(ux, uy, uz):
  assert ux.shape == uy.shape == uz.shape
  size = len(ux)
  energy = {}
  V = reduce(lambda x,y : x*y, ux.shape)
  for x in range(size):
    for y in range(size):
        absk = numpy.sqrt(x*x + y*y)
        egy = abs(ux[x,y]**2 + uy[x,y]**2 + uz[x,y]**2) / V
        if energy.has_key(absk):
          energy[absk] += egy
        else:
          energy.update({absk:egy})
  return numpy.array(energy.keys()), numpy.array(energy.values())

if methods % 19 == 0:
    print "Calculating Energy Spectrums"
    steps = timesslice.keys() #+ timescube.keys()
    steps.sort()
    steps = numpy.array(steps)
    times = timesslice.values() #+ timescube.values()
    times.sort()

    mkdir(outdir+"energyspec")
    mkdir(outdir+"energyspec-data")
    totanal = len(steps)
    resolve = 5
    maxegy = 0
    for i,step in reverse_enumerate(steps[::resolve]):
      rho, rhovx, rhovy, rhovz = numpy.array(load_slice(step,dir,ovars))[:,:,:,0]
      print "Calc+Rndr: %0.2f%%"%((totanal-resolve*i)*100/totanal) + " "*10 + "\r",
      sys.stdout.flush()
      ukx = numpy.fft.fft2(rhovx/rho) 
      uky = numpy.fft.fft2(rhovy/rho)
      ukz = numpy.fft.fft2(rhovz/rho)
      k, energyspectrum = getegyspec(ukx,uky,ukz)
      
      nk = numpy.arange(1,int(ovars['gsize']/2),1)
      ncount = numpy.zeros(int(ovars['gsize']/2))
      negy = numpy.zeros(int(ovars['gsize']/2))
      for j in range(len(k)):
        ni = int(k[j])
        if ni >= 1 and ni <= int(ovars['gsize']/2):
          negy[ni-1] += energyspectrum[j]
          ncount[ni-1] += 1
      negy = negy / ncount
      fit = numpy.polyfit(numpy.log10(nk[int(ovars['gsize']/4):]),numpy.log10(negy[int(ovars['gsize']/4):-1]),2)

      if maxegy == 0:
        maxegy = energyspectrum.max()  
        
      py.clf()
      py.gca().set_xscale('log')
      py.gca().set_yscale('log');
      py.scatter(k,energyspectrum)
      x = numpy.linspace(ovars['lmerge']/2,ovars['gsize']/2.0)
      line = py.plot(x,x**(-fit[0])*negy[int(ovars['lmerge']/2)]/nk[int(ovars['lmerge']/2)]**(-fit[0]),c='r')
      py.legend( (line,), (r'$ \beta = %0.2f$'%fit[0],))
      py.axis([0,ovars['gsize']/2,0,maxegy])
      py.xlim(xmax=ovars['gsize']/2)
      py.axvline(x=ovars['lmerge']/2,c='k')
      py.title("Energy Spectrum at $t = %0.4f t_{merge}$"%(times[i*resolve]/ovars['tmerge']))
      py.xlabel("$log_{10}[k]$")
      py.ylabel("$log_{10}[E(k)]$")
      py.savefig(outdir+"energyspec/analysis-energyspectrum-%08d.png"%step)

      numpy.savez(outdir+"energyspec-data/energyspectrum-k-%08d.npz"%step,k=k,energyspectrum=energyspectrum)
