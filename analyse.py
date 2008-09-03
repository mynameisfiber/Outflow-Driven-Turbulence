#! /usr/bin/env python

import sys
import pylab as py
import numpy

# Get system arguments
try:
	size = int(sys.argv[1])
	start = int(sys.argv[2])
	stop = int(sys.argv[3])
	if(len(sys.argv) == 5):
		step = int(sys.argv[4])
	else:
		step = 1
except Exception, e:
	print "Could not read in parameters"
	exit()

# Now we store the output times
fd = file("outputtimes")
timesvalues, timesname = numpy.transpose([[float(num) for num in line.split()] for line in fd.readlines()])

times = []
kin = []
e = numpy.zeros((size,size,size))

print "start: %d stop: %d step: %d"	% (start,stop,step)
for frame in range(start,stop+1,step):
	# Read file
	filename = "output-%.8d" % frame
	print "Analyzing %s\t\t" % filename,
	sys.stdout.flush()
	try:
		fd = file(filename)
	except Exception, e:
		print "[FAIL]"
		print "Could not open file: %s" % filename
		exit()
	
	# Get the data from the file and shape it into proper arrays
	rho, rhovx, rhovy, rhovz, = numpy.transpose([[float(num) for num in line.split()] for line in fd.readlines()])

	rho = rho.reshape(size,size,size)
	rhovx = rhovx.reshape(size,size,size)
	rhovy = rhovy.reshape(size,size,size)
	rhovz = rhovz.reshape(size,size,size)
	fd.close
	
	# Extract the simulation time of the frame
	time = timesvalues[timesname.searchsorted(frame)]
	times.append(time)

	# And now for the pcolors!
	py.clf()
	py.pcolor(rho[:][:][size/2])
	py.colorbar()
	py.title("Rho @ t = %f" % time)
	py.savefig("images/rho-slice-%.8d.png" % frame)
	
	py.clf()
	py.pcolor(rho.mean(axis=2))
	py.colorbar()
	py.title("Rho @ t = %f" % time)
	py.savefig("images/rho-meanz-%.8d.png" % frame)
	
	# Now we calculate the Kin energy
	kin.append(numpy.mean(1.0/2.0 * (rhovx**2 + rhovy**2 + rhovz**2) / rho))
	
	# And now to update E(k)
	e += numpy.fft.fftn((rhovx**2 + rhovy**2 + rhovz**2) / rho)
	
	# Sigma:
	
	################
	# Old sigma(rho)=v
	################
	# bins = 100
	# rhomin = rho.min()
	# rhomax = rho.max()
	# width = (rhomax - rhomin)/bins
	# if width:
	# 	X = numpy.arange(rhomin,rhomax+width/2,width)
	# 	Y = numpy.zeros(bins+1)
	# 	count = numpy.zeros(bins+1)
	# 	for x in range(len(rho)):
	# 		for y in range(len(rho[0])):
	# 			for z in range(len(rho[0][0])):
	# 				i = int((rho[x,y,z] - rhomin)/width)
	# 				Y[i] += (rhovx[x,y,z]**2 + rhovy[x,y,z]**2 + rhovz[x,y,z]**2)**.5 / rho[x,y,z]
	# 				count[i] += 1
	# 	Y = Y / count
	# 	py.clf()
	# 	a = py.gca()
	# 	a.set_xlim([rhomin,rhomax])
	# 	py.bar(X,numpy.nan_to_num(Y),width=width)
	# 	py.axhline(1)
	# 	py.title("Average velocity of density with %d bins" % bins)
	# 	py.xlabel("rho")
	# 	py.ylabel("Average velocity of cells with rho")
	# 	py.savefig("images/sigma-%.8d.png" % frame)
	
	################
	# New sigma(r)=v
	################
	step = 1
	sigma = numpy.zeros(int(size/step))
	for i in range(1,size+1,step):
		for offset in range(0,i):
			tmp = numpy.zeros(3)
			for x in range(offset,offset+size/i):
				for y in range(offset,offset+size/i):
					for z  in range(offset,offset+size/i):
						tmp[0] += abs(rhovx[x,y,z] / rho[x,y,z])
						tmp[1] += abs(rhovy[x,y,z] / rho[x,y,z])
						tmp[2] += abs(rhovz[x,y,z] / rho[x,y,z])
			tmp *= i*1.0/size
			sigma[i/step-1] += (tmp[0]**2 + tmp[1]**2 + tmp[2]**2)**.5
		if sigma[i/step-1] != 0: sigma[i/step-1] /= i**3
		
	py.clf()
	py.bar(range(1,size+1,step),sigma,width=step)
	py.axhline(1)
	py.title("Average velocity of spheres")
	py.xlabel("r")
	py.ylabel("sigma(r)")
	py.savefig("images/sigma-%.8d.png" % frame)
	print "[DONE]"
	
py.clf()
py.plot(times,kin)
py.title("Mean kinetic energy per timestep")
py.savefig("images/kin.png")
