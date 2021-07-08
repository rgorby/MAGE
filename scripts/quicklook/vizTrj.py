#!/usr/bin/env python
#Make a quick figure of a CHIMP particle trajectory

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
import kaipy.kaiH5 as kh5
import kaipy.chimp.chimph5p as ch5p
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D

#Adds a simple Earth to 3d pics
def addEarth(ax,Re=1):
	N = 100
	T = 0.25 #Transparency
	S = 4 #Stride
	phi = np.linspace(0, 2 * np.pi, N)
	theta = np.linspace(0, np.pi, N)
	xm = Re * np.outer(np.cos(phi), np.sin(theta))
	ym = Re * np.outer(np.sin(phi), np.sin(theta))
	zm = Re * np.outer(np.ones(np.size(phi)), np.cos(theta))

	ax.plot_surface(xm, ym, zm,rstride=S, cstride=S, color='b',alpha=T)

#Forces equal spacing on a 3d plot (surprisingly annoying to do)
def axEqual3d(ax):
	x_limits = ax.get_xlim3d()
	y_limits = ax.get_ylim3d()
	z_limits = ax.get_zlim3d()

	x_range = x_limits[1] - x_limits[0]; x_mean = np.mean(x_limits)
	y_range = y_limits[1] - y_limits[0]; y_mean = np.mean(y_limits)
	z_range = z_limits[1] - z_limits[0]; z_mean = np.mean(z_limits)

	# The plot bounding box is a sphere in the sense of the infinity
	# norm, hence I call half the max range the plot radius.
	plot_radius = 0.5*max([x_range, y_range, z_range])

	ax.set_xlim3d([x_mean - plot_radius, x_mean + plot_radius])
	ax.set_ylim3d([y_mean - plot_radius, y_mean + plot_radius])
	ax.set_zlim3d([z_mean - plot_radius, z_mean + plot_radius])


if __name__ == "__main__":
	#Arg parsing
	id0 = 0
	LW = 0.5
	fs = 16

	MainS = """Plots trajectory of specified CHIMP particle in H5p file
	ID: Particle ID
	Color specified by energy [keV]
	"""
	
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('h5p',metavar='input.h5part',help="Input H5Part file")
	parser.add_argument('-id',type=int,metavar="id",default=id0,help="Display particle trajectory w/ given ID (default: %(default)s)")

	args = parser.parse_args()
	fIn = args.h5p
	id0 = args.id

	print("Reading from %s and looking for particle %d"%(fIn,id0))
	kh5.CheckOrDie(fIn)

	Ntp,nS,nE = ch5p.bndTPs(fIn)
	print("\tFound %d TPs, IDs %d to %d"%(Ntp,nS,nE))

	t,x = ch5p.getH5pid(fIn,"x",id0)
	t,y = ch5p.getH5pid(fIn,"y",id0)
	t,z = ch5p.getH5pid(fIn,"z",id0)
	t,K = ch5p.getH5pid(fIn,"K",id0)

	vMin = 0.0
	vMax = K.max()

	cbLab = "Energy [keV]"
	fig = plt.figure()
	#ax = fig.gca(projection='3d')
	ax = Axes3D(fig)

	addEarth(ax)
	sct3d = ax.scatter(x,y,z,c=K,cmap=plt.get_cmap("cool"),s=20,vmin=vMin,vmax=vMax)
	ax.plot(x,y,z,'k',linewidth=LW)
	plt.colorbar(sct3d,label=cbLab,shrink=0.8)

	ax.set_xlabel("X [Re]")
	ax.set_ylabel("Y [Re]")
	ax.set_zlabel("Z [Re]")
	titS = "Particle trajectory, pID = %d"%(id0)
	ax.set_title(titS,fontsize=fs)
	axEqual3d(ax)
	plt.show()