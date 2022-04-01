import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import struct
import ctypes
import os
import pickle
from numpy import fft as fftnp
from galpy.potential import MWPotential2014
from galpy.actionAngle import actionAngleStaeckel, estimateDeltaStaeckel, actionAngleIsochroneApprox
from galpy.util import bovy_coords as coords
from scipy import stats as stats
import scipy

def loader(filename):
    ffile = open(filename, 'rb')
    t,n,ndim,ng,nd,ns,on=struct.unpack("<diiiiii",ffile.read(32))
    print(t,n,ndim,ng,nd,ns,on)

    catd = {'mass':np.zeros(nd), 'x':np.zeros(nd), 'y':np.zeros(nd),'z':np.zeros(nd),'vx':np.zeros(nd),'vy':np.zeros(nd),'vz':np.zeros(nd),'ID':np.zeros(nd)}
    cats = {'mass':np.zeros(ns), 'x':np.zeros(ns), 'y':np.zeros(ns),'z':np.zeros(ns),'vx':np.zeros(ns),'vy':np.zeros(ns),'vz':np.zeros(ns),'metals':np.zeros(ns), 'tform':np.zeros(ns), 'ID':np.zeros(ns)}

    for i in range(nd):
        mass, x, y, z, vx, vy, vz, IDs = struct.unpack("<fffffffQ", ffile.read(36))
        catd['mass'][i] = mass
        catd['x'][i] = x
        catd['y'][i] = y
        catd['z'][i] = z
        catd['vx'][i] = vx*100.
        catd['vy'][i] = vy*100.
        catd['vz'][i] = vz*100.
        catd['ID'][i] = IDs

    for i in range(ns):
        mass, x, y, z, vx, vy, vz, metals, tform, IDs = struct.unpack("<fffffffffQ", ffile.read(44))
        cats['mass'][i] = mass
        cats['x'][i] = x
        cats['y'][i] = y
        cats['z'][i] = z
        cats['vx'][i] = vx*100.
        cats['vy'][i] = vy*100.
        cats['vz'][i] = vz*100.
        cats['metals'][i] = metals
        cats['tform'][i] = tform
        cats['ID'][i] = IDs
    return(catd,cats,t)

os.system('ls snapshot*-0 > times.txt')
times=np.genfromtxt('times.txt',dtype='str')
#os.system('rm times.txt')
print(times)

ncores=8

ro=8.
vo=220.

orbit=open('bar.txt','w')

barangle=0.
bmbarangle=0.

for i in range(0,510):
    print('Loading '+times[i])
    catd,cats,t=loader(times[i])
    x=cats['x']
    y=cats['y']
    z=cats['z']
    vx=cats['vx']
    vy=cats['vy']
    vz=cats['vz']
    mass=cats['mass']
    idd=cats['ID']

    timefw=t*9.778145/1000.

    for j in range(1,ncores):
        print('Loading '+times[i][:-1]+str(j))
        catd,cats,t=loader(times[i][:-1]+str(j))

        try:
            x=np.hstack((x,cats['x']))
            y=np.hstack((y,cats['y']))
            z=np.hstack((z,cats['z']))
            vx=np.hstack((vx,cats['vx']))
            vy=np.hstack((vy,cats['vy']))
            vz=np.hstack((vz,cats['vz']))
            idd=np.hstack((idd,cats['ID']))
            mass=np.hstack((mass,cats['mass']))
        except:
            print('snapshot__'+times[i]+'-'+str(j),'has no star particles')

    vr,vphi,vzz=coords.rect_to_cyl_vec(vx,vy,vz,x,y,z)
    r,phi,zz=coords.rect_to_cyl(x,y,z)

    barsample=(r<5)*(r>1)

    plt.figure()
    counts, bins, patches=plt.hist(phi[barsample],bins=360,range=[-np.pi,np.pi])
    ff=np.fft.fft(counts-np.mean(counts))
    barangle=-np.angle(ff[2])/2.
    plt.close()

    TU=9.778145/1000.

    if barangle>oldbarangle+np.pi:
        barangle=barangle-np.pi
    if barangle<oldbarangle:
        if (barangle+2*np.pi)>oldbarangle+np.pi:
            barangle=barangle-np.pi

    print("%E %E" %(timefw,barangle), file=orbit, flush=True)
