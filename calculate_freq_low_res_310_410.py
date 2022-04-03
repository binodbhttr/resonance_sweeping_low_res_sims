import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import struct
import ctypes
import os
import sys
import pickle
from galpy.util import coords
import agama
from scipy.stats import binned_statistic
from scipy import stats as stats
import scipy 

datapath="/work2/07428/binod/stampede2/LRBB-IoMW/"
savepath_actions="/work2/07428/binod/stampede2/resonance_sweeping_low_res_sims/actions_angles_freqs/actions/"
savepath_angles="/work2/07428/binod/stampede2/resonance_sweeping_low_res_sims/actions_angles_freqs/angles/"
savepath_freqs="/work2/07428/binod/stampede2/resonance_sweeping_low_res_sims/actions_angles_freqs/freqs/"

times=np.genfromtxt(datapath+'times.txt',dtype='str')
ncores=8
print(len(times),'snapshots on',ncores,'cores')

start=310
end=410


def loader(filename):
    ffile = open(datapath+filename, 'rb')
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


for i in range(start,end,10):
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
    xd=catd['x']
    yd=catd['y']
    zd=catd['z']
    vxd=catd['vx']
    vyd=catd['vy']
    vzd=catd['vz']
    massd=catd['mass']
    iddd=catd['ID']
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
            print('snapshot__'+time[i]+'-'+str(j),'has no star particles')
        try:
            xd=np.hstack((xd,catd['x']))
            yd=np.hstack((yd,catd['y']))
            zd=np.hstack((zd,catd['z']))
            vxd=np.hstack((vxd,catd['vx']))
            vyd=np.hstack((vyd,catd['vy']))
            vzd=np.hstack((vzd,catd['vz']))
            iddd=np.hstack((iddd,catd['ID']))
            massd=np.hstack((massd,catd['mass']))
        except:
            print('snapshot__'+time[i]+'-'+str(j),'has no dark particles')

    massflaglo=5000/2.324876e9 #slightly higher than the min mass of stars and lower than the max mass of stars 
    #min mass of stars for this sim: 4496.5737 max mass of stars for this sim: 13416.486
    massflaghigh= 400000/2.324876e9 # slightly higher than the mass of DM particles 
    # mass of each dark matter particles 348275.70762648

    discindx=(mass<massflaglo)
    bulgeindx=(mass>massflaglo)*(mass<massflaghigh)
    #gets the actions for previous function and saves it
    agama.setUnits(mass=1, length=1, velocity=1)
    load = np.load
    potentialf = agama.Potential
    finder = agama.ActionFinder 

    dm_mass   = massd*2.324876e9  # mass needs to be in solar masses
    disk_mass = mass[discindx]*2.324876e9
    bulg_mass = mass[bulgeindx]*2.324876e9

    bulg_phasespace = np.vstack((x[bulgeindx],y[bulgeindx],z[bulgeindx],vx[bulgeindx],vy[bulgeindx],vz[bulgeindx])).T
    disk_phasespace = np.vstack((x[discindx],y[discindx],z[discindx],vx[discindx],vy[discindx],vz[discindx])).T
    dm_phasespace   = np.vstack((xd,yd,zd,vxd,vyd,vzd)).T

    print(np.shape(bulg_phasespace))

    dark = potentialf(type = "Multipole",particles=(dm_phasespace, dm_mass), symmetry='a', gridsizeR=20, lmax=2)
    disk = potentialf(type = "CylSpline",particles=(disk_phasespace, disk_mass),gridsizeR=20, gridsizeZ=20,
                      mmax=0, Rmin=0.1, symmetry='a',Rmax=70, Zmin=0.02, Zmax=30)            
    bulge = potentialf(type = "Multipole",particles=(bulg_phasespace, bulg_mass),symmetry='a',gridsizeR=20,lmax=2)

    potential = potentialf(dark, disk, bulge) #combining the potentials to get full potential

    w0 = np.zeros((len(disk_phasespace[:,0]), 3)) #creating an array with just the X,Y,Z positions of the Disk

    w0[:,0] = disk_phasespace[:,0]
    w0[:,1] = disk_phasespace[:,1]
    w0[:,2] = disk_phasespace[:,2]

    pot = np.asarray(potential.potential(w0)) #getting the values of the full potential at each point
    force = potential.force(w0) #getting the forces due to the full potential at each point
    #np.save('Pot_initial.npy', pot) #saving for leapfrog code
    #np.save('Force_initial.npy', force) #saving for leapfrog code


    af = finder(potential, interp=False)

    Disk_actions, Disk_angles, Disk_freq = af(disk_phasespace, angles=True) #getting the actions, angles, and freqs

    Jrdisk = Disk_actions[:, 0] #unpacking them all
    Jzdisk = Disk_actions[:, 1]
    Jphidisk = Disk_actions[:, 2]

    Trdisk = Disk_angles[:, 0]
    Tzdisk = Disk_angles[:, 1]
    Tphidisk = Disk_angles[:, 2]

    Ordisk = Disk_freq[:, 0]
    Ozdisk = Disk_freq[:, 1]
    Ophidisk = Disk_freq[:, 2]

    n_counts=len(Jrdisk)
    saved_actions = np.zeros(shape = (4,n_counts))
    saved_angles = np.zeros(shape = (4,n_counts)) #saving them as a np array for ease of storage
    saved_freqs = np.zeros(shape = (4,n_counts)) #saving them as a np array for ease of storage

    saved_actions[0,:] = Jrdisk
    saved_actions[1,:] = Jphidisk
    saved_actions[2,:] = Jzdisk
    saved_actions[3,:] = idd[discindx]

    saved_angles[0,:] = Trdisk
    saved_angles[1,:] = Tphidisk
    saved_angles[2,:] = Tzdisk
    saved_angles[3,:] = idd[discindx]

    saved_freqs[0,:] = Ordisk
    saved_freqs[1,:] = Ophidisk
    saved_freqs[2,:] = Ozdisk
    saved_freqs[3,:] = idd[discindx]

    np.save(savepath_actions+'DiskActions'+str(i)+'.npy', saved_actions)
    np.save(savepath_angles+'DiskAngles'+str(i)+'.npy', saved_angles)
    np.save(savepath_freqs+'DiskFreqs'+str(i)+'.npy', saved_freqs)