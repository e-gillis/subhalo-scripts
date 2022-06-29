#!/usr/local/bin/python3
# coding: utf-8

import numpy as np
from galpy import potential
from galpy.potential import NFWPotential,rtide
from galpy.util import bovy_conversion,bovy_coords
import clustertools as ctools
from galpy.orbit import Orbit,Orbits

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


ro=8.
vo=220.
mo=bovy_conversion.mass_in_msol(ro=ro,vo=vo)

#Parameters from Galacticus

#(1) Cosmological parameters (Planck 2018)
OmegaM = 0.3153
OmegaL = 0.6847
sigma_8 = 0.8111
h0 = 0.6736

#(2) Halo mass defintion

delta_c = 200 #with respect to the critical density at infall redshift.


cosmo = FlatLambdaCDM(H0=h0*100.0,Om0=OmegaM)
t07=cosmo.age(0.7)/u.Gyr
t05=cosmo.age(0.5)/u.Gyr
t0=cosmo.age(0.)/u.Gyr

data_z_0_5 = np.load('M_1.0e13_z_0.5_Mres_1.0e6.npz')

# Hostpotential
hhost=cosmo.H(0.5).value/100
hostpot=NFWPotential(mvir=data_z_0_5['MHost'][0]/1e12,
                     conc=data_z_0_5['cHost'][0],
                     H=hhost*100.0,Om=OmegaM,wrtcrit=True,
                     overdens=delta_c,ro=ro,vo=vo)


# This is where we get our halos
m200,r200,c200=data_z_0_5['M200c'],data_z_0_5['R200c']*1000.0,data_z_0_5['c200c']
zinfall=data_z_0_5['zInfall']
tinfall=cosmo.age(zinfall)/u.Gyr
x,y,z=data_z_0_5['posX']*1000.0,data_z_0_5['posY']*1000.0,data_z_0_5['posZ']*1000.0
vx,vy,vz=data_z_0_5['velX'],data_z_0_5['velY'],data_z_0_5['velZ']

# Conversion to cylindrical coordinates
print("Converting Coordinates")
rad,phi,zed=bovy_coords.rect_to_cyl(x,y,z)
print("Converting Velocities")
vR,vT,vzed=bovy_coords.rect_to_cyl_vec(vx,vy,vz,x,y,z)

np.random.seed(42)
idsubs=(np.arange(len(m200)))
np.random.shuffle(idsubs)
sim_nums=np.arange(1, len(idsubs)+1)

vxvv = np.column_stack([rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi])
orbits = Orbit(vxvv,ro=ro,vo=vo)
f = open("orbital_params.dat", "w")
f.write("# No.  ID       Pericenter Apocenter     Eccen.    Energy\n")


for i in range(len(idsubs)):
    if (i + 1) % 10000 == 0: print(f"Integrating orbit {i+1}")
    indx = idsubs[i]
    os = orbits[indx]
    sim_num = sim_nums[i]
    
    ts=np.linspace(0,10/bovy_conversion.time_in_Gyr(ro=ro,vo=vo),10000)
    os.integrate(ts,hostpot)
    semi=(os.rap()-os.rperi())/2
    f.write((str(sim_num).zfill(6)+' '+str(indx).zfill(6)+"{:12f}"*4+"\n")\
           .format(os.rperi(),os.rap(),os.e(),semi,os.E()))
f.close()