#!/usr/local/bin/python3
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from galpy import potential
from galpy.potential import NFWPotential,HernquistPotential,rtide
from galpy.util import conversion
import scipy.integrate as integrate
from galpy.df import isotropicNFWdf
import clustertools as ctools

# To make directories and scripts
import os
import stat

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

import sympy as sym

# Halo number to start from
nstart = 0
# Halos to extract
nmake = 10
# Pericenter Range
peri_range = [0, 100]

# Concentration Factor
cfactor = 1

ro=8.
vo=220.
mo=conversion.mass_in_msol(ro=ro,vo=vo)

#Check settings for subhalo
def f(x):
    return np.log(1.+x)-x/(1+x)


def get_mfactor(rvir, rt, c):
    # Set constants
    r = sym.Symbol('r')
    rs = rvir/c
    
    # Truncation function from mkhalo
    trunc = 2 / (sym.sech(r/rt) + sym.cosh(r/rt))
    subhalo_rho = trunc / (r/rs * (1 + r/rs)**2)
    # Can drop the 4 pi factor because we take a ratio
    mass_integration = subhalo_rho*r**2
    
    # Set up integral 
    M_factor = sym.Integral(mass_integration, (r, 0, rvir)) /\
               sym.Integral(mass_integration, (r, 0, 10*rt))
    m_eval = M_factor.evalf()
    
    # Normalize according to the truncation function
    m_correction = sym.integrate(trunc*4*sym.pi*r**2, (r, 0, rvir)) /\
                   (4/3 * sym.pi * rvir**3)
    
    return float(m_correction / m_eval)


def get_mkhalo_params(rvir, m200, c, Nacc, fname):
    rt = 1.5 * rvir
    mfactor = get_mfactor(rvir, rt, c)
    
    N_sample = int(1.01 * Nacc * mfactor)
    mass = m200 * mfactor
    rs = rvir / c

    return (fname, N_sample, rs, mass, rt)


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

cosmo.H(0.5)*u.Mpc*u.s/u.km

data_z_0_5 = np.load('M_1.0e13_z_0.5_Mres_1.0e6.npz')
    

hhost=cosmo.H(0.5).value/100


# This is where we get our halos
m200,r200,c200=data_z_0_5['M200c'],data_z_0_5['R200c']*1000.0,data_z_0_5['c200c']
zinfall=data_z_0_5['zInfall']
tinfall=cosmo.age(zinfall)/u.Gyr
x,y,z=data_z_0_5['posX']*1000.0,data_z_0_5['posY']*1000.0,data_z_0_5['posZ']*1000.0
vx,vy,vz=data_z_0_5['velX'],data_z_0_5['velY'],data_z_0_5['velZ']
r=np.sqrt(x*x+y*y+z*z)

c200 = c200 * cfactor

np.random.seed(42)
idsubs=(np.arange(len(m200)))
np.random.shuffle(idsubs)
sim_nums=np.arange(1, len(idsubs)+1)

# Make sure this is the same, but it should be
peri, apo, eccen, energy = np.loadtxt("orbital_params.dat").T[2:]

np.savetxt('idsubs.dat', np.vstack((sim_nums, idsubs)).T, fmt='%06i')
gyr_command_base = 'gyrfalcON in=%s.nemo out=%s.nemo tstop=%f eps=%f kernel=2 '\
                + f'step=%f kmax=%i Nlev=5 fac=0.01 > %s.out_z0%s 2>&1'
mkhalo_base = "mkhalo out=%s.nemo nbody=%i model=NFW r_s=%f M=%f r_t=-%f WD_units=t"

tstep_restart = (t05-t07)/10.  
        
# Add master list of Jobs
joblist = open("jobs_to_run.dat", "a")
finlist = open("finished.dat", "w")
finlist.close()
mass_shuffled = m200[idsubs]

for mass_range in [[1e6, 1e13]]:
    # Cut all our data arrays
    # Cut for m200
    mass_cut = np.logical_and(mass_shuffled >= mass_range[0], mass_shuffled < mass_range[1])
    # Cut for pericenters
    peri_cut = np.logical_and(peri >= peri_range[0], peri < peri_range[1])
    # Now make sure we get the start right 
    index_cut = sim_nums > nstart
    # Combine Everything
    valid_indeces = np.argwhere(np.logical_and(mass_cut, index_cut)).flatten()[:nmake]
    print(valid_indeces)

    for i in valid_indeces:
        indx = idsubs[i]
        sim_num = sim_nums[i]

        # Set Up Subhalo Parameters
        hz = cosmo.H(zinfall[indx]).value/100

        # Make Directory
        dirname = "{:06}_{:06}".format(sim_num, indx)
        os.mkdir(dirname)
        fname = f"{dirname}/{dirname}"
        # Open Job file
        jobfile = open(f"{dirname}/job_commands.sh", "w")
        # Open a Parameters File
        params = open(f"{dirname}/{str(sim_num).zfill(6)}_sim_params.dat", "w")
        param_tup = ("job num", "job id", "start z", "sim start", 
                     "evo time", "tstep", "radius", "x", "y", "z", 
                     "vx", "vy", "vz","mass")
        # Some Fancy Formating to save the information
        params.write("#"+(("%11s "*len(param_tup)) % param_tup) + "\n")
        jobfile.write(r"#!/usr/bin/bash" + "\n"*2)
        jobfile.write("touch running\n")

        # Add job to list
        joblist.write(dirname + "\n")
        
        # Will satisfy Van den Bosch criteria
        epsilon = f(c200[indx])*r200[indx]/(600 * c200[indx]) *\
                  (10/c200[indx])**0.63
        Nacc = int(1000.0*(0.01**(-1./0.8))/0.32)
        
        tstep = 0.05
        tfinal = 2
        zfin = 5
        gyr_outfile = f'{str(indx).zfill(6)}_z00'
        
        # Find optimal timestep
        zmbar = m200[indx]               # In MSun
        rbar  = r200[indx] * 1000        # In pc
        eps_nbody = epsilon / r200[indx] # Unitless

        vbar = 0.06557 * np.sqrt(zmbar / rbar)
        tbar = rbar / (1.023*vbar)

        tmax_nbody = 0.4 * (Nacc/1e5) * (eps_nbody/0.05)
        tmax_gyrs  = tbar/1000 * tmax_nbody
        kmax = int(-np.floor(np.log2(tmax_gyrs)))
        
        print("Using kmax:", kmax)
        
        mkhalo_params = get_mkhalo_params(r200[indx], m200[indx], c200[indx], Nacc, 
                                          str(indx).zfill(6))
        print(mkhalo_params)

        mkhalo_command = mkhalo_base % mkhalo_params
        gyr_command = gyr_command_base % (f'{str(indx).zfill(6)}', gyr_outfile,
                                          tfinal, epsilon, tstep, kmax, 
                                          str(indx).zfill(6), str(zfin)) 

        jobfile.write(mkhalo_command + '\n')
        jobfile.write(gyr_command + '\n')
        params.write((' '+'{:11} '*2+'{:11f} '*(len(param_tup)-3)+'{} '+'\n')\
                     .format(i+1,indx,zinfall[indx],cosmo.age(zinfall[indx])/u.Gyr,tfinal,tstep,
                             r200[indx], 0, 0, 0, 0, 0, 0, m200[indx]))

        # Make sure we have the .dat files we need
        jobfile.write('for f in *z*.nemo; do s2a $f $f".dat"; done\n')

        # Dirty way to keep track of things
        jobfile.write("rm running && touch done" + "\n"*2) 

        # Close files
        jobfile.close()
        params.close()

        # Make job commands executable 
        st = os.stat(f"{dirname}/job_commands.sh")
        os.chmod(f"{dirname}/job_commands.sh", st.st_mode | stat.S_IEXEC)

joblist.close()
