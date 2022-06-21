#!/usr/local/bin/python3

import numpy as np
import sympy as sym

from galpy.potential import NFWPotential,HernquistPotential,rtide
from galpy.util import conversion

# To make directories and scripts
import os
import stat

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--manual", action="store_true", help="Use manual values in setup.py")
parser.add_argument("--nstart", type=int, help="Simulation number to start initializing at")
parser.add_argument("--nmake", type=int, help="Number of simulations to initialize")
parser.add_argument("--mass-bin", type=int, help="Mass bin to sample subhalos from")
parser.add_argument("--baryonic", action="store_true", help="Set to use baryonic potential")
parser.add_argument("--cfactor", default=1, type=int, help="Factor to multiply concentration by")
args = parser.parse_args()


if args.manual:
    # Halo number to start from
    nstart = 500
    # Halos to extract
    nmake = 3
    # Mass range to sample
    mass_range = [1e6, 1e13]
    # Pericenter Range
    peri_range = [0, 100]
    # Use baryonic potential?
    baryonic = True
    # Concentration Factor
    cfactor = 1
else:
    nstart = args.nstart
    nmake = args.nmake
    mass_range = (1e6, 1e7, 1e8, 1e9, 1e10, 1e13)[args.mass_bin-6:args.mass_bin-4]
    peri_range = [0, 100]
    baryonic = args.baryonic
    cfactor = args.cfactor

print(f"Extracting {nmake} halo(s), mass range {np.log10(np.array(mass_range))}, "+\
      f"pericenter range {peri_range}")
print(f"{['Not u', 'U'][baryonic]}sing baryonic potential")


# Galpy Settings
ro = 8.
vo = 220.
mo = conversion.mass_in_msol(ro=ro,vo=vo)


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
t07 = cosmo.age(0.7)/u.Gyr
t05 = cosmo.age(0.5)/u.Gyr
t0 = cosmo.age(0.)/u.Gyr
tstep_restart = (t05-t07)/10.

# cosmo.H(0.5)*u.Mpc*u.s/u.km
data_z_0_5 = np.load('M_1.0e13_z_0.5_Mres_1.0e6.npz')
hhost=cosmo.H(0.5).value/100


# Hostpotential
hhost = cosmo.H(0.5).value / 100
hostpot = NFWPotential(mvir=data_z_0_5['MHost'][0]/1e12,
                       conc=data_z_0_5['cHost'][0],
                       H=hhost*100.0, Om=OmegaM, wrtcrit=True,
                       overdens=delta_c, ro=ro, vo=vo)

if not baryonic:
    # This is easy
    acc_pars = hostpot.nemo_accpars(ro=ro,vo=vo)
    acc_name = hostpot.nemo_accname()

else:
    rs_hern = hostpot.a * ro
    q_hern = 0.06 # will need to check if this is reasonable
    f_hern = 5 # will need to check if this is reasonable
    a_hern = q_hern*rs_hern # Definition
    Ms_nfw = hostpot.mass(a_hern * u.kpc)
    hernquist_norm = 4*f_hern*Ms_nfw
    
    # Units are needed here for Galpy to work
    hernpot = HernquistPotential(amp=2*hernquist_norm*u.Msun, a=a_hern*u.kpc, 
                                 ro=ro, vo=vo)
    
    host_pars, hern_pars = hostpot.nemo_accpars(ro=ro, vo=vo), \
                           hernpot.nemo_accpars(ro=ro, vo=vo)
    host_name, hern_name = hostpot.nemo_accname(), hernpot.nemo_accname()
    
    # Concatenation that gyrfalcON wants
    acc_pars = host_pars+'#'+hern_pars
    acc_name = host_name+'+'+hern_name
    
print("\nHalo Parameters:")
print(acc_pars)
print(acc_name)


# This is where we get our subhalos
m200, r200, c200 = data_z_0_5['M200c'], data_z_0_5['R200c']*1000.0, data_z_0_5['c200c']
zinfall = data_z_0_5['zInfall']
tinfall = cosmo.age(zinfall) / u.Gyr
x, y, z = data_z_0_5['posX']*1000.0, data_z_0_5['posY']*1000.0, data_z_0_5['posZ']*1000.0
vx, vy, vz = data_z_0_5['velX'], data_z_0_5['velY'], data_z_0_5['velZ']
r = np.sqrt(x**2 + y**2 + z**2) 

c200 = c200 * cfactor


# Subhalo Script initialization Routines

# Randomize the galacticus dataset
np.random.seed(42)
idsubs = np.arange(len(m200))
np.random.shuffle(idsubs)
sim_nums=np.arange(1, len(idsubs)+1)

# Load in orbital parameters
peri, apo, eccen, energy = np.loadtxt("orbital_params.dat").T[2:]

# Save idsubs
np.savetxt('idsubs.dat', np.vstack((sim_nums, idsubs)).T, fmt='%06i')

        
# Add master list of Jobs
joblist = open("jobs_to_run.dat", "a")
finlist = open("finished.dat", "w")
finlist.close()

# Cut all our data arrays
mass_shuffled = m200[idsubs]

# Cut for m200
mass_cut = np.logical_and(mass_shuffled >= mass_range[0], mass_shuffled < mass_range[1])
# Cut for pericenters
peri_cut = np.logical_and(peri >= peri_range[0], peri < peri_range[1])
# Now make sure we get the start right 
index_cut = sim_nums > nstart
# Combine Everything
valid_indeces = np.argwhere(np.logical_and(mass_cut, index_cut)).flatten()[:nmake]
print(f"Valid simulation number set: {valid_indeces + 1}\n")


# Commands to format in the job script
mkhalo_base = "mkhalo out=%s.nemo nbody=%i model=NFW r_s=%f M=%f r_t=-%f WD_units=t"
snapshift_base = "snapshift in=%s.nemo out=%s.nemo rshift=%f,%f,%f vshift=%f,%f,%f"
gyr_command_base = "gyrfalcON in=%s.nemo out=%s.nemo tstop=%f eps=%f kernel=2 "\
                + f"step=%f kmax=%i Nlev=1 fac=0.01 accname={acc_name} "\
                + f"accpars={acc_pars} > %s.out_z0%s 2>&1"

for i in valid_indeces:
    indx = idsubs[i]
    sim_num = sim_nums[i]
    
    # Make Directory
    dirname = "{:06}_{:06}".format(sim_num, indx)
    print(f'Initialising subhalo {dirname}:')
    
    try:
        os.mkdir(dirname)
    except FileExistsError:
        print(f"directory {dirname} exixts")
        continue 
    
    # Set Up Subhalo Parameters
    hz = cosmo.H(zinfall[indx]).value/100
    
    print("Subhalo Parameters: Mass log, concentration, radius, infall z")
    print(np.log10(m200[indx]), c200[indx], r200[indx], hz)
    
    
    fname = f"{dirname}/{dirname}"
    
    # Open a Parameters File
    params = open(f"{dirname}/{str(sim_num).zfill(6)}_sim_params.dat", "w")
    param_tup = ("job num", "job id", "start z", "sim start", 
                 "evo time", "tstep", "radius", "x", "y", "z", 
                 "vx", "vy", "vz","mass")
    
    # Some Fancy Formating to save the information
    params.write("#"+(("%11s "*len(param_tup)) % param_tup) + "\n")
    
    # Add job to list
    joblist.write(dirname + "\n")

    # Will satisfy Van den Bosch criteria
    epsilon = f(c200[indx])*r200[indx]/(600 * c200[indx]) *\
              (10/c200[indx])**0.63
    Nacc = int(1000.0*(0.01**(-1./0.8))/0.32)
    
    
    # Open Job file
    jobfile = open(f"{dirname}/job_commands.sh", "w")
    jobfile.write(r"#!/usr/bin/bash" + "\n"*2)
    jobfile.write("touch running\n")
    
    
    # See if it needs restarting or not
    if zinfall[indx] <= 0.7:
        tfinal = t05-cosmo.age(zinfall[indx])/u.Gyr
        gyr_outfile = f'{str(indx).zfill(6)}_z05'
        zfin = 5
        
    else:
        tfinal = t07-cosmo.age(zinfall[indx])/u.Gyr
        gyr_outfile = f'{str(indx).zfill(6)}_z07'
        zfin = 7
    
    tstep=tfinal/10
    
    
    # Find optimal timestep
    eps_nbody = epsilon / r200[indx] # Unitless
    rbar  = r200[indx] * 1000        # In pc
    
    # Copied from clustertools
    vbar = 0.06557 * np.sqrt(m200[indx] / rbar)
    tbar = rbar / (1.023*vbar)

    tmax_nbody = 0.4 * (Nacc/1e5) * (eps_nbody/0.05)
    tmax_gyrs  = tbar/1000 * tmax_nbody
    kmax = -np.floor(np.log2(tmax_gyrs))
    kmax = int(kmax)

    print("Using kmax:", kmax)

    mkhalo_params = get_mkhalo_params(r200[indx], m200[indx], c200[indx], Nacc, 
                                      str(indx).zfill(6)+"_mkhalo")
    mkhalo_command = mkhalo_base % mkhalo_params
    snapshift_command = snapshift_base % (mkhalo_params[0], f'{str(indx).zfill(6)}',
                                          x[indx],  y[indx],  z[indx], 
                                          vx[indx], vy[indx], vz[indx])
    gyr_command = gyr_command_base % (f'{str(indx).zfill(6)}',gyr_outfile,
                                      tfinal,epsilon,tstep,kmax,str(indx).zfill(6), str(zfin))
    
    jobfile.write(mkhalo_command + '\n')
    jobfile.write(snapshift_command + '\n')
    jobfile.write(gyr_command + '\n')
    
    params.write((' '+'{:11} '*2+'{:11f} '*(len(param_tup)-3)+'{} '+'\n')\
                 .format(i+1,indx,zinfall[indx],cosmo.age(zinfall[indx])/u.Gyr,tfinal,tstep,
                         r200[indx],x[indx],y[indx],z[indx],
                         vx[indx],vy[indx],vz[indx],m200[indx]))

    # Add restart commands if needed
    if zinfall[indx] > 0.7:
        tfinal_restart = t05 - cosmo.age(zinfall[indx]) / u.Gyr
        gyr_restart_infile = f'{str(indx).zfill(6)}_restart'
        gyr_restart_outfile = f'{str(indx).zfill(6)}_z05'
        restart_a2s='snaptrim in=%s.nemo out=%s_restart.nemo times=#11' % \
                    (gyr_outfile,str(indx).zfill(6))
        restart_command = gyr_command_base % (gyr_restart_infile, gyr_restart_outfile, 
                                              tfinal_restart, epsilon, tstep_restart, 
                                              kmax,str(indx).zfill(6), str(5))
        jobfile.write(restart_a2s + '\n')
        jobfile.write(restart_command + '\n')
        params.write((' '+'{:11} '*2+'{:11f} '*(len(param_tup)-3)+'{} '+'\n')\
                     .format(i+1,indx,0.7,cosmo.age(0.7)/u.Gyr,tfinal_restart,tstep_restart,
                             r200[indx],x[indx],y[indx],z[indx],vx[indx],vy[indx],vz[indx],
                             m200[indx]))
    
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
    print("Jobfile Written\n")
    
joblist.close()
