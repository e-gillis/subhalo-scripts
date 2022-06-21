#!/usr/local/bin/python3

# Set the right backend
import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import clustertools as ct
import galpy as gal
import astropy.units as u
from glob import glob
import os

from astropy.cosmology import FlatLambdaCDM

#(1) Cosmological parameters (Planck 2018)
OmegaM = 0.3153
OmegaL = 0.6847
sigma_8 = 0.8111
h0 = 0.6736
cosmo = FlatLambdaCDM(H0=h0*100.0,Om0=OmegaM)
t0 = float(cosmo.age(0)/u.Gyr)



def get_rhoprof(cluster, samples, zinfall, log_r=False):
    cluster.to_centre()
    hz = cosmo.H(zinfall).value/100
    rvir = ct.analysis.functions.virial_radius(cluster,
           method="critical_density", 
           Om=OmegaM, H=100*hz)
    
    if log_r:
        r_values = np.logspace(0, rvir, samples + 1)
    else:
        r_values = np.linspace(0, rvir, samples + 1)
    
    m_shell = np.zeros(r_values.shape)[1:]
    
    for i in range(len(m_shell)):
        r_slice = np.logical_and(cluster.r >= r_values[i], cluster.r < r_values[i+1])
        m_shell[i] = np.sum(cluster.m[r_slice])
        
    r_diff = r_values[1:] - r_values[:-1]
    r_mean = np.mean(np.vstack((r_values[1:], r_values[:-1])), axis=0)
    rhoprof = m_shell / (4*np.pi * r_mean**2 * r_diff)
     
    return rhoprof / np.sum(cluster.m[cluster.r < rvir]), r_mean / rvir
    

def rho_array(globstr, samples, xlog=False):
    files = glob(globstr)
    rho_array = np.zeros((len(files), samples))
    r_values = None
    
    for i in range(len(files)):
        print(f"Extracting {files[i]}")
        param_glob = files[i][:files[i].rfind("/")] + "/*params.dat"
        param_name = glob(param_glob)[0]
        params = np.loadtxt(param_name)
        if len(params.shape) > 1: zinfall = params[0][2]
        else:                     zinfall = params[2]
        print(f"Infall Redshift {zinfall}")
        
        cluster = ct.load_cluster(ctype='gyrfalcon', filename=files[i])
        
        rho_array[i], r_values = get_rhoprof(cluster, samples, zinfall)
    
    return rho_array, r_values


def log_scale_data(mean, std):
    log_mean = np.log10(mean)
    log_error = std / (np.abs(mean) * np.log(10))
    
    return log_mean, log_error


def extract_data(directory, redshift, samples=30):
    data_array, r_values = rho_array(directory +\
                                     f"/*_*/z0{redshift}_final_frame.dat", 
                                     samples=samples)
    mean_density = np.mean(data_array, axis=0)
    std_density = np.mean(data_array, axis=0)
    return r_values, mean_density, std_density

# Want to make one plot per mass bin
# Different colors for DM and Baryonic

mass_bins = ["m06", "m07"]

for mass_bin in mass_bins:
    mass_decade = int(mass_bin[-2:])
    print(mass_bin)
    
    bary_dir = "bary/peri100/" + mass_bin
    dm_dir = "dm/peri100/" + mass_bin
    
    print(dm_dir)
    
    bary_c, dm_c = 'r', 'k'
    offset = 0.003
    
    r_values, z05_dm_mean_density, z05_dm_std_density = extract_data(dm_dir, 5)
    r_values, z05_bary_mean_density, z05_bary_std_density = extract_data(bary_dir, 5)
    r_values, z07_dm_mean_density, z07_dm_std_density = extract_data(dm_dir, 7)
    r_values, z07_bary_mean_density, z07_bary_std_density = extract_data(bary_dir, 7)
    
    z07_dm_rho, z07_dm_err = log_scale_data(z07_dm_mean_density, z07_dm_std_density)
    z05_dm_rho, z05_dm_err = log_scale_data(z05_dm_mean_density, z05_dm_std_density)
    z07_bary_rho, z07_bary_err = log_scale_data(z07_bary_mean_density, z07_bary_std_density)
    z05_bary_rho, z05_bary_err = log_scale_data(z05_bary_mean_density, z05_bary_std_density)
    
    np.save("z07_dm_rho.npy", z07_dm_rho)
    np.save("z07_dm_err.npy", z07_dm_err)
    np.save("z05_dm_rho.npy", z05_dm_rho)
    np.save("z05_dm_err.npy", z05_dm_err)
    np.save("z07_bary_rho.npy", z07_bary_rho)
    np.save("z07_bary_err.npy", z07_bary_err)
    np.save("z05_bary_rho.npy", z05_bary_rho)
    np.save("z05_bary_err.npy", z05_bary_err)
    
    
    plt.figure(figsize=(8, 6))
    plt.errorbar(r_values-offset, z07_dm_rho, z07_dm_err,
                 capsize=3, marker='.',
                 label=f"Redshift 0.7 dark matter only",
                 color=dm_c, fmt="--", alpha=0.5)
    plt.errorbar(r_values, z05_dm_rho, z05_dm_err,
                 capsize=3, marker='.',
                 label=f"Redshift 0.5 dark matter only",
                 color=dm_c)
    plt.errorbar(r_values+offset, z07_bary_rho, z07_bary_err, 
                 capsize=3, marker='.',
                 label=f"Redshift 0.7 baryonic",
                 color=bary_c, fmt="--", alpha=0.5)
    plt.errorbar(r_values+2*offset, z05_bary_rho, z05_bary_err,
                 capsize=3, marker='.',
                 label=f"Redshift 0.5 baryonic",
                 color=bary_c)
    plt.xlabel("$r / r_{vir}$")
    plt.ylabel("Normalized Density (kpc$^{-3}$)")
    plt.title(f"$10^{{{mass_decade}}}$ Mass Bin")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{mass_bin}_mean_density_profiles.pdf")
    plt.savefig(f"{mass_bin}_mean_density_profiles.png")
