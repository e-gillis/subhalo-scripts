#!/usr/local/bin/python3

VERSION = "1.1.0"


import matplotlib 
if __name__ == "__main__":
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import patches
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

import numpy as np
import clustertools as ct

import galpy as gal
from galpy.potential import NFWPotential, rtide, HernquistPotential
from galpy.util import coords, conversion
from galpy.orbit import Orbit,Orbits

ro = 8.
vo = 220.
mo = conversion.mass_in_msol(ro=ro,vo=vo)

import astropy.units as u
from glob import glob

from astropy.cosmology import FlatLambdaCDM

import argparse
import pickle

#(1) Cosmological parameters (Planck 2018)
OmegaM = 0.3153
OmegaL = 0.6847
sigma_8 = 0.8111
h0 = 0.6736
cosmo = FlatLambdaCDM(H0=h0*100.0,Om0=OmegaM)
t07 = float(cosmo.age(0.7)/u.Gyr)
t05 = float(cosmo.age(0.5)/u.Gyr)
t0 = float(cosmo.age(0)/u.Gyr)




class Processed_Simulation:
    """A collection of processed clusters with relevant sim information
    
    === Attributes ===
    directory: str
        The directory name of the simulation
    halos: List[Processed_Subhalo]
        A list of Processed Halos (Processed_Halo class) in chronological order
    lookback_times: 1D numpy array
        Array of lookback times in Gyrs corresponding to the halos
    m_ratios: 1D numpy array
        Original subhalo mass fraction enclosed by the infall virial radius
    halo_no: str
        The simulation number
    halo_id: str
        The id corresponding to the subhalo in the galacticus database
    start_sim: float
        The time since the start of the universe corresponding to the subhalo's
        infall in Gyrs
    evo_time: float
        The interation time of the subhalo in Gyrs
    halo_mass: float
        The mass of the subhalo
    pericentre: float
        The pericentre of the orbit
    closest_r: float
        The closest the centre
    """
    
    def __init__(self, directory, sim_type, verbose=False, fix_rvir=False, 
                 rho_bin_num=40):
        # Save if Baryonic or DM
        self.sim_type = sim_type
        # Set the directory of the simulation
        self.directory = directory
        # Make a list of halos to populate
        self.halos = []
        # Saving times like this is also useful, but possibly redundant
        self.lookback_times = []
        # Mass ratio here
        self.m_ratios = []
        # Bound Fraction
        self.bound_fraction = []
        
        # Import Parameters
        param_name = glob(f"{directory}*params.dat")[0]
        params = np.loadtxt(param_name)
        
        fixed_rvir = None
        ignore_idx = None
        
        if len(params.shape) > 1:
            self.halo_no = str(int(params[0, 0])).zfill(6)
            self.halo_id = str(int(params[0, 1])).zfill(6)
            self.start_sim = params[0][3]
            self.evo_time = params[-1][4]
            self.halo_mass = params[0][-1]
            self.zinfall = params[0][2]
            cluster_07_name = glob(f'{directory}*z07.nemo.dat')[0]
            
            if verbose:
                print(f"Opening Cluster {cluster_07_name}")
            
            cluster_07 = ct.load_cluster(ctype='gyrfalcon', origin='galaxy',
                                     units='WDunits',
                                         filename=cluster_07_name)

            for i in range(11):
                cluster_07.to_kpckms()
                self.halos.append(Processed_Halo(cluster_07, self.zinfall, 
                                                 fixed_rvir, ignore_idx, 
                                                 rho_bin_num))
                self.lookback_times.append(self.halos[-1].lookback_time)
                self.m_ratios.append(self.halos[-1].m_encl)
                if fix_rvir:
                    fixed_rvir = self.halos[0].rvir
                ignore_idx = self.halos[-1]._ignore_idx
                if i != 10:
                    cluster_07 = ct.advance_cluster(cluster_07)
                else:
                    cluster_07.to_centre()
                    r_sorts = np.argsort(cluster_07.r)
                    r_cutoff = cluster_07.r[r_sorts[len(cluster_07.r)//100]]
                    kin, pot = ct.analysis.energies(cluster_07, i_d=(cluster_07.r < r_cutoff),
                                                    specific=False)
                    bound_array = (kin + pot) < 0
                    bound_frac = np.sum(bound_array) / np.sum(cluster_07.r < r_cutoff)
                    self.bound_fraction.append(bound_frac)

        else:
            self.start_sim = params[3]
            self.evo_time = params[4]
            self.halo_mass = params[-1]
            self.zinfall = params[2]
            self.halo_no = str(int(params[0])).zfill(6)
            self.halo_id = str(int(params[1])).zfill(6)
        
        cluster_05_name = glob(f'{directory}*z05.nemo.dat')[0]
        
        if verbose:
            print(f"Opening Cluster {cluster_05_name}")
        
        cluster_05 = ct.load_cluster(ctype='gyrfalcon', origin='galaxy',
                                     units='WDunits',
                                     filename=cluster_05_name)
        cluster_05.to_kpckms()
        self.halos.append(Processed_Halo(cluster_05, self.zinfall, fixed_rvir, 
                                         ignore_idx, rho_bin_num))
        self.lookback_times.append(self.halos[-1].lookback_time)
        self.m_ratios.append(self.halos[-1].m_encl)
        
        if fix_rvir:
            fixed_rvir = self.halos[0].rvir
        ignore_idx = self.halos[-1]._ignore_idx
        
        for i in range(10):
            if i == 0:
                if self.halos == []:
                    self.halos.append(Processed_Halo(cluster_05, self.zinfall, 
                                     fixed_rvir, ignore_idx, rho_bin_num))
                    self.lookback_times.append(self.halos[-1].lookback_time)
                    cluster_05 = ct.advance_cluster(cluster_05)
                    cluster_05.to_kpckms()
                    
                    if fixed_rvir: 
                        fixed_rvir = self.halos[0].rvir
                
            self.halos.append(Processed_Halo(cluster_05, self.zinfall, 
                                             fixed_rvir, ignore_idx, 
                                             rho_bin_num))
            self.lookback_times.append(self.halos[-1].lookback_time)
            self.m_ratios.append(self.halos[-1].m_encl)


            ignore_idx = self.halos[-1]._ignore_idx
            
            if i != 9:
                cluster_05 = ct.advance_cluster(cluster_05)
                cluster_05.to_kpckms()
            else:    
                cluster_05.to_centre()
                r_sorts = np.argsort(cluster_05.r)
                r_cutoff = cluster_05.r[r_sorts[len(cluster_05.r)//100]]
                kin, pot = ct.analysis.energies(cluster_05, i_d=(cluster_05.r < r_cutoff),
                                                specific=False)
                bound_array = (kin + pot) < 0
                bound_frac = np.sum(bound_array) / np.sum(cluster_05.r < r_cutoff)
                self.bound_fraction.append(bound_frac)
                    
        self.lookback_times = np.array(self.lookback_times)
        self.m_ratios = np.array(self.m_ratios) / self.halo_mass
        self.bound_fraction = np.array(self.bound_fraction)
        
        for halo in self.halos:
            halo._ignore_idx = np.array([])
        
    
    def __str__(self):
        """Return a string representation of the subhalo
        """
        return f"Subhalo {self.halo_no}-{self.halo_id}"
    
    
    def save(self):
        # Set version
        self.version = VERSION
        
        with open(f"{self.directory[:-1]}_{self.sim_type}.halo", "wb") as f:
            pickle.dump(self, f)
        return None
    
    
    def plot_vprof(self, save_plot=False, logx=True, logy=False, norm_r=True, 
                   norm_v=False):
        virr_plot = plt.figure(figsize=(10, 5))
        virr_axes = virr_plot.add_axes([0,0,0.75,1])
        
        for halo in self.halos:
            time_evolved = t0 - halo.lookback_time - self.start_sim
            time_left = self.evo_time - time_evolved

            ctup = (max(min(time_left/self.evo_time, 1), 0), 0, 
                    max(min(time_evolved/self.evo_time, 1), 0))

            if norm_r:   
                plot_r = halo.r_values / self.halos[0].rvir
                rstr = "$r / r_{vir, infall}$"
            else:        
                plot_r = halo.r_values
                rstr = "$r (kpc)$"
                
            if norm_v:
                plot_v = halo.vprof / self.halos[0].vprof
                vstr = "v / v_{infall}"
                
            else:
                plot_v = halo.vprof
                vstr = "km/s"
            
            virr_axes.plot(plot_r, plot_v, 
                           label=f"{round(halo.lookback_time, 3)} Gyrs",
                           color=ctup, linewidth=0.5)
        
        if logy:
            virr_axes.set_yscale('log')
        if logx:
            virr_axes.set_xscale('log')
        virr_axes.legend(bbox_to_anchor=(1,1), loc="upper left")
        virr_axes.set_title(f"{self.directory[:-1]} Velocity Profile"\
                            .replace("_", "-"))
        virr_axes.set_xlabel(rstr)
        virr_axes.set_ylabel(vstr)
        virr_plot.subplots_adjust(right=0.7)
        if save_plot:
            virr_plot.savefig(f"{self.directory}{self.directory[:-1]}_velprof"+
                              ".pdf", bbox_inches="tight")
            virr_plot.savefig(f"{self.directory}{self.directory[:-1]}_velprof"+\
                              ".png", bbox_inches="tight")
            return None
        plt.show()
 
            
    def plot_rhoprof(self, save_plot=False, norm_rho=True, norm_r=True, 
                     logx=False, logy=True):
        density_plot = plt.figure(figsize=(10, 5))
        density_axes = density_plot.add_axes([0,0,0.75,1])

        for halo in self.halos:
            time_evolved = t0 - halo.lookback_time - self.start_sim
            time_left = self.evo_time - time_evolved

            ctup = (max(min(time_left/self.evo_time, 1), 0), 0, 
                    max(min(time_evolved/self.evo_time, 1), 0))

            if norm_r:   
                plot_r = halo.r_values / self.halos[0].rvir
                rstr = "$r / r_{vir, infall}$"
            else:        
                plot_r = halo.r_values
                rstr = "$r (kpc)$"

            if norm_rho: 
                plot_rho = halo.rhoprof/self.halos[0].rhoprof 
                rhostr = "Normalized Density"
            else:        
                plot_rho = halo.rhoprof
                rhostr = "Density"


            density_axes.plot(plot_r, plot_rho,
                              label=f"{round(halo.lookback_time, 3)} Gyrs",
                              color=ctup, linewidth=0.5)  
        
        if logy:
            density_axes.set_yscale('log')
        if logx:
            density_axes.set_xscale('log')
            
        density_axes.legend(bbox_to_anchor=(1,1), loc="upper left")
        density_axes.set_title(f"{self.halo_no}-{self.halo_id} {self.sim_type}: Density Profile")
        density_axes.set_xlabel(rstr)
        density_axes.set_ylabel(rhostr)
        density_plot.subplots_adjust(right=0.7)
        if save_plot:
            density_plot.savefig(f"{self.directory}{self.directory[:-1]}_rhoprof.pdf",
                               bbox_inches="tight")
            density_plot.savefig(f"{self.directory}{self.directory[:-1]}_rhoprof.png",
                               bbox_inches="tight")
            return None
        plt.show()

   
    def plot_massloss(self, save_plot=False):
        m_array = np.zeros(len(self.halos))
        for i in range(len(self.halos)):
            m_array[i] = self.halos[i].m_encl
        m_ratio_array = m_array / self.halo_mass
         
        fig, ax1 = plt.subplots(figsize=(9, 4))
        ax1.plot(self.lookback_times, m_array)
        ax1.set_xlabel("Lookback Time (Gyrs)")
        ax1.set_ylabel("Mass Enclosed")
        ax2 = ax1.twinx()
        ax2.plot(self.lookback_times, m_ratio_array, alpha=0)
        ax2.set_ylabel("Ratio Enclosed")
        plt.gca().invert_xaxis()
        plt.title((f"{self.directory[:-1]} {self.sim_type}: "+\
                  "Virial Enclosed Mass").replace("_", "-"))
        plt.tight_layout()
        if save_plot:
            plt.savefig(f"{self.directory}{self.directory[:-1]}_m_enclosed.pdf", 
                        bbox_inches="tight")
            plt.savefig(f"{self.directory}{self.directory[:-1]}_m_enclosed.png", 
                        bbox_inches="tight")
            return None
        plt.show()

        
    
    def plot_evo(self, save_plot=False):

        evolution_scatter = plt.figure(figsize=(8.5, 5))
        evolution_axes = evolution_scatter.add_axes([0,0,0.75,1])

        for halo in self.halos:
            time_evolved = t0 - halo.lookback_time - self.start_sim
            time_left = self.evo_time - time_evolved

            ctup = (max(min(time_left/self.evo_time, 1), 0), 0, 
                    max(min(time_evolved/self.evo_time, 1), 0))

            halo_circ = patches.Circle((halo.cx, halo.cy), radius=halo.rvir, 
                                       label=f"{round(halo.lookback_time, 3)} Gyrs",
                                       color=ctup)
            evolution_axes.add_patch(halo_circ)
            evolution_axes.arrow(halo.cx, halo.cy, halo.vx/15, halo.vy/15, 
                                 color=ctup, head_width=5)

            scale = AnchoredSizeBar(evolution_axes.transData,
                              size=1000/15,
                              label="1000 km/s",
                              loc=2,
                              pad=0.5, borderpad=1, sep=5,
                              frameon=True)
            evolution_axes.add_artist(scale)

        evolution_axes.legend(bbox_to_anchor=(1,1), loc="upper left")
        evolution_axes.set_title(f"{self.directory[:-1]} {self.sim_type}: Evolution"\
                                 .replace("_", "-"))
        evolution_axes.axis("equal")
        evolution_axes.set_xlabel("kpc")
        evolution_axes.set_ylabel("kpc")
        evolution_axes.scatter([0], [0], marker='+')
        evolution_scatter.subplots_adjust(right=0.7)

        if save_plot:
            evolution_scatter.savefig(f"{self.directory}{self.directory[:-1]}_evolution.png", 
                                      dpi=500,
                                      bbox_inches="tight")
            return None
        plt.show()
        
        
    def set_orbital_params(self, potentials):
        """Set the pericentre and closest approach of the subhalo. Adds two attributes.

        === Parameters ===
        potentials: list of galpy potentials
            List of potentials to evolve the subhalo in

        Adds pericentre and r_closest attributes
        """
        # Load initial coordinates
        x, y, z = self.halos[0].cx, self.halos[0].cy, self.halos[0].cz
        vx, vy, vz = self.halos[0].vx, self.halos[0].vy, self.halos[0].vz

        # Convert initial coordinates
        rad, phi, zed = coords.rect_to_cyl(x, y, z)
        vR, vT, vzed = coords.rect_to_cyl_vec(vx, vy, vz, x, y, z)

        # Make a coordinate array
        vxvv = np.array([rad/ro, vR/vo, vT/vo, zed/ro, vzed/vo, phi])
        orbit = Orbit(vxvv,ro=ro,vo=vo)

        # Integrate to redshift 0.7
        if len(self.halos) > 10:
            times = np.linspace(0, (t07 - self.start_sim)/conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
            orbit.integrate(times, potentials)
            self.z07_closest = orbit.rperi()

        # Integrate to redshift 0.5
        times = np.linspace(0, (t05 - self.start_sim)/conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
        orbit.integrate(times, potentials)
        self.z05_closest = orbit.rperi()

        # 10 Gyr integration
        times = np.linspace(0, 10/conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
        orbit.integrate(times, potentials)
        self.peri = orbit.rperi()

        

class Processed_Halo:
    """A processed subhalo, imported using clustertools
    
    === Attributes ===
    lookback_time: float
        Lookback time of the halo in Gyrs
    centre: numpy array
        x, y and z position of the subhalo in galacticentric coordinates
    vcentre: numpy array
        Velocity of the cluster centre in km/s
    cx, cy, cz: float
        coordinates of the centre in separate variables
    vx, vy, vz: float
        velocity of the centre in separate variables
    r: float
        Distance from galactic centre in kpc
    rvir: float
        Virial Radius of the subhalo, in kpc
    vmax: float
        Maximum circular velocity of dark matter particles
    rvmax: float
        Madius of maximum circlular velocity of dark matter particles
    m_tot:
        Total Subhalo mass, *including extratidal particles
    m_encl: float
        Mass enclosed by the virial radius in Msun
    r_values: 
        Radii at which the profile characteristics are sampled
    rhoprof:
        Density of the subhalo sampled at r_values
    vprof: 
        velocity profile of the subhalo sampled at r
    """
    
    def __init__(self, cluster, zinfall, fixed_rvir=None, _ignore_idx=None, rho_bin_num=40):
        
        # Define hz from tinfall
        hz = cosmo.H(zinfall).value/100
        tinfall = cosmo.age(zinfall).value
    
        # Get centre, set to cluster
        cx, cy, cz, vx, vy, vz = ct.find_centre(cluster, indx=_ignore_idx)
        
        # Get virial radius
        if fixed_rvir is not None:
            rvir = fixed_rvir
            self.rvir_actual = ct.virial_radius(cluster,
                                                method="critical_density", 
                                                indx=_ignore_idx,
                                                Om=OmegaM, H=100*hz)
        else:
            rvir = ct.virial_radius(cluster,
                   method="critical_density", 
                   indx=_ignore_idx,
                   Om=OmegaM, H=100*hz)
            self.rvir_actual = rvir
            
        cluster_r = ((cluster.x - cx)**2 + (cluster.y - cy)**2 + (cluster.z - cz)**2)**0.5
        _ignore_idx = cluster_r < 3*rvir
                
        r_values = np.linspace(0, rvir, rho_bin_num+1)
        
        # Get density profile
        m_shell = np.zeros(r_values.shape)[1:]

        for i in range(len(m_shell)):
            r_slice = np.logical_and(cluster_r >= r_values[i], cluster_r < r_values[i+1])
            m_shell[i] = np.sum(cluster.m[r_slice])

        r_diff = r_values[1:] - r_values[:-1]
        r_mean = np.mean(np.vstack((r_values[1:], r_values[:-1])), axis=0)
        shell_vol = 4 / 3 * np.pi * (r_values[1:]**3 - r_values[:-1]**3)
        rhoprof = m_shell / shell_vol
        
        # Get velocity profile
        rprof, vprof, rvmax, vmax = ct.analysis.profiles.vcirc_prof(cluster)
        vprof = np.interp(r_mean, rprof, vprof)
        
        # Get mass enclosed
        m_encl = np.sum(cluster.m[cluster_r < rvir])
        m_tot = cluster.mtot # Read Mtot in from params    

        # Set Attributes
        self.lookback_time = t0 - tinfall - cluster.tphys 
        self.centre = np.array((cx, cy, cz))
        self.vcentre = np.array((vx, vy, vz))
        self.cx, self.cy, self.cz = (cx, cy, cz) # Can't hurt having this twice
        self.vx, self.vy, self.vz = (vx, vy, vz) # Can't hurt having this twice
        self.r = np.sqrt(np.sum(self.centre**2))
        self.rvir = rvir
        self.rvmax = rvmax
        self.vmax = vmax
        self.m_tot = m_tot
        
        self.m_encl = m_encl
        self.r_values = r_mean
        self.vprof = vprof
        self.rhoprof = rhoprof
        
        self._ignore_idx = _ignore_idx

        
def load_halo(path):
    with open(path, "rb") as f:
        loaded_halo = pickle.load(f)
    return loaded_halo
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("dir", metavar="directory", 
                        type=str, help="Directory to Process")
    parser.add_argument("sim_type", type=str,
                        help="Type of simulation (Baryonic vs DM)")
    parser.add_argument("--silent", action="store_true", 
                        help="Surpress console outputs")
    args = parser.parse_args()
    
    directory = args.dir
    sim_type = args.sim_type
    verbose = not args.silent
 
    print(f"Extracting {sim_type} type halo from {directory}")
    postprocessed_cluster = Processed_Simulation(directory, sim_type, 
                            verbose=verbose, fix_rvir=True)
    
    delta_c = 200 #with respect to the critical density at infall redshift
    data_z_0_5 = np.load('M_1.0e13_z_0.5_Mres_1.0e6.npz')

    # Hostpotential
    hhost = cosmo.H(0.5).value / 100
    hostpot = NFWPotential(mvir=data_z_0_5['MHost'][0]/1e12,
                           conc=data_z_0_5['cHost'][0],
                           H=hhost*100.0,Om=OmegaM,wrtcrit=True,
                           overdens=delta_c,ro=ro,vo=vo)
    
    if sim_type == 'DM':
        potentials = [hostpot]
        postprocessed_cluster.set_orbital_params(potentials)
        
        
    elif sim_type == 'Baryonic':
        rs_hern = hostpot.a * ro
        q_hern = 0.06 # will need to check if this is reasonable
        f_hern = 5 # will need to check if this is reasonable
        a_hern = q_hern*rs_hern # Definition
        Ms_nfw = hostpot.mass(a_hern * u.kpc)
        hernquist_norm = 4*f_hern*Ms_nfw

        # Units are needed here for Galpy to work
        hernpot = HernquistPotential(amp=2*hernquist_norm*u.Msun, a=a_hern*u.kpc, 
                                     ro=ro, vo=vo)
        
        potentials = [hostpot, hernpot]
        postprocessed_cluster.set_orbital_params(potentials)
        
    
    postprocessed_cluster.save()
    
#    postprocessed_cluster.plot_massloss(save_plot=True)
#    postprocessed_cluster.plot_rhoprof(save_plot=True)
#    postprocessed_cluster.plot_evo(save_plot=True)
    
