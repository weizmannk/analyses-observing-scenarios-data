# This code allow as to plot the density of populations distributions 
# From Farah population (External distribution ) and bayestar_inject , 
# internal distribution. Here Farah data was split split 
# in bbh_astro, farah_bns, farah_nsbh, farah_bbh,  bebore  running,
# the Observing scenarios so here we have 6 polpulations 3 from farah,
# and 3 from bayestar_inject internal distribution.
# We put them together  ie , the 3 polpulations from  bayestar_inject,
# internal distribution together as a same Great Population and 
# the 3 3 polpulations from Frah (External distribution)
# And we plot the density of the ditribution using gaussian_kde
# These plots concern masses and distances.
# mass1, and mass2
# distance and mass1

import os
from tqdm.auto import tqdm
from pathlib import Path
from astropy.table import join, Table, vstack 
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo, z_at_value

import pandas as pd
import numpy as np

from scipy.stats import gaussian_kde
from matplotlib.pyplot import figure
from matplotlib import pyplot as plt
plt.style.use('seaborn-paper')


path_dir = './Farah/runs'


outdir = "./paper_plots"

if not os.path.isdir(outdir):
    os.makedirs(outdir)
    

run_names = run_dirs=  ['O4']
pops =['BNS']

# Neutron star mass limit
ns_max_mass = 3

# Read in a table  all populations, FARAH AND BAYESTAR_INJECT INTERNAL INJECTION PROCESS 
tables = {}
with tqdm(total=len(run_dirs) * len(pops)) as progress:
    for run_name, run_dir in zip(run_names, run_dirs):
        tables[run_name] = {}
        for pop in pops:

            path = Path(f'{path_dir}/{run_dir}/farah')
            injections = Table.read(str(path/'injections.dat'), format='ascii.fast_tab')
            injections.rename_column('simulation_id', 'event_id')

            table = injections
   
            # Split by source frame mass
            z = z_at_value(cosmo.luminosity_distance, table['distance'] * u.Mpc).to_value(u.dimensionless_unscaled)
            zp1 = z + 1

            source_mass1 = table['mass1']/zp1
            source_mass2 = table['mass2']/zp1

            if pop == "BNS":
                table_pop =table[(source_mass1 < ns_max_mass) & (source_mass2 < ns_max_mass)] 
            
            tables[run_name][pop] = {}             
            tables[run_name][pop]['mass1'] = table_pop['mass1']
            tables[run_name][pop]['mass2'] = table_pop['mass2']
            tables[run_name][pop]['distance'] = table_pop['distance']

            del injections, table, table_pop
            progress.update()

            
distributions = ['Farah']
params = ['mass']

with tqdm(total=len(run_names)) as progress:
    for run_name in run_names:
        plt.clf()
        
        # Figure Plot 
        fig, axs = plt.subplots(nrows=1, ncols=2)

        # farah data
        bns_farah   = Table(tables[run_name]['BNS'])
        
        m1    = np.log10(bns_farah['mass1'])
        kde_x = gaussian_kde(m1)  
        
        kde_x_space = np.linspace(m1.min(), m1.max(), 100)
        kde_x_eval = kde_x.evaluate(kde_x_space)
        kde_x_eval /= kde_x_eval.sum()
        log10 = r'$\log_{10}$ (mass1)'


        axs[0].plot(kde_x_space, kde_x_eval, 'k.')
        axs[0].set_xlabel(f'{log10}', fontname="Times New Roman", size=13, fontweight="bold")
        axs[0].text(0.05, 0.95, f"Farah  {run_name}", transform=axs[0].transAxes, color ='blue', va='top', fontname="Times New Roman", size=13, fontweight="bold")
        axs[0].text(0.05, 0.9, f"{pop}", transform=axs[0].transAxes, color ='g', va='top', fontname="Times New Roman", size=13, fontweight="bold")
        axs[0].text(0.7, 0.95, f'1D {log10}', transform=axs[0].transAxes, color ='k', va='top', fontname="Times New Roman", size=13, fontweight="bold")
         
        mass1    = np.log10(bns_farah['mass1'])
        mass2    = np.log10(bns_farah['mass2'])
        xy       = np.vstack([mass1 , mass2])
        
        z     = gaussian_kde(xy)(xy)
        index = z.argsort()
        mass1, mass2, z = mass1[index], mass2[index], z[index]
        
        axs[1].scatter(mass1, mass2, c=z, s=25)
        
        axs[1].set_xlabel(r'$\log_{10}$ (mass1)', fontname="Times New Roman", size=13, fontweight="bold")
        axs[1].set_ylabel(r'$\log_{10}$ (mass2)', fontname="Times New Roman", size=13, fontweight="bold")
        axs[1].text(0.05, 0.95, f'Farah  {run_name}', transform=axs[1].transAxes, color ='blue', va='top', fontname="Times New Roman", size=13, fontweight="bold")
        axs[1].text(0.05, 0.9, f"{pop}", transform=axs[1].transAxes, color ='g', va='top', fontname="Times New Roman", size=13, fontweight="bold")
    
    
    plt.gcf().set_size_inches(12, 6)
    plt.subplots_adjust(left=0.1,
                bottom=0.1,
                right=0.9,
                top=0.9,
                wspace=0.4,
                hspace=0.4)
    fig.tight_layout()
    plt.savefig(f'{outdir}/log10_gaussian_kde_bns_farah_{run_name}.png')
    plt.close()
    progress.update() 
