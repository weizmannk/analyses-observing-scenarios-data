# This code allow as to plot the density of populations distributions 
# From Farah population (External distribution ) and bayestar_inject , 
# internal distribution. Here Farah have run as a one population name "farah",
# and farah cover the try populations (BBH, BNS and NSBH) inside farah pop.
# the Observing scenarios so here we have 4 polpulations 1 from Farah,
# and 3 from bayestar_inject internal distribution.
# We put the 3 polpulations from  bayestar_inject,
# internal distribution together as a same Great Population
# Here we apply the upper-lower limit about 3 in the internal populations 
# And for Farah (External distribution we apply the split in the post run).

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
from scipy.stats import gaussian_kde

import pandas as pd
import numpy as np


from matplotlib.pyplot import figure
from matplotlib import pyplot as plt
import matplotlib
from matplotlib import pyplot as plt
plt.style.use('seaborn-paper')


path_dir = './data/Petrov-Farah-post-runs/runs'

outdir = "./paper_plots"

if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
distributions = ['Petrov', 'Farah']

run_names = run_dirs= ['O3', 'O4', 'O5']
pops      =  ['BNS', 'NSBH', 'BBH'] 

ns_max_mass = 3

# Read in a table  all populations, FARAH AND BAYESTAR_INJECT INTERNAL INJECTION PROCESS 
tables = {}
with tqdm(total=len(run_dirs) * len(pops)) as progress:
    for run_name, run_dir in zip(run_names, run_dirs):
        tables[run_name]={}
        for dist in distributions:
            if dist == 'Petrov':
                tables[run_name][dist] ={}
                for pop in pops: 
                    path = Path(f'{path_dir}/{run_dir}/{pop.lower()}_astro')
                    injections = Table.read(str(path/'injections.dat'), format='ascii')
                    injections.rename_column('simulation_id', 'event_id')
                    
                    tables[run_name][dist][pop] = injections
                    
                    print(f'{run_name} {dist} {pop} : {len(injections)}') 
                    
                    del injections
            else:
                path = Path(f'{path_dir}/{run_dir}/farah')
                injections = Table.read(str(path/'injections.dat'), format='ascii')
                injections.rename_column('simulation_id', 'event_id')
                
                tables[run_name][dist]= injections
                
                print(f'{run_name} {dist} : {len(injections)}')

params = ['mass', 'distance']

with tqdm(total=len(run_names)) as progress:
    for run_name in run_names:
        # Figure Plot 
        plt.clf()
        fig, axs = plt.subplots(nrows=2, ncols=2)
        
        #farah popuation
        Farah = Table(tables[run_name]['Farah'])
        
        # bayestar initial data
        bns_astro   = Table(tables[run_name]['Petrov']['BNS'])
        nsbh_astro  = Table(tables[run_name]['Petrov']['NSBH'])
        bbh_astro   = Table(tables[run_name]['Petrov']['BBH'])
        
         #join all the internal  populations together with vstack, from astropy.table 
        try:
            Petrov =   vstack([bns_astro, nsbh_astro, bbh_astro], join_type='exact')   
            
        except TableMergeError as ex:
            print(ex)    
        else:
            for dist in distributions: 
                if dist == 'Petrov':
                    for param in params:
                        if (dist == 'Petrov') & (param=='mass'):
                            mass1    = np.log10(Petrov['mass1'])
                            mass2    = np.log10(Petrov['mass2'])
                            xy       = np.vstack([mass1 , mass2])
                            
                            z        = gaussian_kde(xy)(xy)
                            index    = z.argsort()
                            mass1, mass2, z = mass1[index], mass2[index], z[index]
                            
                            axs[0, 0].scatter(mass1, mass2, c=z, s=25)
                            
                            #axs[0, 0].set_title(f'{run_name} {dist} ', fontname="Times New Roman", size=13,  fontweight="bold")
                            axs[0, 0].set_xlabel(r'$\log_{10}$ (mass1)', fontname="Times New Roman", size=18, fontweight="bold")
                            axs[0, 0].set_ylabel(r'$\log_{10}$ (mass2)', fontname="Times New Roman", size=18, fontweight="bold")
                            axs[0, 0].text(0.05, 0.95, f'LRR {run_name}', transform=axs[0,0].transAxes,  color='navy', va='top', fontname="Times New Roman", size=18, fontweight="bold")
                            
                        else:
                            distance    = np.log10(Petrov['distance'])
                            mass1       = np.log10(Petrov['mass1'])
                            xy          = np.vstack([distance , mass1])
                                                                 
                            z           = gaussian_kde(xy)(xy)
                            index       = z.argsort()
                            distance, mass2, z = distance[index], mass1[index], z[index]
                                                                 
                            axs[1, 0].scatter(distance, mass1, c=z, s=25)
                                                                 
                            #axs[1, 0].set_title(f'{run_name} {dist} ', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[1, 0].set_xlabel(r'$\log_{10}$ (distance)', fontname="Times New Roman", size=18, fontweight="bold")
                            axs[1, 0].set_ylabel(r'$\log_{10}$ (mass1)', fontname="Times New Roman", size=18, fontweight="bold")
                            axs[1, 0].text(0.05, 0.95, f'LRR {run_name}', transform=axs[1,0].transAxes, color ='navy',  va='top',fontname="Times New Roman", size=18, fontweight="bold")
                            
                else:
                    for param in params:
                        if  (dist == 'Farah') & (param=='mass'):
                            mass1    = np.log10(Farah['mass1'])
                            mass2    = np.log10(Farah['mass2'])
                            xy    = np.vstack([mass1 , mass2])
                            
                            z     = gaussian_kde(xy)(xy)
                            index = z.argsort()
                            mass1, mass2, z = mass1[index], mass2[index], z[index]
                            
                            axs[0, 1].scatter(mass1, mass2, c=z, s=25)
                            
                            #axs[0, 1].set_title(f'{run_name} {dist}', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[0, 1].set_xlabel(r'$\log_{10}$ (mass1)', fontname="Times New Roman", size=18, fontweight="bold")
                            axs[0, 1].set_ylabel(r'$\log_{10}$ (mass2)', fontname="Times New Roman", size=18, fontweight="bold")
                            axs[0, 1].text(0.05, 0.95, f'GWTC-3 {run_name}', transform=axs[0,1].transAxes, color ='navy', va='top', fontname="Times New Roman", size=18, fontweight="bold") 
                        else:
                            distance    = np.log10(Farah['distance'])
                            mass1       = np.log10(Farah['mass1'])
                            xy          = np.vstack([distance , mass1])
                                                                 
                            z     = gaussian_kde(xy)(xy)
                            index = z.argsort()
                            distance, mass2, z = distance[index], mass1[index], z[index]
                                                                 
                            axs[1, 1].scatter(distance, mass1, c=z, s=25)
                                                                 
                            axs[1, 1].set_xlabel(r'$\log_{10}$ (distance)', fontname="Times New Roman", size=18, fontweight="bold")
                            axs[1, 1].set_ylabel(r'$\log_{10}$ (mass1)',    fontname="Times New Roman", size=18, fontweight="bold")
                            axs[1, 1].text(0.05, 0.95, f'GWTC-3 {run_name}', transform=axs[1,1].transAxes, color ='navy',  va='top', fontname="Times New Roman", size=18, fontweight="bold")

        plt.gcf().set_size_inches(12, 12)
        plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
        fig.tight_layout()
        plt.savefig(f'{outdir}/log10_gaussian_kde_distribution_des_masses_{run_name}.png')
        #plt.savefig(f'{outdir}/log10_gaussian_kde_Post-Split_O4_{run_name}.pdf')
        #plt.savefig(f'{outdir}/log10_gaussian_kde_Post-Split_O4_{run_name}.svg')
        plt.close()
        progress.update()
        
