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
    

run_names = run_dirs= ['O3', 'O4', 'O5']
pops      =  ['BNS', 'NSBH', 'BBH', 'Farah']

ns_max_mass = 3

# Read in a table  all populations, FARAH AND BAYESTAR_INJECT INTERNAL INJECTION PROCESS 
tables = {}
with tqdm(total=len(run_dirs) * len(pops)) as progress:
    for run_name, run_dir in zip(run_names, run_dirs):
        tables[run_name] = {}
        for pop in pops:
            #path = Path(f'{path_dir}/{run_dir}/{pop}') 
            if pop.lower() == 'farah':
                #path_dir =  '/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/Leo/runs'
                path = Path(f'{path_dir}/{run_dir}/farah')
            else:
                #path_dir = '/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/New-sim-lower-upper-limit-3/runs'
                path = Path(f'{path_dir}/{run_dir}/{pop.lower()}_astro')
                
            injections = Table.read(str(path/'injections.dat'), format='ascii')
            injections.rename_column('simulation_id', 'event_id')

            table = injections
            
            # Split by source frame mass
            z = z_at_value(cosmo.luminosity_distance, table['distance'] * u.Mpc).to_value(u.dimensionless_unscaled)
            zp1 = z + 1

            tables[run_name][pop] = {}
            source_mass1 = table['mass1']/zp1
            source_mass2 = table['mass2']/zp1
            
            table['mass1'] = source_mass1
            table['mass2'] = source_mass2
            
            # split  farah data in BNS, NSBH , BBH, and put them back together
            if pop.lower() == 'farah':
                BNS   =  table[(source_mass1 < ns_max_mass) & (source_mass2 < ns_max_mass)]
                NSBH  =  table[(source_mass1 >= ns_max_mass) & (source_mass2 < ns_max_mass)]
                BBH   =  table[(source_mass1 >= ns_max_mass) & (source_mass2 >= ns_max_mass)]
                
                #join the all populations masses together with vstack, from astropy.table
                try:
                    farah =   vstack([BNS, NSBH, BBH], join_type='exact')
                
                except TableMergeError as ex:
                    print(ex)    
                else: 
                    tables[run_name][pop]['mass1']    = farah['mass1']
                    tables[run_name][pop]['mass2']    = farah['mass2']
                    tables[run_name][pop]['distance'] = farah['distance']
                    
                    del BNS, NSBH, BBH, farah

            else:
                #tables[run_name][pop]['mass1']    =  source_mass1
                #tables[run_name][pop]['mass2']    = source_mass2
                #tables[run_name][pop]['distance'] = table['distance']
                tables[run_name][pop]['mass1']    =  table['mass1']
                tables[run_name][pop]['mass2']    =  table['mass2']
                tables[run_name][pop]['distance'] =  table['distance']  
                         
            del injections, table, source_mass1, source_mass2, z, zp1
            progress.update()

distributions = ['Petrov', 'Farah']
params = ['mass', 'distance']

with tqdm(total=len(run_names)) as progress:
    for run_name in run_names:
        # Figure Plot 
        plt.clf()
        fig, axs = plt.subplots(nrows=2, ncols=2)
        
        #farah popuation
        Farah = Table(tables[run_name]['Farah'])
        
        # bayestar initial data
        bns_astro   = Table(tables[run_name]['BNS'])
        nsbh_astro  = Table(tables[run_name]['NSBH'])
        bbh_astro   = Table(tables[run_name]['BBH'])
        
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
                            axs[0, 0].set_xlabel(r'$\log_{10}$ (mass1)', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[0, 0].set_ylabel(r'$\log_{10}$ (mass2)', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[0, 0].text(0.05, 0.95, f'{dist} {run_name}', transform=axs[0,0].transAxes,  color='blue', va='top', fontname="Times New Roman", size=13, fontweight="bold")
                            
                        else:
                            distance    = np.log10(Petrov['distance'])
                            mass1       = np.log10(Petrov['mass1'])
                            xy          = np.vstack([distance , mass1])
                                                                 
                            z           = gaussian_kde(xy)(xy)
                            index       = z.argsort()
                            distance, mass2, z = distance[index], mass1[index], z[index]
                                                                 
                            axs[1, 0].scatter(distance, mass1, c=z, s=25)
                                                                 
                            #axs[1, 0].set_title(f'{run_name} {dist} ', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[1, 0].set_xlabel(r'$\log_{10}$ (distance)', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[1, 0].set_ylabel(r'$\log_{10}$ (mass1)', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[1, 0].text(0.05, 0.95, f'{dist} {run_name}', transform=axs[1,0].transAxes, color ='blue',  va='top',fontname="Times New Roman", size=13, fontweight="bold")
                            
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
                            axs[0, 1].set_xlabel(r'$\log_{10}$ (mass1)', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[0, 1].set_ylabel(r'$\log_{10}$ (mass2)', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[0, 1].text(0.05, 0.95, f'{dist} {run_name}', transform=axs[0,1].transAxes, color ='blue', va='top', fontname="Times New Roman", size=13, fontweight="bold") 
                        else:
                            distance    = np.log10(Farah['distance'])
                            mass1       = np.log10(Farah['mass1'])
                            xy    = np.vstack([distance , mass1])
                                                                 
                            z     = gaussian_kde(xy)(xy)
                            index = z.argsort()
                            distance, mass2, z = distance[index], mass1[index], z[index]
                                                                 
                            axs[1, 1].scatter(distance, mass1, c=z, s=25)
                                                                 
                            #axs[1, 1].set_title(f'{run_name}  {dist} ', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[1, 1].set_xlabel(r'$\log_{10}$ (distance)', fontname="Times New Roman", size=13, fontweight="bold")
                            axs[1, 1].set_ylabel(r'$\log_{10}$ (mass1)',    fontname="Times New Roman", size=13, fontweight="bold")
                            axs[1, 1].text(0.05, 0.95, f'{dist} {run_name}', transform=axs[1,1].transAxes, color ='blue',  va='top', fontname="Times New Roman", size=13, fontweight="bold")

        plt.gcf().set_size_inches(12, 12)
        plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
        fig.tight_layout()
        plt.savefig(f'{outdir}/log10_gaussian_kde_Farah_Petrov_{run_name}.png')
        #plt.savefig(f'{outdir}/log10_gaussian_kde_Post-Split_O4_{run_name}.pdf')
        #plt.savefig(f'{outdir}/log10_gaussian_kde_Post-Split_O4_{run_name}.svg')
        plt.close()
        progress.update()
        
