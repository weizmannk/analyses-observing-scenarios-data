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
plt.style.use('seaborn-paper')


path_dir = './outdir/farah_process.csv'

outdir = "./outdir"

if not os.path.isdir(outdir):
    os.makedirs(outdir)

initial = { 'O3':{'BNS': 892762, 'NSBH': 35962, 'BBH': 71276}, 'O4':{'BNS': 892762, 'NSBH': 35962, 'BBH': 71276}, 'O5': {'BNS': 892762, 'NSBH': 35962, 'BBH': 71276}}

detection = {'O3':{'BNS': 460, 'NSBH': 79, 'BBH': 4891}, 'O4':{'BNS': 1004, 'NSBH': 184, 'BBH': 7070}, 'O5':{'BNS': 2003, 'NSBH': 356, 'BBH': 9809}}

resampling = {'O3':{'BNS': 390079, 'NSBH': 49177, 'BBH': 560744} , 'O4':{'BNS': 587016, 'NSBH': 60357, 'BBH': 352627} , 'O5':{'BNS': 768739, 'NSBH': 54642, 'BBH': 176619}
}

steps = [initial, resampling, detection]
run_names = ['O3', 'O4'] #, 'O5']

#Plot Histogram farah and bayestar-intect internal data,
colors = ['r', 'g', 'b']

plt.clf()
# Figure Plot 
fig, axs = plt.subplots(nrows=3, ncols=2) 

for run_name in run_names:

    data = [steps[0][run_name], steps[1][run_name], steps[2][run_name]]

    BNS = [[data[0]['BNS']], [data[1]['BNS']], [data[2]['BNS']]]
    NSBH = [[data[0]['NSBH']], [data[1]['NSBH']], [data[2]['NSBH']]]
    BBH = [[data[0]['BBH']], [data[1]['BBH']], [data[2]['BBH']]] 
   
    if run_name == 'O3':
        for pop in ['BNS', 'NSBH', 'BBH']:
            if pop == 'BNS': 
                #pop
                axs[0, 0].hist(BNS,   histtype='bar',  color=colors, density=True,  label= [f'initial_{pop}', f'resamp_{pop}', f'detection_{pop}'])
                axs[0, 0].legend(prop={'size': 10})
                axs[0, 0].set_title( f'{run_name} - {pop}', fontname="Times New Roman", size=13, fontweight="bold")

            elif pop == 'NSBH':
            #NSBH
                axs[1, 0].hist(NSBH,  histtype='bar', color=colors,  density=True,  label= [f'initial_{pop}', f'resamp_{pop}', f'detection_{pop}'])
                axs[1, 0].set_title( f'{run_name} - {pop}', fontname="Times New Roman", size=13, fontweight="bold")
            else:
            #BBH
                axs[2, 0].hist(BBH,  histtype='bar', color=colors,  density=True,  label= [f'initial_{pop}', f'resamp_{pop}', f'detection_{pop}'])
                axs[2, 0].set_title( f'{run_name} - {pop}', fontname="Times New Roman", size=13, fontweight="bold")
        
    else:
        for pop in ['BNS', 'NSBH', 'BBH']:
            if pop == 'BNS': 
                #pop
                axs[0, 1].hist(BNS,   histtype='bar',  color=colors, density=True,  label= [f'initial_{pop}', f'resamp_{pop}', f'detection_{pop}'])
                axs[0, 1].legend(prop={'size': 10}) 
                axs[0, 1].set_title( f'{run_name} - {pop}', fontname="Times New Roman", size=13, fontweight="bold")

            elif pop == 'NSBH':
            #NSBH
                axs[1, 1].hist(NSBH,  histtype='bar', color=colors,  density=True,  label= [f'initial_{pop}', f'resamp_{pop}', f'detection_{pop}'])
                axs[1, 1].set_title( f'{run_name} - {pop}', fontname="Times New Roman", size=13, fontweight="bold")
                axs[1, 1].set_title( f'{run_name} - {pop}', fontname="Times New Roman", size=13, fontweight="bold")
                
            else:
            #BBH
                axs[2, 1].hist(BBH,  histtype='bar', color=colors,  density=True,  label= [f'initial_{pop}', f'resamp_{pop}', f'detection_{pop}'])
                axs[2, 1].set_title( f'{run_name} - {pop}', fontname="Times New Roman", size=13, fontweight="bold")
            
            
plt.gcf().set_size_inches(12, 12)
plt.subplots_adjust(left=0.1,
        bottom=0.1,
        right=0.9,
        top=0.9,
        wspace=0.4,
        hspace=0.4)
fig.tight_layout()
plt.savefig(f'{outdir}/subpop_number_farah.png')

                