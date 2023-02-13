import os
from tqdm.auto import tqdm

from astropy.table import join, Table, vstack 
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo, z_at_value
from scipy.stats import gaussian_kde

import pandas as pd
import numpy as np

from matplotlib.pyplot import figure
from matplotlib import pyplot as plt
import matplotlib
plt.style.use('seaborn-paper')

data_dir = './data/farah_input/farah.h5'

outdir = "./outdir"


Farah = Table.read(f"{data_dir}")
Farah.sort('mass1')

params = ['mass', 'spin']

plt.clf()
fig, axs = plt.subplots(nrows=2, ncols=1) 
#farah popuation

with tqdm(total=len(params)) as progress:
    for param in params:
        if (param=='mass'):
            mass1    = np.log10(Farah['mass1'][::20])
            mass2    = np.log10(Farah['mass2'][::20])
            xy       = np.vstack([mass1 , mass2])
            
            z        = gaussian_kde(xy)(xy)
            index    = z.argsort()
            mass1, mass2, z = mass1[index], mass2[index], z[index]
            
            axs[0, 0].scatter(mass1, mass2, c=z, s=25)
            
            axs[0, 0].set_title(f'Farah initial {param} distribution', fontname="Times New Roman", size=13,  
                                fontweight="bold")
            axs[0, 0].set_xlabel(r'$\log_{10}$ (mass1)', fontname="Times New Roman", size=13, fontweight="bold")
            axs[0, 0].set_ylabel(r'$\log_{10}$ (mass2)', fontname="Times New Roman", size=13, fontweight="bold")
            
        else:

            spin1z   = np.log10(Farah['spin1z'][::20])
            spin2z   = np.log10(Farah['spin2z'][::20])
            xy    = np.vstack([spin1z , spin2z][::20])
                                                    
            z     = gaussian_kde(xy)(xy)
            index = z.argsort()
            spin1z, spin2z, z = spin1z[index], spin2z[index], z[index]
                                                    
            axs[0, 1].scatter(spin1z, spin2z, c=z, s=25)
                                                    
            #axs[0, 1].set_title(f'{param} ', fontname="Times New Roman", size=13, fontweight="bold")
            axs[0, 1].set_xlabel(r'$\log_{10}$ (spin1z)', f ontname="Times New Roman", size=13,fontweight="bold")
            axs[0, 1].set_ylabel(r'$\log_{10}$ (spin2z)', fontname="Times New Roman", size=13, fontweight="bold")
            progress.update()

plt.gcf().set_size_inches(6, 12)
plt.subplots_adjust(left=0.1,
            bottom=0.1,
            right=0.9,
            top=0.9,
            wspace=0.4,
            hspace=0.4)
fig.tight_layout()
plt.savefig(f'{outdir}/Farah_initial_distribution.png')
#plt.savefig(f'{outdir}/Farah_initial_distribution.pdf')
#plt.savefig(f'{outdir}/Farah_initial_distribution.svg')
plt.close()

       
        
# print the number of each population in Farah distribution
print(f"number of BNS: {len(Farah[(Farah['mass1']<3)])}, number of NSBH:"
f"{len(Farah[(Farah['mass1']>=3) & (Farah['mass2']<3)])}, number of BBH: "
f"{len(Farah[(Farah['mass2']>=3)])}")

