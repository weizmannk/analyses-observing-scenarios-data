import os
from astropy.table import  Table

import pandas as pd
import numpy as np

from matplotlib.pyplot import figure
from matplotlib import pyplot as plt
import matplotlib
plt.style.use('seaborn-paper')

data_dir = './data/farah_input/farah.h5'

outdir = "./outdir"

colors = ['r', 'g', 'b']

Farah = Table.read(f"{data_dir}")
Farah.sort('mass1')

params = ['mass', 'spin']

ns_max_mass = 3

source_mass1 = Farah['mass1']
source_mass2 = Farah['mass2']

BNS   =  Farah[(source_mass1 < ns_max_mass) & (source_mass2 < ns_max_mass)]
NSBH  =  Farah[(source_mass1 >= ns_max_mass) & (source_mass2 < ns_max_mass)]
BBH   =  Farah[(source_mass1 >= ns_max_mass) & (source_mass2 >= ns_max_mass)]

plt.clf()
# Figure Plot 
fig, axs = plt.subplots(nrows=3, ncols=2) 

mass1    = [BNS['mass1'],    NSBH['mass1'], BBH['mass1']]
mass2    = [BNS['mass2'],    NSBH['mass2'], BBH['mass2']]

spin1z    = [BNS['spin1z'],    NSBH['spin1z'], BBH['spin1z']]
spin2z    =   [BNS['spin2z'],    NSBH['spin2z'], BBH['spin2z']]
       
#mass1
axs[0, 0].hist(mass1, bins=30, density=12, histtype='bar', color=colors, label= ['BNS', 'NSBH', 'BBH'])
axs[0, 1].hist(mass1, bins=300, density =12, histtype='step' , color=colors,   label= ['BNS', 'NSBH', 'BBH'], linewidth=2)
axs[0, 0].legend(prop={'size': 10})
axs[0, 0].set_title(r' (mass1)' + f' Farah', fontname="Times New Roman", size=13, fontweight="bold")
#axs[0, 1].legend(prop={'size': 10})
axs[0, 1].set_title(r'(mass1)' + f' Farah', fontname="Times New Roman", size=13, fontweight="bold")

#mass2
axs[1, 0].hist(mass2, bins=30, density=12, histtype='bar', color=colors,  label= ['BNS', 'NSBH', 'BBH'])
axs[1, 1].hist(mass2, bins=300, density =2, histtype='step' , color=colors,   label= ['BNS', 'NSBH', 'BBH'], linewidth=2)
#axs[1, 0].legend(prop={'size': 10})
axs[1, 0].set_title(r' (mass2)' + f' Farah', fontname="Times New Roman", size=13, fontweight="bold")
axs[1, 1].legend(prop={'size': 10})
axs[1, 1].set_title(r'(mass2)' + f' Farah', fontname="Times New Roman", size=13, fontweight="bold") 

#spinz1
axs[2, 0].hist(spin1z, bins=30, density=12, histtype='bar', color=colors,  label= ['BNS', 'NSBH', 'BBH'])

axs[2, 1].hist(spin1z, bins=300, density =12, histtype='step' , color=colors,   label= ['BNS', 'NSBH', 'BBH'], linewidth=2)
#axs[2, 0].legend(prop={'size': 10})
axs[2, 0].set_title(r' (spin1z)' + f'  Farah', fontname="Times New Roman", size=13, fontweight="bold")
#axs[2, 1].legend(prop={'size': 10})
axs[2, 1].set_title(r'(spin1z)' + f'  Farah', fontname="Times New Roman", size=13, fontweight="bold")

plt.gcf().set_size_inches(12, 12)
plt.subplots_adjust(left=0.1,
        bottom=0.1,
        right=0.9,
        top=0.9,
        wspace=0.4,
        hspace=0.4)
fig.tight_layout()
plt.savefig(f'{outdir}/hist_Farah.png')


                        
