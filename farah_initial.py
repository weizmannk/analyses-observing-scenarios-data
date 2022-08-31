import os
from astropy.table import  Table

import pandas as pd
import numpy as np

from matplotlib.pyplot import figure
from matplotlib import pyplot as plt
import matplotlib
plt.style.use('seaborn-paper')

data_dir = './data/Pretrov_Farah_input'

outdir = "./outdir"

colors = ['r', 'g']


Farah = Table.read(f"{data_dir}/farah.h5")
Farah.sort('mass1')
Petrov = Table.read(f"{data_dir}/petrov.h5")
Farah.sort('mass1')

params = ['mass', 'spin']

ns_max_mass = 3

F_source_mass1 = Farah['mass1']
F_source_mass2 = Farah['mass2']

P_source_mass1 = Petrov['mass1']
P_source_mass2 = Petrov['mass2']

BNS   = [Petrov[(P_source_mass1 < ns_max_mass) & (P_source_mass2 < ns_max_mass)], Farah[(F_source_mass1 < ns_max_mass) & (F_source_mass2 < ns_max_mass)]]

NSBH  =  [Petrov[(P_source_mass1 >= ns_max_mass) & (P_source_mass2 < ns_max_mass)], Farah[(F_source_mass1 >= ns_max_mass) & (F_source_mass2 < ns_max_mass)]]
BBH   = [Petrov[(P_source_mass1 >= ns_max_mass) & (P_source_mass2 >= ns_max_mass)], Farah[(F_source_mass1 >= ns_max_mass) & (F_source_mass2 >= ns_max_mass)]] 

tables = {"BNS": BNS, 
         "NSBH": NSBH, 
         "BBH": BBH
        }
         
plt.clf()
# Figure Plot 
fig, axs = plt.subplots(nrows=3, ncols=2)

for pop in ['BNS', 'NSBH', 'BBH']:
        mass1 = [tables[pop][0]['mass1'], tables[pop][1]['mass1']]
        mass2 = [tables[pop][0]['mass2'], tables[pop][1]['mass2']] 

        if pop == 'BNS':

                #mass1
                axs[0, 0].hist(mass1, bins=200, density =12, histtype='step' , color=colors,   label= [f'Petrov_{pop}', f'Farah_{pop}'], linewidth=1.3)
                axs[0, 0].set_title(r'(mass1)'+ f' {pop}', fontname="Times New Roman", size=13, fontweight="bold")
                axs[0, 0].legend(prop={'size': 10})
          
                #mass2
                axs[0, 1].hist(mass2, bins=200, density =2, histtype='step' , color=colors,   label= [f'Petrov_{pop}', f'Farah_{pop}'], linewidth=1.3)
                axs[0, 1].legend(prop={'size': 10})
                axs[0, 1].set_title(r'(mass2)'+ f' {pop}', fontname="Times New Roman", size=13, fontweight="bold") 

        elif pop == 'NSBH':
   
                #mass1
                axs[1, 0].hist(mass1, bins=200, density =12, histtype='step' , color=colors,   label= [f'Petrov_{pop}', f'Farah_{pop}'], linewidth=1.3)
                axs[1, 0].legend(prop={'size': 10})
                axs[1, 0].set_title(r'(mass1)'+ f' {pop}', fontname="Times New Roman", size=13, fontweight="bold")
 
                #mass2
                axs[1, 1].hist(mass2, bins=200, density =2, histtype='step' , color=colors,   label= [f'Petrov_{pop}', f'Farah_{pop}'], linewidth=1.3)
                axs[1, 1].legend(prop={'size': 10})
                axs[1, 1].set_title(r'(mass2)'+ f' {pop}', fontname="Times New Roman", size=13, fontweight="bold") 

        else:
                   
                #mass1
                axs[2, 0].hist(mass1, bins=200, density =12, histtype='step' , color=colors,   label= [f'Petrov_{pop}', f'Farah_{pop}'], linewidth=1.3)
                axs[2, 0].legend(prop={'size': 10})
                axs[2, 0].set_title(r'(mass1)'+ f' {pop}', fontname="Times New Roman", size=13, fontweight="bold")

                #mass2
                axs[2, 1].hist(mass2, bins=200, density =2, histtype='step' , color=colors,   label= [f'Petrov_{pop}', f'Farah_{pop}'], linewidth=1.3)
                axs[2, 1].legend(prop={'size': 10})
                axs[2, 1].set_title(r'(mass2)'+ f' {pop}', fontname="Times New Roman", size=13, fontweight="bold") 


plt.gcf().set_size_inches(12, 12)
plt.subplots_adjust(left=0.1,
        bottom=0.1,
        right=0.9,
        top=0.9,
        wspace=0.4,
        hspace=0.4)
fig.tight_layout()
plt.savefig(f'{outdir}/hist_Petrov_Farah_initial.png')