import os
from tqdm.auto import tqdm
from pathlib import Path
import shutil
from astropy.table import join, Table, Column, vstack
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo, z_at_value
import numpy as np

from scipy.stats import gaussian_kde
from matplotlib.pyplot import figure
from matplotlib import pyplot as plt
plt.style.use('seaborn-paper')

path_dir = '/home/wkiendrebeogo/Projets/LVK-collaboration/GihubStock/analyses-observing-scenarios-data/data/Petrov-Farah-post-runs' #./Farah/'

# the distribution flders
distribution = ['Farah']

pops = ['BNS']
run_names = run_dirs=  ['O5']

draw_number = 10


outdir = '/home/wkiendrebeogo/Projets/LVK-collaboration/GihubStock/analyses-observing-scenarios-data/paper_plot'
if not os.path.isdir(outdir):
    os.makedirs(outdir)

for dist in distribution:
    if dist == 'Farah':

        # For splitting into BNS, NSBH, and BBH populations
        ns_max_mass = 3

        for run_name, run_dir in zip(tqdm(run_names), run_dirs):

            path = Path(f'{path_dir}/runs/{run_dir}/farah')
            injections = Table.read(str(path/'injections.dat'), format='ascii.fast_tab')

            table = injections

            # Split by source frame mass
            z = z_at_value(cosmo.luminosity_distance, table['distance'] * u.Mpc).to_value(u.dimensionless_unscaled)
            zp1 = z + 1

            source_mass1 = table['mass1']/zp1
            source_mass2 = table['mass2']/zp1

            print("===============================================================")
            print(f'The number of subpopulation in {run_name} is : ')

            for pop in pops:
                pop_dir = Path(f"{outdir}/{run_dir}/{pop.lower()}_farah")
                if not os.path.isdir(pop_dir):
                    os.makedirs(pop_dir)

                if pop == 'BNS':
                    data = table[(source_mass1 < ns_max_mass) & (source_mass2 < ns_max_mass)]
                elif pop == 'NSBH':
                    data= table[(source_mass1 >= ns_max_mass) & (source_mass2 < ns_max_mass)]
                else:
                    if pop== 'BBH':
                        data = table[(source_mass1 >= ns_max_mass) & (source_mass2 >= ns_max_mass)]

                
                mass_tot = Column(data['mass1']+data['mass2'], name = 'massTot')
                data.add_column(mass_tot)

                ## Figure gaussian_kde distribution for BNS 
                
                plt.clf()
                
                # Figure Plot 
                fig, axs = plt.subplots(nrows=1, ncols=1)

                # farah data
                bns_farah   = data
                
                massTot    = np.log10(bns_farah['massTot'])
                distance   = np.log10(bns_farah['distance'])
                
                xy    = np.vstack([distance , massTot])
                                                        
                z     = gaussian_kde(xy)(xy)
                index = z.argsort()
                distance, massTot, z = distance[index], massTot[index], z[index]
                                                        
                axs.scatter(distance, massTot, c=z, s=25)
                
                axs.set_xlabel(r'$\log_{10}$ (distance)', fontname="Times New Roman", size=13, fontweight="bold")
                axs.set_ylabel(r'$\log_{10}$ (Total mass)',    fontname="Times New Roman", size=13, fontweight="bold")
      
                axs.text(0.05, 0.95, f'Farah  {run_name}', transform=axs.transAxes, color ='blue', va='top', fontname="Times New Roman", size=13, fontweight="bold")
                axs.text(0.05, 0.9, f"{pop}", transform=axs.transAxes, color ='g', va='top', fontname="Times New Roman", size=13, fontweight="bold")
                
                plt.gcf().set_size_inches(6, 6)
                plt.subplots_adjust(left=0.1,
                            bottom=0.1,
                            right=0.9,
                            top=0.9,
                            wspace=0.4,
                            hspace=0.4)
                fig.tight_layout()
                plt.savefig(f'{outdir}/log10_gaussian_kde_bns_farah_mass_distance_{run_name}.png')
                plt.close()
                
                
                print(f'{pop} {len(data)} ; ')
                
                print(" ") 

                del injections, table, source_mass1, source_mass2, z, zp1,  data, mass_tot, bns_farah 
                print("***************************************************************")
