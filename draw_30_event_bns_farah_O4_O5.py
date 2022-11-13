import os
from tqdm.auto import tqdm
from pathlib import Path
import shutil
from astropy.table import join, Table, vstack
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo, z_at_value
import numpy as np


path_dir = '/home/wkiendrebeogo/Projets/LVK-collaboration/GihubStock/analyses-observing-scenarios-data/data/Petrov-Farah-post-runs' #./Farah/'

# the distribution flders
distribution = ['Farah']

pops = ['BNS']
run_names = run_dirs=  ['O4', 'O5']


draw_number = 30

for dist in distribution:
    outdir = f'outdir/drawn_{draw_number}_CBC/runs'
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

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
            
            table['mass1'] = source_mass1
            table['mass2'] = source_mass2


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

            
            
                data.sort('distance')
                
                #Remove the masses > 2.43 for EOS
                NS_EOS_limit = 2.43
                 
                data = data[data['mass1']< NS_EOS_limit]
                
                data[:30].write(Path(f"{pop_dir}/{run_name}_real_mass_injections.dat"), format='ascii.tab', overwrite=True)
                
                ### detecteur masses with remoe real mass > 2.43
                
                # Split by source frame mass
                
                data_detector = data
                z1 = z_at_value(cosmo.luminosity_distance, data_detector['distance'] * u.Mpc).to_value(u.dimensionless_unscaled)
                zp2 = z1 + 1 

                data_detector['mass1']  = data_detector['mass1']*zp2
                data_detector['mass2']  = data_detector['mass2']*zp2
                
                data_detector[:30].write(Path(f"{pop_dir}/{run_name}_detector_mass_injections.dat"), format='ascii.tab', overwrite=True) 
                print(f'{pop} real masses {len(data)} and detector masses {len(data_detector)}  ; ') 
                 
            del injections, table, source_mass1, source_mass2, z, zp1, zp2, z1, data,  data_detector 
               
            
            
               
               
               
               
               
               
               
               
               
               
               
               
               
               
               
               
"""               
                ## draw an 10 random index for GW+EM stuff
                np.random.seed(150914) 
                Data_draw= np.random.choice(data, draw_number)
                
                ## re-put draw data into astropy table        
                Data = Table(Data_draw)
                Data.write(Path(f"{pop_dir}/michael_injections.dat"), format='ascii.tab', overwrite=True)
                
                ### draw using weizmann and Abby method
                ## We draw  split data in 10 groups of 100 events then draw one 
                ## in each group 
                  
                data.sort('mass1')               
                
                event =[] 
                indice = 0
                while indice  < len(data):
                    if indice == 900:
                        indice +=4                  
                    draw = np.random.choice(data[indice:int(indice+100)], 1) 
                    event.append(draw) 
                    
                    indice += 100
                
                EVENT = vstack([Table(event[i]) for i in range(len(event))], join_type= 'exact')
                
                EVENT.write(Path(f"{pop_dir}/abby_injections.dat"), format='ascii.tab', overwrite=True) 
                
                print(f'{pop} {len(data)} ; ')
                
                print(" ")
                
                print(f"we drawn {draw_number} {pop} from {len(data)} in {run_name} ")
                print(" ")
                print(("draw 10 envent  from 1004"))
                print(" ") 
                print(Data)
                
                print('===============================================================')
                print('Create 10 groups of 100 events then draw ')
                print(" ")
                print(EVENT) 

            del injections, table, source_mass1, source_mass2, z, zp1,  data, Data, EVENT
            print("***************************************************************")

"""
