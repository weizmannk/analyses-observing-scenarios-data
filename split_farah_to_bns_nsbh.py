import os
from tqdm.auto import tqdm
from pathlib import Path
from astropy.table import join, Table
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo, z_at_value

import pandas as pd
import numpy as np


path_dir = '/home/wkiendrebeogo/Projets/LVK-collaboration/GitLabStock/observing-scenarios-simulations/runs'

outdir = './data/upper-lower-limit/runs'


pops = ['BNS', 'NSBH', 'BBH']
run_names = run_dirs=  ['O3', 'O4', 'O5']

# For splitting into BNS, NSBH, and BBH populations
ns_max_mass = 3

for run_name, run_dir in zip(tqdm(run_names), run_dirs):
           
    path = Path(f'{path_dir}/{run_dir}/farah')
    injections = Table.read(str(path/'injections.dat'), format='ascii.fast_tab')
     
    table = injections
    
    # Split by source frame mass
    z = z_at_value(cosmo.luminosity_distance, table['distance'] * u.Mpc).to_value(u.dimensionless_unscaled)
    zp1 = z + 1
    

    source_mass1 = table['mass1']/zp1
    source_mass2 = table['mass2']/zp1
    
    for pop in pops:
        pop_dir = Path(f"{outdir}/{run_dir}/farah_{pop.lower()}")
        if not os.path.isdir(pop_dir):
            os.makedirs(pop_dir)
            
        if pop == 'BNS': 
            data = table[(source_mass1 < ns_max_mass) & (source_mass2 < ns_max_mass)] 
        elif pop == 'NSBH':
            data= table[(source_mass1 >= ns_max_mass) & (source_mass2 < ns_max_mass)]
        else:
           data = table[(source_mass1 >= ns_max_mass) & (source_mass2 >= ns_max_mass)]
        
        data.write(Path(f"{pop_dir}/injections.dat"), format='ascii.tab', overwrite=True)
         
    del injections, table, source_mass1, source_mass2, z, zp1,  data
