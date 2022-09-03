import os
from tqdm.auto import tqdm
from pathlib import Path
from astropy.table import  Table 
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo, z_at_value

import pandas as pd
import numpy as np


# Read Ligo xml file generate by ligo.skymap
path_dir_resamp = "./data/inital-distribution-in-bayestar-inject_xml"
#path_dir_detect = "./data/Post-Split-UpperLowerLimit/runs"
path_dir_detect   =  "./data/runs"

outdir = "./outdir"

if not os.path.isdir(outdir):
    os.makedirs(outdir)
 
ns_max_mass  = 3.0

run_names =  ['O3', 'O4', 'O5']
steps = ['initial', 'resampling', 'detection']
 
tables = {}
with tqdm(total=len(steps)*len(run_names)) as progress:
    for step in steps:
         
        tables[step] = {}  
        for run_name in run_names:
            if step == 'initial':
                
                path = Path('./farah/farah.h5')
                            
                injections = Table.read(path)
                table = injections
                
                source_mass1 = table['mass1']
                source_mass2 = table['mass2'] 
            
            else: 
                if  step  == 'resampling':
                    path = Path(f'{path_dir_resamp}/injections_xml_{run_name}.h5')
                    injections = Table.read(path) 
                else:
                    path = Path(f'{path_dir_detect}/{run_name}/farah/injections.dat')
                    injections = Table.read(path, format ="ascii")
                    injections.rename_column('simulation_id', 'event_id')
                
                table = injections 
            
                # Split by source frame mass
                z = z_at_value(cosmo.luminosity_distance, table['distance'] * u.Mpc).to_value(u.dimensionless_unscaled)
                zp1 = z + 1
                
                source_mass1 = table['mass1']/zp1
                source_mass2 = table['mass2']/zp1

            
            table['mass1'] = source_mass1
            table['mass2'] = source_mass2
            
            tables[step][run_name] = {} 
            tables[step][run_name]['BNS']   =  len(table[(source_mass1 < ns_max_mass)  & (source_mass2 < ns_max_mass)])
            tables[step][run_name]['NSBH']  =  len(table[(source_mass1 >= ns_max_mass) & (source_mass2 < ns_max_mass)])
            tables[step][run_name]['BBH']   =  len(table[(source_mass1 >= ns_max_mass) & (source_mass2 >= ns_max_mass)])

            del injections, table, source_mass1, source_mass2
            progress.update()
        

df =pd.DataFrame(data=tables) 
df.to_csv(f"{outdir}/farah_process.csv")


with open(f'{outdir}/farah_process.tex', 'w') as f:
    f.write(df.to_latex(
        longtable=True, 
        caption=' Farah process, number of each subpopulation from inital sample  \\ \
        to event detected via reampling process (xml file).', 
        label='tab:farah-distribution')
    )
